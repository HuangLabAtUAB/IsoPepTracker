#===============================================================================
# DATA PROCESSING FUNCTIONS
#===============================================================================

# Performance optimization: set data.table options for faster processing
if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::setDTthreads(0)  # Use all available cores
}

# Pre-process data for the app - simplified

#===============================================================================
# rMATS-SPECIFIC DATA LOADING FUNCTIONS
#===============================================================================

#' Load transcript exons and CDS data from rMATS GTF file
#' 
#' @param rmats_gtf_file Path to rMATS GTF file
#' @param transcript_ids Vector of transcript IDs to load
#' @return List with success, exons, and cds data
load_rmats_transcript_exons <- function(rmats_gtf_file, transcript_ids) {
  if (!file.exists(rmats_gtf_file)) {
    cat("ERROR: rMATS GTF file not found:", rmats_gtf_file, "\n")
    return(list(success = FALSE, message = paste("rMATS GTF file not found:", rmats_gtf_file)))
  }
  
  cat("🔍 Loading rMATS GTF data from:", rmats_gtf_file, "\n")
  cat("📋 Requested transcript IDs:", paste(transcript_ids, collapse = ", "), "\n")
  
  tryCatch({
    # Import all features from the rMATS GTF
    cat("📥 Importing exon features from rMATS GTF...\n")
    exons <- import(rmats_gtf_file, format = "gtf", feature.type = "exon")
    cat("✅ Found", length(exons), "total exon features\n")
    
    # Import CDS features
    cat("📥 Importing CDS features from rMATS GTF...\n")
    cds <- import(rmats_gtf_file, format = "gtf", feature.type = "CDS")
    cat("✅ Found", length(cds), "total CDS features\n")
    
    # Check for transcript_id attribute
    if (!"transcript_id" %in% names(mcols(exons))) {
      cat("ERROR: rMATS GTF file does not contain transcript_id attribute for exons\n")
      return(list(success = FALSE, message = "rMATS GTF file does not contain transcript_id attribute for exons"))
    }
    
    if (length(cds) > 0 && !"transcript_id" %in% names(mcols(cds))) {
      cat("WARNING: rMATS GTF file does not contain transcript_id attribute for CDS features\n")
      cds <- GRanges() # Empty GRanges if no transcript_id
    }
    
    # Debug: Show available transcript IDs in GTF
    available_exon_ids <- unique(exons$transcript_id)
    cat("📋 Available exon transcript IDs in GTF:", paste(head(available_exon_ids, 10), collapse = ", "))
    if (length(available_exon_ids) > 10) cat(" (showing first 10)")
    cat("\n")
    
    if (length(cds) > 0) {
      available_cds_ids <- unique(cds$transcript_id)
      cat("📋 Available CDS transcript IDs in GTF:", paste(head(available_cds_ids, 10), collapse = ", "))
      if (length(available_cds_ids) > 10) cat(" (showing first 10)")
      cat("\n")
    }
    
    # Extract exons and CDS for each requested transcript
    result <- list(exons = list(), cds = list())
    successful_matches <- 0
    
    for (tx_id in transcript_ids) {
      cat("\n🔍 Processing transcript:", tx_id, "\n")
      
      # rMATS transcripts have quoted gene IDs in GTF like: "ENSG00000099721.15".inclusion
      # Extract the gene ID and suffix to build the correct format
      if (grepl("\\.(inclusion|exclusion)$", tx_id)) {
        # Split transcript ID into gene part and suffix
        parts <- strsplit(tx_id, "\\.")[[1]]
        suffix <- parts[length(parts)]  # "inclusion" or "exclusion"
        gene_part <- paste(parts[-length(parts)], collapse = ".")  # Everything before the suffix
        
        # Build the quoted format that matches the GTF (rtracklayer removes outer quotes)
        gtf_format <- paste0('"', gene_part, '".', suffix)
        
        cat("🔄 Converting rMATS transcript ID:", tx_id, "→", gtf_format, "\n")
        
        # Debug: Show exact comparison
        cat("🔍 Debug matching:", paste0("'", gtf_format, "'"), "against available IDs\n")
        available_matches <- as.character(exons$transcript_id)
        cat("🔍 First few available transcript IDs:", paste(head(available_matches, 3), collapse = ", "), "\n")
        cat("🔍 Exact match test:", gtf_format %in% available_matches, "\n")
        
        # Try exact match with robust GRanges filtering
        tx_exons <- exons[which(as.character(exons$transcript_id) == gtf_format)]
        tx_cds <- cds[which(as.character(cds$transcript_id) == gtf_format)]
        
        # CRITICAL DEBUG: Show filtering results
        cat("🚨 FILTERING RESULT:", length(tx_exons), "exons found for format:", gtf_format, "\n")
        
        # If no exact match, try alternative formats
        if (length(tx_exons) == 0) {
          cat("⚠️  No exact match found, trying alternative formats...\n")
          
          # Try without quotes: ENSG00000099721.15.inclusion
          alt_format1 <- paste0(gene_part, ".", suffix)
          tx_exons <- exons[which(as.character(exons$transcript_id) == alt_format1)]
          tx_cds <- cds[which(as.character(cds$transcript_id) == alt_format1)]
          
          if (length(tx_exons) > 0) {
            cat("✅ Found match with format:", alt_format1, "\n")
          } else {
            # Try with different quote format: 'ENSG00000099721.15'.inclusion
            alt_format2 <- paste0("'", gene_part, "'.", suffix)
            tx_exons <- exons[which(as.character(exons$transcript_id) == alt_format2)]
            tx_cds <- cds[which(as.character(cds$transcript_id) == alt_format2)]
            
            if (length(tx_exons) > 0) {
              cat("✅ Found match with format:", alt_format2, "\n")
            } else {
              cat("❌ No match found with any format\n")
              cat("Available transcript IDs that contain '", suffix, "':\n")
              matching_ids <- available_exon_ids[grepl(suffix, available_exon_ids)]
              if (length(matching_ids) > 0) {
                cat("  ", paste(head(matching_ids, 5), collapse = "\n  "), "\n")
              } else {
                cat("  None found\n")
              }
            }
          }
        }
      } else {
        # Fallback to direct match for non-rMATS transcripts
        cat("🔄 Direct matching for non-rMATS transcript:", tx_id, "\n")
        tx_exons <- exons[which(as.character(exons$transcript_id) == tx_id)]
        tx_cds <- cds[which(as.character(cds$transcript_id) == tx_id)]
      }
      
      # If we found exons, add them to the result
      if (length(tx_exons) > 0) {
        # Sort exons by position
        tx_exons <- tx_exons[order(start(tx_exons))]
        result$exons[[tx_id]] <- tx_exons
        successful_matches <- successful_matches + 1
        
        cat("✅ Found", length(tx_exons), "exons for transcript", tx_id)
        cat(" (", paste(start(tx_exons), "-", end(tx_exons), collapse = ", "), ")\n")
      } else {
        cat("❌ No exons found for transcript", tx_id, "\n")
      }
      
      # Add CDS if found
      if (length(tx_cds) > 0) {
        # Sort CDS by position
        tx_cds <- tx_cds[order(start(tx_cds))]
        result$cds[[tx_id]] <- tx_cds
        
        cat("✅ Found", length(tx_cds), "CDS regions for transcript", tx_id)
        cat(" (", paste(start(tx_cds), "-", end(tx_cds), collapse = ", "), ")\n")
      } else {
        cat("ℹ️  No CDS found for transcript", tx_id, "\n")
      }
    }
    
    cat("\n📊 Summary: Successfully loaded structure data for", successful_matches, "out of", length(transcript_ids), "transcripts\n")
    
    # Check if we successfully loaded any transcript data
    if (successful_matches == 0) {
      cat("❌ FAILURE: No transcript structure data could be loaded from rMATS GTF\n")
      return(list(
        success = FALSE,
        message = paste("No exon/CDS data found for any of the requested transcripts:", paste(transcript_ids, collapse = ", ")),
        exons = list(),
        cds = list()
      ))
    }
    
    cat("✅ SUCCESS: rMATS structure data loading completed\n")
    return(list(
      success = TRUE,
      exons = result$exons,
      cds = result$cds,
      message = paste("Successfully loaded structure data for", successful_matches, "out of", length(transcript_ids), "transcripts from rMATS GTF"),
      loaded_transcripts = successful_matches,
      requested_transcripts = length(transcript_ids)
    ))
    
  }, error = function(e) {
    cat("❌ ERROR in load_rmats_transcript_exons:", e$message, "\n")
    return(list(
      success = FALSE, 
      message = paste("Error loading rMATS GTF:", e$message),
      exons = list(),
      cds = list()
    ))
  })
}
prepare_data <- function(peptides, as_database) {
  # Extract unique genes
  genes <- unique(peptides$geneID)
  gene_symbols <- unique(peptides$geneSymbol)
  gene_lookup <- setNames(gene_symbols, genes)
  
  # Create a map of transcripts per gene
  transcripts_per_gene <- split(peptides$txID, peptides$geneID)
  
  # Create AS event mappings
  as_events <- as_database[, c("refTx", "asTx", "geneID", "eventID", "AS_type", 
                              "AS_range", "AS_range1", "AS_range2")]
  
  # Extract peptides by protease type
  proteases <- c("trp", "chymo", "aspn", "lysc", "lysn", "gluc")
  
  # Return as a list
  return(list(
    genes = genes,
    gene_symbols = gene_symbols,
    gene_lookup = gene_lookup,
    transcripts_per_gene = transcripts_per_gene,
    as_events = as_events,
    proteases = proteases,
    original_peptides = peptides
  ))
}

# Function to load gene details from GTF
load_gene_details <- function(gene_id) {
  if (!file.exists(gtf_file)) {
    return(list(success = FALSE, message = "GTF file not found"))
  }
  
  # Extract gene features to find chromosome and location
  gene_features <- import(gtf_file, format = "gtf", feature.type = "gene")
  
  # Try with and without version number
  gene_id_no_version <- gsub("\\.[0-9]+$", "", gene_id)
  gene_matches <- gene_features[gene_features$gene_id == gene_id | 
                               gene_features$gene_id == gene_id_no_version]
  
  if (length(gene_matches) == 0) {
    return(list(success = FALSE, 
                message = paste("Gene", gene_id, "not found in GTF file")))
  }
  
  # Get chromosome and position
  gene_match <- gene_matches[1]
  chr <- as.character(seqnames(gene_match))
  start_pos <- start(gene_match)
  end_pos <- end(gene_match)
  
  # Return gene details
  return(list(
    success = TRUE,
    gene_id = gene_match$gene_id,
    chromosome = chr,
    start = start_pos,
    end = end_pos
  ))
}


# Function to load exons and CDS for specific transcripts
load_transcript_exons <- function(gene_details, transcript_ids) {
  if (!gene_details$success) {
    return(list(success = FALSE, message = gene_details$message))
  }
  
  # Only load exons and CDS in the gene region (with some padding)
  padding <- 10000
  gene_region <- GRanges(
    seqnames = gene_details$chromosome,
    ranges = IRanges(
      start = max(1, gene_details$start - padding),
      end = gene_details$end + padding
    )
  )
  
  # Import exons
  exons <- import(gtf_file, format = "gtf", feature.type = "exon", which = gene_region)
  
  # Import CDS features
  cds <- import(gtf_file, format = "gtf", feature.type = "CDS", which = gene_region)
  
  # Make sure we have transcript IDs
  if (!"transcript_id" %in% names(mcols(exons))) {
    return(list(success = FALSE, message = "GTF file does not contain transcript_id attribute for exons"))
  }
  
  if (length(cds) > 0 && !"transcript_id" %in% names(mcols(cds))) {
    message("WARNING: GTF file does not contain transcript_id attribute for CDS features")
    cds <- GRanges() # Empty GRanges if no transcript_id
  }
  
  # Extract exons and CDS for each requested transcript
  result <- list(exons = list(), cds = list())
  
  for (tx_id in transcript_ids) {
    # Try exact match with version
    tx_exons <- exons[exons$transcript_id == tx_id]
    tx_cds <- cds[cds$transcript_id == tx_id]
    
    # If that fails, try without version
    if (length(tx_exons) == 0) {
      tx_id_no_version <- gsub("\\.[0-9]+$", "", tx_id)
      tx_exons <- exons[grep(paste0("^", tx_id_no_version, "(\\..*)?$"), exons$transcript_id)]
      tx_cds <- cds[grep(paste0("^", tx_id_no_version, "(\\..*)?$"), cds$transcript_id)]
    }
    
    # If we found exons, add them to the result
    if (length(tx_exons) > 0) {
      # Sort exons by position
      tx_exons <- tx_exons[order(start(tx_exons))]
      result$exons[[tx_id]] <- tx_exons
      
      # Add CDS if found
      if (length(tx_cds) > 0) {
        tx_cds <- tx_cds[order(start(tx_cds))]
        result$cds[[tx_id]] <- tx_cds
      }
    }
  }
  
  # Return results
  if (length(result$exons) == 0) {
    return(list(
      success = FALSE,
      message = "No exons found for any transcripts"
    ))
  }
  
  return(list(
    success = TRUE,
    exons = result$exons,
    cds = result$cds
  ))
}

#===============================================================================
# DECISION TREE BASED ISOFORM DETECTABILITY SCORING
#===============================================================================

# Calculate decision tree-based detectability scores for all isoforms
calculate_isoform_detectability_scores <- function(all_peptides_df, highlight_peptides, 
                                                  all_transcripts, protease, miscleavage_type) {
  
  # Initialize results
  scores_list <- list()
  detailed_breakdown <- list()
  protease_recommendations <- list()
  
  # Process each transcript
  for (tx in all_transcripts) {
    tx_peptides <- all_peptides_df[all_peptides_df$transcript == tx, ]
    
    if (nrow(tx_peptides) == 0) {
      # No peptides found
      scores_list[[tx]] <- list(
        transcript = tx,
        tier = "Tier 4",
        confidence = "Very Low",
        score = 0,
        reasoning = "No peptides found"
      )
      next
    }
    
    # Calculate peptide specificity for this transcript
    unique_peptides <- 0
    subset_peptides <- 0
    widely_shared_peptides <- 0
    universal_peptides <- 0
    
    for (pep in tx_peptides$peptide) {
      # Count other transcripts with this peptide
      other_tx <- unique(all_peptides_df$transcript[
        all_peptides_df$peptide == pep & all_peptides_df$transcript != tx
      ])
      other_count <- length(other_tx)
      total_isoforms <- length(all_transcripts)
      
      if (other_count == 0) {
        unique_peptides <- unique_peptides + 1
      } else if (other_count <= 2) {
        subset_peptides <- subset_peptides + 1
      } else if (other_count < (total_isoforms - 2)) {
        widely_shared_peptides <- widely_shared_peptides + 1
      } else {
        universal_peptides <- universal_peptides + 1
      }
    }
    
    total_peptides <- nrow(tx_peptides)
    
    # Calculate peptide quality metrics
    optimal_length_peptides <- sum(nchar(tx_peptides$peptide) >= 7 & nchar(tx_peptides$peptide) <= 25)
    quality_ratio <- optimal_length_peptides / total_peptides
    
    # Calculate coverage (simplified - based on peptide span)
    if (nrow(tx_peptides) > 0) {
      coverage_span <- max(tx_peptides$end) - min(tx_peptides$start)
      # Estimate transcript length (this is simplified)
      estimated_length <- max(tx_peptides$end) - min(tx_peptides$start) + 1000
      coverage_ratio <- min(1.0, coverage_span / estimated_length)
    } else {
      coverage_ratio <- 0
    }
    
    # Decision Tree Logic
    tier <- "Tier 4"
    confidence <- "Very Low"
    reasoning <- ""
    score <- 0
    
    # Node 1: Unique Peptide Availability
    if (unique_peptides >= 3) {
      # High Confidence Path
      if (quality_ratio >= 0.7 && coverage_ratio >= 0.5) {
        tier <- "Tier 1"
        confidence <- "Very High"
        score <- 90 + min(10, unique_peptides)
        reasoning <- "Multiple unique peptides + excellent quality + high coverage"
      } else if (quality_ratio >= 0.5 || coverage_ratio >= 0.3) {
        tier <- "Tier 2"
        confidence <- "High"
        score <- 75 + min(10, unique_peptides)
        reasoning <- "Multiple unique peptides + good supporting metrics"
      } else {
        tier <- "Tier 2"
        confidence <- "Medium-High"
        score <- 65 + min(10, unique_peptides)
        reasoning <- "Multiple unique peptides but quality/coverage concerns"
      }
    } else if (unique_peptides >= 1) {
      # Medium Confidence Path
      unique_ratio <- (unique_peptides + subset_peptides) / total_peptides
      if (unique_ratio >= 0.6 && quality_ratio >= 0.6) {
        tier <- "Tier 2"
        confidence <- "Medium-High"
        score <- 60 + (unique_peptides * 5)
        reasoning <- "Some unique peptides + high specificity ratio + good quality"
      } else if (unique_ratio >= 0.3 || coverage_ratio >= 0.4) {
        tier <- "Tier 3"
        confidence <- "Medium"
        score <- 45 + (unique_peptides * 5)
        reasoning <- "Some unique peptides with moderate supporting metrics"
      } else {
        tier <- "Tier 3"
        confidence <- "Low-Medium"
        score <- 35 + (unique_peptides * 5)
        reasoning <- "Limited unique peptides with weak supporting metrics"
      }
    } else {
      # Low Confidence Path - No unique peptides
      subset_ratio <- subset_peptides / total_peptides
      if (subset_ratio >= 0.5 && quality_ratio >= 0.7 && coverage_ratio >= 0.4) {
        tier <- "Tier 3"
        confidence <- "Low-Medium"
        score <- 30
        reasoning <- "No unique peptides but excellent subset specificity"
      } else if (subset_ratio >= 0.3 && (quality_ratio >= 0.5 || coverage_ratio >= 0.3)) {
        tier <- "Tier 4"
        confidence <- "Low"
        score <- 20
        reasoning <- "No unique peptides, limited subset specificity"
      } else {
        tier <- "Tier 4"
        confidence <- "Very Low"
        score <- 10
        reasoning <- "Poor detectability across all metrics"
      }
    }
    
    # Store results
    scores_list[[tx]] <- list(
      transcript = tx,
      tier = tier,
      confidence = confidence,
      score = score,
      reasoning = reasoning
    )
    
    # Detailed breakdown
    detailed_breakdown[[tx]] <- list(
      transcript = tx,
      total_peptides = total_peptides,
      unique_peptides = unique_peptides,
      subset_peptides = subset_peptides,
      widely_shared = widely_shared_peptides,
      universal = universal_peptides,
      quality_ratio = round(quality_ratio, 3),
      coverage_ratio = round(coverage_ratio, 3),
      unique_ratio = round((unique_peptides + subset_peptides) / total_peptides, 3)
    )
    
    # Protease recommendations (simplified)
    recommendation <- "Current protease suitable"
    if (score < 50) {
      if (protease == "trp") {
        recommendation <- "Consider Chymotrypsin or LysC for different cleavage patterns"
      } else if (protease == "chymo") {
        recommendation <- "Consider Trypsin or AspN for different specificity"
      } else {
        recommendation <- "Try Trypsin (most common) or evaluate peptide length distribution"
      }
    }
    
    protease_recommendations[[tx]] <- list(
      transcript = tx,
      current_protease = protease,
      recommendation = recommendation,
      score = score
    )
  }
  
  return(list(
    scores = scores_list,
    detailed = detailed_breakdown,
    protease_rec = protease_recommendations
  ))
}