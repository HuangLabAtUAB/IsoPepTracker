# rmats_novel_peptide_generator.R
# Self-contained peptide generator adapted from novel_peptide_generator.R for rMATS events
# Handles all 5 rMATS event types: SE, RI, MXE, A3SS, A5SS

library(stringr)
library(data.table)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(Biostrings)
library(IRanges)
library(cleaver)
library(seqinr)

#===============================================================================
# CONFIGURATION
#===============================================================================

# All enzymes must be processed for complete output structure
enzymes <- c("trp", "chymo", "aspn", "lysc", "lysn", "gluc")

# Enzyme to cleaver mapping
enzyme_map <- list(
  trp = "trypsin",
  chymo = "chymotrypsin-high",
  aspn = "asp-n endopeptidase",
  lysc = "lysc",
  lysn = "lysn",
  gluc = "glutamyl endopeptidase"
)

#===============================================================================
# RMATS-SPECIFIC PROTEIN EXTRACTION
#===============================================================================

# Extract protein sequences from rMATS translation results
extract_rmats_proteins <- function(event_coords, translation_results, phase_results, event_type) {
  cat("Extracting rMATS protein sequences for event type:", event_type, "\n")
  
  # Get gene information - clean up gene_id if it has extra quotes
  gene_symbol <- event_coords$gene_symbol
  gene_id <- gsub('"', '', event_coords$gene_id)  # Remove any extra quotes
  
  # Initialize protein data structure
  protein_data_list <- list()
  
  # Process inclusion isoform
  if (!is.null(translation_results$inclusion_protein) && nchar(translation_results$inclusion_protein) > 0) {
    inclusion_data <- data.table(
      proteinID = paste0(gene_symbol, ".inclusion"),
      txID = paste0(gene_symbol, ".inclusion"), 
      geneID = gene_id,
      geneSymbol = gene_symbol,
      numAA = as.character(nchar(translation_results$inclusion_protein)),
      seq = translation_results$inclusion_protein
    )
    protein_data_list[["inclusion"]] <- inclusion_data
  }
  
  # Process exclusion isoform
  if (!is.null(translation_results$exclusion_protein) && nchar(translation_results$exclusion_protein) > 0) {
    exclusion_data <- data.table(
      proteinID = paste0(gene_symbol, ".exclusion"),
      txID = paste0(gene_symbol, ".exclusion"),
      geneID = gene_id, 
      geneSymbol = gene_symbol,
      numAA = as.character(nchar(translation_results$exclusion_protein)),
      seq = translation_results$exclusion_protein
    )
    protein_data_list[["exclusion"]] <- exclusion_data
  }
  
  # Combine all proteins
  if (length(protein_data_list) > 0) {
    protein_data <- rbindlist(protein_data_list)
  } else {
    stop("No valid protein sequences found for rMATS event")
  }
  
  # Remove stop codons (*) from sequences if present
  protein_data[, seq := gsub("\\*$", "", seq)]
  
  cat("Generated", nrow(protein_data), "protein sequences\n")
  cat("Sequences:", paste(protein_data$proteinID, collapse = ", "), "\n")
  
  return(protein_data)
}

#===============================================================================
# GENOMIC COORDINATE MAPPING FOR RMATS
#===============================================================================

# Create GTF structures from rMATS phase results
create_rmats_gtf_from_phase <- function(phase_results, event_type) {
  cat("Creating GTF structures from rMATS phase results\n")
  
  gtf_list <- list()
  
  # Process inclusion GTF
  if (!is.null(phase_results$inclusion_cds_gtf) && nrow(phase_results$inclusion_cds_gtf) > 0) {
    inclusion_gtf <- as.data.frame(phase_results$inclusion_cds_gtf)
    # Extract gene_id from attributes if needed
    if (!("gene_id" %in% colnames(inclusion_gtf))) {
      # Extract gene_id and clean it properly
      gene_id_match <- regmatches(inclusion_gtf$attributes, 
                                 regexec('gene_id "([^"]+)"', inclusion_gtf$attributes))
      if (length(gene_id_match[[1]]) > 1) {
        inclusion_gtf$gene_id <- gene_id_match[[1]][2]  # Extract the captured group
      } else {
        inclusion_gtf$gene_id <- event_coords$gene_id %||% "unknown"
      }
    }
    inclusion_gtf$transcript_id <- paste0(inclusion_gtf$gene_id[1], ".inclusion")
    inclusion_gtf$type <- "CDS"
    # Ensure consistent data types
    inclusion_gtf$phase <- as.numeric(inclusion_gtf$frame)
    gtf_list[["inclusion"]] <- inclusion_gtf
  }
  
  # Process exclusion GTF  
  if (!is.null(phase_results$exclusion_cds_gtf) && nrow(phase_results$exclusion_cds_gtf) > 0) {
    exclusion_gtf <- as.data.frame(phase_results$exclusion_cds_gtf)
    # Extract gene_id from attributes if needed
    if (!("gene_id" %in% colnames(exclusion_gtf))) {
      # Extract gene_id and clean it properly
      gene_id_match <- regmatches(exclusion_gtf$attributes, 
                                 regexec('gene_id "([^"]+)"', exclusion_gtf$attributes))
      if (length(gene_id_match[[1]]) > 1) {
        exclusion_gtf$gene_id <- gene_id_match[[1]][2]  # Extract the captured group
      } else {
        exclusion_gtf$gene_id <- event_coords$gene_id %||% "unknown"
      }
    }
    exclusion_gtf$transcript_id <- paste0(exclusion_gtf$gene_id[1], ".exclusion")
    exclusion_gtf$type <- "CDS"
    # Ensure consistent data types
    exclusion_gtf$phase <- as.numeric(exclusion_gtf$frame)
    gtf_list[["exclusion"]] <- exclusion_gtf
  }
  
  if (length(gtf_list) == 0) {
    stop("No valid GTF structures could be created from phase results")
  }
  
  # Combine GTF entries using do.call to avoid data.table issues
  combined_gtf <- do.call(rbind, gtf_list)
  
  cat("Created GTF with", nrow(combined_gtf), "CDS entries\n")
  return(combined_gtf)
}

# Convert rMATS GTF to GRanges for genomic mapping
convert_rmats_gtf_to_granges <- function(rmats_gtf) {
  cat("Converting rMATS GTF to GRanges\n")
  
  # Create GRanges object
  gtf_gr <- GRanges(
    seqnames = paste0("chr", rmats_gtf$seqname),
    ranges = IRanges(start = rmats_gtf$start, end = rmats_gtf$end),
    strand = rmats_gtf$strand,
    type = rmats_gtf$type,
    transcript_id = rmats_gtf$transcript_id,
    gene_id = rmats_gtf$gene_id,
    phase = rmats_gtf$phase
  )
  
  cat("Created GRanges with", length(gtf_gr), "entries\n")
  return(gtf_gr)
}

#===============================================================================
# MAIN RMATS GENERATOR FUNCTION
#===============================================================================

generate_rmats_peptide_data <- function(event_coords, translation_results, phase_results, event_type, 
                                       missedCleavages = 0, minLength = 7, maxLength = 60) {
  cat("=== STARTING RMATS PEPTIDE DATA GENERATION ===\n")
  cat("Event type:", event_type, "\n")
  cat("Gene:", event_coords$gene_symbol, "\n\n")
  
  # Step 1: Extract protein sequences from rMATS results
  cat("Step 1: Extracting protein sequences...\n")
  protein_seqs <- extract_rmats_proteins(event_coords, translation_results, phase_results, event_type)
  
  # Ensure no factors in protein_seqs
  for (col in names(protein_seqs)) {
    if (is.factor(protein_seqs[[col]])) {
      protein_seqs[[col]] <- as.character(protein_seqs[[col]])
    }
  }
  
  # Step 2: Process ALL enzymes to match expected output structure
  cat("\nStep 2: Processing enzyme digestion...\n")
  for (enzyme in enzymes) {
    cat("Processing enzyme:", enzyme, "\n")
    cleaver_enzyme <- enzyme_map[[enzyme]]
    
    # Digest proteins with the enzyme
    protein_seqs[, paste0(enzyme, "Peps") := sapply(seq, function(x) {
      tryCatch({
        cleaver::cleave(
          x, enzym = cleaver_enzyme,
          missedCleavages = missedCleavages,
          custom = NULL, unique = TRUE
        )
      }, error = function(e) {
        warning("Error digesting protein with ", enzyme, ": ", e$message)
        return(character(0))
      })
    })]
    
    # Filter by length
    protein_seqs[, paste0(enzyme, "Peps") := sapply(get(paste0(enzyme, "Peps")), function(x) {
      if (length(x) == 0) return(character(0))
      aaNo <- sapply(x, nchar)
      x[aaNo >= minLength & aaNo <= maxLength]
    })]
  }
  
  # Step 3: Calculate peptide positions for ALL enzymes
  cat("\nStep 3: Calculating peptide positions...\n")
  
  # Function to get peptide positions (same as novel_peptide_generator.R)
  get_peptide_positions <- function(txID_i, peptides) {
    protein_seq <- protein_seqs[txID == txID_i, seq]
    
    if (length(protein_seq) == 0 || is.null(peptides) || length(peptides) == 0 || all(is.na(peptides))) {
      return(NA)
    }
    
    protein_seq <- as.character(protein_seq[[1]])
    peptide_positions_list <- lapply(peptides, function(pep) {
      match_positions <- tryCatch({
        str_locate_all(protein_seq, fixed(pep))[[1]]
      }, error = function(e) {
        warning("Error locating peptide positions: ", e$message)
        return(matrix(nrow=0, ncol=2, dimnames=list(NULL, c("start", "end"))))
      })
      
      if (nrow(match_positions) == 0) {
        return(NULL)
      }
      
      data.frame(
        peptide = pep,
        aa_start = match_positions[, "start"],
        aa_end = match_positions[, "end"],
        stringsAsFactors = FALSE
      )
    })
    
    peptide_positions_df <- do.call(rbind, peptide_positions_list)
    
    if (!is.null(peptide_positions_df) && nrow(peptide_positions_df) > 0) {
      return(peptide_positions_df)
    } else {
      return(NA)
    }
  }
  
  # Apply get_peptide_positions to ALL enzyme peptide types
  for (enzyme in enzymes) {
    enzyme_pep <- paste0(enzyme, "Peps")
    peptide_col <- paste0(enzyme_pep, "_positions")
    
    protein_seqs[, (peptide_col) := lapply(.I, function(i) {
      txID_i <- txID[i]
      peptides_i <- get(enzyme_pep)[[i]]
      
      if (is.null(peptides_i) || length(peptides_i) == 0 || all(is.na(peptides_i))) {
        return(NA)
      }
      
      positions_df <- tryCatch({
        result <- get_peptide_positions(txID_i, peptides_i)
        # Ensure no factors in the result
        if (is.data.frame(result)) {
          for (col in names(result)) {
            if (is.factor(result[[col]])) {
              result[[col]] <- as.character(result[[col]])
            }
          }
        }
        result
      }, error = function(e) {
        warning("Error calculating peptide positions for ", txID_i, ": ", e$message)
        return(NA)
      })
      
      return(positions_df)
    })]
  }
  
  # Step 4: Create GTF structures from phase results for genomic mapping
  cat("\nStep 4: Creating genomic coordinate mapping...\n")
  rmats_gtf <- create_rmats_gtf_from_phase(phase_results, event_type)
  gtf_gr <- convert_rmats_gtf_to_granges(rmats_gtf)
  
  # Extract CDS information for mapping
  cds_gr <- gtf_gr[gtf_gr$type == "CDS"]
  phases <- mcols(cds_gr)$phase
  cds_by_tx <- split(cds_gr, cds_gr$transcript_id)
  
  # Step 5: Functions for genomic mapping (same as novel_peptide_generator.R)
  
  # Validate strand consistency
  validate_strand_consistency <- function(cds_tx) {
    strands <- unique(as.character(strand(cds_tx)))
    if (length(strands) != 1) {
      stop("Inconsistent strand information across CDS entries for the transcript.")
    }
    return(strands[1])
  }
  
  # Retrieve and sort CDS information
  retrieve_sorted_cds <- function(txID, cds_by_tx) {
    cds_tx <- cds_by_tx[[txID]]
    if (is.null(cds_tx)) {
      warning(paste("Transcript ID", txID, "not found in CDS entries."))
      return(NULL)
    }
    
    # Validate strand consistency
    strand_tx <- validate_strand_consistency(cds_tx)
    
    # Sort CDS based on strand
    if (strand_tx == "+") {
      cds_tx <- cds_tx[order(start(cds_tx))]
    } else {
      cds_tx <- cds_tx[order(start(cds_tx), decreasing = TRUE)]
    }
    
    return(list(cds_tx = cds_tx, strand_tx = strand_tx))
  }
  
  # Calculate cumulative CDS positions
  calculate_cumulative_cds <- function(cds_tx, phases, strand_tx) {
    cumulative_cds <- data.frame(
      exon_idx = integer(),
      cds_start_nt = integer(),
      cds_end_nt = integer(),
      genomic_start = integer(),
      genomic_end = integer(),
      stringsAsFactors = FALSE
    )
    
    cumulative_length <- 0
    for (i in seq_along(cds_tx)) {
      exon <- cds_tx[i]
      exon_phase <- phases[i]
      exon_length <- width(exon)
      
      if (i == 1) {
        cds_start_nt <- 1 + exon_phase
      } else {
        cds_start_nt <- cumulative_length + 1 + exon_phase
      }
      
      cds_end_nt <- cds_start_nt + (exon_length - exon_phase) - 1
      
      if (strand_tx == "+") {
        genomic_start <- start(exon) + exon_phase
        genomic_end <- end(exon)
      } else {
        genomic_start <- end(exon) - exon_phase
        genomic_end <- start(exon)
      }
      
      cumulative_cds <- rbind(cumulative_cds, data.frame(
        exon_idx = i,
        cds_start_nt = cds_start_nt,
        cds_end_nt = cds_end_nt,
        genomic_start = genomic_start,
        genomic_end = genomic_end,
        stringsAsFactors = FALSE
      ))
      
      cumulative_length <- cds_end_nt
    }
    
    return(cumulative_cds)
  }
  
  # Map CDS nucleotide to genomic position
  map_cds_nt_to_genomic <- function(cds_cumulative, nt_pos, strand_tx) {
    exon_row <- cds_cumulative[cds_cumulative$cds_start_nt <= nt_pos & cds_cumulative$cds_end_nt >= nt_pos, ]
    
    if (nrow(exon_row) == 0) {
      return(NA)
    }
    
    offset <- nt_pos - exon_row$cds_start_nt
    
    if (strand_tx == "+") {
      genomic_pos <- exon_row$genomic_start + offset
    } else {
      genomic_pos <- exon_row$genomic_start - offset
    }
    
    return(genomic_pos)
  }
  
  # Map peptide to genomic coordinates
  map_peptide_to_genomic <- function(peptide, cds_cumulative, strand_tx, chrom) {
    peptide_name <- peptide$peptide
    aa_start <- peptide$aa_start
    aa_end <- peptide$aa_end
    
    cds_start_pos <- (aa_start - 1) * 3 + 1
    cds_end_pos <- aa_end * 3
    
    # Find exons spanned by peptide
    exons_spanned <- which(
      (cds_cumulative$cds_start_nt <= cds_end_pos) &
        (cds_cumulative$cds_end_nt >= cds_start_pos)
    )
    
    if (length(exons_spanned) == 0) {
      warning(paste("Peptide", peptide_name, "has nucleotide positions outside CDS range."))
      return(NULL)
    }
    
    peptide_gr_list <- list()
    
    for (exon_idx in exons_spanned) {
      exon_info <- cds_cumulative[exon_idx, ]
      
      overlap_start <- max(cds_start_pos, exon_info$cds_start_nt)
      overlap_end <- min(cds_end_pos, exon_info$cds_end_nt)
      
      genomic_start <- map_cds_nt_to_genomic(cds_cumulative, overlap_start, strand_tx)
      genomic_end <- map_cds_nt_to_genomic(cds_cumulative, overlap_end, strand_tx)
      
      if (is.na(genomic_start) || is.na(genomic_end)) {
        warning(paste("Peptide", peptide_name, "has mapping issues in exon", exon_idx))
        next
      }
      
      ir_start <- min(genomic_start, genomic_end)
      ir_end <- max(genomic_start, genomic_end)
      
      # Create GRanges object for this peptide segment
      partial_gr <- GRanges(
        seqnames = chrom,
        ranges = IRanges(start = ir_start, end = ir_end),
        strand = strand_tx,
        peptide = peptide_name
      )
      
      peptide_gr_list[[length(peptide_gr_list) + 1]] <- partial_gr
    }
    
    if (length(peptide_gr_list) == 0) {
      return(NULL)
    }
    
    peptide_gr <- do.call(c, peptide_gr_list)
    return(peptide_gr)
  }
  
  # Step 6: Initialize mapped_ranges columns for ALL enzymes
  for (enzyme in enzymes) {
    mapped_ranges_column_name <- paste0(enzyme, "Peps_mapped_ranges")
    protein_seqs[[mapped_ranges_column_name]] <- vector("list", nrow(protein_seqs))
  }
  
  # Step 7: Map peptides to genomic coordinates for ALL enzymes
  cat("\nStep 5: Mapping peptides to genomic coordinates...\n")
  for (entry_idx in 1:nrow(protein_seqs)) {
    entry <- protein_seqs[entry_idx, ]
    txID_i <- entry$txID
    
    # Retrieve sorted CDS and strand information
    cds_info <- tryCatch(
      retrieve_sorted_cds(txID_i, cds_by_tx),
      error = function(e) {
        warning(paste("Skipping transcript", txID_i, "due to error:", e$message))
        return(NULL)
      }
    )
    
    if (is.null(cds_info)) {
      next
    }
    
    cds_tx <- cds_info$cds_tx
    strand_tx <- cds_info$strand_tx
    chrom <- as.character(seqnames(cds_tx))[1]
    
    # Extract phase information for sorted CDS
    phases_tx <- mcols(cds_tx)$phase
    
    # Calculate cumulative CDS positions considering phase
    cds_cumulative <- calculate_cumulative_cds(cds_tx, phases_tx, strand_tx)
    
    # Loop over ALL enzymes
    for (enzyme in enzymes) {
      # Construct the column name for the peptide positions
      positions_column_name <- paste0(enzyme, "Peps_positions")
      
      # Check if the peptide positions column exists and is not empty
      if (!is.null(entry[[positions_column_name]]) && length(entry[[positions_column_name]]) > 0) {
        # Retrieve peptide positions data frame
        peptides_df <- as.data.frame(entry[[positions_column_name]][[1]])
        
        if (is.data.frame(peptides_df) && nrow(peptides_df) > 0) {
          # Initialize list to store GRanges for this enzyme
          peptide_gr_list <- list()
          
          # Iterate through each peptide
          for (peptide_idx in 1:nrow(peptides_df)) {
            peptide <- as.list(peptides_df[peptide_idx, ])
            
            # Map peptide to genomic ranges
            peptide_gr <- tryCatch({
              map_peptide_to_genomic(peptide, cds_cumulative, strand_tx, chrom)
            }, error = function(e) {
              warning("Error mapping peptide to genomic coordinates: ", e$message)
              return(NULL)
            })
            
            if (!is.null(peptide_gr)) {
              # Tag each peptide GRanges with transcript ID for correct faceting
              mcols(peptide_gr)$txID <- txID_i
              peptide_gr_list[[length(peptide_gr_list) + 1]] <- peptide_gr
            }
          }
          
          # Aggregate mapped peptides for this transcript and enzyme
          if (length(peptide_gr_list) > 0) {
            all_peptides_genomic_gr <- do.call(c, peptide_gr_list)
            # Store in the protein_seqs data frame
            mapped_ranges_column_name <- paste0(enzyme, "Peps_mapped_ranges")
            protein_seqs[[mapped_ranges_column_name]][[entry_idx]] <- all_peptides_genomic_gr
          }
        }
      }
    }
  }
  
  # Step 8: Ensure exact column order matching ENSG*.rds structure (24 columns)
  # CRITICAL: This must match the exact order from the working novel_peptide_generator.R
  expected_columns <- c(
    "proteinID", "txID", "geneID", "geneSymbol", "numAA", "seq",
    "trpPeps", "chymoPeps", "aspnPeps", "lyscPeps", "lysnPeps", "glucPeps",
    "trpPeps_positions", "chymoPeps_positions", "aspnPeps_positions", 
    "lyscPeps_positions", "lysnPeps_positions", "glucPeps_positions",
    "trpPeps_mapped_ranges", "aspnPeps_mapped_ranges", "chymoPeps_mapped_ranges",
    "glucPeps_mapped_ranges", "lysnPeps_mapped_ranges", "lyscPeps_mapped_ranges"
  )
  
  # Reorder columns to match expected structure
  protein_seqs <- protein_seqs[, ..expected_columns]
  
  # Convert to data.frame (not data.table) to match expected format
  return_data <- as.data.frame(protein_seqs)
  
  # Explicitly ensure all character columns stay as character
  for (col in names(return_data)) {
    if (is.character(return_data[[col]])) {
      return_data[[col]] <- as.character(return_data[[col]])
    }
  }
  
  cat("\n=== RMATS PEPTIDE DATA GENERATION COMPLETED ===\n")
  cat("Generated", nrow(return_data), "rows with", ncol(return_data), "columns\n")
  cat("Column names:", paste(names(return_data), collapse = ", "), "\n")
  
  # Event-type specific validation
  cat("\nEvent type specific validation for", event_type, ":\n")
  if (event_type == "RI") {
    cat("- RI events may have short exclusion proteins due to early stop codons\n")
  } else if (event_type == "SE") {
    cat("- SE events should show clear inclusion/exclusion differences\n") 
  } else if (event_type == "MXE") {
    cat("- MXE events have mutually exclusive exon differences\n")
  } else if (event_type %in% c("A3SS", "A5SS")) {
    cat("- Alternative splice site events may show subtle protein differences\n")
  }
  
  return(return_data)
}

#===============================================================================
# WRAPPER FUNCTION FOR APP INTEGRATION  
#===============================================================================

# Wrapper function that matches the signature expected by the app
create_rmats_rds_dataframe_new <- function(event_coords, translation_results, phase_results = NULL, 
                                          event_type = "SE", missedCleavages = 0, minLength = 7, maxLength = 60) {
  
  if (is.null(phase_results)) {
    stop("phase_results parameter is required for rMATS RDS generation")
  }
  
  cat("=== CREATING RMATS RDS DATAFRAME ===\n")
  cat("Using new self-contained generator\n")
  
  # Generate the peptide data using our new self-contained function
  peptide_data <- generate_rmats_peptide_data(
    event_coords = event_coords,
    translation_results = translation_results, 
    phase_results = phase_results,
    event_type = event_type,
    missedCleavages = missedCleavages,
    minLength = minLength,
    maxLength = maxLength
  )
  
  cat("Successfully created RDS dataframe with", nrow(peptide_data), "rows\n")
  return(peptide_data)
}

cat("=== RMATS NOVEL PEPTIDE GENERATOR LOADED ===\n")
cat("Available functions:\n")
cat("- generate_rmats_peptide_data(): Main generator function\n") 
cat("- create_rmats_rds_dataframe_new(): App integration wrapper\n")
cat("Ready to process all rMATS event types: SE, RI, MXE, A3SS, A5SS\n")