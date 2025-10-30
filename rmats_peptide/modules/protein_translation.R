# Protein Translation Module (Step 6)
# Comprehensive analysis with case handling for different biological scenarios

# Step 6: Comprehensive Protein Translation Analysis
analyze_protein_translation <- function(step5_results) {
  
  cat("=== STEP 6: COMPREHENSIVE PROTEIN TRANSLATION ANALYSIS ===\n")
  cat("Event Type:", step5_results$event_type, "\n")
  cat("Gene ID:", step5_results$gene_id, "\n")
  
  # Check BSgenome availability first
  if (!check_bsgenome_availability()) {
    return(create_error_result("BSgenome.Hsapiens.UCSC.hg38 package not available"))
  }
  
  # Determine case type based on Step 5 results
  case_type <- determine_case_type(step5_results)
  cat("Translation Case:", case_type, "\n")
  
  # Initialize result structure
  result <- list(
    event_type = step5_results$event_type,
    gene_id = step5_results$gene_id,
    case_type = case_type,
    inclusion_protein = NULL,
    exclusion_protein = NULL,
    inclusion_analysis = NULL,
    exclusion_analysis = NULL,
    functional_consequence = NULL,
    nmd_predictions = list(inclusion_nmd = FALSE, exclusion_nmd = FALSE),
    warnings = c(),
    summary = "",
    success = TRUE
  )
  
  # Process based on case type
  if (case_type == "both_translatable") {
    result <- process_both_translatable(step5_results, result)
  } else if (case_type == "inclusion_only") {
    result <- process_inclusion_only(step5_results, result)
  } else if (case_type == "exclusion_only") {
    result <- process_exclusion_only(step5_results, result)
  } else if (case_type == "neither_translatable") {
    result <- process_neither_translatable(step5_results, result)
  }
  
  # Generate comprehensive summary
  result$summary <- generate_biological_summary(result)
  
  cat("\n=== STEP 6 RESULTS ===\n")
  cat("Case Type:", result$case_type, "\n")
  cat("Functional Consequence:", result$functional_consequence, "\n")
  cat("Warnings:", paste(result$warnings, collapse = ", "), "\n")
  cat("Summary:", result$summary, "\n")
  
  return(result)
}

# Determine translation case type based on Step 5 results
determine_case_type <- function(step5_results) {
  
  inclusion_translatable <- step5_results$inclusion_translatable %||% FALSE
  exclusion_translatable <- step5_results$exclusion_translatable %||% FALSE
  
  if (inclusion_translatable && exclusion_translatable) {
    return("both_translatable")
  } else if (inclusion_translatable && !exclusion_translatable) {
    return("inclusion_only")
  } else if (!inclusion_translatable && exclusion_translatable) {
    return("exclusion_only") 
  } else {
    return("neither_translatable")
  }
}

# Check if BSgenome package is available
check_bsgenome_availability <- function() {
  tryCatch({
    # Check if BSgenome is installed
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      cat("BSgenome.Hsapiens.UCSC.hg38 not installed. Checking for alternative genome sources...\n")
      
      # Check if we have local genome file
      if (file.exists("rmats_peptide/GRCh38.fa") || file.exists("reference/GRCh38.fa")) {
        cat("Local genome file found - will use for sequence extraction\n")
        library(Biostrings, quietly = TRUE)
        library(rtracklayer, quietly = TRUE)
        return(TRUE)
      }
      
      cat("WARNING: BSgenome.Hsapiens.UCSC.hg38 not available and no local genome found\n")
      cat("To install BSgenome, run: BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')\n")
      return(FALSE)
    }
    
    library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)
    library(rtracklayer, quietly = TRUE)
    library(Biostrings, quietly = TRUE)
    return(TRUE)
    
  }, error = function(e) {
    cat("ERROR: Required packages not available:", e$message, "\n")
    return(FALSE)
  })
}

# Process Case 1: Both isoforms translatable
process_both_translatable <- function(step5_results, result) {
  cat("\n--- Processing Both Translatable Case ---\n")
  
  # Translate inclusion isoform
  inclusion_result <- translate_single_isoform(step5_results$inclusion_cds_gtf, "inclusion")
  result$inclusion_protein <- inclusion_result$protein
  result$inclusion_analysis <- inclusion_result$analysis
  
  # Translate exclusion isoform
  exclusion_result <- translate_single_isoform(step5_results$exclusion_cds_gtf, "exclusion")
  result$exclusion_protein <- exclusion_result$protein
  result$exclusion_analysis <- exclusion_result$analysis
  
  # Collect warnings
  result$warnings <- c(result$warnings, inclusion_result$warnings, exclusion_result$warnings)
  
  # NMD predictions
  result$nmd_predictions$inclusion_nmd <- inclusion_result$analysis$nmd_candidate
  result$nmd_predictions$exclusion_nmd <- exclusion_result$analysis$nmd_candidate
  
  # Compare proteins for functional consequence
  result$functional_consequence <- compare_proteins(inclusion_result, exclusion_result)
  
  return(result)
}

# Process Case 2: Only inclusion translatable
process_inclusion_only <- function(step5_results, result) {
  cat("\n--- Processing Inclusion Only Case ---\n")
  
  # Translate inclusion isoform
  inclusion_result <- translate_single_isoform(step5_results$inclusion_cds_gtf, "inclusion")
  result$inclusion_protein <- inclusion_result$protein
  result$inclusion_analysis <- inclusion_result$analysis
  result$warnings <- c(result$warnings, inclusion_result$warnings)
  result$nmd_predictions$inclusion_nmd <- inclusion_result$analysis$nmd_candidate
  
  # Exclusion isoform analysis
  result$exclusion_analysis <- list(
    translatable = FALSE,
    reason = step5_results$exclusion_reason %||% "First exon not in CDS",
    length = 0,
    nmd_candidate = FALSE
  )
  
  result$functional_consequence <- "inclusion_coding_exclusion_noncoding"
  
  return(result)
}

# Process Case 3: Only exclusion translatable
process_exclusion_only <- function(step5_results, result) {
  cat("\n--- Processing Exclusion Only Case ---\n")
  
  # Translate exclusion isoform
  exclusion_result <- translate_single_isoform(step5_results$exclusion_cds_gtf, "exclusion")
  result$exclusion_protein <- exclusion_result$protein
  result$exclusion_analysis <- exclusion_result$analysis
  result$warnings <- c(result$warnings, exclusion_result$warnings)
  result$nmd_predictions$exclusion_nmd <- exclusion_result$analysis$nmd_candidate
  
  # Inclusion isoform analysis
  result$inclusion_analysis <- list(
    translatable = FALSE,
    reason = step5_results$inclusion_reason %||% "First exon not in CDS",
    length = 0,
    nmd_candidate = FALSE
  )
  
  result$functional_consequence <- "exclusion_coding_inclusion_noncoding"
  
  return(result)
}

# Process Case 4: Neither translatable
process_neither_translatable <- function(step5_results, result) {
  cat("\n--- Processing Neither Translatable Case ---\n")
  
  result$inclusion_analysis <- list(
    translatable = FALSE,
    reason = step5_results$inclusion_reason %||% "First exon not in CDS",
    length = 0,
    nmd_candidate = FALSE
  )
  
  result$exclusion_analysis <- list(
    translatable = FALSE,
    reason = step5_results$exclusion_reason %||% "First exon not in CDS", 
    length = 0,
    nmd_candidate = FALSE
  )
  
  result$functional_consequence <- "noncoding_alternative_splicing"
  
  return(result)
}

# Translate a single isoform using BSgenome with PHASE-AWARE extraction
translate_single_isoform <- function(cds_gtf, isoform_type) {
  
  cat("  Translating", isoform_type, "isoform...\n")
  
  if (is.null(cds_gtf) || nrow(cds_gtf) == 0) {
    return(list(
      protein = NULL,
      analysis = list(translatable = FALSE, reason = "No CDS GTF", length = 0, nmd_candidate = FALSE),
      warnings = c()
    ))
  }
  
  tryCatch({
    # Extract phase information and coordinates
    strand_info <- cds_gtf$strand[1]
    phases <- cds_gtf$frame  # These are now ACTUAL CDS phases from Step 5
    
    cat("    Strand:", strand_info, "\n")
    cat("    CDS phases:", paste(phases, collapse = ", "), "\n")
    
    # Convert GTF to GRanges and sort by biological order (5' to 3')
    cds_ranges <- GRanges(
      seqnames = cds_gtf$seqname,
      ranges = IRanges(start = cds_gtf$start, end = cds_gtf$end),
      strand = cds_gtf$strand
    )
    
    # Sort by biological order (maintain the minus strand fix)
    if (strand_info == "-") {
      # For minus strand: use REVERSE genomic order to get correct 5' to 3' biological sequence
      order_indices <- order(start(cds_ranges), decreasing = TRUE)
      cds_ranges <- cds_ranges[order_indices]
      phases <- phases[order_indices]  # Reorder phases to match
      cat("    Minus strand: sorted", length(cds_ranges), "exons in reverse genomic order (5' to 3' biological)\n")
    } else {
      # For plus strand: use ascending genomic order (5' to 3' = left to right)
      order_indices <- order(start(cds_ranges))
      cds_ranges <- cds_ranges[order_indices]
      phases <- phases[order_indices]  # Reorder phases to match
      cat("    Plus strand: sorted", length(cds_ranges), "exons by genomic coordinates\n")
    }
    
    cat("    Extracting", length(cds_ranges), "CDS exons with correct phase interpretation\n")
    
    # CRITICAL: Use biologically correct CDS sequence building 
    cds_result <- build_correct_cds_sequence(cds_ranges, phases, strand_info)
    full_cds <- cds_result$sequence
    warnings <- cds_result$warnings
    
    cat("    Final CDS length:", length(full_cds), "bp\n")
    cat("    Phase applied:", cds_result$phase_applied, "\n")
    cat("    Ready for translation:", length(full_cds) > 0, "\n")
    
    # Translate to protein (accept biological reality - don't force clean results)
    if (length(full_cds) == 0) {
      warnings <- c(warnings, "empty_sequence_after_phase")
      cat("    WARNING: Empty sequence after phase application\n")
      protein <- ""
      protein_length <- 0
    } else {
      # Translate and handle stop codons properly - stop at first stop codon
      # Debug the CDS object before translation
      cat("    DEBUG: CDS class:", class(full_cds), "length:", length(full_cds), "\n")
      if (length(full_cds) < 1) {
        cat("    ERROR: Empty CDS sequence\n")
        full_protein <- ""
      } else {
        # Check if full_cds is a proper DNAString object
        if (!inherits(full_cds, "DNAString") && !inherits(full_cds, "DNAStringSet")) {
          cat("    ERROR: CDS is not a DNAString object, attempting conversion\n")
          full_cds <- DNAString(as.character(full_cds))
        }
        
        tryCatch({
          full_protein <- as.character(Biostrings::translate(full_cds))
        }, error = function(e) {
          cat("    ERROR in translation:", e$message, "\n")
          cat("    CDS content preview:", substr(as.character(full_cds), 1, 60), "...\n")
          warnings <<- c(warnings, "translation_error")
          full_protein <<- ""
        })
      }
      
      # Find first stop codon and truncate translation there (biologically accurate)
      first_stop <- regexpr("\\*", full_protein)
      if (first_stop > 0) {
        protein <- substr(full_protein, 1, first_stop - 1)  # Exclude the stop codon
        cat("    NOTE: Translation stopped at first stop codon (position", first_stop, ")\n")
      } else {
        protein <- full_protein  # No stop codons found
      }
      protein_length <- nchar(protein)
      
      if (length(full_cds) %% 3 != 0) {
        cat("    NOTE: Translating sequence not divisible by 3 - may have incomplete final codon\n")
      }
    }
    
    cat("    Protein length:", protein_length, "amino acids\n")
    
    # Analyze stop codons for NMD prediction
    nmd_analysis <- analyze_stop_codons(protein)
    
    if (nmd_analysis$premature_stops > 0) {
      warnings <- c(warnings, "premature_stop_codons")
      cat("    WARNING:", nmd_analysis$premature_stops, "premature stop codon(s) detected\n")
    }
    
    # Create analysis result
    analysis <- list(
      translatable = TRUE,
      reason = "Successfully translated",
      length = protein_length,
      cds_length = length(full_cds),
      nmd_candidate = nmd_analysis$nmd_candidate,
      first_stop_position = nmd_analysis$first_stop_position,
      total_stops = nmd_analysis$total_stops,
      premature_stops = nmd_analysis$premature_stops
    )
    
    return(list(
      protein = protein,
      analysis = analysis,
      warnings = warnings
    ))
    
  }, error = function(e) {
    cat("    ERROR in translation:", e$message, "\n")
    return(list(
      protein = NULL,
      analysis = list(translatable = FALSE, reason = paste("Translation error:", e$message), length = 0, nmd_candidate = FALSE),
      warnings = c("translation_error")
    ))
  })
}

# Analyze stop codons for NMD prediction
analyze_stop_codons <- function(protein) {
  
  # Find all stop codon positions
  stop_positions <- gregexpr("\\*", protein)[[1]]
  
  if (stop_positions[1] == -1) {
    # No stop codons found
    return(list(
      nmd_candidate = FALSE,
      first_stop_position = NA,
      total_stops = 0,
      premature_stops = 0
    ))
  }
  
  total_stops <- length(stop_positions)
  first_stop <- stop_positions[1]
  protein_length <- nchar(protein)
  
  # Consider premature if stop codon is not at the very end
  # More sophisticated NMD rules could be implemented here
  premature_stops <- sum(stop_positions < protein_length)
  
  nmd_candidate <- (first_stop < protein_length) || (premature_stops > 1)
  
  return(list(
    nmd_candidate = nmd_candidate,
    first_stop_position = first_stop,
    total_stops = total_stops,
    premature_stops = premature_stops
  ))
}

# Compare proteins to determine functional consequence
compare_proteins <- function(inclusion_result, exclusion_result) {
  
  inc_protein <- inclusion_result$protein
  exc_protein <- exclusion_result$protein
  
  # Both failed translation
  if (is.null(inc_protein) && is.null(exc_protein)) {
    return("both_translation_failed")
  }
  
  # One failed translation
  if (is.null(inc_protein) || is.null(exc_protein)) {
    return("one_translation_failed")
  }
  
  # Both successful - compare sequences
  if (identical(inc_protein, exc_protein)) {
    return("identical_proteins")
  }
  
  # Different proteins - analyze type of difference
  inc_length <- nchar(inc_protein)
  exc_length <- nchar(exc_protein)
  
  inc_nmd <- inclusion_result$analysis$nmd_candidate
  exc_nmd <- exclusion_result$analysis$nmd_candidate
  
  # NMD differences
  if (inc_nmd != exc_nmd) {
    return("nmd_difference")
  }
  
  # Length differences
  if (inc_length != exc_length) {
    return("length_difference")
  }
  
  # Same length but different sequence
  return("sequence_difference")
}

# Generate comprehensive biological summary
generate_biological_summary <- function(result) {
  
  case_type <- result$case_type
  consequence <- result$functional_consequence
  
  if (case_type == "both_translatable") {
    if (consequence == "identical_proteins") {
      return("Alternative splicing produces identical proteins - no functional impact.")
    } else if (consequence == "nmd_difference") {
      return("Alternative splicing affects NMD susceptibility - potential regulatory impact.")
    } else if (consequence == "length_difference") {
      inc_len <- result$inclusion_analysis$length
      exc_len <- result$exclusion_analysis$length
      diff <- abs(inc_len - exc_len)
      return(paste0("Proteins differ by ", diff, " amino acids (", inc_len, " vs ", exc_len, ") - functional domains may be affected."))
    } else if (consequence == "sequence_difference") {
      return("Alternative splicing produces different protein sequences - functional impact likely.")
    }
  } else if (case_type == "inclusion_only") {
    return("Inclusion isoform produces protein, exclusion isoform is non-coding - major functional impact.")
  } else if (case_type == "exclusion_only") {
    return("Exclusion isoform produces protein, inclusion isoform is non-coding - major functional impact.")
  } else if (case_type == "neither_translatable") {
    return("Alternative splicing occurs in non-coding region - likely regulatory or structural impact.")
  }
  
  return("Unable to determine functional consequence.")
}

# Create error result structure
create_error_result <- function(error_message) {
  return(list(
    case_type = "error",
    inclusion_protein = NULL,
    exclusion_protein = NULL,
    functional_consequence = "analysis_failed",
    nmd_predictions = list(inclusion_nmd = FALSE, exclusion_nmd = FALSE),
    warnings = c("analysis_error"),
    summary = paste0("Analysis failed: ", error_message),
    success = FALSE
  ))
}

# CORRECT: Properly handle GTF phase for multi-exon CDS sequences
build_correct_cds_sequence <- function(cds_ranges, phases, strand_info) {
  
  cat("    Building CDS sequence with proper phase handling...\n")
  warnings <- c()
  
  # Extract all sequences using BSgenome or local genome
  sequence_set <- tryCatch({
    # Try BSgenome first if available
    if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      getSeq(BSgenome.Hsapiens.UCSC.hg38, cds_ranges)
    } else {
      # Fallback to local genome file
      genome_file <- if (file.exists("rmats_peptide/GRCh38.fa")) {
        "rmats_peptide/GRCh38.fa"
      } else if (file.exists("reference/GRCh38.fa")) {
        "reference/GRCh38.fa"
      } else {
        stop("No genome source available for sequence extraction")
      }
      
      # Read local genome
      genome <- readDNAStringSet(genome_file)
      names(genome) <- gsub(" .*", "", names(genome))  # Clean chromosome names
      
      # Extract sequences from local genome
      getSeq(genome, cds_ranges)
    }
  }, error = function(e) {
    cat("    ERROR: Failed to extract sequences:", e$message, "\n")
    stop("Cannot proceed without genome sequence")
  })
  
  # CRITICAL FIX: getSeq returns DNAStringSet, extract individual DNAStrings
  sequences <- list()
  for (i in seq_along(sequence_set)) {
    sequences[[i]] <- sequence_set[[i]]
  }
  
  cat("    Extracted", length(sequences), "exons\n")
  for (i in seq_along(sequences)) {
    cat("      Exon", i, ": length =", length(sequences[[i]]), "bp, phase =", phases[i], "\n")
  }
  
  # GTF Phase specification:
  # - Phase indicates how many bases from the start of the exon need to be skipped to reach the next codon boundary
  # - Phase 0: Start of exon aligns with codon boundary (no skip)
  # - Phase 1: Skip 1 base from start to reach codon boundary  
  # - Phase 2: Skip 2 bases from start to reach codon boundary
  
  # For the FIRST exon: phase tells us where the CDS starts within the exon
  first_exon <- sequences[[1]]
  first_phase <- phases[1]
  
  # Apply phase correction to first exon
  if (first_phase > 0 && length(first_exon) > first_phase) {
    first_exon <- subseq(first_exon, start = first_phase + 1)
    cat("    First exon: applied phase", first_phase, "correction (removed", first_phase, "bases from start)\n")
  } else {
    cat("    First exon: no phase correction needed (phase", first_phase, ")\n")
  }
  
  # For SUBSEQUENT exons: handle phase transitions correctly
  processed_sequences <- list(first_exon)
  
  if (length(sequences) > 1) {
    for (i in 2:length(sequences)) {
      current_exon <- sequences[[i]]
      current_phase <- phases[i]
      
      # Calculate cumulative length so far to check frame alignment
      cumulative_length <- sum(sapply(processed_sequences, length))
      expected_phase <- cumulative_length %% 3
      
      cat("    Exon", i, ": phase =", current_phase, ", expected based on cumulative =", expected_phase, "\n")
      
      # For alternative splicing: document phase mismatches but preserve original coordinates
      # The GTF coordinates should be trusted for alternative splicing analysis
      if (current_phase != expected_phase) {
        cat("    NOTE: Exon", i, "phase mismatch (got", current_phase, ", expected", expected_phase, ")\n")
        cat("    This is expected in alternative splicing when exons are skipped\n")
        warnings <- c(warnings, paste0("phase_mismatch_exon_", i))
      } else {
        cat("    Exon", i, ": phase matches expected\n")
      }
      
      # Always use the original exon sequence - trust the CDS coordinates
      processed_sequences[[i]] <- current_exon
    }
  }
  
  # Concatenate all processed sequences
  complete_cds <- do.call(c, processed_sequences)
  final_length <- length(complete_cds)
  
  cat("    Complete CDS sequence length:", final_length, "bp\n")
  cat("    Applied phase correction:", first_phase, "to first exon\n")
  
  # Check if divisible by 3 (proper CDS should be)
  if (final_length %% 3 != 0) {
    remainder <- final_length %% 3
    cat("    NOTE: Sequence length not divisible by 3 (remainder:", remainder, ")\n")
    cat("    This may indicate incomplete CDS boundaries in alternative splicing\n")
    warnings <- c(warnings, paste0("not_divisible_by_3_remainder_", remainder))
  } else {
    cat("    Sequence length is divisible by 3 - perfect for translation\n")
  }
  
  return(list(
    sequence = complete_cds,
    warnings = warnings,
    original_length = sum(sapply(sequences, length)),
    final_length = final_length,
    phase_applied = first_phase
  ))
}

# Helper function for null coalescing
`%||%` <- function(lhs, rhs) {
  if (is.null(lhs) || (is.logical(lhs) && !lhs)) rhs else lhs
}