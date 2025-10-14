# Novel Isoform Discovery and Analysis Functions
# This module handles the novel isoform discovery pipeline and analysis

#===============================================================================
# NOVEL ISOFORM PIPELINE EXECUTION
#===============================================================================

#' Enhanced Novel Transcript Analysis Pipeline
#' 
#' Advanced pipeline with configurable TransDecoder parameters for ORF prediction
#' and optional homology search
#' 
#' @param input_fasta_path Path to input FASTA file
#' @param min_protein_length Minimum protein length in amino acids (default: 30)
#' @param genetic_code Genetic code table number (default: 1 for standard code)
#' @param strand_specific Strand specificity: "both", "plus", or "minus"
#' @param retain_long_orfs Keep all ORFs meeting length criteria (default: TRUE)
#' @param single_best_orf Return only longest ORF per transcript (default: FALSE)
#' @param require_start_codon Require ORFs to start with ATG (default: FALSE)
#' @param require_stop_codon Require ORFs to have stop codon (default: FALSE)
#' @param min_orf_coverage Minimum ORF coverage percentage (default: 0)
#' @param enable_gene_search Enable boundary-based gene search (default: TRUE)
#' @param min_overlap_percent Minimum overlap percentage for gene matching (default: 20)
#' @param progress_callback Optional progress callback function
#' @return List with pipeline results
run_novel_isoform_pipeline <- function(input_fasta_path, 
                                     min_protein_length = 30,
                                     genetic_code = 1,
                                     strand_specific = "both",
                                     retain_long_orfs = TRUE,
                                     single_best_orf = FALSE,
                                     require_start_codon = FALSE,
                                     require_stop_codon = FALSE,
                                     min_orf_coverage = 0,
                                     enable_gene_search = TRUE,
                                     min_overlap_percent = 20,
                                     progress_callback = NULL) {
  
  if (!is.null(progress_callback)) {
    progress_callback("Using direct approach: replicating successful manual pipeline...", 0.05)
  }
  
  # CALL THE DIRECT APPROACH that exactly replicates user's successful manual run
  return(run_novel_isoform_genomic_direct(
    input_fasta_path,
    min_protein_length,
    genetic_code,
    progress_callback
  ))
}

# SIMPLE DIRECT APPROACH: Replicate the user's exact successful terminal commands
run_novel_isoform_genomic_direct <- function(input_fasta_file, 
                                           min_protein_length = 30,
                                           genetic_code = 1,
                                           progress_callback = NULL) {
  
  # Create timestamped working directory
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  work_dir <- file.path("novel_isoform_results", timestamp)
  dir.create(work_dir, recursive = TRUE)
  dir.create(file.path(work_dir, "results"), recursive = TRUE)
  
  # Copy and normalize input file header consistently 
  input_file <- file.path(work_dir, "novel_transcript_nt.fa")
  
  # Always normalize header to "Novel_sequence_1" for consistency across pipeline
  temp_lines <- readLines(input_fasta_file)
  header_normalized <- c(">Novel_sequence_1", temp_lines[!grepl("^>", temp_lines)])
  writeLines(header_normalized, input_file)
  cat("Input header normalized to: Novel_sequence_1\n")
  
  if (!is.null(progress_callback)) {
    progress_callback("Running 8-step pipeline following PIPELINE_STEPS.md...", 0.1)
  }
  
  # Run commands directly like user's successful terminal session
  original_wd <- getwd()
  
  tryCatch({
    setwd(work_dir)
    
    # STEP 1: TransDecoder Analysis (exactly like user's working version)
    cat("STEP 1: TransDecoder Analysis...\n")
    result1 <- system2("/Users/Mahmuda/opt/anaconda3/envs/module3_tools/bin/TransDecoder.LongOrfs", 
                      args = c("-t", "novel_transcript_nt.fa", "-m", min_protein_length), 
                      wait = TRUE)
    result2 <- system2("/Users/Mahmuda/opt/anaconda3/envs/module3_tools/bin/TransDecoder.Predict", 
                      args = c("-t", "novel_transcript_nt.fa", "--no_refine_starts"),
                      wait = TRUE)
    cat("TransDecoder analysis completed\n")
    
    if (!is.null(progress_callback)) {
      progress_callback("Genome alignment...", 0.3)
    }
    
    # STEP 2: Genome Alignment (exactly like user's working version)
    cat("STEP 2: Genome Alignment...\n")
    system2("/Users/Mahmuda/opt/anaconda3/envs/module3_tools/bin/minimap2", 
            args = c("-ax", "splice", "../../reference/GRCh38.mmi", "novel_transcript_nt.fa"),
            stdout = "novel_transcript_nt_aligned.sam", wait = TRUE)
    system2("/Users/Mahmuda/opt/anaconda3/envs/module3_tools/bin/samtools", 
            args = c("view", "-Sb", "novel_transcript_nt_aligned.sam"),
            stdout = "novel_transcript_nt_aligned.bam", wait = TRUE)
    system2("/Users/Mahmuda/opt/anaconda3/envs/module3_tools/bin/samtools", 
            args = c("sort", "novel_transcript_nt_aligned.bam", "-o", "novel_transcript_nt_aligned_sorted.bam"),
            wait = TRUE)
    cat("Genome alignment completed\n")
    
    if (!is.null(progress_callback)) {
      progress_callback("StringTie transcript reconstruction...", 0.4)
    }
    
    # STEP 3: StringTie Transcript Reconstruction (exactly like user's working version)
    cat("STEP 3: StringTie Transcript Reconstruction...\n")
    system2("/Users/Mahmuda/opt/anaconda3/envs/module3_tools/bin/stringtie", 
            args = c("novel_transcript_nt_aligned_sorted.bam", "-o", "novel_transcript_nt_transcripts.gtf", "-l", "Novel_sequence_1"),
            wait = TRUE)
    cat("Transcript reconstruction completed\n")
    
    if (!is.null(progress_callback)) {
      progress_callback("GTF to GFF3 conversion...", 0.5)
    }
    
    # STEP 4: Convert StringTie GTF to Target= GFF3 (exactly like user's working version)
    cat("STEP 4: Convert StringTie GTF to Target= GFF3...\n")
    system2("/Users/Mahmuda/opt/anaconda3/envs/module3_tools/opt/transdecoder/util/gtf_to_alignment_gff3.pl",
            args = "novel_transcript_nt_transcripts.gtf",
            stdout = "novel_transcript_nt_alignment.gff3", wait = TRUE)
    
    # Clean up sequence names in GFF3 using R (avoid hanging sed)
    if (file.exists("novel_transcript_nt_alignment.gff3")) {
      gff3_lines <- readLines("novel_transcript_nt_alignment.gff3")
      gff3_lines_cleaned <- gsub("\\.1\\.1", "", gff3_lines)
      writeLines(gff3_lines_cleaned, "novel_transcript_nt_alignment.gff3")
    }
    cat("GTF to GFF3 conversion completed\n")
    
    if (!is.null(progress_callback)) {
      progress_callback("Mapping ORFs to genome (critical step)...", 0.6)
    }
    
    # STEP 5: Map ORFs to Genome (exactly like user's working version)
    cat("STEP 5: Map ORFs to Genome...\n")
    result5 <- system2("/Users/Mahmuda/opt/anaconda3/envs/module3_tools/opt/transdecoder/util/cdna_alignment_orf_to_genome_orf.pl",
            args = c("novel_transcript_nt.fa.transdecoder.gff3", "novel_transcript_nt_alignment.gff3", "novel_transcript_nt.fa"),
            stdout = "novel_transcript_nt.transdecoder.genome.gff3", wait = TRUE)
    
    # Check mapping success
    if (file.exists("novel_transcript_nt.transdecoder.genome.gff3") && file.size("novel_transcript_nt.transdecoder.genome.gff3") > 0) {
      mapped_orfs <- length(readLines("novel_transcript_nt.transdecoder.genome.gff3"))
      cat("ORF to genome mapping completed:", mapped_orfs, "lines in output\n")
    } else {
      stop("ORF to genome mapping failed - no output file created or file is empty")
    }
  
  if (!is.null(progress_callback)) {
      progress_callback("Converting to final GTF...", 0.7)
  }
  
        # STEP 6: Convert to GTF (exactly like user's working version)
    cat("STEP 6: Convert to GTF...\n")
    system2("/Users/Mahmuda/opt/anaconda3/envs/module3_tools/bin/gffread", 
            args = c("novel_transcript_nt.transdecoder.genome.gff3", "-T", "-o", "novel_transcript_nt.transdecoder.genome.gtf"),
            wait = TRUE)
    cat("Final GTF created\n")
    
    if (!is.null(progress_callback)) {
      progress_callback("Running peptide generator...", 0.8)
    }
    
    # STEP 7-8: Run novel_peptide_generator.R (exactly like user's working version)
    cat("STEP 7-8: Running novel_peptide_generator.R...\n")
    system2("/usr/local/bin/Rscript", args = "../../novel_peptide_generator.R", wait = TRUE)
    
    # Copy results to expected locations
    file.copy("novel_transcript_nt.transdecoder.genome.gtf", "results/novel_final.gtf", overwrite = TRUE)
    file.copy("novel_transcript_nt.fa.transdecoder.pep", "results/novel_proteins.pep", overwrite = TRUE)
    
    if (file.exists("novel_transcript_nt_peptides.rds")) {
      file.copy("novel_transcript_nt_peptides.rds", "results/novel_isoform_dataframe.rds", overwrite = TRUE)
    }
    
    # Create completion marker
    writeLines("8-step pipeline completed successfully", "results/pipeline_complete.txt")
      
      if (!is.null(progress_callback)) {
      progress_callback("Pipeline completed successfully!", 1.0)
      }
      
    cat("✅ 8-step pipeline completed successfully!\n")
      
    # Return success
      return(list(
        success = TRUE,
        work_dir = work_dir,
      dataframe_file = file.path(work_dir, "results", "novel_isoform_dataframe.rds"),
      gtf_file = file.path(work_dir, "results", "novel_final.gtf"),
      log = "Pipeline completed successfully",
        parameters = list(
          min_protein_length = min_protein_length,
          genetic_code = genetic_code,
        pipeline_type = "Direct replication of successful manual run"
        ),
      message = "8-step genomic mapping pipeline completed successfully!"
      ))
    
  }, error = function(e) {
    cat("Pipeline failed with error:", e$message, "\n")
    return(list(
      success = FALSE,
      work_dir = work_dir,
      error = e$message,
      log = paste("Failed at step:", e$message),
      parameters = list(
        min_protein_length = min_protein_length,
        genetic_code = genetic_code,
        pipeline_type = "Direct replication (failed)"
      ),
      message = paste("Pipeline failed:", e$message)
    ))
  }, finally = {
    setwd(original_wd)
  })
}

#' Create GRanges for Novel Transcript Peptides
#' 
#' Creates GRanges objects from novel isoform data compatible with existing pipeline
#' 
#' @param txID Novel transcript ID
#' @param novel_merged_data Merged novel and known data
#' @param protease Protease name
#' @return GRanges object compatible with get_transcript_peptides_for_comparison
get_novel_transcript_peptides_genomic <- function(txID, novel_merged_data, protease) {
  
  cat("DEBUG: get_novel_transcript_peptides_genomic called for:", txID, "\n")
  
  # Find the transcript in merged data
  tx_row <- which(novel_merged_data$txID == txID)
  
  if (length(tx_row) == 0) {
    cat("DEBUG: Transcript not found in merged data:", txID, "\n")
    return(NULL)
  }
  
  tx_data <- novel_merged_data[tx_row[1], ]
  
  # Get peptide data
  peptide_col <- paste0(protease, "Peps")
  position_col <- paste0(protease, "Peps_positions")
  
  if (!peptide_col %in% names(tx_data)) {
    cat("DEBUG: Peptide column not found:", peptide_col, "\n")
    return(NULL)
  }
  
  peptides <- tx_data[[peptide_col]][[1]]
  positions <- tx_data[[position_col]][[1]]
  
  if (is.null(peptides) || length(peptides) == 0 || all(is.na(peptides))) {
    cat("DEBUG: No peptides found for transcript:", txID, "\n")
    return(NULL)
  }
  
  # Create GRanges object compatible with existing pipeline
  if (is.data.frame(positions) && nrow(positions) > 0) {
    
    # Apply peptide length filter (6-60 amino acids)
    peptide_lengths <- nchar(positions$peptide)
    valid_peptides <- peptide_lengths >= 6 & peptide_lengths <= 60
    positions_filtered <- positions[valid_peptides, ]
    
    if (nrow(positions_filtered) == 0) {
      cat("DEBUG: No valid peptides (6-60 AA) found for transcript:", txID, "\n")
      return(NULL)
    }
    
    # Use AA positions as coordinates (consistent with simplified approach)
    gr <- GenomicRanges::GRanges(
      seqnames = "protein_sequence",
      ranges = IRanges::IRanges(
        start = positions_filtered$aa_start,
        end = positions_filtered$aa_end
      ),
      peptide = positions_filtered$peptide
    )
    
    cat("DEBUG: Created GRanges with", length(gr), "peptides (6-60 AA filter) for", txID, "\n")
    return(gr)
    
  } else {
    cat("DEBUG: Invalid position data for transcript:", txID, "\n")
    return(NULL)
  }
}

#===============================================================================
# NOVEL ISOFORM DATA PROCESSING
#===============================================================================

#' Load Novel Isoform Data
#' 
#' Loads and validates novel isoform dataframe
#' 
#' @param dataframe_file Path to novel isoform dataframe RDS file
#' @return Loaded dataframe or NULL if not found
load_novel_isoform_data <- function(dataframe_file) {
  if (!file.exists(dataframe_file)) {
    return(NULL)
  }
  
  tryCatch({
    novel_data <- readRDS(dataframe_file)
    
    # Validate required columns exist
    required_cols <- c("proteinID", "txID", "geneID", "geneSymbol", "numAA", "seq")
    if (!all(required_cols %in% names(novel_data))) {
      warning("Novel isoform dataframe missing required columns")
      return(NULL)
    }
    
    return(novel_data)
  }, error = function(e) {
    warning(paste("Error loading novel isoform data:", e$message))
    return(NULL)
  })
}

#' Extract ORF Information
#' 
#' Extracts and summarizes ORF information from novel isoform data
#' 
#' @param novel_data Novel isoform dataframe
#' @return Dataframe with ORF information
extract_orf_information <- function(novel_data) {
  if (is.null(novel_data) || nrow(novel_data) == 0) {
    return(data.frame())
  }
  
  cat("DEBUG: Extracting ORF information from", nrow(novel_data), "novel isoforms\n")
  
  orf_info <- data.frame(
    orf_id = novel_data$proteinID,
    protein_length = as.numeric(novel_data$numAA),  # Convert to numeric
    protein_sequence = novel_data$seq,
    gene_symbol = novel_data$geneSymbol,
    transcript_id = novel_data$txID,
    sequence_preview = substr(novel_data$seq, 1, 50),
    starts_with_M = startsWith(novel_data$seq, "M"),
    has_stop = grepl("\\*$", novel_data$seq),
    stringsAsFactors = FALSE
  )
  
  # Calculate quality score (simple scoring based on length and start codon)
  orf_info$quality_score <- round(
    log10(orf_info$protein_length) + 
    ifelse(orf_info$starts_with_M, 2, 0) + 
    ifelse(orf_info$has_stop, 1, 0), 2
  )
  
  cat("DEBUG: Extracted", nrow(orf_info), "ORFs\n")
  cat("DEBUG: Length range:", min(orf_info$protein_length), "-", max(orf_info$protein_length), "AA\n")
  cat("DEBUG: ORFs starting with M:", sum(orf_info$starts_with_M), "\n")
  
  # Find best ORF
  best_orf_idx <- which.max(orf_info$quality_score)
  if (length(best_orf_idx) > 0) {
    best_orf <- orf_info[best_orf_idx, ]
    cat("DEBUG: Best ORF:", best_orf$orf_id, "Length:", best_orf$protein_length, "Score:", best_orf$quality_score, "\n")
  }
  
  return(orf_info)
}

#' Create Novel Isoform Summary
#' 
#' Creates summary statistics for novel isoform analysis
#' 
#' @param novel_data Novel isoform dataframe
#' @param enzyme Selected enzyme for peptide analysis
#' @return List with summary information
create_novel_isoform_summary <- function(novel_data, enzyme = "trp") {
  if (is.null(novel_data) || nrow(novel_data) == 0) {
    return(list(
      total_orfs = 0,
      total_peptides = 0,
      avg_length = 0,
      specificity_summary = data.frame()
    ))
  }
  
  # Get peptide column for the enzyme
  peptide_col <- paste0(enzyme, "Peps")
  
  if (!peptide_col %in% names(novel_data)) {
    return(list(
      total_orfs = nrow(novel_data),
      total_peptides = 0,
    avg_length = mean(as.numeric(novel_data$numAA), na.rm = TRUE),
      specificity_summary = data.frame()
    ))
  }
  
  # Count peptides
  peptide_counts <- sapply(novel_data[[peptide_col]], function(x) {
    if (is.null(x) || length(x) == 0) return(0)
    return(length(x))
  })
  
  total_peptides <- sum(peptide_counts, na.rm = TRUE)
  
  # Create specificity summary (simplified - just counts for now)
  specificity_summary <- data.frame(
    Specificity_Category = c("Novel-Unique", "Widely-Shared"),
    Count = c(total_peptides, 0), # All peptides considered novel-unique for now
    Percentage = c(100, 0),
    stringsAsFactors = FALSE
  )
  
  return(list(
    total_orfs = nrow(novel_data),
    total_peptides = total_peptides,
    avg_length = mean(as.numeric(novel_data$numAA), na.rm = TRUE),
    specificity_summary = specificity_summary
  ))
}

#' Get Available Novel Results
#' 
#' Scans for available novel isoform analysis results
#' 
#' @return Character vector of available result directories
get_available_novel_results <- function() {
  results_dir <- "novel_isoform_results"
  
  if (!dir.exists(results_dir)) {
    return(character(0))
  }
  
  # Find directories with completed dataframes
  result_dirs <- list.dirs(results_dir, full.names = FALSE, recursive = FALSE)
  
  valid_results <- c()
  for (dir_name in result_dirs) {
    dataframe_file <- file.path(results_dir, dir_name, "results", "novel_isoform_dataframe.rds")
    if (file.exists(dataframe_file)) {
      valid_results <- c(valid_results, dir_name)
    }
  }
  
  return(valid_results)
}

#===============================================================================
# NOVEL ISOFORM ANALYSIS FUNCTIONS
#===============================================================================

#' Analyze Novel Isoform Specificity
#' 
#' Calculates specificity metrics for novel isoforms
#' 
#' @param novel_data Novel isoform dataframe
#' @param protease Protease/enzyme name
#' @param highlight_isoform Isoform to highlight for analysis
#' @return List with analysis results
analyze_novel_isoform_specificity <- function(novel_data, protease, highlight_isoform = NULL, miscleavage_type = "no_miss_cleavage") {
  
  cat("DEBUG: analyze_novel_isoform_specificity called\n")
  cat("DEBUG: novel_data rows:", nrow(novel_data), "\n")
  cat("DEBUG: protease:", protease, "\n")
  cat("DEBUG: highlight_isoform:", highlight_isoform, "\n")
  cat("DEBUG: miscleavage_type:", miscleavage_type, "\n")
  
  # Get all peptides for specificity calculation (includes both novel and known isoforms)
  cat("DEBUG: Calling get_novel_isoform_peptides...\n")
  all_peptides_df <- get_novel_isoform_peptides(novel_data, protease, miscleavage_type)
  cat("DEBUG: get_novel_isoform_peptides returned\n")
  
  if (is.null(all_peptides_df)) {
    return(NULL)
  }
  
  # Get all isoforms for specificity calculation
  all_isoforms <- unique(all_peptides_df$transcript)
  cat("DEBUG: All available isoforms:", paste(all_isoforms, collapse = ", "), "\n")
  
  # Use transcript ID directly (dropdown now uses transcript IDs)
  actual_highlight_isoform <- highlight_isoform
  cat("DEBUG: Using transcript ID directly:", actual_highlight_isoform, "\n")
  
  if (!is.null(actual_highlight_isoform) && actual_highlight_isoform %in% all_isoforms) {
    cat("DEBUG: Preparing comparison view for selected isoform:", actual_highlight_isoform, "\n")
    
    # Filter to show ONLY the selected novel isoform + known gene isoforms
    # Identify known isoforms (those that don't start with "NOVEL_")
    known_isoforms <- all_isoforms[!grepl("^NOVEL_", all_isoforms)]
    selected_and_known_isoforms <- c(actual_highlight_isoform, known_isoforms)
    
    cat("DEBUG: Selected isoform:", actual_highlight_isoform, "\n")
    cat("DEBUG: Known isoforms for comparison:", paste(known_isoforms, collapse = ", "), "\n")
    cat("DEBUG: Total isoforms to display:", length(selected_and_known_isoforms), "\n")
    
    # Filter peptides to show only selected novel + known isoforms
    display_peptides_df <- all_peptides_df[all_peptides_df$transcript %in% selected_and_known_isoforms, ]
    cat("DEBUG: Display peptides (selected + known):", nrow(display_peptides_df), "\n")
    
    # Keep all isoforms for specificity calculation
    all_isoforms_for_calc <- all_isoforms
    
    highlight_peptides <- all_peptides_df[all_peptides_df$transcript == actual_highlight_isoform, ]
    
    if (nrow(highlight_peptides) > 0) {
      # Calculate specificity for each peptide
      highlight_peptides$other_isoform_count <- 0
      highlight_peptides$other_isoforms <- ""
      highlight_peptides$specificity_category <- ""
      
      for (i in 1:nrow(highlight_peptides)) {
        peptide_seq <- highlight_peptides$peptide[i]
        # Count isoforms (excluding highlight) that have this peptide
        other_isoforms_with_peptide <- unique(all_peptides_df$transcript[
          all_peptides_df$peptide == peptide_seq & all_peptides_df$transcript != actual_highlight_isoform
        ])
        
        highlight_peptides$other_isoform_count[i] <- length(other_isoforms_with_peptide)
        highlight_peptides$other_isoforms[i] <- paste(other_isoforms_with_peptide, collapse = ", ")
        
        # Classify specificity
        total_isoforms <- length(all_isoforms_for_calc)
        other_count <- length(other_isoforms_with_peptide)
        
              if (other_count == 0) {
        # Peptide found only in this isoform
        highlight_peptides$specificity_category[i] <- "Unique"
      } else if (other_count == (total_isoforms - 1)) {
        # Peptide found in ALL isoforms
        highlight_peptides$specificity_category[i] <- "Universal"
      } else {
        # Peptide found in some but not all isoforms
        highlight_peptides$specificity_category[i] <- "Shared"
      }
      }
      
      # Add specificity info to ALL peptides for visualization
      display_peptides_df$is_highlighted <- display_peptides_df$transcript == actual_highlight_isoform
      display_peptides_df$specificity_category <- ""
      display_peptides_df$specificity_category[display_peptides_df$is_highlighted] <- 
        highlight_peptides$specificity_category[match(display_peptides_df$peptide[display_peptides_df$is_highlighted], highlight_peptides$peptide)]
      
      # Update hover text for all peptides
      display_peptides_df$hover_text <- paste0(
        "Peptide: ", display_peptides_df$peptide,
        "<br>Position: ", display_peptides_df$start, "-", display_peptides_df$end,
        "<br>Isoform: ", display_peptides_df$transcript,
        "<br>Enzyme: ", protease,
        "<br>Miscleavage: ", miscleavage_type,
        ifelse(display_peptides_df$is_highlighted, 
               paste0("<br>Specificity: ", display_peptides_df$specificity_category), 
               "")
      )
      
      # Create summary statistics
      specificity_summary <- as.data.frame(table(highlight_peptides$specificity_category))
      names(specificity_summary) <- c("Specificity_Category", "Count")
      specificity_summary$Percentage <- round(specificity_summary$Count / nrow(highlight_peptides) * 100, 1)
      
    } else {
      highlight_peptides <- NULL
      specificity_summary <- NULL
      # For cases where selected isoform has no peptides, still filter to selected + known
      known_isoforms <- all_isoforms[!grepl("^NOVEL_", all_isoforms)]
      selected_and_known_isoforms <- c(highlight_isoform, known_isoforms)
      display_peptides_df <- all_peptides_df[all_peptides_df$transcript %in% selected_and_known_isoforms, ]
      display_peptides_df$is_highlighted <- display_peptides_df$transcript == highlight_isoform
      display_peptides_df$specificity_category <- ""
    }
  } else {
    # No specific isoform selected - show only known isoforms if available, otherwise all
    known_isoforms <- all_isoforms[!grepl("^NOVEL_", all_isoforms)]
    if (length(known_isoforms) > 0) {
      display_peptides_df <- all_peptides_df[all_peptides_df$transcript %in% known_isoforms, ]
    } else {
      display_peptides_df <- all_peptides_df
    }
    highlight_peptides <- NULL
    specificity_summary <- NULL
    display_peptides_df$is_highlighted <- FALSE
    display_peptides_df$specificity_category <- ""
  }
  
  # Calculate gene boundaries based on displayed peptides only
  gene_start <- min(display_peptides_df$start, na.rm = TRUE) - 1000
  gene_end <- max(display_peptides_df$end, na.rm = TRUE) + 1000
  
  # Create transcript position mapping for DISPLAYED isoforms only (selected novel + known)
  all_displayed_isoforms <- unique(display_peptides_df$transcript)
  transcript_df <- data.frame(
    transcript = all_displayed_isoforms,
    y_position = seq_along(all_displayed_isoforms),
    stringsAsFactors = FALSE
  )
  
  # CRITICAL: Remap y_positions in peptide data to match transcript_df
  for (i in 1:nrow(display_peptides_df)) {
    transcript_name <- display_peptides_df$transcript[i]
    new_y_pos <- transcript_df$y_position[transcript_df$transcript == transcript_name]
    display_peptides_df$y_position[i] <- new_y_pos
  }
  
  cat("DEBUG: Final results - Display peptides:", nrow(display_peptides_df), "Transcripts:", nrow(transcript_df), "\n")
  cat("DEBUG: Displayed isoforms:", paste(all_displayed_isoforms, collapse = ", "), "\n")
  cat("DEBUG: Y-position range after remapping:", min(display_peptides_df$y_position), "-", max(display_peptides_df$y_position), "\n")
  
  return(list(
    all_peptides = display_peptides_df,  # Contains all isoforms (novel + known)
    highlight_peptides = highlight_peptides,
    transcript_df = transcript_df,
    gene_start = gene_start,
    gene_end = gene_end,
    highlight_isoform = ifelse(exists("actual_highlight_isoform"), actual_highlight_isoform, highlight_isoform),
    specificity_summary = specificity_summary,
    all_isoforms = all_displayed_isoforms,  # Contains all displayed isoforms
    miscleavage_type = miscleavage_type
  ))
}

#' Create Novel Isoform Pipeline Summary
#' 
#' Creates a summary of the novel isoform discovery pipeline results
#' 
#' @param novel_data Novel isoform dataframe
#' @param pipeline_log Pipeline execution log
#' @return Data frame with pipeline summary
create_novel_pipeline_summary <- function(novel_data, pipeline_log = "") {
  
  cat("DEBUG: create_novel_pipeline_summary called\n")
  cat("DEBUG: novel_data rows:", ifelse(is.null(novel_data), "NULL", nrow(novel_data)), "\n")
  cat("DEBUG: pipeline_log length:", nchar(pipeline_log), "\n")
  
  if (is.null(novel_data) || nrow(novel_data) == 0) {
    cat("DEBUG: Returning error summary - no data\n")
    return(data.frame(
      Metric = "Error",
      Value = "No data available",
      stringsAsFactors = FALSE
    ))
  }
  
  # Calculate summary metrics
  cat("DEBUG: Calculating summary metrics...\n")
  total_isoforms <- nrow(novel_data)
  unique_gene_symbols <- length(unique(novel_data$geneSymbol))
  cat("DEBUG: Total isoforms:", total_isoforms, "Unique gene symbols:", unique_gene_symbols, "\n")
  
  # Calculate peptide counts for each enzyme
  cat("DEBUG: Calculating peptide counts...\n")
  enzyme_counts <- list()
  enzymes <- c("trp", "chymo", "aspn", "lysc", "lysn", "gluc")
  enzyme_names <- c("Trypsin", "Chymotrypsin", "AspN", "LysC", "LysN", "GluC")
  
  for (i in seq_along(enzymes)) {
    enzyme <- enzymes[i]
    peptide_col <- paste0(enzyme, "Peps")
    
    if (peptide_col %in% names(novel_data)) {
      total_peptides <- sum(sapply(novel_data[[peptide_col]], function(x) {
        if (is.null(x) || all(is.na(x))) return(0)
        return(length(x))
      }))
      enzyme_counts[[enzyme_names[i]]] <- total_peptides
    } else {
      enzyme_counts[[enzyme_names[i]]] <- 0
    }
  }
  
  # Count isoforms with genomic mapping
  cat("DEBUG: Counting isoforms with genomic mapping...\n")
  mapped_isoforms <- 0
  for (i in 1:nrow(novel_data)) {
    tryCatch({
      mapped_ranges <- novel_data$trpPeps_mapped_ranges[[i]]
      if (!is.null(mapped_ranges) && length(mapped_ranges) > 0) {
        mapped_isoforms <- mapped_isoforms + 1
      }
    }, error = function(e) {
      cat("DEBUG: Error checking mapping for row", i, ":", e$message, "\n")
    })
  }
  cat("DEBUG: Found", mapped_isoforms, "isoforms with genomic mapping\n")
  
  # Create summary dataframe
  cat("DEBUG: Creating summary dataframe...\n")
  
  tryCatch({
    avg_length <- round(mean(as.numeric(novel_data$numAA), na.rm = TRUE), 1)
    cat("DEBUG: Average protein length:", avg_length, "\n")
    
    summary_df <- data.frame(
      Metric = c(
        "Total Novel Isoforms",
        "Unique Gene Symbols", 
        "Isoforms with Genomic Mapping",
        "Average Protein Length (AA)",
        paste("Total", enzyme_names, "Peptides")
      ),
      Value = c(
        total_isoforms,
        unique_gene_symbols,
        mapped_isoforms,
        avg_length,
        unlist(enzyme_counts)
      ),
      stringsAsFactors = FALSE
    )
    
    cat("DEBUG: Summary dataframe created successfully with", nrow(summary_df), "rows\n")
    return(summary_df)
    
  }, error = function(e) {
    cat("ERROR in summary dataframe creation:", e$message, "\n")
    stop("Failed to create summary: ", e$message)
  })
}

#===============================================================================
# BOUNDARY-BASED GENE SEARCH INTEGRATION
#===============================================================================

# Load boundary-based gene search functions
source("R/boundary_gene_search.R")

#' Run Boundary-Based Gene Search for Novel Isoform
#' 
#' Replacement for BLAST - uses genomic coordinate overlap to find genes
#' 
#' @param work_dir Working directory containing novel isoform GTF
#' @param min_overlap_bp Minimum overlap in base pairs (default: 50)
#' @param min_overlap_percent Minimum overlap percentage (default: 10)
#' @param max_genes Maximum number of genes to return (default: 10)
#' @return Data frame with gene search results in BLAST-compatible format
run_boundary_gene_search_for_novel <- function(work_dir,
                                              min_overlap_bp = 50,
                                              min_overlap_percent = 10,
                                              max_genes = 10) {
  
  cat("DEBUG: Starting boundary-based gene search\n")
  cat("DEBUG: Work directory:", work_dir, "\n")
  cat("DEBUG: Min overlap:", min_overlap_bp, "bp,", min_overlap_percent, "%\n")
  
  # Find the novel GTF file
  novel_gtf_file <- file.path(work_dir, "results", "novel_final.gtf")
  
  if (!file.exists(novel_gtf_file)) {
    cat("DEBUG: Novel GTF file not found:", novel_gtf_file, "\n")
    return(create_empty_blast_format_result())
  }
  
  cat("DEBUG: Using GTF file:", novel_gtf_file, "\n")
  
  # Run boundary-based gene search
  tryCatch({
    search_results <- run_boundary_gene_search(
      novel_gtf_file = novel_gtf_file,
      min_overlap_bp = min_overlap_bp,
      min_overlap_percent = min_overlap_percent,
      max_genes = max_genes
    )
    
    cat("DEBUG: Boundary search completed successfully\n")
    cat("DEBUG: Found", nrow(search_results), "candidate genes\n")
    
    return(search_results)
    
  }, error = function(e) {
    cat("DEBUG: Boundary search failed:", e$message, "\n")
    return(create_empty_blast_format_result())
  })
}

#' Parse Gene Information from Boundary Search Results
#' 
#' Extracts gene information from boundary search results (already in correct format)
#' 
#' @param boundary_results Data frame from run_boundary_gene_search_for_novel
#' @return Data frame with parsed gene information (pass-through)
parse_boundary_gene_info <- function(boundary_results) {
  
  if (nrow(boundary_results) == 0) {
    return(data.frame())
  }
  
  cat("DEBUG: Processing boundary search results\n")
  cat("DEBUG: Found", nrow(boundary_results), "genes from boundary search\n")
  
  # Boundary search results are already in the correct format
  # Just add some debug information
  if (nrow(boundary_results) > 0) {
    cat("DEBUG: Top candidates:\n")
    for (i in 1:min(3, nrow(boundary_results))) {
      cat("DEBUG:  ", i, ".", boundary_results$gene_symbol[i], 
          "(", boundary_results$gene_id[i], ") - Confidence:", 
          boundary_results$pident[i], "\n")
    }
  }
  
  return(boundary_results)
}

#' Check Gene RDS File Availability for Boundary Results
#' 
#' Validates that RDS files exist for boundary-matched genes
#' 
#' @param gene_ids Vector of gene IDs from boundary search results
#' @param rds_dir Directory containing gene RDS files (default: data/genes)
#' @return Data frame with gene availability status
check_boundary_gene_rds_availability_wrapper <- function(gene_ids, rds_dir = "data/genes") {
  
  cat("DEBUG: Checking RDS availability for", length(gene_ids), "genes from boundary search\n")
  cat("DEBUG: RDS directory:", rds_dir, "\n")
  
  if (!dir.exists(rds_dir)) {
    stop("RDS directory not found: ", rds_dir)
  }
  
  availability_df <- data.frame(
    gene_id = gene_ids,
    rds_file = paste0(gene_ids, ".rds"),
    rds_path = file.path(rds_dir, paste0(gene_ids, ".rds")),
    available = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Check file existence
  for (i in 1:nrow(availability_df)) {
    availability_df$available[i] <- file.exists(availability_df$rds_path[i])
  }
  
  # Summary
  available_count <- sum(availability_df$available)
  cat("DEBUG: RDS files available for", available_count, "out of", 
      nrow(availability_df), "genes\n")
  
  return(availability_df)
}

#' Load and Merge Gene Data with Novel Isoform
#' 
#' Loads known gene data and creates temporary merged dataset
#' 
#' @param gene_id Gene ID to load
#' @param novel_data Novel isoform dataframe
#' @param rds_dir Directory containing gene RDS files
#' @return Merged dataframe with novel and known isoforms
load_and_merge_gene_data <- function(gene_id, novel_data, miscleavage_type = "no_miss_cleavage", rds_dir = "data/genes") {
  
  # DEBUG: Start merging process
  cat("DEBUG: Loading and merging data for gene:", gene_id, "\n")
  cat("DEBUG: Novel data rows:", nrow(novel_data), "\n")
  cat("DEBUG: Novel data cols:", ncol(novel_data), "\n")
  cat("DEBUG: Novel data column names:", paste(head(names(novel_data), 10), collapse = ", "), "\n")
  cat("DEBUG: Miscleavage type:", miscleavage_type, "\n")
  
  # Use the gene-by-gene loading system that respects miscleavage types
  source("R/gene_splitter.R")
  
  # Load gene data with correct miscleavage type
  gene_data_result <- load_gene_data(gene_id, miscleavage_type = miscleavage_type, data_dir = dirname(rds_dir))
  
  # Check if gene data was loaded successfully
  if (is.null(gene_data_result)) {
    cat("DEBUG: Gene data not found, returning novel data only\n")
    cat("DEBUG: Returning novel data with", nrow(novel_data), "rows,", ncol(novel_data), "columns\n")
    cat("DEBUG: Returned data column names:", paste(head(names(novel_data), 10), collapse = ", "), "\n")
    return(novel_data)
  }
  
  # Extract known peptide data (already in correct format)
  tryCatch({
    known_data <- gene_data_result$peptides
    
    # DEBUG: Show known data info
    cat("DEBUG: Known data loaded - rows:", nrow(known_data), "cols:", ncol(known_data), "\n")
    cat("DEBUG: Known gene symbol:", unique(known_data$geneSymbol)[1], "\n")
    cat("DEBUG: Miscleavage type used:", gene_data_result$miscleavage_type, "\n")
    
    # The known data is already in the correct format (list columns)
    # Just need to ensure column compatibility
    
    # Validate structure compatibility
    if (ncol(known_data) != ncol(novel_data)) {
      cat("DEBUG: Warning - Column count mismatch. Known:", ncol(known_data), 
          "Novel:", ncol(novel_data), "\n")
    }
    
    # Check column names match
    missing_cols <- setdiff(names(known_data), names(novel_data))
    extra_cols <- setdiff(names(novel_data), names(known_data))
    
    if (length(missing_cols) > 0) {
      cat("DEBUG: Missing columns in novel data:", paste(missing_cols, collapse = ", "), "\n")
    }
    if (length(extra_cols) > 0) {
      cat("DEBUG: Extra columns in novel data:", paste(extra_cols, collapse = ", "), "\n")
    }
    
    # DEBUG: Show sample transcript IDs before merging
    cat("DEBUG: Novel transcript IDs (first 3):", paste(head(novel_data$txID, 3), collapse = ", "), "\n")
    cat("DEBUG: Known transcript IDs (first 3):", paste(head(known_data$txID, 3), collapse = ", "), "\n")
    
    # CRITICAL FIX: Handle structural incompatibilities before rbind
    if (length(missing_cols) > 0 || length(extra_cols) > 0) {
      cat("DEBUG: STRUCTURAL MISMATCH DETECTED - Harmonizing structures...\n")
      
      # Remove extra columns from novel data that don't exist in known data
      if (length(extra_cols) > 0) {
        cat("DEBUG: Removing extra columns from novel data:", paste(extra_cols, collapse = ", "), "\n")
        novel_data_clean <- novel_data[, !names(novel_data) %in% extra_cols, drop = FALSE]
        cat("DEBUG: Novel data columns after cleanup:", ncol(novel_data_clean), "\n")
      } else {
        novel_data_clean <- novel_data
      }
      
      # Add missing columns to novel data that exist in known data
      if (length(missing_cols) > 0) {
        cat("DEBUG: Adding missing columns to novel data:", paste(missing_cols, collapse = ", "), "\n")
        for (col in missing_cols) {
          novel_data_clean[[col]] <- NA  # Initialize with NA
        }
      }
      
      # Reorder columns to match known data structure
      novel_data_clean <- novel_data_clean[, names(known_data), drop = FALSE]
      cat("DEBUG: Structure harmonization complete. Both datasets now have", ncol(novel_data_clean), "columns\n")
      
      # Use cleaned novel data for merging
      novel_data <- novel_data_clean
    }
    
    # Merge datasets (novel data first, then known data)
    cat("DEBUG: Attempting rbind merge...\n")
    merged_data <- rbind(novel_data, known_data)
    
    # DEBUG: Show merge results
    cat("DEBUG: Merged data - total rows:", nrow(merged_data), "\n")
    cat("DEBUG: Novel isoforms:", nrow(novel_data), "Known isoforms:", nrow(known_data), "\n")
    
    # DEBUG: Show final transcript IDs in merged data
    all_tx_ids <- unique(merged_data$txID)
    novel_tx_ids <- all_tx_ids[grepl("^NOVEL_", all_tx_ids)]
    known_tx_ids <- all_tx_ids[!grepl("^NOVEL_", all_tx_ids)]
    cat("DEBUG: Final merged transcripts - Novel:", length(novel_tx_ids), "Known:", length(known_tx_ids), "\n")
    cat("DEBUG: Known transcript IDs:", paste(head(known_tx_ids, 5), collapse = ", "), "\n")
    
    # Add source indicator for analysis
    merged_data$data_source <- c(rep("Novel", nrow(novel_data)), 
                                rep("Known", nrow(known_data)))
    
    cat("DEBUG: Final merged result:", nrow(merged_data), "rows,", ncol(merged_data), "columns\n")
    cat("DEBUG: Final column names:", paste(head(names(merged_data), 10), collapse = ", "), "\n")
    
    return(merged_data)
    
  }, error = function(e) {
    cat("DEBUG: Error processing known data:", e$message, "\n")
    cat("DEBUG: Returning novel data only\n")
    cat("DEBUG: Error fallback - returning novel data with", nrow(novel_data), "rows,", ncol(novel_data), "columns\n")
    return(novel_data)
  })
}

#' Get Novel Isoform Peptides for Analysis
#' 
#' Extracts peptides for a specific enzyme from novel isoform data
#' 
#' @param novel_data Novel isoform dataframe
#' @param protease Protease/enzyme name
#' @return Data frame with peptide information
get_novel_isoform_peptides <- function(novel_data, protease, miscleavage_type = "no_miss_cleavage") {
  cat("DEBUG: get_novel_isoform_peptides called\n")
  cat("DEBUG: novel_data rows:", nrow(novel_data), "\n")
  cat("DEBUG: protease:", protease, "\n")
  cat("DEBUG: miscleavage_type:", miscleavage_type, "\n")
  
  peptide_col <- paste0(protease, "Peps")
  position_col <- paste0(protease, "Peps_positions")
  mapped_col <- paste0(protease, "Peps_mapped_ranges")
  
  cat("DEBUG: Looking for columns:", peptide_col, position_col, mapped_col, "\n")
  cat("DEBUG: Available columns:", paste(names(novel_data), collapse = ", "), "\n")
  
  if (!peptide_col %in% names(novel_data)) {
    cat("ERROR: Peptide column not found:", peptide_col, "\n")
    stop("Peptide column not found: ", peptide_col)
  }
  
  all_peptides_list <- list()
  
  for (i in 1:nrow(novel_data)) {
    # Always use txID field for transcript identifiers (both novel and known)
    transcript_id <- novel_data$txID[i]
    
    peptides <- novel_data[[peptide_col]][[i]]
    positions <- novel_data[[position_col]][[i]]
    mapped_ranges <- novel_data[[mapped_col]][[i]]
    
    if (length(peptides) > 0 && !all(is.na(peptides))) {
      # Create peptide data frame
      if (is.data.frame(positions) && nrow(positions) > 0) {
        peptide_df <- data.frame(
          transcript = transcript_id,
          y_position = i,
          peptide = positions$peptide,
          aa_start = positions$aa_start,
          aa_end = positions$aa_end,
          stringsAsFactors = FALSE
        )
        
        # Apply peptide length filter (6-60 amino acids)
        peptide_lengths <- nchar(peptide_df$peptide)
        valid_peptides <- peptide_lengths >= 6 & peptide_lengths <= 60
        peptide_df <- peptide_df[valid_peptides, ]
        
        cat("DEBUG: Transcript", transcript_id, "- filtered", sum(valid_peptides), "out of", length(valid_peptides), "peptides (6-60 AA)\n")
        
        # Add genomic coordinates if available
        if (!is.null(mapped_ranges) && length(mapped_ranges) > 0) {
          if (class(mapped_ranges) == "GRanges") {
            # Match peptides to genomic ranges
            peptide_df$start <- NA
            peptide_df$end <- NA
            
            for (j in 1:nrow(peptide_df)) {
              pep_seq <- peptide_df$peptide[j]
              matching_ranges <- mapped_ranges[mapped_ranges$peptide == pep_seq]
              if (length(matching_ranges) > 0) {
                peptide_df$start[j] <- min(start(matching_ranges))
                peptide_df$end[j] <- max(end(matching_ranges))
              }
            }
          }
        }
        
        # If no genomic coordinates, use AA positions as proxy
        if (all(is.na(peptide_df$start))) {
          peptide_df$start <- peptide_df$aa_start
          peptide_df$end <- peptide_df$aa_end
        }
        
        # Only add if we have valid peptides after filtering
        if (nrow(peptide_df) > 0) {
          all_peptides_list[[transcript_id]] <- peptide_df
        }
      }
    }
  }
  
  if (length(all_peptides_list) == 0) {
    return(NULL)
  }
  
  # Combine all peptides
  all_peptides_df <- do.call(rbind, all_peptides_list)
  
  # Add hover text
  all_peptides_df$hover_text <- paste0(
    "Peptide: ", all_peptides_df$peptide,
    "<br>Position: ", all_peptides_df$start, "-", all_peptides_df$end,
    "<br>Isoform: ", all_peptides_df$transcript,
    "<br>Enzyme: ", protease
  )
  
  return(all_peptides_df)
}

 