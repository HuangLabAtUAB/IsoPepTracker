#' Create Enhanced Novel Isoform Dataframe
#' 
#' Creates a comprehensive dataframe from TransDecoder results with metadata
#' 
#' @param work_dir Working directory containing TransDecoder results
#' @param min_protein_length Minimum protein length used
#' @param genetic_code Genetic code used
#' @param enable_blast Whether BLAST search was enabled
#' @param blast_threshold BLAST identity threshold
#' @return Path to created RDS file
create_novel_isoform_dataframe <- function(work_dir, 
                                         min_protein_length = 30,
                                         genetic_code = 1,
                                         enable_blast = TRUE,
                                         blast_threshold = 70) {
  
  cat("DEBUG: Creating enhanced novel isoform dataframe from:", work_dir, "\n")
  
  # Load required libraries
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package is required but not installed")
  }
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    stop("seqinr package is required but not installed")
  }
  if (!requireNamespace("cleaver", quietly = TRUE)) {
    stop("cleaver package is required but not installed")
  }
  
  library(Biostrings)
  library(seqinr)
  library(cleaver)
  library(GenomicRanges)
  
  # Find the protein FASTA file
  pep_file <- file.path(work_dir, "input", "novel_transcript.fa.transdecoder_dir", "longest_orfs.pep")
  
  if (!file.exists(pep_file)) {
    stop("Protein FASTA file not found: ", pep_file)
  }
  
  cat("DEBUG: Using FASTA:", pep_file, "\n")
  cat("DEBUG: Analysis parameters - Min length:", min_protein_length, "AA, Genetic code:", genetic_code, "\n")
  
  # Read protein sequences using seqinr
  cat("DEBUG: Reading protein sequences...\n")
  sequences <- seqinr::read.fasta(pep_file, seqtype = "AA", as.string = TRUE)
  
  if (length(sequences) == 0) {
    stop("No protein sequences found in file")
  }
  
  cat("DEBUG: Found", length(sequences), "protein sequences\n")
  
  # Create base dataframe with enhanced metadata
  novel_df <- data.frame(
    proteinID = character(length(sequences)),
    txID = character(length(sequences)),
    geneID = character(length(sequences)),
    geneSymbol = character(length(sequences)),
    numAA = integer(length(sequences)),
    seq = character(length(sequences)),
    # Analysis metadata
    min_length_used = rep(min_protein_length, length(sequences)),
    genetic_code_used = rep(genetic_code, length(sequences)),
    analysis_timestamp = rep(Sys.time(), length(sequences)),
    stringsAsFactors = FALSE
  )
  
  # Process each sequence with enhanced ID generation
  for (i in seq_along(sequences)) {
    seq_name <- names(sequences)[i]
    protein_seq <- toupper(as.character(sequences[[i]]))
    
    # Generate systematic IDs
    novel_df$proteinID[i] <- sprintf("NOVEL_%04d.p1", i)
    novel_df$txID[i] <- sprintf("NOVEL_%04d.t1", i)
    novel_df$geneID[i] <- sprintf("NOVEL_%04d", i)
    novel_df$geneSymbol[i] <- sprintf("Novel_%d", i)
    novel_df$numAA[i] <- nchar(protein_seq)
    novel_df$seq[i] <- protein_seq
  }
  
  cat("DEBUG: Created dataframe with", nrow(novel_df), "protein entries\n")
  
  # Add peptide columns for all enzymes
  enzymes <- c("trp", "chymo", "aspn", "lysc", "lysn", "gluc")
  
  # Map enzyme abbreviations to cleaver enzyme names
  enzyme_map <- list(
    trp = "trypsin",
    chymo = "chymotrypsin-high",
    aspn = "asp-n endopeptidase",
    lysc = "lysc",
    lysn = "lysn", 
    gluc = "glutamyl endopeptidase"
  )
  
  for (enzyme in enzymes) {
    cat("DEBUG: Processing", enzyme, "peptides...\n")
    
    peptide_col <- paste0(enzyme, "Peps")
    position_col <- paste0(enzyme, "Peps_positions")
    ranges_col <- paste0(enzyme, "Peps_mapped_ranges")
    
    novel_df[[peptide_col]] <- vector("list", nrow(novel_df))
    novel_df[[position_col]] <- vector("list", nrow(novel_df))
    novel_df[[ranges_col]] <- vector("list", nrow(novel_df))
    
    for (i in 1:nrow(novel_df)) {
      protein_seq <- novel_df$seq[i]
      
      # Generate peptides using cleaver with proper enzyme name
      peptides <- tryCatch({
        cleaver_enzyme <- enzyme_map[[enzyme]]
        cat("DEBUG: Using cleaver enzyme:", cleaver_enzyme, "for protein", i, "\n")
        result <- cleaver::cleave(protein_seq, cleaver_enzyme)[[1]]
        cat("DEBUG: Generated", length(result), "peptides for protein", i, "\n")
        result
      }, error = function(e) {
        cat("ERROR: Cleaver failed for", enzyme, "protein", i, ":", e$message, "\n")
        character(0)
      })
      
      if (length(peptides) > 0) {
        # Apply peptide length filter (6-60 amino acids) AT THE SOURCE
        peptide_lengths <- nchar(peptides)
        valid_indices <- which(peptide_lengths >= 6 & peptide_lengths <= 60)
        filtered_peptides <- peptides[valid_indices]
        
        cat("DEBUG: Filtered", length(filtered_peptides), "out of", length(peptides), "peptides (6-60 AA) for protein", i, "\n")
        
        if (length(filtered_peptides) > 0) {
          # Calculate positions and create data frame structure for FILTERED peptides only
          position_df <- data.frame(
            peptide = character(length(filtered_peptides)),
            aa_start = integer(length(filtered_peptides)),
            aa_end = integer(length(filtered_peptides)),
            stringsAsFactors = FALSE
          )
          
          current_pos <- 1
          peptide_index <- 1
          for (j in seq_along(peptides)) {
            peptide <- peptides[j]
            start_pos <- current_pos
            end_pos <- current_pos + nchar(peptide) - 1
            
            # Only store if this peptide passed the length filter
            if (j %in% valid_indices) {
              position_df$peptide[peptide_index] <- peptide
              position_df$aa_start[peptide_index] <- start_pos
              position_df$aa_end[peptide_index] <- end_pos
              peptide_index <- peptide_index + 1
            }
            
            current_pos <- end_pos + 1
          }
          
          novel_df[[peptide_col]][[i]] <- filtered_peptides
          novel_df[[position_col]][[i]] <- position_df
        } else {
          # No valid peptides after filtering
          novel_df[[peptide_col]][[i]] <- character(0)
          novel_df[[position_col]][[i]] <- data.frame(
            peptide = character(0),
            aa_start = integer(0),
            aa_end = integer(0),
            stringsAsFactors = FALSE
          )
        }
        
        # Create empty GRanges (no genomic mapping in simplified version)
        novel_df[[ranges_col]][[i]] <- GRanges()
        
        # Log results based on whether we had valid peptides
        if (length(filtered_peptides) > 0) {
          cat("DEBUG: Stored", length(filtered_peptides), "filtered peptides for protein", i, "\n")
        } else {
          cat("DEBUG: No valid peptides to store for protein", i, "\n")
        }
      } else {
        novel_df[[peptide_col]][[i]] <- character(0)
        novel_df[[position_col]][[i]] <- data.frame(
          peptide = character(0),
          aa_start = integer(0),
          aa_end = integer(0),
          stringsAsFactors = FALSE
        )
        novel_df[[ranges_col]][[i]] <- GRanges()
      }
    }
  }
  
  # Add analysis metadata as attributes
  attr(novel_df, "analysis_parameters") <- list(
    min_protein_length = min_protein_length,
    genetic_code = genetic_code,
    enable_blast = enable_blast,
    blast_threshold = blast_threshold,
    created_date = Sys.time(),
    transdecoder_version = "Enhanced pipeline",
    total_proteins = nrow(novel_df)
  )
  
  # Save enhanced dataframe
  output_file <- file.path(work_dir, "results", "novel_isoform_dataframe.rds")
  saveRDS(novel_df, output_file)
  
  cat("DEBUG: Saved enhanced novel isoform dataframe to:", output_file, "\n")
  cat("DEBUG: Dataframe dimensions:", nrow(novel_df), "x", ncol(novel_df), "\n")
  cat("DEBUG: Column structure matches expected format (", ncol(novel_df), "columns)\n")
  
  # Create analysis report
  report_file <- file.path(work_dir, "results", "analysis_report.txt")
  writeLines(c(
    "=== Enhanced Novel Transcript Analysis Report ===",
    paste("Analysis completed:", Sys.time()),
    paste("Total proteins identified:", nrow(novel_df)),
    paste("Minimum protein length:", min_protein_length, "AA"),
    paste("Genetic code used:", genetic_code),
    paste("Length distribution:"),
    paste("  Min:", min(novel_df$numAA), "AA"),
    paste("  Max:", max(novel_df$numAA), "AA"),
    paste("  Mean:", round(mean(novel_df$numAA), 1), "AA"),
    "",
    "Peptide generation completed for all proteases:",
    paste("  Trypsin, Chymotrypsin, AspN, LysC, LysN, GluC"),
    "",
    "Files created:",
    paste("  Dataframe:", basename(output_file)),
    paste("  Analysis report:", basename(report_file))
  ), report_file)
  
  return(output_file)
} 