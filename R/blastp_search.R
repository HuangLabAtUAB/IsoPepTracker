# BLASTP Search Integration
# This module handles BLASTP-based peptide searching against the protein database

#===============================================================================
# BLASTP SEARCH FUNCTIONS
#===============================================================================

#' Run BLASTP Search for Peptide Query
#' 
#' Executes BLASTP search against the gencode protein database
#' 
#' @param peptide_query Peptide sequence to search
#' @param evalue E-value threshold for BLAST hits (default: 10)
#' @param max_target_seqs Maximum number of target sequences to return (default: 500)
#' @param identity_threshold Minimum identity percentage (default: 30)
#' @param progress_callback Optional progress callback function
#' @return List containing BLAST results and metadata
run_blastp_peptide_search <- function(peptide_query,
                                    evalue = 10,
                                    max_target_seqs = 500,
                                    identity_threshold = 0,
                                    progress_callback = NULL) {
  
  # Input validation
  if (is.null(peptide_query) || nchar(peptide_query) == 0) {
    return(list(
      success = FALSE,
      error = "Empty peptide query provided",
      results = data.frame()
    ))
  }
  
  # Clean peptide sequence (remove non-amino acid characters)
  peptide_clean <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", toupper(peptide_query))
  
  # Debug: Show what peptide we're actually searching for
  cat("DEBUG: Input peptide query:\n")
  cat("  - Original:", peptide_query, "\n")
  cat("  - Cleaned:", peptide_clean, "\n")
  cat("  - Length:", nchar(peptide_clean), "amino acids\n")
  
  if (nchar(peptide_clean) < 3) {
    return(list(
      success = FALSE,
      error = "Peptide query too short (minimum 3 amino acids)",
      results = data.frame()
    ))
  }
  
  if (!is.null(progress_callback)) {
    progress_callback("Preparing BLASTP search...", 0.1)
  }
  
  # Create temporary files
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S_%f")
  temp_dir <- file.path(tempdir(), paste0("blastp_search_", timestamp))
  dir.create(temp_dir, recursive = TRUE)
  
  query_file <- file.path(temp_dir, "query.fasta")
  output_file <- file.path(temp_dir, "blast_results.txt")
  
  tryCatch({
    # Write query sequence to FASTA file
    query_fasta <- paste0(">peptide_query\n", peptide_clean, "\n")
    writeLines(query_fasta, query_file)
    
    # Debug: Verify query file contents
    cat("DEBUG: Query FASTA file contents:\n")
    cat(paste(readLines(query_file), collapse = "\n"), "\n")
    
    if (!is.null(progress_callback)) {
      progress_callback("Running BLASTP search...", 0.3)
    }
    
    # Path to BLAST binary and database - using absolute paths like Novel tab
    blastp_path <- "/Users/Mahmuda/opt/anaconda3/envs/module3_tools/bin/blastp"
    
    # Get absolute path to database (same pattern as Novel tab's reference handling)
    app_dir <- getwd()
    db_path <- file.path(app_dir, "reference", "gencode_v38_pc")
    
    # Debug information
    cat("DEBUG: BLASTP Configuration\n")
    cat("  - App working directory:", app_dir, "\n")
    cat("  - BLASTP binary:", blastp_path, "\n")
    cat("  - Database path:", db_path, "\n")
    cat("  - Query file:", query_file, "\n")
    cat("  - Output file:", output_file, "\n")
    cat("  - Temp directory:", temp_dir, "\n")
    
    # Validate paths before execution
    if (!file.exists(blastp_path)) {
      return(list(
        success = FALSE,
        error = paste("BLASTP binary not found at:", blastp_path),
        results = data.frame()
      ))
    }
    
    # Check database files exist
    required_db_files <- c(
      paste0(db_path, ".phr"),
      paste0(db_path, ".pin"), 
      paste0(db_path, ".psq")
    )
    
    missing_files <- c()
    for (db_file in required_db_files) {
      if (!file.exists(db_file)) {
        missing_files <- c(missing_files, db_file)
      }
    }
    
    if (length(missing_files) > 0) {
      return(list(
        success = FALSE,
        error = paste("Missing BLAST database files:", paste(missing_files, collapse = ", ")),
        results = data.frame()
      ))
    }
    
    # Build BLASTP command with ultra-permissive parameters for testing
    blastp_args <- c(
      "-query", query_file,
      "-db", db_path,
      "-out", output_file,
      "-outfmt", "6",
      "-evalue", "10",  # Very permissive E-value
      "-max_target_seqs", as.character(max_target_seqs),
      "-word_size", "2",  # Small word size for short peptides
      "-threshold", "11",  # Lower threshold for better sensitivity
      "-task", "blastp-short",  # Optimized for short sequences
      "-comp_based_stats", "0",  # Disable composition-based statistics
      "-soft_masking", "false"  # Disable soft masking
    )
    
    # Debug: show full command
    full_command <- paste(blastp_path, paste(blastp_args, collapse = " "))
    cat("DEBUG: Full BLASTP command:", full_command, "\n")
    
    # Execute BLASTP following Novel tab pattern
    blast_result <- system2(
      command = blastp_path, 
      args = blastp_args, 
      wait = TRUE,
      stdout = TRUE,
      stderr = TRUE
    )
    
    cat("DEBUG: BLASTP system2() result:\n")
    if (is.character(blast_result)) {
      cat("  - STDOUT/STDERR output:", paste(blast_result, collapse = "\n"), "\n")
    } else {
      cat("  - Exit code:", blast_result, "\n")
    }
    cat("  - Output file exists:", file.exists(output_file), "\n")
    if (file.exists(output_file)) {
      cat("  - Output file size:", file.size(output_file), "bytes\n")
    }
    
    # Check BLASTP execution results
    exit_code <- 0
    stderr_output <- ""
    
    if (is.character(blast_result)) {
      stderr_output <- paste(blast_result, collapse = "\n")
    } else {
      exit_code <- blast_result
    }
    
    # Check for actual system errors (non-zero exit code or missing output file)
    if (exit_code != 0 || !file.exists(output_file)) {
      error_msg <- paste(
        "BLASTP execution failed.",
        "Command:", full_command,
        "Working dir:", getwd(),
        "Exit code:", exit_code,
        if (nchar(stderr_output) > 0) paste("STDERR:", stderr_output) else "",
        sep = "\n"
      )
      
      return(list(
        success = FALSE,
        error = error_msg,
        results = data.frame()
      ))
    }
    
    # If output file exists but is empty, this means no matches found (not an error)
    if (file.size(output_file) == 0) {
      cat("DEBUG: No matches found by BLAST - output file is empty (0 bytes)\n")
      
      return(list(
        success = TRUE,
        query = peptide_clean,
        total_hits = 0,
        results = data.frame(
          query_id = character(0),
          subject_id = character(0),
          gene_id = character(0),
          gene_symbol = character(0),
          transcript_id = character(0),
          protein_id = character(0),
          identity_percent = numeric(0),
          alignment_length = numeric(0),
          evalue = numeric(0),
          bit_score = numeric(0),
          stringsAsFactors = FALSE
        ),
        message = paste("No matches found for peptide", peptide_clean, "in the protein database."),
        search_parameters = list(
          evalue = evalue,
          max_target_seqs = max_target_seqs,
          identity_threshold = identity_threshold
        )
      ))
    }
    
    if (!is.null(progress_callback)) {
      progress_callback("Processing BLAST results...", 0.7)
    }
    
    # Debug: Show raw BLAST output before parsing
    if (file.exists(output_file) && file.size(output_file) > 0) {
      blast_raw_lines <- readLines(output_file)
      cat("DEBUG: Raw BLAST output (first 10 lines):\n")
      cat(paste(head(blast_raw_lines, 10), collapse = "\n"), "\n")
      cat("DEBUG: Total raw BLAST lines:", length(blast_raw_lines), "\n")
    }
    
    # Parse BLAST output
    blast_results <- parse_blastp_output(output_file, identity_threshold, peptide_clean)
    
    if (!is.null(progress_callback)) {
      progress_callback("Mapping to gene information...", 0.9)
    }
    
    # Map protein IDs to gene information
    enhanced_results <- map_proteins_to_genes(blast_results)
    
    if (!is.null(progress_callback)) {
      progress_callback("Search completed!", 1.0)
    }
    
    return(list(
      success = TRUE,
      query = peptide_clean,
      total_hits = nrow(enhanced_results),
      results = enhanced_results,
      search_parameters = list(
        evalue = evalue,
        max_target_seqs = max_target_seqs,
        identity_threshold = identity_threshold
      )
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      error = paste("BLASTP search failed:", e$message),
      results = data.frame()
    ))
  }, finally = {
    # Clean up temporary files
    if (dir.exists(temp_dir)) {
      unlink(temp_dir, recursive = TRUE)
    }
  })
}

#' Parse BLASTP Output
#' 
#' Parses tabular BLASTP output into a data frame
#' 
#' @param output_file Path to BLAST output file
#' @param identity_threshold Minimum identity percentage to include
#' @return Data frame with parsed BLAST results
parse_blastp_output <- function(output_file, identity_threshold = 70, query_peptide = "") {
  
  if (!file.exists(output_file) || file.size(output_file) == 0) {
    return(data.frame(
      query_id = character(0),
      subject_id = character(0),
      identity_percent = numeric(0),
      alignment_length = numeric(0),
      mismatches = numeric(0),
      gap_opens = numeric(0),
      query_start = numeric(0),
      query_end = numeric(0),
      subject_start = numeric(0),
      subject_end = numeric(0),
      evalue = numeric(0),
      bit_score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Read BLAST results
  blast_lines <- readLines(output_file)
  
  if (length(blast_lines) == 0) {
    return(data.frame(
      query_id = character(0),
      subject_id = character(0),
      identity_percent = numeric(0),
      alignment_length = numeric(0),
      mismatches = numeric(0),
      gap_opens = numeric(0),
      query_start = numeric(0),
      query_end = numeric(0),
      subject_start = numeric(0),
      subject_end = numeric(0),
      evalue = numeric(0),
      bit_score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Parse each line (default BLAST tabular format has 12 columns)
  results_list <- list()
  
  for (i in seq_along(blast_lines)) {
    line_parts <- strsplit(blast_lines[i], "\t")[[1]]
    
    if (length(line_parts) >= 12) {
      identity_pct <- as.numeric(line_parts[3])
      
      # Apply identity threshold filter
      if (identity_pct >= identity_threshold) {
        results_list[[i]] <- data.frame(
          query_id = line_parts[1],
          subject_id = line_parts[2],
          identity_percent = identity_pct,
          alignment_length = as.numeric(line_parts[4]),
          mismatches = as.numeric(line_parts[5]),
          gap_opens = as.numeric(line_parts[6]),
          query_start = as.numeric(line_parts[7]),
          query_end = as.numeric(line_parts[8]),
          subject_start = as.numeric(line_parts[9]),
          subject_end = as.numeric(line_parts[10]),
          evalue = as.numeric(line_parts[11]),
          bit_score = as.numeric(line_parts[12]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(results_list) == 0) {
    return(data.frame(
      query_id = character(0),
      subject_id = character(0),
      identity_percent = numeric(0),
      alignment_length = numeric(0),
      mismatches = numeric(0),
      gap_opens = numeric(0),
      query_start = numeric(0),
      query_end = numeric(0),
      subject_start = numeric(0),
      subject_end = numeric(0),
      evalue = numeric(0),
      bit_score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Combine results
  results_df <- do.call(rbind, results_list)
  
  # Sort by bit score (descending) and E-value (ascending)
  results_df <- results_df[order(-results_df$bit_score, results_df$evalue), ]
  
  return(results_df)
}

#' Map Protein IDs to Gene Information
#' 
#' Maps BLAST protein hits to gene and transcript information
#' 
#' @param blast_results Data frame with BLAST results
#' @return Enhanced data frame with gene information
map_proteins_to_genes <- function(blast_results) {
  
  if (nrow(blast_results) == 0) {
    return(data.frame(
      query_id = character(0),
      subject_id = character(0),
      gene_id = character(0),
      gene_symbol = character(0),
      transcript_id = character(0),
      protein_id = character(0),
      identity_percent = numeric(0),
      alignment_length = numeric(0),
      evalue = numeric(0),
      bit_score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Extract protein/transcript IDs from subject IDs
  # Format is typically: ENST00000123456.1|ENSP00000123456.1|ENSG00000123456.1|OTTHUMT00000123456.1|SYMBOL-001|7|
  
  enhanced_results <- blast_results
  enhanced_results$gene_id <- NA
  enhanced_results$gene_symbol <- NA
  enhanced_results$transcript_id <- NA
  enhanced_results$protein_id <- NA
  
  for (i in 1:nrow(blast_results)) {
    subject_id <- blast_results$subject_id[i]
    
    # Parse the complex protein ID format
    id_parts <- strsplit(subject_id, "\\|")[[1]]
    
    if (length(id_parts) >= 7) {
      # Extract IDs - Format: ENSP|ENST|ENSG|OTTHUMG|OTTHUMT|TRANSCRIPT-NAME|GENE_SYMBOL|LENGTH
      protein_id <- id_parts[1]      # ENSP (protein ID)
      transcript_id <- id_parts[2]   # ENST (transcript ID)
      gene_id <- id_parts[3]         # ENSG (gene ID)
      gene_symbol <- id_parts[7]     # Gene symbol (OR4F5)
      
      enhanced_results$transcript_id[i] <- transcript_id
      enhanced_results$protein_id[i] <- protein_id
      enhanced_results$gene_id[i] <- gene_id
      enhanced_results$gene_symbol[i] <- gene_symbol
    } else {
      # Fallback parsing for simpler formats
      enhanced_results$protein_id[i] <- subject_id
      enhanced_results$transcript_id[i] <- subject_id
      enhanced_results$gene_id[i] <- "Unknown"
      enhanced_results$gene_symbol[i] <- "Unknown"
    }
  }
  
  # Remove hits without proper gene information
  enhanced_results <- enhanced_results[!is.na(enhanced_results$gene_id) & 
                                     enhanced_results$gene_id != "Unknown", ]
  
  # Select and reorder columns for final output
  final_results <- enhanced_results[, c(
    "query_id", "subject_id", "gene_id", "gene_symbol", "transcript_id", "protein_id",
    "identity_percent", "alignment_length", "evalue", "bit_score"
  )]
  
  return(final_results)
}

#' Validate BLAST Database
#' 
#' Checks if BLAST protein database files exist and are accessible
#' 
#' @param db_path Path to BLAST database (without extension)
#' @return List with validation results
validate_blast_database <- function(db_path = NULL) {
  
  # Use same path resolution as main BLASTP function
  if (is.null(db_path)) {
    app_dir <- getwd()
    db_path <- file.path(app_dir, "reference", "gencode_v38_pc")
  }
  
  required_extensions <- c(".phr", ".pin", ".psq")
  existing_files <- c()
  missing_files <- c()
  
  # Debug information
  cat("DEBUG: Database validation\n")
  cat("  - App working directory:", getwd(), "\n")
  cat("  - Database base path:", db_path, "\n")
  
  for (ext in required_extensions) {
    file_path <- paste0(db_path, ext)
    cat("  - Checking file:", file_path, "- exists:", file.exists(file_path), "\n")
    
    if (file.exists(file_path)) {
      existing_files <- c(existing_files, file_path)
    } else {
      missing_files <- c(missing_files, file_path)
    }
  }
  
  is_valid <- length(missing_files) == 0
  
  return(list(
    valid = is_valid,
    db_path = db_path,
    existing_files = existing_files,
    missing_files = missing_files,
    message = if (is_valid) {
      paste("BLAST database is valid at:", db_path)
    } else {
      paste("Missing BLAST database files:", paste(missing_files, collapse = ", "))
    }
  ))
}

#' Get BLASTP Search Summary
#' 
#' Creates a summary of BLASTP search results
#' 
#' @param blast_search_result Result from run_blastp_peptide_search
#' @return Character string with search summary
get_blastp_search_summary <- function(blast_search_result) {
  
  if (!blast_search_result$success) {
    return(paste("Search failed:", blast_search_result$error))
  }
  
  total_hits <- blast_search_result$total_hits
  unique_genes <- length(unique(blast_search_result$results$gene_id))
  unique_proteins <- length(unique(blast_search_result$results$protein_id))
  
  if (total_hits == 0) {
    return("No matches found")
  }
  
  # Get statistics from results
  results <- blast_search_result$results
  best_hit_identity <- max(results$identity_percent, na.rm = TRUE)
  best_hit_evalue <- min(results$evalue, na.rm = TRUE)
  
  summary_text <- paste0(
    "Found ", total_hits, " BLASTP hits in ", unique_genes, " genes (",
    unique_proteins, " proteins). ",
    "Best hit: ", round(best_hit_identity, 1), "% identity (E=", 
    format(best_hit_evalue, scientific = TRUE, digits = 2), ")"
  )
  
  return(summary_text)
}

#' Cross-reference BLAST Results with Peptide Databases
#' 
#' Enhances BLAST results by cross-referencing with internal peptide databases
#' to provide additional context and information
#' 
#' @param blast_results Data frame with BLAST results from map_proteins_to_genes
#' @param peptide_query Original peptide query string
#' @return Enhanced data frame with additional database information
cross_reference_blast_with_peptide_databases <- function(blast_results, peptide_query) {
  
  # Input validation
  if (is.null(blast_results) || nrow(blast_results) == 0) {
    return(data.frame(
      query_id = character(0),
      subject_id = character(0),
      gene_id = character(0),
      gene_symbol = character(0),
      transcript_id = character(0),
      protein_id = character(0),
      identity_percent = numeric(0),
      alignment_length = numeric(0),
      evalue = numeric(0),
      bit_score = numeric(0),
      database_source = character(0),
      peptide_context = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Start with the original BLAST results
  enhanced_results <- blast_results
  
  # Add database source information
  enhanced_results$database_source <- "GENCODE_v38_Proteins"
  
  # Add peptide context information
  enhanced_results$peptide_context <- paste0(
    "Query: ", peptide_query, 
    " (Length: ", nchar(peptide_query), " AA)"
  )
  
  # Try to cross-reference with internal peptide databases if they exist
  tryCatch({
    # Check if peptide databases are available in the global environment
    if (exists("peptides") && is.data.frame(get("peptides"))) {
      internal_peptides <- get("peptides")
      
      # Try to find matches in internal database by gene ID
      for (i in 1:nrow(enhanced_results)) {
        gene_id <- enhanced_results$gene_id[i]
        
        # Look for peptides from this gene in the internal database
        gene_peptides <- internal_peptides[internal_peptides$geneID == gene_id, ]
        
        if (nrow(gene_peptides) > 0) {
          # Count how many peptides we have for this gene
          peptide_count <- nrow(gene_peptides)
          transcript_count <- length(unique(gene_peptides$txID))
          
          enhanced_results$peptide_context[i] <- paste0(
            enhanced_results$peptide_context[i],
            " | Internal DB: ", peptide_count, " peptides from ", 
            transcript_count, " transcripts"
          )
        }
      }
      
      cat("DEBUG: Cross-referenced BLAST results with internal peptide database\n")
      cat("  - Internal database contains", nrow(internal_peptides), "peptides\n")
      
    } else {
      cat("DEBUG: No internal peptide database found for cross-referencing\n")
    }
    
    # Check for AS database as well
    if (exists("as_database") && is.data.frame(get("as_database"))) {
      as_db <- get("as_database")
      
      # Try to find AS events for the genes found
      for (i in 1:nrow(enhanced_results)) {
        gene_id <- enhanced_results$gene_id[i]
        
        # Look for AS events involving this gene
        gene_as_events <- as_db[as_db$geneID == gene_id, ]
        
        if (nrow(gene_as_events) > 0) {
          as_event_count <- nrow(gene_as_events)
          enhanced_results$peptide_context[i] <- paste0(
            enhanced_results$peptide_context[i],
            " | AS Events: ", as_event_count
          )
        }
      }
      
      cat("DEBUG: Cross-referenced with AS events database\n")
      cat("  - AS database contains", nrow(as_db), "events\n")
    }
    
  }, error = function(e) {
    cat("DEBUG: Error during cross-referencing:", e$message, "\n")
    # Continue without cross-referencing if there's an error
  })
  
  # Add match quality assessment
  enhanced_results$match_quality <- ifelse(
    enhanced_results$identity_percent >= 100, "Perfect",
    ifelse(enhanced_results$identity_percent >= 90, "Excellent",
           ifelse(enhanced_results$identity_percent >= 80, "Good",
                  ifelse(enhanced_results$identity_percent >= 70, "Fair", "Poor")))
  )
  
  # Sort by identity percentage (descending) and E-value (ascending)
  enhanced_results <- enhanced_results[order(-enhanced_results$identity_percent, enhanced_results$evalue), ]
  
  cat("DEBUG: Enhanced", nrow(enhanced_results), "BLAST results with database cross-referencing\n")
  
  return(enhanced_results)
}