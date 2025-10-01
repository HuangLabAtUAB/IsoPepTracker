#===============================================================================
# GENE-SPECIFIC RDS FILE SPLITTER
#===============================================================================

#' Split large peptide and AS event RDS files into gene-specific files
#'
#' @param no_miss_peptides_file Path to the no miscleavage peptides RDS file
#' @param upto_two_miss_peptides_file Path to the up to two miscleavage peptides RDS file
#' @param as_events_file Path to the large AS events RDS file
#' @param output_dir Directory where gene-specific files will be stored
#' @param create_index Whether to create an index file (default: TRUE)
#' @return Invisibly returns a list with information about the created files
#'
split_rds_by_gene <- function(no_miss_peptides_file, upto_two_miss_peptides_file, as_events_file, output_dir = "data", create_index = TRUE) {
  # Load the large files
  message("Loading no miscleavage peptides file: ", no_miss_peptides_file)
  no_miss_peptides <- readRDS(no_miss_peptides_file)
  
  message("Loading up to two miscleavage peptides file: ", upto_two_miss_peptides_file)
  upto_two_miss_peptides <- readRDS(upto_two_miss_peptides_file)
  
  message("Loading AS events file: ", as_events_file)
  as_database <- readRDS(as_events_file)
  
  # Create output directories if they don't exist
  no_miss_dir <- file.path(output_dir, "genes", "no_miss_cleavage")
  upto_two_miss_dir <- file.path(output_dir, "genes", "upto_two_misscleavage")
  as_events_dir <- file.path(output_dir, "as_events")
  index_dir <- file.path(output_dir, "index")
  
  for (dir in c(no_miss_dir, upto_two_miss_dir, as_events_dir, index_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
  
  # Get unique gene IDs from both datasets
  all_gene_ids <- unique(c(no_miss_peptides$geneID, upto_two_miss_peptides$geneID))
  message("Found ", length(all_gene_ids), " unique genes across both miscleavage datasets")
  
  # Initialize index data frame
  index_data <- data.frame(
    geneID = character(),
    geneSymbol = character(),
    as_events_file = character(),
    num_transcripts = integer(),
    num_as_events = integer(),
    has_no_miss_cleavage = logical(),
    has_upto_two_misscleavage = logical(),
    num_peptides_no_miss = integer(),
    num_peptides_upto_two_miss = integer(),
    stringsAsFactors = FALSE
  )
  
  # Process each gene
  for (i in seq_along(all_gene_ids)) {
    gene_id <- all_gene_ids[i]
    
    # Show progress
    if (i %% 100 == 0 || i == length(all_gene_ids)) {
      message("Processing gene ", i, " of ", length(all_gene_ids), " (", round(i/length(all_gene_ids)*100), "%)")
    }
    
    # Extract gene-specific data
    gene_no_miss <- no_miss_peptides[no_miss_peptides$geneID == gene_id, ]
    gene_upto_two_miss <- upto_two_miss_peptides[upto_two_miss_peptides$geneID == gene_id, ]
    gene_as_events <- as_database[as_database$geneID == gene_id, ]
    
    # Determine gene symbol (try from either dataset)
    gene_symbol <- NA
    if (nrow(gene_no_miss) > 0) {
      gene_symbol <- unique(gene_no_miss$geneSymbol)[1]
    } else if (nrow(gene_upto_two_miss) > 0) {
      gene_symbol <- unique(gene_upto_two_miss$geneSymbol)[1]
    }
    if (is.na(gene_symbol)) gene_symbol <- gene_id
    
    # Create file paths
    no_miss_file <- file.path(no_miss_dir, paste0(gene_id, ".rds"))
    upto_two_miss_file <- file.path(upto_two_miss_dir, paste0(gene_id, ".rds"))
    as_events_file <- file.path(as_events_dir, paste0(gene_id, ".rds"))
    
    # Track availability and counts
    has_no_miss <- nrow(gene_no_miss) > 0
    has_upto_two_miss <- nrow(gene_upto_two_miss) > 0
    num_peptides_no_miss <- nrow(gene_no_miss)
    num_peptides_upto_two_miss <- nrow(gene_upto_two_miss)
    
    # Save gene-specific files (only if data exists)
    if (has_no_miss) {
      saveRDS(gene_no_miss, no_miss_file)
    }
    if (has_upto_two_miss) {
      saveRDS(gene_upto_two_miss, upto_two_miss_file)
    }
    
    # AS events are shared - save once per gene
    saveRDS(gene_as_events, as_events_file)
    
    # Get transcript count from whichever dataset has data
    num_transcripts <- 0
    if (has_no_miss) {
      num_transcripts <- length(unique(gene_no_miss$txID))
    } else if (has_upto_two_miss) {
      num_transcripts <- length(unique(gene_upto_two_miss$txID))
    }
    
    # Add to index
    index_data <- rbind(index_data, data.frame(
      geneID = gene_id,
      geneSymbol = gene_symbol,
      as_events_file = as_events_file,
      num_transcripts = num_transcripts,
      num_as_events = nrow(gene_as_events),
      has_no_miss_cleavage = has_no_miss,
      has_upto_two_misscleavage = has_upto_two_miss,
      num_peptides_no_miss = num_peptides_no_miss,
      num_peptides_upto_two_miss = num_peptides_upto_two_miss,
      stringsAsFactors = FALSE
    ))
  }
  
  # Save index file
  if (create_index) {
    index_file <- file.path(index_dir, "gene_index.rds")
    saveRDS(index_data, index_file)
    message("Created index file: ", index_file)
    message("Summary:")
    message("  - Genes with no miscleavage data: ", sum(index_data$has_no_miss_cleavage))
    message("  - Genes with up to two miscleavage data: ", sum(index_data$has_upto_two_misscleavage))
    message("  - Genes with both: ", sum(index_data$has_no_miss_cleavage & index_data$has_upto_two_misscleavage))
  }
  
  # Return information about the created files
  invisible(list(
    num_genes = length(all_gene_ids),
    index_data = index_data
  ))
}

#' Find gene file with flexible version handling
#'
#' @param gene_id Gene ID (with or without version number)
#' @param miscleavage_type Type of miscleavage data to load
#' @param data_dir Base directory for gene data
#' @return The actual filename found, or NULL if not found
#'
find_gene_file_with_version <- function(gene_id, miscleavage_type = "no_miss_cleavage", data_dir = "data") {
  # First, try exact match (preserves existing behavior)
  exact_file <- file.path(data_dir, "genes", miscleavage_type, paste0(gene_id, ".rds"))
  if (file.exists(exact_file)) {
    return(exact_file)
  }
  
  # If exact match fails, try version-agnostic search
  # Strip version number if present
  base_gene_id <- sub("\\.[0-9]+$", "", gene_id)
  
  # Search directory for files matching the base pattern
  search_dir <- file.path(data_dir, "genes", miscleavage_type)
  if (!dir.exists(search_dir)) {
    return(NULL)
  }
  
  # Look for files that match the base gene ID with any version
  pattern <- paste0("^", base_gene_id, "\\.[0-9]+\\.rds$")
  all_files <- list.files(search_dir, pattern = pattern, full.names = TRUE)
  
  if (length(all_files) > 0) {
    # Return the first match (or could implement logic to pick latest version)
    cat("DEBUG: Found gene file with version:", basename(all_files[1]), "for requested gene:", gene_id, "\n")
    return(all_files[1])
  }
  
  # No file found
  return(NULL)
}

#' Load gene-specific data with flexible version handling
#'
#' @param gene_id Gene ID to load (with or without version number)
#' @param miscleavage_type Type of miscleavage data to load ("no_miss_cleavage" or "upto_two_misscleavage")
#' @param data_dir Base directory for gene data (default: "data")
#' @return List containing peptides and AS events for the gene
#'
load_gene_data <- function(gene_id, miscleavage_type = "no_miss_cleavage", data_dir = "data") {
  # Validate miscleavage_type
  if (!miscleavage_type %in% c("no_miss_cleavage", "upto_two_misscleavage")) {
    stop("miscleavage_type must be either 'no_miss_cleavage' or 'upto_two_misscleavage'")
  }
  
  # Use flexible file lookup for peptide file
  peptide_file <- find_gene_file_with_version(gene_id, miscleavage_type, data_dir)
  if (is.null(peptide_file)) {
    warning("Peptide file not found for gene ", gene_id, " with miscleavage type ", miscleavage_type)
    return(NULL)
  }
  
  # For AS events, use flexible lookup based on the found peptide file's gene ID
  # Extract the actual gene ID with version from the found peptide file
  actual_gene_id <- sub("\\.rds$", "", basename(peptide_file))
  as_events_file <- file.path(data_dir, "as_events", paste0(actual_gene_id, ".rds"))
  
  # Load peptides
  gene_peptides <- readRDS(peptide_file)
  
  # Load AS events if available (shared between miscleavage types)
  gene_as_events <- NULL
  if (file.exists(as_events_file)) {
    gene_as_events <- readRDS(as_events_file)
  }
  
  # Return as list
  return(list(
    peptides = gene_peptides,
    as_events = gene_as_events,
    miscleavage_type = miscleavage_type
  ))
}

#' Load gene index
#'
#' @param data_dir Base directory for gene data (default: "data")
#' @return Data frame containing gene index information
#'
load_gene_index <- function(data_dir = "data") {
  index_file <- file.path(data_dir, "index", "gene_index.rds")
  
  if (!file.exists(index_file)) {
    stop("Gene index file not found: ", index_file)
  }
  
  readRDS(index_file)
}

#' Run the gene splitting process for both miscleavage types
#'
#' @param database_dir Directory containing the new miscleavage database files
#' @param output_dir Output directory for split files (default: "data")
#'
run_miscleavage_split <- function(database_dir = "database", output_dir = "data") {
  # Define file paths
  no_miss_file <- file.path(database_dir, "PEPTIDE_GENOMIC_LOCATION_no_misscleavage.RDS")
  upto_two_miss_file <- file.path(database_dir, "PEPTIDE_GENOMIC_LOCATION_added_misscleavage.RDS")
  
  # Check if files exist
  if (!file.exists(no_miss_file)) {
    stop("No miscleavage file not found: ", no_miss_file)
  }
  if (!file.exists(upto_two_miss_file)) {
    stop("Up to two miscleavage file not found: ", upto_two_miss_file)
  }
  
  # AS events file path
  as_events_file <- "AS_filtered.rds"
  
  message("Starting gene splitting process for miscleavage databases...")
  result <- split_rds_by_gene(
    no_miss_peptides_file = no_miss_file,
    upto_two_miss_peptides_file = upto_two_miss_file,
    as_events_file = as_events_file,
    output_dir = output_dir
  )
  
  message("Gene splitting completed successfully!")
  return(result)
}

#' Search for peptides across all genes and isoforms
#'
#' @param peptide_query Peptide sequence to search for
#' @param protease Protease to search in ("trp", "chymo", "aspn", "lysc", "lysn", "gluc")
#' @param miscleavage_type Miscleavage type ("no_miss_cleavage" or "upto_two_misscleavage")
#' @param exact_match Whether to search for exact matches only (default: FALSE)
#' @param search_dir Directory containing search databases (default: "data/search")
#' @return Data frame with search results
#'
search_peptides <- function(peptide_query, protease, miscleavage_type, exact_match = FALSE, search_dir = "data/search") {
  # Validate inputs
  valid_proteases <- c("trp", "chymo", "aspn", "lysc", "lysn", "gluc")
  valid_miscleavage <- c("no_miss_cleavage", "upto_two_misscleavage")
  
  if (!protease %in% valid_proteases) {
    stop("Invalid protease. Must be one of: ", paste(valid_proteases, collapse = ", "))
  }
  
  if (!miscleavage_type %in% valid_miscleavage) {
    stop("Invalid miscleavage_type. Must be one of: ", paste(valid_miscleavage, collapse = ", "))
  }
  
  # Convert miscleavage type to file naming convention
  miscleavage_suffix <- ifelse(miscleavage_type == "no_miss_cleavage", "no_miss", "upto2miss")
  
  # Load appropriate search database
  db_file <- file.path(search_dir, paste0("peptide_search_", miscleavage_suffix, "_", protease, ".rds"))
  
  if (!file.exists(db_file)) {
    stop("Search database not found: ", db_file, "\nPlease run create_peptide_search_databases_fast() first.")
  }
  
  search_db <- readRDS(db_file)
  
  # Search logic
  if (exact_match) {
    results <- search_db[search_db$peptide == peptide_query, ]
  } else {
    results <- search_db[grepl(peptide_query, search_db$peptide, ignore.case = TRUE), ]
  }
  
  # Add search metadata and gene symbols
  if (nrow(results) > 0) {
    results$protease_used <- protease
    results$miscleavage_type_used <- miscleavage_type
    results$search_query <- peptide_query
    
    # Add gene symbols from total database if available
    tryCatch({
      if (file.exists("total_database.rds")) {
        total_db <- readRDS("total_database.rds")
        gene_mapping <- unique(total_db[, c("geneID", "geneSymbol")])
        results <- merge(results, gene_mapping, by = "geneID", all.x = TRUE, suffixes = c("", "_db"))
        # Replace NA geneSymbol with the one from database
        if ("geneSymbol_db" %in% names(results)) {
          results$geneSymbol <- ifelse(is.na(results$geneSymbol), results$geneSymbol_db, results$geneSymbol)
          results$geneSymbol_db <- NULL
        }
      }
    }, error = function(e) {
      message("Warning: Could not add gene symbols: ", e$message)
    })
  }
  
  return(results)
}

#' Parse FASTA header into geneID and transcriptID
#'
#' @param header FASTA header string
#' @return List with transcriptID and geneID
#'
parse_header <- function(header) {
  parts <- strsplit(sub("^>", "", header), "\\|")[[1]]
  list(
    transcriptID = parts[2],
    geneID = parts[3]
  )
}

#' Create peptide search databases from GENCODE FASTA (FAST VERSION)
#'
#' @param fasta_file Path to GENCODE protein FASTA file
#' @param output_dir Directory where search databases will be stored (default: "data/search")
#' @return Information about created search databases
#'
create_peptide_search_databases_fast <- function(fasta_file, output_dir = "data/search") {
  # Load required libraries
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    stop("Package 'seqinr' is required but not installed.")
  }
  if (!requireNamespace("cleaver", quietly = TRUE)) {
    stop("Package 'cleaver' is required but not installed.")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  
  library(seqinr)
  library(cleaver)
  library(data.table)
  
  # Check if FASTA file exists
  if (!file.exists(fasta_file)) {
    stop("FASTA file not found: ", fasta_file)
  }
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  message("Loading FASTA file: ", fasta_file)
  fasta_seqs <- read.fasta(fasta_file, seqtype = "AA", as.string = TRUE)
  message("Loaded ", length(fasta_seqs), " protein sequences")
  
  # Enzyme mapping
  enzymes <- list(
    trp = "trypsin",
    chymo = "chymotrypsin-high", 
    aspn = "asp-n endopeptidase",
    lysc = "lysc",
    lysn = "lysn",
    gluc = "glutamyl endopeptidase"
  )
  
  # Miscleavage settings
  miscleavage_settings <- list(
    no_miss = 0,
    upto2miss = 2
  )
  
  results_summary <- list()
  
  # Process each combination of enzyme and miscleavage
  for (miss_name in names(miscleavage_settings)) {
    miss_value <- miscleavage_settings[[miss_name]]
    
    message("\n=== Processing ", miss_name, " (missedCleavages = ", miss_value, ") ===")
    
    for (enzyme_name in names(enzymes)) {
      enzyme_string <- enzymes[[enzyme_name]]
      
      message("Processing ", enzyme_name, " (", enzyme_string, ")...")
      
      # Create peptide data for this enzyme/miscleavage combination
      peptide_data <- rbindlist(lapply(names(fasta_seqs), function(hdr) {
        seq_string <- as.character(fasta_seqs[[hdr]])
        hdr_info <- parse_header(hdr)
        
        # Skip sequences that don't parse properly
        if (is.null(hdr_info$geneID) || is.null(hdr_info$transcriptID)) {
          return(NULL)
        }
        
        # Cleave protein sequence
        peptides_raw <- tryCatch({
          unlist(cleave(seq_string, 
                       enzym = enzyme_string,
                       missedCleavages = miss_value, 
                       unique = TRUE))
        }, error = function(e) {
          message("Warning: Failed to cleave sequence for ", hdr_info$transcriptID)
          return(character(0))
        })
        
        if (length(peptides_raw) == 0) return(NULL)
        
        # Filter peptides by length (6-60 amino acids)
        peptides <- peptides_raw[nchar(peptides_raw) > 6 & nchar(peptides_raw) < 61]
        
        if (length(peptides) == 0) return(NULL)
        
        # Return data.table
        data.table(
          geneID = hdr_info$geneID,
          geneSymbol = NA,  # We'll need to add this from a mapping file if needed
          txID = hdr_info$transcriptID,
          peptide = peptides
        )
      }))
      
      # Remove any NULL results
      if (is.null(peptide_data) || nrow(peptide_data) == 0) {
        message("  No peptides generated for ", enzyme_name, " with ", miss_name)
        next
      }
      
      # Save database
      output_file <- file.path(output_dir, paste0("peptide_search_", miss_name, "_", enzyme_name, ".rds"))
      saveRDS(peptide_data, output_file)
      
      # Store results summary
      results_summary[[paste0(miss_name, "_", enzyme_name)]] <- list(
        file = output_file,
        rows = nrow(peptide_data),
        size_mb = round(file.size(output_file) / 1024 / 1024, 2)
      )
      
      message("  Saved ", nrow(peptide_data), " peptide entries (", 
              round(file.size(output_file) / 1024 / 1024, 1), " MB)")
    }
  }
  
  # Print summary
  message("\n=== SUMMARY ===")
  total_files <- length(results_summary)
  total_rows <- sum(sapply(results_summary, function(x) x$rows))
  total_size <- sum(sapply(results_summary, function(x) x$size_mb))
  
  message("Created ", total_files, " search database files")
  message("Total peptide entries: ", format(total_rows, big.mark = ","))
  message("Total size: ", round(total_size, 1), " MB")
  message("Files saved in: ", output_dir)
  
  return(invisible(results_summary))
} 