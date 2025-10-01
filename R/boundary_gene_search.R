#===============================================================================
# BOUNDARY-BASED GENE SEARCH FOR NOVEL ISOFORMS
# Replacement for BLAST-based gene matching using genomic coordinates
#===============================================================================

# Load required functions
source("R/create_gene_boundary_database.R")
source("R/gene_boundary_matcher.R")

#' Run Boundary-Based Gene Search for Novel Isoform
#' 
#' Replaces BLAST search with boundary-based gene matching using genomic coordinates
#' 
#' @param novel_gtf_file Path to novel isoform GTF file
#' @param work_dir Working directory containing novel isoform results
#' @param min_overlap_bp Minimum overlap in base pairs (default: 50)
#' @param min_overlap_percent Minimum overlap percentage (default: 20)
#' @param max_genes Maximum number of genes to return (default: 10)
#' @return Data frame with matched genes in BLAST-compatible format
#' @export
run_boundary_gene_search <- function(novel_gtf_file = NULL,
                                    work_dir = NULL,
                                    min_overlap_bp = 50,
                                    min_overlap_percent = 10,
                                    max_genes = 10) {
  
  cat("=== Running Boundary-Based Gene Search ===\n")
  cat("Minimum overlap:", min_overlap_bp, "bp,", min_overlap_percent, "%\n")
  cat("Max genes:", max_genes, "\n")
  
  # Determine GTF file path
  if (is.null(novel_gtf_file) && !is.null(work_dir)) {
    novel_gtf_file <- file.path(work_dir, "results", "novel_final.gtf")
  }
  
  if (is.null(novel_gtf_file) || !file.exists(novel_gtf_file)) {
    stop("Novel GTF file not found: ", novel_gtf_file)
  }
  
  cat("Novel GTF file:", novel_gtf_file, "\n")
  
  # Load gene boundary database (create if doesn't exist)
  boundary_db_file <- "data/gene_boundaries.rds"
  
  if (!file.exists(boundary_db_file)) {
    cat("Creating gene boundary database...\n")
    boundary_db <- create_gene_boundary_database(
      gtf_file = "reference/gencode.v38.annotation.gtf",
      output_file = boundary_db_file,
      gene_types = c("protein_coding")
    )
  } else {
    cat("Loading existing gene boundary database...\n")
    boundary_db <- load_gene_boundary_database(boundary_db_file)
  }
  
  # Extract coordinates from novel GTF
  cat("Extracting coordinates from novel GTF...\n")
  novel_coordinates <- extract_coordinates_from_gtf(novel_gtf_file)
  
  if (nrow(novel_coordinates) == 0) {
    warning("No coordinates extracted from GTF file")
    return(create_empty_blast_format_result())
  }
  
  cat("Extracted", nrow(novel_coordinates), "coordinate ranges\n")
  
  # Search for overlapping genes
  cat("Searching for overlapping genes...\n")
  gene_matches <- search_genes_for_novel_isoform(
    boundary_db = boundary_db,
    novel_coordinates = novel_coordinates,
    search_method = "multiple",
    min_overlap_bp = min_overlap_bp,
    min_overlap_percent = min_overlap_percent,
    max_genes = max_genes
  )
  
  if (nrow(gene_matches) == 0) {
    cat("No overlapping genes found\n")
    return(create_empty_blast_format_result())
  }
  
  # Convert to BLAST-compatible format for downstream compatibility
  blast_format_results <- convert_to_blast_format(gene_matches)
  
  cat("Gene search completed successfully\n")
  cat("Found", nrow(blast_format_results), "candidate genes\n")
  
  return(blast_format_results)
}

#' Convert Boundary Search Results to BLAST Format
#' 
#' Converts boundary search results to BLAST-compatible format for downstream processing
#' 
#' @param gene_matches Gene matches from boundary search
#' @return Data frame in BLAST results format
convert_to_blast_format <- function(gene_matches) {
  
  if (nrow(gene_matches) == 0) {
    return(create_empty_blast_format_result())
  }
  
  # Create BLAST-compatible results format
  blast_results <- data.frame(
    qseqid = rep("novel_protein", nrow(gene_matches)),
    sseqid = paste0(gene_matches$gene_id, "|", gene_matches$gene_symbol),
    pident = pmin(gene_matches$max_confidence_score, 100),  # Use confidence as identity
    length = gene_matches$total_overlap_bp,
    mismatch = 0,
    gapopen = 0,
    qstart = 1,
    qend = 100,  # Placeholder values
    sstart = 1,
    send = 100,
    evalue = 1e-10,  # Dummy e-value (boundary search doesn't use e-values)
    bitscore = gene_matches$max_confidence_score * 10,  # Scale confidence to bitscore range
    stitle = paste0(gene_matches$gene_id, "|", gene_matches$gene_symbol, "|", 
                   gene_matches$gene_type, "|", gene_matches$chr, ":", 
                   gene_matches$gene_start, "-", gene_matches$gene_end),
    stringsAsFactors = FALSE
  )
  
  # Add parsed gene information columns (compatible with parse_blast_gene_info output)
  blast_results$gene_id <- gene_matches$gene_id
  blast_results$gene_symbol <- gene_matches$gene_symbol
  blast_results$transcript_id <- gene_matches$gene_id  # Use gene_id as transcript_id
  blast_results$protein_id <- gene_matches$gene_id
  blast_results$transcript_id_simple <- gsub("\\.[0-9]+$", "", gene_matches$gene_id)
  
  return(blast_results)
}

#' Create Empty BLAST Format Result
#' 
#' Creates an empty data frame in BLAST results format
#' 
#' @return Empty data frame with BLAST columns
create_empty_blast_format_result <- function() {
  
  empty_result <- data.frame(
    qseqid = character(0),
    sseqid = character(0),
    pident = numeric(0),
    length = numeric(0),
    mismatch = numeric(0),
    gapopen = numeric(0),
    qstart = numeric(0),
    qend = numeric(0),
    sstart = numeric(0),
    send = numeric(0),
    evalue = numeric(0),
    bitscore = numeric(0),
    stitle = character(0),
    gene_id = character(0),
    gene_symbol = character(0),
    transcript_id = character(0),
    protein_id = character(0),
    transcript_id_simple = character(0),
    stringsAsFactors = FALSE
  )
  
  return(empty_result)
}

#' Extract Novel Isoform Genomic Summary
#' 
#' Extracts genomic summary information from novel isoform GTF
#' 
#' @param novel_gtf_file Path to novel GTF file
#' @return List with genomic summary information
#' @export
extract_novel_genomic_summary <- function(novel_gtf_file) {
  
  if (!file.exists(novel_gtf_file)) {
    return(list(
      total_exons = 0,
      chromosomes = character(0),
      genomic_span = 0,
      coordinates = data.frame()
    ))
  }
  
  coordinates <- extract_coordinates_from_gtf(novel_gtf_file)
  
  if (nrow(coordinates) == 0) {
    return(list(
      total_exons = 0,
      chromosomes = character(0),
      genomic_span = 0,
      coordinates = data.frame()
    ))
  }
  
  # Calculate genomic span
  genomic_span <- max(coordinates$end) - min(coordinates$start) + 1
  
  # Get unique chromosomes
  chromosomes <- unique(coordinates$chr)
  
  summary_info <- list(
    total_exons = nrow(coordinates),
    chromosomes = chromosomes,
    genomic_span = genomic_span,
    coordinates = coordinates,
    genomic_range = paste0(chromosomes[1], ":", min(coordinates$start), "-", max(coordinates$end))
  )
  
  return(summary_info)
}

#' Check Gene RDS Availability for Boundary Results
#' 
#' Checks if RDS files are available for boundary search results
#' 
#' @param boundary_results Boundary search results in BLAST format
#' @param rds_dir Directory containing gene RDS files
#' @return Data frame with gene availability status
#' @export
check_boundary_gene_rds_availability <- function(boundary_results, rds_dir = "data/genes") {
  
  if (nrow(boundary_results) == 0) {
    return(data.frame())
  }
  
  gene_ids <- boundary_results$gene_id
  
  availability_df <- data.frame(
    gene_id = gene_ids,
    gene_symbol = boundary_results$gene_symbol,
    rds_file = paste0(gene_ids, ".rds"),
    rds_path = file.path(rds_dir, paste0(gene_ids, ".rds")),
    available = FALSE,
    confidence_score = boundary_results$pident,  # Use pident as confidence score
    stringsAsFactors = FALSE
  )
  
  # Check file existence
  for (i in 1:nrow(availability_df)) {
    availability_df$available[i] <- file.exists(availability_df$rds_path[i])
  }
  
  # Sort by availability and confidence score
  availability_df <- availability_df[order(-availability_df$available, -availability_df$confidence_score), ]
  
  available_count <- sum(availability_df$available)
  cat("RDS files available for", available_count, "out of", nrow(availability_df), "genes\n")
  
  return(availability_df)
}

#' Test Boundary Gene Search
#' 
#' Tests the boundary-based gene search functionality
#' 
#' @export
test_boundary_gene_search <- function() {
  
  cat("=== Testing Boundary Gene Search ===\n")
  
  # Check if we have a test GTF file from novel isoform results
  test_dirs <- list.dirs("novel_isoform_results", recursive = FALSE)
  
  if (length(test_dirs) == 0) {
    cat("No novel isoform results found - skipping test\n")
    return(invisible(NULL))
  }
  
  # Use the most recent result
  latest_dir <- test_dirs[length(test_dirs)]
  test_gtf <- file.path(latest_dir, "results", "novel_final.gtf")
  
  if (!file.exists(test_gtf)) {
    cat("No test GTF file found - skipping test\n")
    return(invisible(NULL))
  }
  
  cat("Testing with GTF file:", test_gtf, "\n")
  
  # Test coordinate extraction
  cat("Testing coordinate extraction...\n")
  genomic_summary <- extract_novel_genomic_summary(test_gtf)
  
  if (genomic_summary$total_exons > 0) {
    cat("✅ Coordinate extraction successful\n")
    cat("   Exons found:", genomic_summary$total_exons, "\n")
    cat("   Chromosomes:", paste(genomic_summary$chromosomes, collapse = ", "), "\n")
    cat("   Genomic range:", genomic_summary$genomic_range, "\n")
  } else {
    cat("❌ Coordinate extraction failed\n")
    return(invisible(NULL))
  }
  
  # Test boundary gene search
  cat("Testing boundary gene search...\n")
  
  tryCatch({
    search_results <- run_boundary_gene_search(
      novel_gtf_file = test_gtf,
      min_overlap_bp = 500,
      min_overlap_percent = 10,
      max_genes = 5
    )
    
    if (nrow(search_results) > 0) {
      cat("✅ Boundary gene search successful\n")
      cat("   Genes found:", nrow(search_results), "\n")
      cat("   Top candidate:", search_results$gene_symbol[1], 
          "(", search_results$gene_id[1], ")\n")
      cat("   Confidence:", search_results$pident[1], "\n")
      
      # Test RDS availability
      availability <- check_boundary_gene_rds_availability(search_results)
      
      if (nrow(availability) > 0) {
        cat("✅ RDS availability check completed\n")
        available_genes <- sum(availability$available)
        cat("   Available genes:", available_genes, "out of", nrow(availability), "\n")
      }
      
    } else {
      cat("⚠️  No genes found (this may be normal for some novel sequences)\n")
    }
    
    cat("Boundary gene search test completed successfully\n")
    return(search_results)
    
  }, error = function(e) {
    cat("❌ Boundary gene search test failed:", e$message, "\n")
    return(NULL)
  })
}