#===============================================================================
# GENE BOUNDARY DATABASE GENERATOR
# Parse GENCODE GTF to create gene boundary database for fast interval searches
#===============================================================================

#' Create Gene Boundary Database from GENCODE GTF
#' 
#' Parses GENCODE GTF file to extract gene boundaries and create searchable database
#' 
#' @param gtf_file Path to GENCODE GTF annotation file
#' @param output_file Path to save the gene boundary database RDS file
#' @param gene_types Vector of gene types to include (default: protein_coding genes)
#' @return Gene boundary database as data.frame
#' @export
create_gene_boundary_database <- function(gtf_file = "reference/gencode.v38.annotation.gtf",
                                         output_file = "data/gene_boundaries.rds",
                                         gene_types = c("protein_coding")) {
  
  cat("Creating gene boundary database from GTF file...\n")
  cat("GTF file:", gtf_file, "\n")
  cat("Output file:", output_file, "\n")
  
  if (!file.exists(gtf_file)) {
    stop("GTF file not found: ", gtf_file)
  }
  
  # Load required libraries
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required but not installed")
  }
  library(data.table)
  
  cat("Reading GTF file (this may take a moment)...\n")
  
  # Read GTF file with data.table for speed
  gtf_data <- fread(gtf_file, sep = "\t", header = FALSE, 
                    col.names = c("seqname", "source", "feature", "start", "end", 
                                "score", "strand", "frame", "attributes"),
                    quote = "", stringsAsFactors = FALSE)
  
  cat("GTF file loaded:", nrow(gtf_data), "total features\n")
  
  # Filter for gene features only
  gene_features <- gtf_data[feature == "gene"]
  cat("Gene features found:", nrow(gene_features), "\n")
  
  if (nrow(gene_features) == 0) {
    stop("No gene features found in GTF file")
  }
  
  # Parse gene attributes to extract gene_id, gene_name, and gene_type
  cat("Parsing gene attributes...\n")
  
  gene_boundaries <- data.frame(
    chr = gene_features$seqname,
    start = gene_features$start,
    end = gene_features$end,
    strand = gene_features$strand,
    gene_id = character(nrow(gene_features)),
    gene_symbol = character(nrow(gene_features)),
    gene_type = character(nrow(gene_features)),
    width = gene_features$end - gene_features$start + 1,
    stringsAsFactors = FALSE
  )
  
  # Parse attributes string to extract gene information
  for (i in 1:nrow(gene_features)) {
    attributes <- gene_features$attributes[i]
    
    # Extract gene_id
    gene_id_match <- regmatches(attributes, regexpr('gene_id "([^"]+)"', attributes))
    if (length(gene_id_match) > 0) {
      gene_boundaries$gene_id[i] <- gsub('gene_id "([^"]+)"', '\\1', gene_id_match)
    }
    
    # Extract gene_name (gene symbol)
    gene_name_match <- regmatches(attributes, regexpr('gene_name "([^"]+)"', attributes))
    if (length(gene_name_match) > 0) {
      gene_boundaries$gene_symbol[i] <- gsub('gene_name "([^"]+)"', '\\1', gene_name_match)
    } else {
      # Fallback to gene_id if no gene_name
      gene_boundaries$gene_symbol[i] <- gene_boundaries$gene_id[i]
    }
    
    # Extract gene_type
    gene_type_match <- regmatches(attributes, regexpr('gene_type "([^"]+)"', attributes))
    if (length(gene_type_match) > 0) {
      gene_boundaries$gene_type[i] <- gsub('gene_type "([^"]+)"', '\\1', gene_type_match)
    }
    
    # Progress indicator
    if (i %% 5000 == 0) {
      cat("Processed", i, "of", nrow(gene_features), "genes\n")
    }
  }
  
  cat("Gene attribute parsing completed\n")
  
  # Filter by gene types if specified
  if (!is.null(gene_types) && length(gene_types) > 0) {
    before_filter <- nrow(gene_boundaries)
    gene_boundaries <- gene_boundaries[gene_boundaries$gene_type %in% gene_types, ]
    cat("Filtered by gene types", paste(gene_types, collapse = ", "), "\n")
    cat("Genes before filter:", before_filter, "After filter:", nrow(gene_boundaries), "\n")
  }
  
  # Remove genes with missing essential information
  complete_genes <- complete.cases(gene_boundaries[, c("chr", "start", "end", "gene_id")])
  gene_boundaries <- gene_boundaries[complete_genes, ]
  
  # Sort by chromosome and position for efficient searching
  gene_boundaries <- gene_boundaries[order(gene_boundaries$chr, gene_boundaries$start), ]
  
  # Add unique index for fast lookup
  gene_boundaries$boundary_id <- 1:nrow(gene_boundaries)
  
  # Create chromosome index for fast filtering
  chr_index <- split(1:nrow(gene_boundaries), gene_boundaries$chr)
  
  # Summary statistics
  cat("\n=== Gene Boundary Database Summary ===\n")
  cat("Total genes in database:", nrow(gene_boundaries), "\n")
  cat("Chromosomes covered:", length(unique(gene_boundaries$chr)), "\n")
  cat("Gene types included:", paste(unique(gene_boundaries$gene_type), collapse = ", "), "\n")
  cat("Gene width range:", min(gene_boundaries$width), "-", max(gene_boundaries$width), "bp\n")
  
  # Show chromosome distribution
  chr_counts <- table(gene_boundaries$chr)
  cat("Top chromosomes by gene count:\n")
  print(head(sort(chr_counts, decreasing = TRUE), 10))
  
  # Create final database structure
  boundary_database <- list(
    genes = gene_boundaries,
    chr_index = chr_index,
    metadata = list(
      created_date = Sys.time(),
      gtf_file = gtf_file,
      gene_types_included = gene_types,
      total_genes = nrow(gene_boundaries),
      version = "1.0"
    )
  )
  
  # Save database
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(boundary_database, output_file)
  
  cat("Gene boundary database saved to:", output_file, "\n")
  cat("Database creation completed successfully!\n")
  
  return(boundary_database)
}

#' Load Gene Boundary Database
#' 
#' Loads the pre-computed gene boundary database
#' 
#' @param database_file Path to gene boundary database RDS file
#' @return Gene boundary database
#' @export
load_gene_boundary_database <- function(database_file = "data/gene_boundaries.rds") {
  
  if (!file.exists(database_file)) {
    stop("Gene boundary database not found: ", database_file, 
         "\nPlease run create_gene_boundary_database() first.")
  }
  
  cat("Loading gene boundary database from:", database_file, "\n")
  boundary_db <- readRDS(database_file)
  
  cat("Database loaded successfully\n")
  cat("Total genes:", boundary_db$metadata$total_genes, "\n")
  cat("Created:", boundary_db$metadata$created_date, "\n")
  
  return(boundary_db)
}

#' Get Gene Boundaries by Chromosome
#' 
#' Efficiently retrieves gene boundaries for a specific chromosome
#' 
#' @param boundary_db Gene boundary database
#' @param chromosome Chromosome name (e.g., "chr1", "chrX")
#' @return Subset of gene boundaries for the specified chromosome
#' @export
get_genes_by_chromosome <- function(boundary_db, chromosome) {
  
  if (!chromosome %in% names(boundary_db$chr_index)) {
    return(data.frame())
  }
  
  gene_indices <- boundary_db$chr_index[[chromosome]]
  return(boundary_db$genes[gene_indices, ])
}

#' Test Gene Boundary Database Creation
#' 
#' Creates a test database and validates functionality
#' 
#' @export
test_gene_boundary_database <- function() {
  
  cat("=== Testing Gene Boundary Database Creation ===\n")
  
  # Test database creation
  gtf_file <- "reference/gencode.v38.annotation.gtf"
  test_output <- "data/test_gene_boundaries.rds"
  
  if (!file.exists(gtf_file)) {
    cat("GTF file not found - skipping test\n")
    return(invisible(NULL))
  }
  
  # Create database with all gene types for testing
  test_db <- create_gene_boundary_database(
    gtf_file = gtf_file,
    output_file = test_output,
    gene_types = NULL  # Include all gene types for comprehensive testing
  )
  
  # Test loading
  loaded_db <- load_gene_boundary_database(test_output)
  
  # Test chromosome filtering
  chr1_genes <- get_genes_by_chromosome(loaded_db, "chr1")
  cat("chr1 genes found:", nrow(chr1_genes), "\n")
  
  # Test specific gene lookup
  if (nrow(chr1_genes) > 0) {
    sample_gene <- chr1_genes[1, ]
    cat("Sample gene from chr1:\n")
    cat("  Gene ID:", sample_gene$gene_id, "\n")
    cat("  Symbol:", sample_gene$gene_symbol, "\n")
    cat("  Position:", sample_gene$start, "-", sample_gene$end, "\n")
    cat("  Type:", sample_gene$gene_type, "\n")
  }
  
  cat("✅ Gene boundary database test completed successfully!\n")
  
  # Cleanup test file
  if (file.exists(test_output)) {
    unlink(test_output)
  }
  
  return(loaded_db)
}