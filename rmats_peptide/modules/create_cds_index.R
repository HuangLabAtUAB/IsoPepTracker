# Create CDS Index RDS File
# Single row per CDS exon with all essential information for fast searching

library(data.table)
library(readr)

# Create CDS index from GTF file - one row per CDS exon
create_cds_index_rds <- function(gtf_file_path, output_rds_path = NULL) {
  
  if (is.null(output_rds_path)) {
    output_rds_path <- paste0(dirname(gtf_file_path), "/cds_index.RDS")
  }
  
  cat("Creating CDS index from GTF:", gtf_file_path, "\n")
  cat("Output RDS:", output_rds_path, "\n")
  
  if (!file.exists(gtf_file_path)) {
    stop("GTF file not found: ", gtf_file_path)
  }
  
  # Read GTF file
  cat("Reading GTF file...\n")
  gtf_data <- fread(gtf_file_path, sep = "\t", header = FALSE, skip = "#",
                   col.names = c("seqname", "source", "feature", "start", "end", 
                                "score", "strand", "frame", "attributes"))
  
  cat("Total GTF records:", nrow(gtf_data), "\n")
  
  # Filter for CDS features only (these have phase information)
  cds_data <- gtf_data[feature == "CDS"]
  cat("CDS records found:", nrow(cds_data), "\n")
  
  if (nrow(cds_data) == 0) {
    stop("No CDS records found in GTF file. Cannot create CDS index.")
  }
  
  # Parse attributes to extract essential information
  cat("Parsing GTF attributes...\n")
  cds_data[, gene_id := extract_gtf_attribute(attributes, "gene_id")]
  cds_data[, transcript_id := extract_gtf_attribute(attributes, "transcript_id")]
  cds_data[, gene_name := extract_gtf_attribute(attributes, "gene_name")]
  cds_data[, exon_number := as.integer(extract_gtf_attribute_unquoted(attributes, "exon_number"))]
  
  # Use frame as phase (0, 1, or 2)
  cds_data[, phase := ifelse(frame == ".", NA, as.integer(frame))]
  
  # Create final CDS index - ONE ROW PER CDS EXON
  cds_index <- cds_data[, .(
    # Coordinates for exact matching
    chromosome = seqname,
    start = start,
    end = end,
    strand = strand,
    
    # Gene/transcript information
    gene_id = gene_id,
    gene_name = gene_name,
    transcript_id = transcript_id,
    exon_number = exon_number,
    
    # Phase information (critical for translation)
    phase = phase,
    
    # Additional useful info
    cds_length = end - start + 1,
    
    # Source info
    source = source
  )]
  
  # Remove rows with missing essential information
  cds_index <- cds_index[!is.na(gene_id) & !is.na(transcript_id) & !is.na(phase)]
  
  cat("Final CDS index rows:", nrow(cds_index), "\n")
  cat("Genes in index:", length(unique(cds_index$gene_id)), "\n")
  cat("Transcripts in index:", length(unique(cds_index$transcript_id)), "\n")
  
  # Add metadata
  attr(cds_index, "source_gtf") <- gtf_file_path
  attr(cds_index, "created_at") <- Sys.time()
  attr(cds_index, "version") <- "1.0"
  attr(cds_index, "description") <- "CDS exon index for rMATS functional analysis"
  
  # Save as RDS
  cat("Saving CDS index to:", output_rds_path, "\n")
  saveRDS(cds_index, output_rds_path)
  
  # Print summary
  cat("\n=== CDS INDEX CREATED SUCCESSFULLY ===\n")
  cat("Output file:", output_rds_path, "\n")
  cat("Index size:", round(file.size(output_rds_path) / 1024^2, 2), "MB\n")
  cat("CDS exons indexed:", nrow(cds_index), "\n")
  cat("Unique genes:", length(unique(cds_index$gene_id)), "\n")
  cat("Unique transcripts:", length(unique(cds_index$transcript_id)), "\n")
  
  return(output_rds_path)
}

# Extract attribute value from GTF attributes string (quoted)
extract_gtf_attribute <- function(attributes, attr_name) {
  pattern <- paste0(attr_name, '\\s+"([^"]+)"')
  matches <- regmatches(attributes, regexpr(pattern, attributes, perl = TRUE))
  result <- gsub(paste0(attr_name, '\\s+"([^"]+)"'), "\\1", matches, perl = TRUE)
  result[result == "" | length(result) == 0] <- NA
  return(result)
}

# Extract attribute value from GTF attributes string (unquoted numbers)
extract_gtf_attribute_unquoted <- function(attributes, attr_name) {
  pattern <- paste0(attr_name, '\\s+(\\d+)')
  matches <- regmatches(attributes, regexpr(pattern, attributes, perl = TRUE))
  result <- gsub(paste0(attr_name, '\\s+(\\d+)'), "\\1", matches, perl = TRUE)
  result[result == "" | length(result) == 0] <- NA
  return(result)
}


# Load the CDS index RDS file
load_cds_index <- function(rds_path) {
  
  if (!file.exists(rds_path)) {
    stop("CDS index RDS file not found: ", rds_path)
  }
  
  cat("Loading CDS index from:", rds_path, "\n")
  cds_index <- readRDS(rds_path)
  
  cat("CDS index loaded:\n")
  cat("- CDS exons:", nrow(cds_index), "\n")
  cat("- Genes:", length(unique(cds_index$gene_id)), "\n")
  cat("- Created:", attr(cds_index, "created_at"), "\n")
  
  return(cds_index)
}

# Search CDS index for exact coordinate matches
search_cds_index <- function(cds_index, search_chromosome, search_start, search_end, search_strand, search_gene_id = NULL) {
  
  # Filter by gene if specified (much faster)
  if (!is.null(search_gene_id)) {
    search_data <- cds_index[gene_id == search_gene_id]  # Use data.table syntax
  } else {
    search_data <- cds_index
  }
  
  # Exact coordinate matching
  matches <- search_data[
    chromosome == search_chromosome & 
    start == search_start & 
    end == search_end & 
    strand == search_strand
  ]
  
  if (nrow(matches) == 0) {
    return(list(
      found = FALSE,
      matches = data.table(),
      reason = paste("No exact CDS match:", search_chromosome, search_start, search_end, search_strand)
    ))
  }
  
  return(list(
    found = TRUE,
    matches = matches,
    count = nrow(matches)
  ))
}

# Create example script to build CDS index
create_cds_index_script <- function(gtf_path = "gencode.v38.annotation.gtf") {
  
  script_content <- paste0('
# Build CDS Index Script
# Run this once to create CDS index from your GTF file

source("modules/create_cds_index.R")

# Specify your GTF file
gtf_file <- "', gtf_path, '"

# Create CDS index (one row per CDS exon)
index_path <- create_cds_index_rds(gtf_file)

cat("\\nCDS index created!\\n")
cat("To use in analysis:\\n")
cat("cds_index <- load_cds_index(\\"", basename(index_path), "\\")\\n")
')
  
  writeLines(script_content, "build_cds_index.R")
  cat("CDS index builder script created: build_cds_index.R\n")
  cat("Run this script once to build your CDS index\n")
}