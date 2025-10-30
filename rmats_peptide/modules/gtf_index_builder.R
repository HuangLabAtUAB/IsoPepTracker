# GTF Index Builder - Creates compact searchable index from large GTF files
# This avoids loading the entire GTF file during analysis

library(data.table)
library(readr)

# Build compact GTF index from full GTF file
build_gtf_index_from_file <- function(gtf_file_path, index_output_path = NULL) {
  
  if (is.null(index_output_path)) {
    index_output_path <- paste0(dirname(gtf_file_path), "/", 
                               gsub("\\.(gtf|gff).*$", "_index.RDS", basename(gtf_file_path)))
  }
  
  cat("Building GTF index from:", gtf_file_path, "\n")
  cat("Output index will be saved to:", index_output_path, "\n")
  
  # Check if GTF file exists
  if (!file.exists(gtf_file_path)) {
    stop("GTF file not found: ", gtf_file_path)
  }
  
  # Read GTF file efficiently
  cat("Reading GTF file...\n")
  
  # GTF columns: seqname, source, feature, start, end, score, strand, frame, attributes
  gtf_data <- fread(gtf_file_path, sep = "\t", header = FALSE, comment.char = "#",
                   col.names = c("seqname", "source", "feature", "start", "end", 
                                "score", "strand", "frame", "attributes"))
  
  cat("Total GTF records:", nrow(gtf_data), "\n")
  
  # Filter for exons only (we only need exons for splicing analysis)
  exon_data <- gtf_data[feature == "exon"]
  cat("Exon records:", nrow(exon_data), "\n")
  
  # Parse attributes to extract gene_id, transcript_id, exon_id
  cat("Parsing GTF attributes...\n")
  exon_data[, gene_id := extract_gtf_attribute(attributes, "gene_id")]
  exon_data[, transcript_id := extract_gtf_attribute(attributes, "transcript_id")]
  exon_data[, exon_id := extract_gtf_attribute(attributes, "exon_id")]
  
  # Extract phase information (if available)
  exon_data[, phase_start := ifelse(frame == ".", NA, as.integer(frame))]
  exon_data[, phase_end := calculate_end_phase(start, end, phase_start)]
  
  # Mark which exons are in CDS (have phase information)
  exon_data[, in_cds := !is.na(phase_start)]
  
  # Create compact index structure
  cat("Building compact index structure...\n")
  
  # Index by gene_id for fast gene-specific searches
  by_gene_index <- split(exon_data, exon_data$gene_id)
  
  # Create coordinate lookup table (more compact than full data)
  coord_index <- exon_data[, .(
    seqname, start, end, strand, 
    gene_id, transcript_id, exon_id,
    phase_start, phase_end, in_cds
  )]
  
  # Create final index structure
  gtf_index <- list(
    # Gene-specific index (most important for our use case)
    by_gene = by_gene_index,
    
    # Coordinate lookup table
    coord_table = coord_index,
    
    # Summary statistics
    stats = list(
      total_exons = nrow(exon_data),
      genes_indexed = length(unique(exon_data$gene_id)),
      cds_exons = sum(exon_data$in_cds, na.rm = TRUE),
      chromosomes = unique(exon_data$seqname)
    ),
    
    # Metadata
    source_file = gtf_file_path,
    index_built_at = Sys.time(),
    index_version = "1.0"
  )
  
  # Save index to disk
  cat("Saving index to:", index_output_path, "\n")
  saveRDS(gtf_index, index_output_path)
  
  # Print summary
  cat("\n=== GTF INDEX BUILT SUCCESSFULLY ===\n")
  cat("Source GTF:", gtf_file_path, "\n")
  cat("Index file:", index_output_path, "\n")
  cat("Total exons indexed:", gtf_index$stats$total_exons, "\n")
  cat("Genes indexed:", gtf_index$stats$genes_indexed, "\n")
  cat("CDS exons:", gtf_index$stats$cds_exons, "\n")
  cat("Chromosomes:", length(gtf_index$stats$chromosomes), "\n")
  cat("Index file size:", round(file.size(index_output_path) / 1024^2, 2), "MB\n")
  
  return(index_output_path)
}

# Load pre-built GTF index from disk
load_gtf_index <- function(index_file_path) {
  
  if (!file.exists(index_file_path)) {
    stop("GTF index file not found: ", index_file_path, 
         "\nPlease build the index first using build_gtf_index_from_file()")
  }
  
  cat("Loading GTF index from:", index_file_path, "\n")
  gtf_index <- readRDS(index_file_path)
  
  cat("Index loaded successfully:\n")
  cat("- Source GTF:", gtf_index$source_file, "\n")
  cat("- Built at:", gtf_index$index_built_at, "\n")
  cat("- Genes indexed:", gtf_index$stats$genes_indexed, "\n")
  cat("- Total exons:", gtf_index$stats$total_exons, "\n")
  
  return(gtf_index)
}

# Extract specific attribute from GTF attributes column
extract_gtf_attribute <- function(attributes, attr_name) {
  pattern <- paste0(attr_name, '\\s+"([^"]+)"')
  matches <- regmatches(attributes, regexpr(pattern, attributes, perl = TRUE))
  result <- gsub(paste0(attr_name, '\\s+"([^"]+)"'), "\\1", matches, perl = TRUE)
  result[result == ""] <- NA
  return(result)
}

# Calculate end phase from start phase and exon length
calculate_end_phase <- function(start, end, start_phase) {
  if (any(is.na(start_phase))) {
    return(rep(NA, length(start)))
  }
  
  exon_length <- end - start + 1
  end_phase <- (start_phase + exon_length) %% 3
  return(end_phase)
}

# Search the pre-built index (much faster than searching full GTF)
search_gtf_index <- function(gtf_index, chromosome, start, end, strand, gene_id = NULL) {
  
  # Use gene-specific index if gene_id provided (fastest)
  if (!is.null(gene_id) && gene_id %in% names(gtf_index$by_gene)) {
    search_data <- gtf_index$by_gene[[gene_id]]
  } else {
    # Fall back to coordinate table search
    search_data <- gtf_index$coord_table
  }
  
  # Exact coordinate matching
  if (is.data.table(search_data)) {
    matches <- search_data[
      seqname == chromosome & 
      start == !!start & 
      end == !!end & 
      strand == !!strand
    ]
  } else {
    # Handle data.frame
    matches <- search_data[
      search_data$seqname == chromosome &
      search_data$start == start &
      search_data$end == end &
      search_data$strand == strand,
    ]
  }
  
  if (nrow(matches) == 0) {
    return(list(
      found = FALSE,
      matches = data.frame(),
      reason = paste("No exact match for", chromosome, start, end, strand)
    ))
  }
  
  return(list(
    found = TRUE,
    matches = as.data.frame(matches),
    count = nrow(matches)
  ))
}

# Create example index builder script
create_index_builder_script <- function(gtf_file_path, output_script = "build_gtf_index.R") {
  
  script_content <- paste0('
# GTF Index Builder Script
# Run this once to create a compact index from your GTF file

source("modules/gtf_index_builder.R")

# Specify your GTF file path
gtf_file <- "', gtf_file_path, '"

# Build the index (this may take several minutes for large GTF files)
index_path <- build_gtf_index_from_file(gtf_file)

cat("GTF index created successfully!\\n")
cat("Index saved to:", index_path, "\\n")
cat("\\nTo use the index in your analysis:")
cat("gtf_index <- load_gtf_index(\\"", index_path, "\\")\\n")
')
  
  writeLines(script_content, output_script)
  cat("Index builder script created:", output_script, "\n")
  cat("Run this script once to build your GTF index\n")
  
  return(output_script)
}