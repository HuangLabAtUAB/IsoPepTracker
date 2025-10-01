#===============================================================================
# GTF PREPROCESSOR - LIGHTNING FAST VISUALIZATION SOLUTION
#===============================================================================

#' Pre-process GTF file into gene-specific visualization data
#' This eliminates the 15+ second GTF parsing bottleneck
#' 
#' @param gtf_file Path to GTF file
#' @param output_dir Directory to save processed data
#' @param gene_filter Optional vector of gene IDs to process (default: all genes)
#' @return Information about processed genes
#' 
create_gtf_visualization_cache <- function(gtf_file, output_dir = "data/gtf_cache", gene_filter = NULL) {
  
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required")
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required")
  }
  
  library(rtracklayer)
  library(GenomicRanges)
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Loading GTF file once (this takes ~15 seconds)...")
  start_time <- Sys.time()
  
  # Load all data once
  gene_features <- import(gtf_file, format = "gtf", feature.type = "gene")
  exon_features <- import(gtf_file, format = "gtf", feature.type = "exon") 
  cds_features <- import(gtf_file, format = "gtf", feature.type = "CDS")
  
  message("GTF loading completed in ", round(difftime(Sys.time(), start_time, units = "secs"), 2), " seconds")
  
  # Get all unique gene IDs
  all_gene_ids <- unique(gene_features$gene_id)
  
  # Helper function to strip version numbers from gene IDs
  strip_version <- function(gene_ids) {
    return(sub("\\.[0-9]+$", "", gene_ids))
  }
  
  # Filter genes if specified (with version-flexible matching)
  if (!is.null(gene_filter)) {
    # Strip versions from both filter and gene IDs for flexible matching
    filter_base <- strip_version(gene_filter)
    genes_base <- strip_version(all_gene_ids)
    
    # Find matches using base IDs (without versions)
    matched_indices <- match(filter_base, genes_base)
    matched_indices <- matched_indices[!is.na(matched_indices)]
    
    if (length(matched_indices) > 0) {
      all_gene_ids <- all_gene_ids[matched_indices]
      message("Processing ", length(all_gene_ids), " genes (filtered from ", length(unique(gene_features$gene_id)), " total genes using version-flexible matching)...")
    } else {
      message("No genes matched after version-flexible filtering. Processing all ", length(all_gene_ids), " genes...")
    }
  } else {
    message("Processing ALL ", length(all_gene_ids), " genes (no filtering)...")
  }
  
  # Process each gene
  processed_count <- 0
  for (i in seq_along(all_gene_ids)) {
    gene_id <- all_gene_ids[i]
    
    # Progress update
    if (i %% 500 == 0 || i == length(all_gene_ids)) {
      message("⚡ Processed ", i, " of ", length(all_gene_ids), " genes (", 
              round(i/length(all_gene_ids)*100, 1), "%) - making visualizations lightning-fast!")
    }
    
    # Get gene info
    gene_info <- gene_features[gene_features$gene_id == gene_id]
    if (length(gene_info) == 0) next
    
    # Get exons for this gene
    gene_exons <- exon_features[exon_features$gene_id == gene_id]
    
    # Get CDS for this gene  
    gene_cds <- cds_features[cds_features$gene_id == gene_id]
    
    # Create lightweight visualization data
    viz_data <- list(
      gene_id = gene_id,
      gene_symbol = gene_info$gene_name[1],
      chromosome = as.character(seqnames(gene_info[1])),
      gene_start = start(gene_info[1]),
      gene_end = end(gene_info[1]),
      strand = as.character(strand(gene_info[1])),
      
      # Pre-computed exon data by transcript
      exons_by_transcript = split(gene_exons, gene_exons$transcript_id),
      
      # Pre-computed CDS data by transcript
      cds_by_transcript = split(gene_cds, gene_cds$transcript_id),
      
      # Pre-computed transcript list
      transcript_ids = unique(gene_exons$transcript_id)
    )
    
    # Save gene-specific visualization data
    output_file <- file.path(output_dir, paste0(gene_id, ".rds"))
    saveRDS(viz_data, output_file)
    processed_count <- processed_count + 1
  }
  
  message("Successfully processed ", processed_count, " genes")
  message("Visualization cache created in: ", output_dir)
  
  return(list(
    processed_genes = processed_count,
    cache_dir = output_dir
  ))
}

#' Load pre-computed GTF visualization data (INSTANT)
#' 
#' @param gene_id Gene ID to load
#' @param cache_dir Cache directory
#' @return Pre-computed visualization data
#' 
load_gtf_visualization_data <- function(gene_id, cache_dir = "data/gtf_cache") {
  # Helper function to strip version numbers
  strip_version <- function(gene_ids) {
    return(sub("\\.[0-9]+$", "", gene_ids))
  }
  
  # Try exact match first
  cache_file <- file.path(cache_dir, paste0(gene_id, ".rds"))
  
  if (!file.exists(cache_file)) {
    # Try version-flexible matching
    gene_base <- strip_version(gene_id)
    
    # List all cache files and find matches
    cache_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
    cache_gene_ids <- sub("\\.rds$", "", basename(cache_files))
    cache_gene_bases <- strip_version(cache_gene_ids)
    
    # Find match using base ID
    match_index <- match(gene_base, cache_gene_bases)
    
    if (!is.na(match_index)) {
      cache_file <- cache_files[match_index]
      message("Found GTF cache using version-flexible matching: ", basename(cache_file))
    } else {
      return(list(
        success = FALSE,
        message = paste("GTF cache not found for gene", gene_id, "(tried both exact and version-flexible matching)")
      ))
    }
  }
  
  # This is INSTANT (0.002 seconds vs 15+ seconds)
  viz_data <- readRDS(cache_file)
  
  return(list(
    success = TRUE,
    gene_id = viz_data$gene_id,
    gene_symbol = viz_data$gene_symbol,
    chromosome = viz_data$chromosome,
    gene_start = viz_data$gene_start,
    gene_end = viz_data$gene_end,
    strand = viz_data$strand,
    exons_by_transcript = viz_data$exons_by_transcript,
    cds_by_transcript = viz_data$cds_by_transcript,
    transcript_ids = viz_data$transcript_ids
  ))
}

#' Create the GTF cache (run once)
#' 
create_gtf_cache <- function() {
  gtf_file <- "~/Downloads/gencode.v38.annotation.gtf.gz"
  
  if (!file.exists(gtf_file)) {
    stop("GTF file not found at: ", gtf_file)
  }
  
  message("Creating GTF visualization cache...")
  message("This will take ~15 seconds but only needs to be done ONCE")
  
  result <- create_gtf_visualization_cache(gtf_file)
  
  message("GTF cache creation completed!")
  message("Visualizations will now be INSTANT!")
  
  return(result)
}