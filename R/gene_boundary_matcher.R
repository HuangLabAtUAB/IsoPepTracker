#===============================================================================
# GENE BOUNDARY MATCHER
# Boundary-based gene selection using genomic coordinate interval intersection
#===============================================================================

#' Find Genes by Genomic Boundary Overlap
#' 
#' Identifies genes that overlap with given genomic coordinates using interval intersection
#' 
#' @param boundary_db Gene boundary database (from load_gene_boundary_database)
#' @param chr Chromosome name (e.g., "chr1", "chrX") 
#' @param start Start coordinate
#' @param end End coordinate
#' @param min_overlap_bp Minimum overlap in base pairs (default: 1)
#' @param min_overlap_percent Minimum overlap percentage (default: 0)
#' @param max_genes Maximum number of genes to return (default: 10)
#' @return Data.frame with overlapping genes and overlap statistics
#' @export
find_genes_by_boundary_overlap <- function(boundary_db, chr, start, end,
                                          min_overlap_bp = 1,
                                          min_overlap_percent = 0,
                                          max_genes = 10) {
  
  cat("Searching for genes overlapping:", chr, ":", start, "-", end, "\n")
  
  # Validate inputs
  if (is.na(chr) || is.na(start) || is.na(end)) {
    warning("Invalid coordinates provided")
    return(data.frame())
  }
  
  if (start > end) {
    warning("Start coordinate is greater than end coordinate")
    return(data.frame())
  }
  
  # Get genes on the specified chromosome
  chr_genes <- get_genes_by_chromosome(boundary_db, chr)
  
  if (nrow(chr_genes) == 0) {
    cat("No genes found on chromosome:", chr, "\n")
    return(data.frame())
  }
  
  cat("Checking", nrow(chr_genes), "genes on", chr, "\n")
  
  # Calculate overlaps using interval intersection
  overlapping_genes <- data.frame()
  
  for (i in 1:nrow(chr_genes)) {
    gene <- chr_genes[i, ]
    
    # Calculate overlap using interval intersection
    overlap_start <- max(start, gene$start)
    overlap_end <- min(end, gene$end)
    
    # Check if there's actual overlap
    if (overlap_start <= overlap_end) {
      overlap_bp <- overlap_end - overlap_start + 1
      query_length <- end - start + 1
      gene_length <- gene$end - gene$start + 1
      
      # Calculate overlap percentages
      overlap_percent_query <- (overlap_bp / query_length) * 100
      overlap_percent_gene <- (overlap_bp / gene_length) * 100
      
      # Apply filters
      if (overlap_bp >= min_overlap_bp && 
          overlap_percent_query >= min_overlap_percent) {
        
        overlapping_genes <- rbind(overlapping_genes, data.frame(
          gene_id = gene$gene_id,
          gene_symbol = gene$gene_symbol,
          gene_type = gene$gene_type,
          chr = gene$chr,
          gene_start = gene$start,
          gene_end = gene$end,
          gene_strand = gene$strand,
          gene_width = gene$width,
          overlap_start = overlap_start,
          overlap_end = overlap_end,
          overlap_bp = overlap_bp,
          overlap_percent_query = round(overlap_percent_query, 2),
          overlap_percent_gene = round(overlap_percent_gene, 2),
          confidence_score = calculate_overlap_confidence(overlap_percent_query, overlap_percent_gene, overlap_bp),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  if (nrow(overlapping_genes) == 0) {
    cat("No overlapping genes found\n")
    return(data.frame())
  }
  
  # Sort by confidence score (descending)
  overlapping_genes <- overlapping_genes[order(-overlapping_genes$confidence_score), ]
  
  # Limit results
  if (nrow(overlapping_genes) > max_genes) {
    overlapping_genes <- overlapping_genes[1:max_genes, ]
  }
  
  cat("Found", nrow(overlapping_genes), "overlapping genes\n")
  
  return(overlapping_genes)
}

#' Calculate Overlap Confidence Score
#' 
#' Calculates a confidence score based on overlap metrics
#' 
#' @param overlap_percent_query Percentage of query region overlapped
#' @param overlap_percent_gene Percentage of gene overlapped  
#' @param overlap_bp Overlap in base pairs
#' @return Confidence score (0-100)
calculate_overlap_confidence <- function(overlap_percent_query, overlap_percent_gene, overlap_bp) {
  
  # Base score from query overlap percentage (heavily weighted)
  base_score <- overlap_percent_query * 0.6
  
  # Bonus for gene overlap percentage
  gene_overlap_bonus <- overlap_percent_gene * 0.2
  
  # Bonus for substantial overlap in base pairs
  bp_bonus <- pmin(log10(overlap_bp + 1) * 5, 20)  # Cap at 20 points
  
  # Calculate final score
  confidence <- base_score + gene_overlap_bonus + bp_bonus
  
  # Cap at 100
  confidence <- pmin(confidence, 100)
  
  return(round(confidence, 2))
}

#' Find Genes by Multiple Coordinate Ranges
#' 
#' Finds genes that overlap with any of multiple coordinate ranges (e.g., exons)
#' 
#' @param boundary_db Gene boundary database
#' @param coordinate_ranges Data.frame with chr, start, end columns
#' @param min_overlap_bp Minimum overlap in base pairs
#' @param min_overlap_percent Minimum overlap percentage
#' @param max_genes Maximum number of genes to return
#' @return Data.frame with overlapping genes
#' @export
find_genes_by_multiple_ranges <- function(boundary_db, coordinate_ranges,
                                         min_overlap_bp = 1,
                                         min_overlap_percent = 0,
                                         max_genes = 10) {
  
  if (nrow(coordinate_ranges) == 0) {
    return(data.frame())
  }
  
  cat("Searching for genes overlapping", nrow(coordinate_ranges), "coordinate ranges\n")
  
  all_overlapping_genes <- list()
  
  for (i in 1:nrow(coordinate_ranges)) {
    range_data <- coordinate_ranges[i, ]
    
    overlapping <- find_genes_by_boundary_overlap(
      boundary_db = boundary_db,
      chr = range_data$chr,
      start = range_data$start,
      end = range_data$end,
      min_overlap_bp = min_overlap_bp,
      min_overlap_percent = min_overlap_percent,
      max_genes = max_genes
    )
    
    if (nrow(overlapping) > 0) {
      overlapping$range_index <- i
      all_overlapping_genes[[i]] <- overlapping
    }
  }
  
  if (length(all_overlapping_genes) == 0) {
    return(data.frame())
  }
  
  # Combine all results
  combined_results <- do.call(rbind, all_overlapping_genes)
  
  # Aggregate by gene (sum overlaps, take max confidence)
  aggregated_results <- aggregate(
    cbind(overlap_bp, confidence_score) ~ gene_id + gene_symbol + gene_type + 
          chr + gene_start + gene_end + gene_strand + gene_width,
    data = combined_results,
    FUN = function(x) c(sum = sum(x), max = max(x))
  )
  
  # Flatten the aggregated columns
  final_results <- data.frame(
    gene_id = aggregated_results$gene_id,
    gene_symbol = aggregated_results$gene_symbol,
    gene_type = aggregated_results$gene_type,
    chr = aggregated_results$chr,
    gene_start = aggregated_results$gene_start,
    gene_end = aggregated_results$gene_end,
    gene_strand = aggregated_results$gene_strand,
    gene_width = aggregated_results$gene_width,
    total_overlap_bp = aggregated_results$overlap_bp[, "sum"],
    max_confidence_score = aggregated_results$confidence_score[, "max"],
    num_ranges_overlapped = as.numeric(table(combined_results$gene_id)[aggregated_results$gene_id]),
    stringsAsFactors = FALSE
  )
  
  # Sort by max confidence score
  final_results <- final_results[order(-final_results$max_confidence_score), ]
  
  # Limit results
  if (nrow(final_results) > max_genes) {
    final_results <- final_results[1:max_genes, ]
  }
  
  cat("Found", nrow(final_results), "genes overlapping multiple ranges\n")
  
  return(final_results)
}

#' Extract Genomic Coordinates from GTF File
#' 
#' Extracts genomic coordinates from a novel isoform GTF file
#' 
#' @param gtf_file Path to GTF file
#' @return Data.frame with genomic coordinates
#' @export
extract_coordinates_from_gtf <- function(gtf_file) {
  
  if (!file.exists(gtf_file)) {
    warning("GTF file not found: ", gtf_file)
    return(data.frame())
  }
  
  cat("Extracting coordinates from GTF file:", gtf_file, "\n")
  
  # Read GTF file
  gtf_lines <- readLines(gtf_file)
  
  # Filter for exon features
  exon_lines <- gtf_lines[grepl("\\sexon\\s", gtf_lines)]
  
  if (length(exon_lines) == 0) {
    # Try with transcript features if no exons
    exon_lines <- gtf_lines[grepl("\\stranscript\\s", gtf_lines)]
  }
  
  if (length(exon_lines) == 0) {
    warning("No exon or transcript features found in GTF file")
    return(data.frame())
  }
  
  coordinates <- data.frame()
  
  for (line in exon_lines) {
    parts <- strsplit(line, "\\t")[[1]]
    
    if (length(parts) >= 5) {
      coordinates <- rbind(coordinates, data.frame(
        chr = parts[1],
        start = as.numeric(parts[4]),
        end = as.numeric(parts[5]),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Remove duplicates and sort
  coordinates <- unique(coordinates)
  coordinates <- coordinates[order(coordinates$chr, coordinates$start), ]
  
  cat("Extracted", nrow(coordinates), "coordinate ranges\n")
  
  return(coordinates)
}

#' Search Genes for Novel Isoform
#' 
#' Main function to search for genes overlapping a novel isoform
#' 
#' @param boundary_db Gene boundary database
#' @param novel_coordinates Novel isoform coordinates (chr, start, end) or GTF file path
#' @param search_method "single" for single range, "multiple" for multiple ranges
#' @param min_overlap_bp Minimum overlap in base pairs
#' @param min_overlap_percent Minimum overlap percentage
#' @param max_genes Maximum number of genes to return
#' @return Data.frame with candidate genes
#' @export
search_genes_for_novel_isoform <- function(boundary_db, novel_coordinates,
                                          search_method = "single",
                                          min_overlap_bp = 100,
                                          min_overlap_percent = 10,
                                          max_genes = 10) {
  
  cat("=== Searching Genes for Novel Isoform ===\n")
  cat("Search method:", search_method, "\n")
  cat("Minimum overlap:", min_overlap_bp, "bp,", min_overlap_percent, "%\n")
  cat("Max genes:", max_genes, "\n")
  
  # Handle different input types
  if (is.character(novel_coordinates) && length(novel_coordinates) == 1) {
    # Assume it's a GTF file path
    cat("Extracting coordinates from GTF file\n")
    coordinates <- extract_coordinates_from_gtf(novel_coordinates)
    search_method <- "multiple"
  } else if (is.data.frame(novel_coordinates)) {
    coordinates <- novel_coordinates
  } else {
    stop("Invalid novel_coordinates input. Provide data.frame or GTF file path.")
  }
  
  if (nrow(coordinates) == 0) {
    warning("No valid coordinates provided")
    return(data.frame())
  }
  
  # Perform search based on method
  if (search_method == "single" && nrow(coordinates) == 1) {
    results <- find_genes_by_boundary_overlap(
      boundary_db = boundary_db,
      chr = coordinates$chr[1],
      start = coordinates$start[1],
      end = coordinates$end[1],
      min_overlap_bp = min_overlap_bp,
      min_overlap_percent = min_overlap_percent,
      max_genes = max_genes
    )
  } else {
    results <- find_genes_by_multiple_ranges(
      boundary_db = boundary_db,
      coordinate_ranges = coordinates,
      min_overlap_bp = min_overlap_bp,
      min_overlap_percent = min_overlap_percent,
      max_genes = max_genes
    )
  }
  
  if (nrow(results) > 0) {
    cat("Gene search completed successfully\n")
    cat("Top candidate genes:\n")
    for (i in 1:min(3, nrow(results))) {
      gene <- results[i, ]
      cat(sprintf("  %d. %s (%s) - Confidence: %.1f\n", 
                  i, gene$gene_symbol, gene$gene_id, 
                  ifelse("confidence_score" %in% names(gene), 
                         gene$confidence_score, 
                         gene$max_confidence_score)))
    }
  } else {
    cat("No candidate genes found\n")
  }
  
  return(results)
}