# Coordinate System Utilities for rMATS Analysis
# Critical functions for handling BED-style coordinates from rMATS

# CORRECTED: Convert rMATS coordinates to 1-based genomic coordinates
# Critical: Different coordinate systems in rMATS:
# - *_0base coordinates: 0-based inclusive (need +1)
# - *ES coordinates: 1-based exclusive (already correct for genomic)
# - *EE, *End coordinates: 1-based exclusive (already correct for genomic)

rmats_0base_to_genomic <- function(start_0base) {
  if (any(start_0base < 0, na.rm = TRUE)) {
    warning("Negative coordinates detected - check input data")
  }
  return(start_0base + 1)
}

# ES/EE coordinates are already 1-based (splice site coordinates)
rmats_1base_to_genomic <- function(coord_1base) {
  if (any(coord_1base <= 0, na.rm = TRUE)) {
    warning("Invalid coordinates detected - check input data")
  }
  return(coord_1base)
}

# Alias functions for clarity (referenced in other functions)
rmats_start_to_genomic <- function(start_coord) {
  return(rmats_0base_to_genomic(start_coord))
}

rmats_end_to_genomic <- function(end_coord) {
  return(rmats_1base_to_genomic(end_coord))
}

# Comprehensive coordinate validation
validate_coordinates <- function(start_genomic, end_genomic, allow_na = FALSE) {
  if (!allow_na && (any(is.na(start_genomic)) || any(is.na(end_genomic)))) {
    stop("NA coordinates detected - cannot proceed with analysis")
  }
  
  # Remove NA values for validation
  valid_indices <- !is.na(start_genomic) & !is.na(end_genomic)
  start_valid <- start_genomic[valid_indices]
  end_valid <- end_genomic[valid_indices]
  
  if (any(start_valid >= end_valid)) {
    invalid_indices <- which(start_genomic >= end_genomic)
    stop(paste("Invalid coordinates: start >= end at positions:", 
               paste(invalid_indices, collapse = ", ")))
  }
  
  if (any(start_valid <= 0) || any(end_valid <= 0)) {
    stop("Coordinates must be positive (1-based genomic coordinates)")
  }
  
  return(TRUE)
}

# Convert full rMATS coordinate set for an event
convert_rmats_coordinates <- function(coords_list) {
  result <- list()
  
  for (coord_name in names(coords_list)) {
    if (grepl("Start.*0base", coord_name) || grepl("ES$", coord_name)) {
      # These are 0-based starts, convert to 1-based
      result[[coord_name]] <- rmats_start_to_genomic(coords_list[[coord_name]])
    } else if (grepl("End$", coord_name) || grepl("EE$", coord_name)) {
      # These are already 1-based ends
      result[[coord_name]] <- rmats_end_to_genomic(coords_list[[coord_name]])
    } else {
      # Keep other coordinates as-is but validate
      result[[coord_name]] <- coords_list[[coord_name]]
    }
  }
  
  return(result)
}

# Build exon coordinate list from rMATS event
build_exon_coordinates <- function(start_coords, end_coords, strand = "+") {
  if (length(start_coords) != length(end_coords)) {
    stop("Start and end coordinate vectors must have same length")
  }
  
  # Convert coordinates
  genomic_starts <- rmats_start_to_genomic(start_coords)
  genomic_ends <- rmats_end_to_genomic(end_coords)
  
  # Validate each exon
  validate_coordinates(genomic_starts, genomic_ends)
  
  # Create exon list
  exons <- list()
  for (i in seq_along(genomic_starts)) {
    exons[[i]] <- list(
      start = genomic_starts[i],
      end = genomic_ends[i],
      length = genomic_ends[i] - genomic_starts[i] + 1
    )
  }
  
  # Sort exons by genomic position (handle strand later)
  exon_order <- order(sapply(exons, function(x) x$start))
  exons <- exons[exon_order]
  
  return(exons)
}

# Calculate total isoform length from exon list
calculate_isoform_length <- function(exon_list) {
  if (length(exon_list) == 0) return(0)
  
  total_length <- sum(sapply(exon_list, function(exon) exon$length))
  return(total_length)
}

# Check for overlapping exons (quality control)
check_exon_overlaps <- function(exon_list) {
  if (length(exon_list) <= 1) return(FALSE)
  
  for (i in 1:(length(exon_list) - 1)) {
    current_end <- exon_list[[i]]$end
    next_start <- exon_list[[i + 1]]$start
    
    if (current_end >= next_start) {
      warning(paste("Overlapping exons detected between positions", i, "and", i + 1))
      return(TRUE)
    }
  }
  
  return(FALSE)
}

# Format coordinates for display
format_coordinates <- function(start, end, chromosome = NULL, strand = NULL) {
  coord_str <- paste0(start, "-", end)
  
  if (!is.null(chromosome)) {
    coord_str <- paste0(chromosome, ":", coord_str)
  }
  
  if (!is.null(strand)) {
    coord_str <- paste0(coord_str, " (", strand, ")")
  }
  
  return(coord_str)
}

# Create coordinate summary for an event
summarize_event_coordinates <- function(event_data, event_type) {
  summary <- list(
    event_type = event_type,
    chromosome = event_data$chr,
    strand = event_data$strand,
    gene_id = event_data$GeneID,
    gene_symbol = event_data$geneSymbol
  )
  
  # Add event-specific coordinate information
  if (event_type == "SE") {
    summary$skipped_exon <- format_coordinates(
      rmats_start_to_genomic(event_data$exonStart_0base),
      event_data$exonEnd,
      event_data$chr,
      event_data$strand
    )
    summary$upstream_exon <- format_coordinates(
      rmats_start_to_genomic(event_data$upstreamES),
      event_data$upstreamEE,
      event_data$chr,
      event_data$strand
    )
    summary$downstream_exon <- format_coordinates(
      rmats_start_to_genomic(event_data$downstreamES),
      event_data$downstreamEE,
      event_data$chr,
      event_data$strand
    )
  }
  
  return(summary)
}