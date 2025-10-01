# rMATS Event Parser Module
# Handles coordinate-precise parsing of all rMATS event types
# Critical: Properly handles BED-style coordinate conversion

source("R/coordinate_utils.R")

# Main function to detect event type from rMATS data
detect_event_type <- function(rmats_data) {
  colnames_data <- colnames(rmats_data)
  
  if ("exonStart_0base" %in% colnames_data && "upstreamES" %in% colnames_data) {
    return("SE")
  } else if ("longExonStart_0base" %in% colnames_data && "shortES" %in% colnames_data) {
    return("A3SS")
  } else if ("shortExonEnd" %in% colnames_data && "longExonEnd" %in% colnames_data) {
    return("A5SS")
  } else if (("1stExonStart_0base" %in% colnames_data && "2ndExonStart_0base" %in% colnames_data) ||
             ("X1stExonStart_0base" %in% colnames_data && "X2ndExonStart_0base" %in% colnames_data)) {
    return("MXE")
  } else if ("riExonStart_0base" %in% colnames_data) {
    return("RI")
  } else {
    stop("Unknown rMATS event type - check column names")
  }
}

# Validate rMATS data structure
validate_rmats_data <- function(rmats_data, event_type = NULL) {
  if (is.null(event_type)) {
    event_type <- detect_event_type(rmats_data)
  }
  
  # Special handling for MXE which can have X-prefixed columns
  if (event_type == "MXE") {
    # Check if we have either original or X-prefixed columns
    original_cols <- c("1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd")
    x_prefixed_cols <- c("X1stExonStart_0base", "X1stExonEnd", "X2ndExonStart_0base", "X2ndExonEnd")
    
    has_original <- all(original_cols %in% colnames(rmats_data))
    has_x_prefixed <- all(x_prefixed_cols %in% colnames(rmats_data))
    
    if (!has_original && !has_x_prefixed) {
      stop(paste("Missing required MXE columns. Need either:",
                 "\n  Original format:", paste(original_cols, collapse = ", "),
                 "\n  Or X-prefixed format:", paste(x_prefixed_cols, collapse = ", ")))
    }
    
    # Check other required columns (non-MXE specific)
    other_required <- c("ID", "GeneID", "geneSymbol", "chr", "strand", 
                       "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
                       "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
                       "IncFormLen", "SkipFormLen", "PValue", "FDR")
    missing_other <- setdiff(other_required, colnames(rmats_data))
    
    if (length(missing_other) > 0) {
      stop(paste("Missing required columns for", event_type, "events:", 
                 paste(missing_other, collapse = ", ")))
    }
    
  } else {
    # Standard validation for non-MXE events
    required_cols <- get_required_columns(event_type)
    missing_cols <- setdiff(required_cols, colnames(rmats_data))
    
    if (length(missing_cols) > 0) {
      stop(paste("Missing required columns for", event_type, "events:", 
                 paste(missing_cols, collapse = ", ")))
    }
  }
  
  # Check for coordinate consistency
  validate_coordinate_consistency(rmats_data, event_type)
  
  return(TRUE)
}

# Get required columns for each event type
get_required_columns <- function(event_type) {
  # Base columns that appear in all rMATS files
  base_cols <- c("ID", "GeneID", "geneSymbol", "chr", "strand")
  
  # Event-specific coordinate columns
  event_specific <- switch(event_type,
    "SE" = c("exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", 
             "downstreamES", "downstreamEE"),
    "A3SS" = c("longExonStart_0base", "longExonEnd", "shortES", "shortEE", 
               "flankingES", "flankingEE"),
    "A5SS" = c("longExonStart_0base", "longExonEnd", "shortES", "shortEE", 
               "flankingES", "flankingEE"),
    "MXE" = c("1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd",
              "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"),
    "RI" = c("riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", 
             "downstreamES", "downstreamEE"),
    stop("Unknown event type")
  )
  
  # Standard rMATS statistical columns (these appear in all files)
  stats_cols <- c("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
                  "IncFormLen", "SkipFormLen", "PValue", "FDR", 
                  "IncLevel1", "IncLevel2", "IncLevelDifference")
  
  return(c(base_cols, event_specific, stats_cols))
}

# Validate coordinate consistency within events
validate_coordinate_consistency <- function(rmats_data, event_type) {
  for (i in 1:nrow(rmats_data)) {
    row_data <- rmats_data[i, ]
    
    tryCatch({
      coords <- extract_event_coordinates(row_data, event_type)
      
      # Check that all coordinates are positive
      all_coords <- unlist(coords[grepl("_coords$", names(coords))])
      if (any(all_coords <= 0, na.rm = TRUE)) {
        warning(paste("Non-positive coordinates found in row", i))
      }
      
    }, error = function(e) {
      warning(paste("Coordinate validation failed for row", i, ":", e$message))
    })
  }
  
  return(TRUE)
}

# Extract and convert coordinates for a single event
extract_event_coordinates <- function(event_row, event_type) {
  if (is.null(event_type)) {
    event_type <- detect_event_type(data.frame(t(event_row)))
  }
  
  result <- list(
    event_type = event_type,
    gene_id = as.character(event_row$GeneID),
    gene_symbol = as.character(event_row$geneSymbol),
    chromosome = as.character(event_row$chr),
    strand = as.character(event_row$strand)
  )
  
  # Extract event-specific coordinates
  if (event_type == "SE") {
    result <- extract_SE_coordinates(event_row, result)
  } else if (event_type == "A3SS") {
    result <- extract_A3SS_coordinates(event_row, result)
  } else if (event_type == "A5SS") {
    result <- extract_A5SS_coordinates(event_row, result)
  } else if (event_type == "MXE") {
    result <- extract_MXE_coordinates(event_row, result)
  } else if (event_type == "RI") {
    result <- extract_RI_coordinates(event_row, result)
  } else {
    stop(paste("Unsupported event type:", event_type))
  }
  
  return(result)
}

# SE (Skipped Exon) coordinate extraction
extract_SE_coordinates <- function(event_row, result) {
  # CORRECTED: Handle coordinate systems properly
  # exonStart_0base: 0-based, needs +1
  skipped_start <- rmats_0base_to_genomic(event_row$exonStart_0base)
  # exonEnd: 1-based, already correct
  skipped_end <- rmats_1base_to_genomic(event_row$exonEnd)
  
  # CRITICAL FIX: ES coordinates are 0-based, EE coordinates are 1-based
  upstream_start <- rmats_0base_to_genomic(event_row$upstreamES)
  upstream_end <- rmats_1base_to_genomic(event_row$upstreamEE)
  
  downstream_start <- rmats_0base_to_genomic(event_row$downstreamES)
  downstream_end <- rmats_1base_to_genomic(event_row$downstreamEE)
  
  # Validate coordinates
  validate_coordinates(skipped_start, skipped_end)
  validate_coordinates(upstream_start, upstream_end)
  validate_coordinates(downstream_start, downstream_end)
  
  # Build exon coordinates (biological names, not positional)
  result$skipped_exon_coords <- list(start = skipped_start, end = skipped_end)
  result$upstream_exon_coords <- list(start = upstream_start, end = upstream_end)
  result$downstream_exon_coords <- list(start = downstream_start, end = downstream_end)
  
  # Build isoform structures with strand-aware ordering
  result$inclusion_isoform <- build_SE_inclusion_isoform(result)
  result$exclusion_isoform <- build_SE_exclusion_isoform(result)
  
  return(result)
}

# A3SS (Alternative 3' Splice Site) coordinate extraction
extract_A3SS_coordinates <- function(event_row, result) {
  # CORRECTED: longExonStart_0base is 0-based, longExonEnd is 1-based
  long_start <- rmats_0base_to_genomic(event_row$longExonStart_0base)
  long_end <- rmats_1base_to_genomic(event_row$longExonEnd)
  
  # CRITICAL FIX: ES coordinates are 0-based, EE coordinates are 1-based
  short_start <- rmats_0base_to_genomic(event_row$shortES)
  short_end <- rmats_1base_to_genomic(event_row$shortEE)
  
  flanking_start <- rmats_0base_to_genomic(event_row$flankingES)
  flanking_end <- rmats_1base_to_genomic(event_row$flankingEE)
  
  # Validate coordinates
  validate_coordinates(long_start, long_end)
  validate_coordinates(short_start, short_end)
  validate_coordinates(flanking_start, flanking_end)
  
  result$long_exon_coords <- list(start = long_start, end = long_end)
  result$short_exon_coords <- list(start = short_start, end = short_end)
  result$flanking_exon_coords <- list(start = flanking_start, end = flanking_end)
  
  # Build isoform structures
  result$inclusion_isoform <- build_A3SS_inclusion_isoform(result)
  result$exclusion_isoform <- build_A3SS_exclusion_isoform(result)
  
  return(result)
}

# A5SS (Alternative 5' Splice Site) coordinate extraction
extract_A5SS_coordinates <- function(event_row, result) {
  # CORRECTED: longExonStart_0base is 0-based, longExonEnd is 1-based
  long_start <- rmats_0base_to_genomic(event_row$longExonStart_0base)
  long_end <- rmats_1base_to_genomic(event_row$longExonEnd)
  
  # CRITICAL FIX: ES coordinates are 0-based, EE coordinates are 1-based
  short_start <- rmats_0base_to_genomic(event_row$shortES)
  short_end <- rmats_1base_to_genomic(event_row$shortEE)
  
  flanking_start <- rmats_0base_to_genomic(event_row$flankingES)
  flanking_end <- rmats_1base_to_genomic(event_row$flankingEE)
  
  # Validate coordinates
  validate_coordinates(long_start, long_end)
  validate_coordinates(short_start, short_end)
  validate_coordinates(flanking_start, flanking_end)
  
  result$long_exon_coords <- list(start = long_start, end = long_end)
  result$short_exon_coords <- list(start = short_start, end = short_end)
  result$flanking_exon_coords <- list(start = flanking_start, end = flanking_end)
  
  # Build isoform structures
  result$inclusion_isoform <- build_A5SS_inclusion_isoform(result)
  result$exclusion_isoform <- build_A5SS_exclusion_isoform(result)
  
  return(result)
}

# MXE (Mutually Exclusive Exons) coordinate extraction
extract_MXE_coordinates <- function(event_row, result) {
  # CORRECTED: 0base coordinates need conversion, ES/EE are 1-based
  # Handle both original column names and X-prefixed names (from read.table)
  first_start <- rmats_0base_to_genomic(
    if ("1stExonStart_0base" %in% names(event_row)) event_row$`1stExonStart_0base` 
    else event_row$X1stExonStart_0base
  )
  first_end <- rmats_1base_to_genomic(
    if ("1stExonEnd" %in% names(event_row)) event_row$`1stExonEnd`
    else event_row$X1stExonEnd
  )
  
  second_start <- rmats_0base_to_genomic(
    if ("2ndExonStart_0base" %in% names(event_row)) event_row$`2ndExonStart_0base`
    else event_row$X2ndExonStart_0base
  )
  second_end <- rmats_1base_to_genomic(
    if ("2ndExonEnd" %in% names(event_row)) event_row$`2ndExonEnd`
    else event_row$X2ndExonEnd
  )
  
  upstream_start <- rmats_0base_to_genomic(event_row$upstreamES)
  upstream_end <- rmats_1base_to_genomic(event_row$upstreamEE)
  
  downstream_start <- rmats_0base_to_genomic(event_row$downstreamES)
  downstream_end <- rmats_1base_to_genomic(event_row$downstreamEE)
  
  # Validate all coordinates
  validate_coordinates(first_start, first_end)
  validate_coordinates(second_start, second_end)
  validate_coordinates(upstream_start, upstream_end)
  validate_coordinates(downstream_start, downstream_end)
  
  result$first_exon_coords <- list(start = first_start, end = first_end)
  result$second_exon_coords <- list(start = second_start, end = second_end)
  result$upstream_exon_coords <- list(start = upstream_start, end = upstream_end)
  result$downstream_exon_coords <- list(start = downstream_start, end = downstream_end)
  
  # Build isoform structures
  result$inclusion_isoform <- build_MXE_inclusion_isoform(result)
  result$exclusion_isoform <- build_MXE_exclusion_isoform(result)
  
  return(result)
}

# RI (Retained Intron) coordinate extraction
extract_RI_coordinates <- function(event_row, result) {
  # CORRECTED: 0base coordinates need conversion, ES/EE are 1-based
  ri_start <- rmats_0base_to_genomic(event_row$riExonStart_0base)
  ri_end <- rmats_1base_to_genomic(event_row$riExonEnd)
  
  upstream_start <- rmats_0base_to_genomic(event_row$upstreamES)
  upstream_end <- rmats_1base_to_genomic(event_row$upstreamEE)
  
  downstream_start <- rmats_0base_to_genomic(event_row$downstreamES)
  downstream_end <- rmats_1base_to_genomic(event_row$downstreamEE)
  
  # Validate coordinates
  validate_coordinates(ri_start, ri_end)
  validate_coordinates(upstream_start, upstream_end)
  validate_coordinates(downstream_start, downstream_end)
  
  result$retained_intron_coords <- list(start = ri_start, end = ri_end)
  result$upstream_exon_coords <- list(start = upstream_start, end = upstream_end)
  result$downstream_exon_coords <- list(start = downstream_start, end = downstream_end)
  
  # Build isoform structures
  result$inclusion_isoform <- build_RI_inclusion_isoform(result)
  result$exclusion_isoform <- build_RI_exclusion_isoform(result)
  
  return(result)
}

# Isoform building functions for each event type

build_SE_inclusion_isoform <- function(coords) {
  # CORRECTED: Build inclusion isoform in biological 5' to 3' order
  # For translation, we need biological order, not genomic order
  
  if (coords$strand == "+") {
    # Plus strand: biological = genomic order
    # upstream (5') -> skipped -> downstream (3')
    exons <- list(
      coords$upstream_exon_coords,
      coords$skipped_exon_coords,
      coords$downstream_exon_coords
    )
  } else {
    # Minus strand: biological order is REVERSE of genomic order
    # downstream (5') -> skipped -> upstream (3')
    exons <- list(
      coords$downstream_exon_coords,  # 5' end
      coords$skipped_exon_coords,     # middle
      coords$upstream_exon_coords     # 3' end
    )
  }
  
  # For minus strand, we need to keep biological order for proper translation
  return(exons)
}

build_SE_exclusion_isoform <- function(coords) {
  # CORRECTED: Build exclusion isoform in biological 5' to 3' order
  
  if (coords$strand == "+") {
    # Plus strand: biological = genomic order
    # upstream (5') -> downstream (3') (skip middle)
    exons <- list(
      coords$upstream_exon_coords,
      coords$downstream_exon_coords
    )
  } else {
    # Minus strand: biological order is REVERSE of genomic order
    # downstream (5') -> upstream (3') (skip middle)
    exons <- list(
      coords$downstream_exon_coords,  # 5' end
      coords$upstream_exon_coords     # 3' end
    )
  }
  
  # Keep biological order for proper translation
  return(exons)
}

build_A3SS_inclusion_isoform <- function(coords) {
  # A3SS inclusion uses long exon (includes more of the alternative 3' end)
  exons <- list(
    coords$flanking_exon_coords,
    coords$long_exon_coords
  )
  return(order_exons_by_biological_position(exons, coords$strand))
}

build_A3SS_exclusion_isoform <- function(coords) {
  # A3SS exclusion uses short exon (less of the alternative 3' end)
  exons <- list(
    coords$flanking_exon_coords,
    coords$short_exon_coords
  )
  return(order_exons_by_biological_position(exons, coords$strand))
}

build_A5SS_inclusion_isoform <- function(coords) {
  # A5SS inclusion uses long exon (includes more of the alternative 5' end)
  exons <- list(
    coords$long_exon_coords,
    coords$flanking_exon_coords
  )
  return(order_exons_by_biological_position(exons, coords$strand))
}

build_A5SS_exclusion_isoform <- function(coords) {
  # A5SS exclusion uses short exon (less of the alternative 5' end)
  exons <- list(
    coords$short_exon_coords,
    coords$flanking_exon_coords
  )
  return(order_exons_by_biological_position(exons, coords$strand))
}

build_MXE_inclusion_isoform <- function(coords) {
  # MXE inclusion: upstream -> first exon -> downstream
  exons <- list(
    coords$upstream_exon_coords,
    coords$first_exon_coords,
    coords$downstream_exon_coords
  )
  return(order_exons_by_biological_position(exons, coords$strand))
}

build_MXE_exclusion_isoform <- function(coords) {
  # MXE exclusion: upstream -> second exon -> downstream
  exons <- list(
    coords$upstream_exon_coords,
    coords$second_exon_coords,
    coords$downstream_exon_coords
  )
  return(order_exons_by_biological_position(exons, coords$strand))
}

build_RI_inclusion_isoform <- function(coords) {
  # RI inclusion: The intron is retained and becomes part of the mature mRNA
  # For downstream functional analysis, we need THREE separate components:
  # 1. Upstream exon
  # 2. Retained intron (treated as an exon for translation)
  # 3. Downstream exon
  # This allows proper phase analysis and translation of each component
  
  if (coords$strand == "+") {
    # Plus strand: biological 5' to 3' order
    # upstream -> retained_intron -> downstream
    exons <- list(
      coords$upstream_exon_coords,     # 5' exon
      coords$retained_intron_coords,   # Retained intron (as exon)
      coords$downstream_exon_coords    # 3' exon
    )
  } else {
    # Minus strand: biological 5' to 3' order (reversed from genomic)
    # downstream -> retained_intron -> upstream
    exons <- list(
      coords$downstream_exon_coords,   # 5' end (biological)
      coords$retained_intron_coords,   # Retained intron (as exon)
      coords$upstream_exon_coords      # 3' end (biological)
    )
  }
  
  return(exons)
}

build_RI_exclusion_isoform <- function(coords) {
  # RI exclusion: intron is spliced out - only upstream and downstream
  
  if (coords$strand == "+") {
    # Plus strand: biological = genomic order
    # upstream (5') -> downstream (3')
    exons <- list(
      coords$upstream_exon_coords,
      coords$downstream_exon_coords
    )
  } else {
    # Minus strand: biological order is REVERSE of genomic order
    # downstream (5') -> upstream (3')
    exons <- list(
      coords$downstream_exon_coords,  # 5' end
      coords$upstream_exon_coords     # 3' end
    )
  }
  
  return(exons)
}

# Helper function to order exons by genomic position
order_exons_by_genomic_position <- function(exon_list) {
  # Always order by genomic start position (ascending)
  # This gives us the correct order for sequence extraction
  exon_order <- order(sapply(exon_list, function(x) x$start))
  return(exon_list[exon_order])
}

# Helper function to order exons by biological position (5' to 3')
order_exons_by_biological_position <- function(exon_list, strand) {
  if (length(exon_list) <= 1) {
    return(exon_list)
  }
  
  # CRITICAL FIX: Ensure sapply returns a numeric vector, not a list
  start_positions <- unlist(sapply(exon_list, function(x) x$start))
  
  if (strand == "+") {
    # Plus strand: biological = genomic order (5' to 3' = left to right)
    exon_order <- order(start_positions)
    return(exon_list[exon_order])
  } else {
    # Minus strand: biological order is reverse of genomic order (5' to 3' = right to left)
    exon_order <- order(start_positions, decreasing = TRUE)
    return(exon_list[exon_order])
  }
}

# Summary function for event coordinates
summarize_event_coordinates <- function(event_coords) {
  summary <- list(
    event_type = event_coords$event_type,
    gene_symbol = event_coords$gene_symbol,
    chromosome = event_coords$chromosome,
    strand = event_coords$strand,
    inclusion_exons = length(event_coords$inclusion_isoform),
    exclusion_exons = length(event_coords$exclusion_isoform)
  )
  
  # Add event-specific coordinate summaries
  if (event_coords$event_type == "SE") {
    summary$skipped_exon <- format_coordinates(
      event_coords$skipped_exon_coords$start,
      event_coords$skipped_exon_coords$end,
      event_coords$chromosome,
      event_coords$strand
    )
  }
  
  return(summary)
}

# Quality control function
check_event_quality <- function(event_coords) {
  quality <- list(
    coordinate_valid = TRUE,
    exon_order_valid = TRUE,
    strand_consistent = TRUE,
    warnings = c()
  )
  
  tryCatch({
    # Check for overlapping exons in isoforms
    inc_overlaps <- check_exon_overlaps(event_coords$inclusion_isoform)
    exc_overlaps <- check_exon_overlaps(event_coords$exclusion_isoform)
    
    if (inc_overlaps || exc_overlaps) {
      quality$exon_order_valid <- FALSE
      quality$warnings <- c(quality$warnings, "Overlapping exons detected")
    }
    
  }, error = function(e) {
    quality$coordinate_valid <- FALSE
    quality$warnings <- c(quality$warnings, paste("Validation error:", e$message))
  })
  
  return(quality)
}

# Build GTF-like structures for isoforms (Step 2)
build_gtf_structures <- function(event_coords, event_type) {
  gene_id <- event_coords$gene_id
  gene_symbol <- event_coords$gene_symbol
  chromosome <- event_coords$chromosome
  strand <- event_coords$strand
  
  # Generate synthetic transcript IDs for inclusion and exclusion isoforms
  inclusion_transcript_id <- paste0(gene_id, ".inclusion")
  exclusion_transcript_id <- paste0(gene_id, ".exclusion")
  
  # Build GTF-like records for inclusion isoform
  # CRITICAL FIX: Handle missing coordinates by filtering out invalid exons
  valid_inclusion_exons <- Filter(function(x) !is.null(x$start) && !is.null(x$end) && length(x$start) > 0, 
                                  event_coords$inclusion_isoform)
  
  if (length(valid_inclusion_exons) == 0) {
    stop("No valid exons found in inclusion isoform - cannot build GTF structure")
  }
  
  inclusion_gtf <- data.frame(
    seqname = chromosome,
    source = "rmats_prediction", 
    feature = "exon",
    start = unlist(sapply(valid_inclusion_exons, function(x) x$start)),
    end = unlist(sapply(valid_inclusion_exons, function(x) x$end)),
    score = ".",
    strand = strand,
    frame = ".",
    attributes = paste0(
      "gene_id \"", gene_id, "\"; ",
      "transcript_id \"", inclusion_transcript_id, "\"; ",
      "gene_name \"", gene_symbol, "\"; ",
      "exon_number \"", seq_along(valid_inclusion_exons), "\"; ",
      "isoform_type \"inclusion\";"
    ),
    stringsAsFactors = FALSE
  )
  
  # Build GTF-like records for exclusion isoform
  # CRITICAL FIX: Handle missing coordinates by filtering out invalid exons
  valid_exclusion_exons <- Filter(function(x) !is.null(x$start) && !is.null(x$end) && length(x$start) > 0, 
                                  event_coords$exclusion_isoform)
  
  if (length(valid_exclusion_exons) == 0) {
    stop("No valid exons found in exclusion isoform - cannot build GTF structure")
  }
  
  exclusion_gtf <- data.frame(
    seqname = chromosome,
    source = "rmats_prediction",
    feature = "exon", 
    start = unlist(sapply(valid_exclusion_exons, function(x) x$start)),
    end = unlist(sapply(valid_exclusion_exons, function(x) x$end)),
    score = ".",
    strand = strand,
    frame = ".",
    attributes = paste0(
      "gene_id \"", gene_id, "\"; ",
      "transcript_id \"", exclusion_transcript_id, "\"; ",
      "gene_name \"", gene_symbol, "\"; ",
      "exon_number \"", seq_along(valid_exclusion_exons), "\"; ",
      "isoform_type \"exclusion\";"
    ),
    stringsAsFactors = FALSE
  )
  
  # Return structured result
  list(
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    chromosome = chromosome,
    strand = strand,
    event_type = event_type,
    inclusion_transcript_id = inclusion_transcript_id,
    exclusion_transcript_id = exclusion_transcript_id,
    inclusion_gtf = inclusion_gtf,
    exclusion_gtf = exclusion_gtf,
    combined_gtf = rbind(inclusion_gtf, exclusion_gtf),
    # Keep original coordinates for compatibility
    inclusion_isoform = event_coords$inclusion_isoform,
    exclusion_isoform = event_coords$exclusion_isoform
  )
}