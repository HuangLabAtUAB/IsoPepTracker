#===============================================================================
# COORDINATE COMPRESSION MODULE
# Handles intron compression for improved visualization clarity
#===============================================================================

# Fixed size for compressed introns (in display units)
COMPRESSED_INTRON_SIZE <- 50

#' Create compression mapping for genomic coordinates
#' 
#' @param exons_by_transcript List of GRanges objects for each transcript
#' @return List containing compression functions and mapping data
create_compression_map <- function(exons_by_transcript) {
  # Collect all exons from all transcripts
  all_exons <- unlist(GRangesList(exons_by_transcript))
  
  # Reduce overlapping ranges to create superset
  superset <- reduce(all_exons)
  
  # Sort superset by position
  superset <- sort(superset)
  
  # Create compression mapping data frame
  compressed_coords <- data.frame(
    original_start = start(superset),
    original_end = end(superset),
    stringsAsFactors = FALSE
  )
  
  # Calculate compressed coordinates
  current_pos <- start(superset)[1]  # Start from first exon's position
  compressed_coords$compressed_start <- NA
  compressed_coords$compressed_end <- NA
  
  for(i in seq_len(nrow(compressed_coords))) {
    if(i == 1) {
      # First exon starts at its original position
      compressed_coords$compressed_start[i] <- current_pos
      current_pos <- current_pos + (compressed_coords$original_end[i] - compressed_coords$original_start[i])
      compressed_coords$compressed_end[i] <- current_pos
    } else {
      # Add fixed intron size
      current_pos <- current_pos + COMPRESSED_INTRON_SIZE
      compressed_coords$compressed_start[i] <- current_pos
      # Add exon length
      current_pos <- current_pos + (compressed_coords$original_end[i] - compressed_coords$original_start[i])
      compressed_coords$compressed_end[i] <- current_pos
    }
  }
  
  # Function to compress a single genomic position
  compress_position <- function(pos) {
    # If position is before first exon, return as is
    if(pos <= compressed_coords$original_start[1]) {
      return(pos)
    }
    
    # If position is after last exon, adjust by total compression
    if(pos >= compressed_coords$original_end[nrow(compressed_coords)]) {
      total_compression <- (compressed_coords$original_end[nrow(compressed_coords)] - 
                          compressed_coords$original_start[1]) -
                         (compressed_coords$compressed_end[nrow(compressed_coords)] - 
                          compressed_coords$compressed_start[1])
      return(pos - total_compression)
    }
    
    # Find which interval contains this position
    for(i in seq_len(nrow(compressed_coords))) {
      # If position is in current exon
      if(pos >= compressed_coords$original_start[i] && 
         pos <= compressed_coords$original_end[i]) {
        offset <- pos - compressed_coords$original_start[i]
        return(compressed_coords$compressed_start[i] + offset)
      }
      
      # If position is in intron before next exon
      if(i < nrow(compressed_coords) && 
         pos > compressed_coords$original_end[i] && 
         pos < compressed_coords$original_start[i + 1]) {
        # Linear interpolation within intron
        intron_pos <- (pos - compressed_coords$original_end[i]) / 
                     (compressed_coords$original_start[i + 1] - compressed_coords$original_end[i])
        return(compressed_coords$compressed_end[i] + 
               intron_pos * COMPRESSED_INTRON_SIZE)
      }
    }
    return(pos)  # Fallback
  }
  
  # Function to decompress a display position back to genomic
  decompress_position <- function(display_pos) {
    # Reverse mapping - find which compressed interval contains this position
    for(i in seq_len(nrow(compressed_coords))) {
      if(display_pos >= compressed_coords$compressed_start[i] && 
         display_pos <= compressed_coords$compressed_end[i]) {
        # Position is in a compressed exon
        offset <- display_pos - compressed_coords$compressed_start[i]
        return(compressed_coords$original_start[i] + offset)
      }
      
      # Check if in compressed intron
      if(i < nrow(compressed_coords) &&
         display_pos > compressed_coords$compressed_end[i] &&
         display_pos < compressed_coords$compressed_start[i + 1]) {
        # In compressed intron - interpolate back
        intron_progress <- (display_pos - compressed_coords$compressed_end[i]) / COMPRESSED_INTRON_SIZE
        original_intron_start <- compressed_coords$original_end[i]
        original_intron_end <- compressed_coords$original_start[i + 1]
        original_intron_length <- original_intron_end - original_intron_start
        return(original_intron_start + intron_progress * original_intron_length)
      }
    }
    return(display_pos)  # Fallback
  }
  
  # Return compression map with all needed components
  return(list(
    compress = compress_position,
    decompress = decompress_position,
    coords = compressed_coords,
    superset = superset
  ))
}

#' Create intron marker data for zigzag visualization between actual transcript exons
#' 
#' @param exons_by_transcript List of GRanges for each transcript
#' @param compression_map Compression mapping from create_compression_map
#' @param y_position Y position for the markers (default 0 for relative positioning)
#' @return Data frame with zigzag line coordinates
create_intron_markers_for_transcripts <- function(exons_by_transcript, compression_map, y_position = 0) {
  intron_markers <- data.frame(
    x = numeric(),
    y = numeric(),
    group = character(),
    transcript = character(),
    stringsAsFactors = FALSE
  )
  
  # For each transcript, find gaps between consecutive exons
  for(tx_name in names(exons_by_transcript)) {
    tx_exons <- exons_by_transcript[[tx_name]]
    tx_exons <- sort(tx_exons)  # Sort by position
    
    # Skip if transcript has less than 2 exons (no introns to mark)
    if(length(tx_exons) < 2) next
    
    # Find gaps between consecutive exons in this transcript
    for(i in 2:length(tx_exons)) {
      # Get the gap between this exon and the previous one
      prev_exon_end <- end(tx_exons[i-1])
      curr_exon_start <- start(tx_exons[i])
      
      # Calculate intron size
      intron_size <- curr_exon_start - prev_exon_end
      
      # Put zigzag lines between ALL consecutive exons (no thresholds!)
      # Compress the coordinates
      compressed_prev_end <- compression_map$compress(prev_exon_end)
      compressed_curr_start <- compression_map$compress(curr_exon_start)
      
      # Create straight line for intron gap (no zigzag)
      x_points <- c(compressed_prev_end, compressed_curr_start)
      # Use straight line at the y_position
      y_points <- c(y_position, y_position)
      
      group_id <- paste(tx_name, i, sep = "_")
      
      intron_markers <- rbind(intron_markers, data.frame(
        x = x_points,
        y = y_points,
        group = group_id,
        transcript = tx_name,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(intron_markers)
}

# Keep the old function for backward compatibility
create_intron_markers <- function(compressed_coords, y_position = 0) {
  # This is now a simplified version - may not be needed anymore
  return(data.frame(x = numeric(), y = numeric(), group = integer(), stringsAsFactors = FALSE))
}

#' Create intron length labels for compressed regions
#' 
#' @param compressed_coords Data frame with compressed coordinate mapping
#' @param y_position Y position for the labels
#' @return Data frame with label positions and text
create_intron_labels <- function(compressed_coords, y_position = 1) {
  intron_labels <- data.frame(
    x = numeric(),
    y = numeric(),
    label = character(),
    stringsAsFactors = FALSE
  )
  
  for(i in 2:nrow(compressed_coords)) {
    # Calculate actual intron size
    actual_size <- compressed_coords$original_start[i] - compressed_coords$original_end[i-1]
    
    # Only label significant introns
    if(actual_size > 5000) {
      # Position label in middle of compressed intron
      label_x <- (compressed_coords$compressed_end[i-1] + compressed_coords$compressed_start[i]) / 2
      
      # Format size nicely
      if(actual_size > 1000000) {
        size_label <- paste0(round(actual_size/1000000, 1), "Mb")
      } else if(actual_size > 1000) {
        size_label <- paste0(round(actual_size/1000, 1), "kb")
      } else {
        size_label <- paste0(actual_size, "bp")
      }
      
      intron_labels <- rbind(intron_labels, data.frame(
        x = label_x,
        y = y_position - 0.5,  # Position below the line
        label = size_label,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(intron_labels)
}

#' Compress a data frame of genomic coordinates
#' 
#' @param coord_df Data frame with start and end columns
#' @param compression_map Compression map from create_compression_map
#' @return Data frame with compressed coordinates
compress_coordinates_df <- function(coord_df, compression_map) {
  if(nrow(coord_df) == 0) return(coord_df)
  
  coord_df$start <- sapply(coord_df$start, compression_map$compress)
  coord_df$end <- sapply(coord_df$end, compression_map$compress)
  
  return(coord_df)
}

#' Check if compression would be beneficial for a gene
#' 
#' @param gene_start Start position of gene
#' @param gene_end End position of gene
#' @param exon_count Number of exons
#' @return Logical indicating if compression is recommended
should_use_compression <- function(gene_start, gene_end, exon_count = NULL) {
  gene_span <- gene_end - gene_start
  
  # Recommend compression for genes > 50kb
  # or genes with many exons (likely to have long introns)
  if(gene_span > 50000) {
    return(TRUE)
  }
  
  if(!is.null(exon_count) && exon_count > 10 && gene_span > 20000) {
    return(TRUE)
  }
  
  return(FALSE)
}

# REMOVED: create_compressed_axis_breaks function
# This function was causing coordinate mismatch issues by overcomplicating the axis creation.
# Now using simple approach: real GTF coordinates as labels positioned at compressed coordinates.