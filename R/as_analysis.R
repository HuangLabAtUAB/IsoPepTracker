#===============================================================================
# ALTERNATIVE SPLICING FUNCTIONS
#===============================================================================

# Function to create superset of all exons for a gene
create_exon_superset <- function(exons_result) {
  if (!exons_result$success) return(NULL)
  
  # Collect all exons from all transcripts
  all_exons <- unlist(GRangesList(exons_result$exons))
  
  # Reduce overlapping ranges to create superset
  superset <- reduce(all_exons)
  
  return(superset)
}

# Function to create AS visualization and event table
create_as_view <- function(exons_result, gene_id, gene_symbol, as_database) {
  # Get superset of exons
  superset <- create_exon_superset(exons_result)
  if (is.null(superset)) return(NULL)
  
  # Get AS events for this gene
  gene_events <- as_database[as_database$geneID == gene_id,]
  
  # Sort superset by position
  superset <- sort(superset)
  
  # Extract chromosome and strand information
  chromosome <- as.character(seqnames(superset)[1])
  strand_char <- as.character(strand(superset)[1])
  if (strand_char == "+" || strand_char == "-") {
    strand_direction <- strand_char
  } else {
    strand_direction <- "unknown"
  }
  
  # Create strand direction display (plots are always shown 5' to 3')
  strand_display <- ""
  
  # Create compression mapping
  INTRON_SIZE <- 50  # Fixed size for compressed introns
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
      current_pos <- current_pos + INTRON_SIZE
      compressed_coords$compressed_start[i] <- current_pos
      # Add exon length
      current_pos <- current_pos + (compressed_coords$original_end[i] - compressed_coords$original_start[i])
      compressed_coords$compressed_end[i] <- current_pos
    }
  }
  
  # Function to compress a genomic position
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
               intron_pos * INTRON_SIZE)
      }
    }
    return(pos)  # Fallback
  }
  
  # Create superset data frame with compressed coordinates
  superset_df <- data.frame(
    start = compressed_coords$compressed_start,
    end = compressed_coords$compressed_end,
    type = "superset",
    stringsAsFactors = FALSE
  )
  
  # Process AS events
  event_regions <- list()
  
  for(event_id in unique(gene_events$eventID)) {
    event <- gene_events[gene_events$eventID == event_id,][1,]
    
    # Extract AS ranges
    as_range <- NULL
    as_range1 <- NULL
    as_range2 <- NULL
    
    if("AS_range" %in% names(gene_events) && !is.null(event$AS_range[[1]])) {
      as_range <- event$AS_range[[1]]
    }
    if("AS_range1" %in% names(gene_events) && !is.null(event$AS_range1[[1]])) {
      as_range1 <- event$AS_range1[[1]]
    }
    if("AS_range2" %in% names(gene_events) && !is.null(event$AS_range2[[1]])) {
      as_range2 <- event$AS_range2[[1]]
    }
    
    i <- which(unique(gene_events$eventID) == event_id)
    
    # Process ranges with compression
    process_range <- function(range, type_label, event_type, hover_prefix) {
      if(!is.null(range) && length(range) > 0) {
        compressed_start <- compress_position(start(range))
        compressed_end <- compress_position(end(range))
        
        data.frame(
          start = compressed_start,
          end = compressed_end,
          type = type_label,
          event_type = event_type,
          event_id = event$eventID,
          y_pos = i + 1.5,
          hover_text = clean_hover_text(paste0(
            hover_prefix, "<br>",
            event$eventID, "<br>",
            "Original position: ", start(range), "-", end(range)
          )),
          stringsAsFactors = FALSE
        )
      }
    }
    
    # Add processed ranges to event_regions
    if(event$AS_type == "SE") {
      event_regions[[length(event_regions) + 1]] <- 
        process_range(as_range, "SE", "SE", "SE: Skipped Exon")
    } else if(event$AS_type == "RI") {
      event_regions[[length(event_regions) + 1]] <- 
        process_range(as_range, "RI", "RI", "RI: Retained Intron")
    } else if(event$AS_type == "MX") {
      event_regions[[length(event_regions) + 1]] <- 
        process_range(as_range1, "MXE", "MX", "MXE: Mutually Exclusive Exon 1")
      event_regions[[length(event_regions) + 1]] <- 
        process_range(as_range2, "MXE", "MX", "MXE: Mutually Exclusive Exon 2")
    } else if(event$AS_type == "A3") {
      event_regions[[length(event_regions) + 1]] <- 
        process_range(as_range, "A3SS", "A3", "A3SS: Alternative 3' Splice Site")
    } else if(event$AS_type == "A5") {
      event_regions[[length(event_regions) + 1]] <- 
        process_range(as_range, "A5SS", "A5", "A5SS: Alternative 5' Splice Site")
    }
  }
  
  # Create the plot
  if(length(event_regions) > 0) {
    event_df <- do.call(rbind, event_regions)
    
    # Create clean display table
    events_table <- data.frame(
      EventID = unique(event_df$event_id),
      Type = sapply(unique(event_df$event_id), function(id) {
        # find the event_type for this event ID
        unique(event_df$event_type[event_df$event_id == id])
      }),
      RefTranscripts = sapply(unique(event_df$event_id), function(id) {
        paste(unique(gene_events$refTx[gene_events$eventID == id]), collapse = ", ")
      }),
      AltTranscripts = sapply(unique(event_df$event_id), function(id) {
        paste(unique(gene_events$asTx[gene_events$eventID == id]), collapse = ", ")
      }),
      stringsAsFactors = FALSE
    )
    
    # Source compression functions if not already loaded
    if (!exists("create_intron_markers_for_transcripts")) {
      source("R/coordinate_compression.R")
    }
    
    # Create a mock exons_by_transcript for the superset view
    # This represents the merged exon structure
    mock_exons_by_transcript <- list(
      "superset" = superset
    )
    
    # Create compression map for intron markers
    compression_map <- list(
      compress = compress_position,
      coords = compressed_coords
    )
    
    # Use the proper function to create intron markers at y_position = 1
    intron_markers <- create_intron_markers_for_transcripts(
      mock_exons_by_transcript, 
      compression_map, 
      y_position = 1
    )
    
    p <- ggplot() +
      # Add intron markers
      geom_line(data = intron_markers,
               aes(x = x, y = y, group = group),
               color = "grey50", size = 0.5) +
      
      # Add superset exons
      geom_rect(data = superset_df,
               aes(xmin = start, xmax = end,
                   ymin = 0.7, ymax = 1.3),
               fill = "grey70", color = "black") +
      
      # Add AS event regions
      geom_rect(data = event_df,
               aes(xmin = start, xmax = end,
                   ymin = y_pos - 0.2, ymax = y_pos + 0.2,
                   fill = type),
               color = "black", alpha = 0.8) +
      
      scale_fill_manual(
        values = c(
          "SE" = "#FF0000",
          "RI" = "#00FF00",
          "A3SS" = "#235e24",
          "A5SS" = "#FF00FF",
          "MXE" = "#0000FF"
        ),
        labels = c(
          "SE" = "SE (Skipped Exon)",
          "RI" = "RI (Retained Intron)",
          "A3SS" = "A3SS (Alt 3' Splice Site)",
          "A5SS" = "A5SS (Alt 5' Splice Site)",
          "MXE" = "MXE (Mutually Exclusive Exon)"
        ),
        name = "AS Event Types"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        legend.position = "bottom"
      ) +
      labs(
        title = paste0("Alternative Splicing Events - ", 
                      gene_symbol, " (", gene_id, ")"),
        x = paste0("Genomic Position (chromosome ", chromosome, ")"),
        y = ""
      )
    
    # Use ACTUAL raw GTF coordinates at each exon position
    axis_breaks <- compressed_coords$compressed_start  # Where each exon appears on plot
    axis_labels <- as.character(compressed_coords$original_start)  # Raw GTF coordinates
    
    # Apply axis with actual GTF coordinates
    p <- p + scale_x_continuous(breaks = axis_breaks, labels = axis_labels)
    
    # Convert to plotly with proper event registration
    p_plotly <- ggplotly(p, tooltip = "text", source = "as_plot") %>% 
      layout(hoverlabel = list(bgcolor = "white", font = list(size = 10)))
    
    # Make sure to register the click event
    p_plotly <- event_register(p_plotly, "plotly_click")
    
    return(list(
      plot = p_plotly,
      events = events_table,
      event_regions = event_regions
    ))
  } else {
    p <- ggplot() +
      geom_rect(data = superset_df,
               aes(xmin = start, xmax = end,
                   ymin = 0.7, ymax = 1.3),
               fill = "grey70", color = "black") +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()
      ) +
      labs(
        title = paste0("Gene Structure - ", gene_symbol, " (", gene_id, ")"),
        subtitle = "No alternative splicing events found",
        x = paste0("Genomic Position (chromosome ", chromosome, ")"),
        y = ""
      )
    
    # Use ACTUAL raw GTF coordinates at each exon position
    axis_breaks <- compressed_coords$compressed_start  # Where each exon appears on plot
    axis_labels <- as.character(compressed_coords$original_start)  # Raw GTF coordinates
    
    # Apply axis with actual GTF coordinates
    p <- p + scale_x_continuous(breaks = axis_breaks, labels = axis_labels)
    
    # Convert to plotly
    p_plotly <- ggplotly(p)
    
    return(list(
      plot = p_plotly,
      events = data.frame(
        EventID = character(0),
        Type = character(0),
        RefTranscripts = character(0),
        AltTranscripts = character(0),
        stringsAsFactors = FALSE
      ),
      event_regions = list()
    ))
  }
}

# Function to get peptides for a transcript (used in AS comparison)
get_transcript_peptides_for_comparison <- function(txID, processed_data, protease) {
  
  # Check if this is a novel transcript (any transcript not in the known processed_data)
  is_novel <- FALSE
  
  # Try to find in known data first - try multiple approaches
  known_result <- tryCatch({
    get_transcript_peptides_genomic(txID, processed_data, protease)
  }, error = function(e) {
    NULL
  })
  
  # If get_transcript_peptides_genomic failed, try direct access to processed_data
  if (is.null(known_result) && !is.null(processed_data$original_peptides)) {
    known_result <- tryCatch({
      # Direct lookup in original_peptides
      tx_matches <- processed_data$original_peptides$txID == txID
      if (any(tx_matches)) {
        tx_data <- processed_data$original_peptides[tx_matches, ]
        mapped_ranges_col <- paste0(protease, "Peps_mapped_ranges")
        if (mapped_ranges_col %in% names(tx_data) && !is.null(tx_data[[mapped_ranges_col]][[1]])) {
          tx_data[[mapped_ranges_col]][[1]]  # Return the GRanges object
        } else {
          NULL
        }
      } else {
        NULL
      }
    }, error = function(e) {
      NULL
    })
  }
  
  # Only consider it novel if we truly can't find it in canonical data
  if (is.null(known_result)) {
    # Final check: is this actually a canonical transcript ID pattern?
    is_canonical <- grepl("^ENS[A-Z]*T[0-9]+", txID) # ENSEMBL transcript pattern
    if (is_canonical) {
      cat("WARNING: Canonical transcript", txID, "not found in processed data\n")
      return(NULL)  # Return NULL instead of processing as novel
    }
    is_novel <- TRUE
    cat("DEBUG: Processing novel transcript:", txID, "\n")
  }
  
  if (is_novel) {
    # Novel transcripts - use the 24-column RDS structure with _mapped_ranges
    # Check if novel data is available in the server environment
    if (exists("novel_pipeline_results", envir = .GlobalEnv)) {
      novel_results <- get("novel_pipeline_results", envir = .GlobalEnv)()
      
      if (!is.null(novel_results) && novel_results$success && !is.null(novel_results$dataframe_file) && file.exists(novel_results$dataframe_file)) {
        
        tryCatch({
          # Load the 24-column novel RDS
          novel_data <- readRDS(novel_results$dataframe_file)
          
          # Find the transcript
          tx_row <- which(novel_data$txID == txID)
          
          if (length(tx_row) == 0) {
            cat("DEBUG: Novel transcript not found in RDS:", txID, "\n")
            return(NULL)
          }
          
          tx_data <- novel_data[tx_row[1], ]
          
          # Get mapped_ranges for the protease (contains GRanges with genomic coordinates)
          mapped_ranges_col <- paste0(protease, "Peps_mapped_ranges")
          
          if (!mapped_ranges_col %in% names(tx_data)) {
            cat("DEBUG: Mapped ranges column not found:", mapped_ranges_col, "\n")
            return(NULL)
          }
          
          mapped_ranges <- tx_data[[mapped_ranges_col]][[1]]
          
          if (is.null(mapped_ranges) || length(mapped_ranges) == 0) {
            cat("DEBUG: No mapped ranges found for novel transcript:", txID, "\n")
            return(NULL)
          }
          
          # Combine all GRanges objects into one (same format as known genes)
          if (is.list(mapped_ranges) && length(mapped_ranges) > 0) {
            # Check if first element is a GRanges
            if (inherits(mapped_ranges[[1]], "GRanges")) {
              combined_gr <- do.call(c, mapped_ranges)
              cat("DEBUG: Created combined GRanges with", length(combined_gr), "peptides for novel transcript", txID, "\n")
              return(combined_gr)
            } else {
              cat("DEBUG: Mapped ranges do not contain GRanges objects for:", txID, "\n")
              return(NULL)
            }
          } else {
            cat("DEBUG: Invalid mapped ranges structure for novel transcript:", txID, "\n")
            return(NULL)
          }
          
        }, error = function(e) {
          cat("ERROR processing novel transcript", txID, ":", e$message, "\n")
          return(NULL)
        })
        
      } else {
        cat("DEBUG: Novel pipeline results not available for transcript:", txID, "\n")
      }
    } else {
      cat("DEBUG: novel_pipeline_results not found in global environment\n")
    }
    
    # Check for rMATS transcripts if novel didn't work
    if (exists("rmats_pipeline_results", envir = .GlobalEnv)) {
      rmats_results <- get("rmats_pipeline_results", envir = .GlobalEnv)()
      
      if (!is.null(rmats_results) && rmats_results$success && !is.null(rmats_results$dataframe_file) && file.exists(rmats_results$dataframe_file)) {
        
        tryCatch({
          # Load the 24-column rMATS RDS
          rmats_data <- readRDS(rmats_results$dataframe_file)
          
          # Find the transcript
          tx_row <- which(rmats_data$txID == txID)
          
          if (length(tx_row) == 0) {
            cat("DEBUG: rMATS transcript not found in RDS:", txID, "\n")
            return(NULL)
          }
          
          tx_data <- rmats_data[tx_row[1], ]
          
          # Get mapped_ranges for the protease (contains GRanges with genomic coordinates)
          mapped_ranges_col <- paste0(protease, "Peps_mapped_ranges")
          
          if (!mapped_ranges_col %in% names(tx_data)) {
            cat("DEBUG: Mapped ranges column not found for rMATS:", mapped_ranges_col, "\n")
            return(NULL)
          }
          
          mapped_ranges <- tx_data[[mapped_ranges_col]][[1]]
          
          if (is.null(mapped_ranges) || length(mapped_ranges) == 0) {
            cat("DEBUG: No mapped ranges found for rMATS transcript:", txID, "\n")
            return(NULL)
          }
          
          # Handle GRanges structure (same format as known genes)
          if (is.list(mapped_ranges) && length(mapped_ranges) > 0) {
            # Check if first element is a GRanges
            if (inherits(mapped_ranges[[1]], "GRanges")) {
              # If single GRanges object in list, extract it directly
              if (length(mapped_ranges) == 1) {
                combined_gr <- mapped_ranges[[1]]
              } else {
                # Multiple GRanges objects, combine them
                combined_gr <- do.call(c, mapped_ranges)
              }
              cat("DEBUG: Created combined GRanges with", length(combined_gr), "peptides for rMATS transcript", txID, "\n")
              return(combined_gr)
            } else {
              cat("DEBUG: Mapped ranges do not contain GRanges objects for rMATS:", txID, "\n")
              return(NULL)
            }
          } else if (inherits(mapped_ranges, "GRanges")) {
            # Direct GRanges object
            cat("DEBUG: Using direct GRanges with", length(mapped_ranges), "peptides for rMATS transcript", txID, "\n")
            return(mapped_ranges)
          } else {
            cat("DEBUG: Invalid mapped ranges structure for rMATS transcript:", txID, "\n")
            return(NULL)
          }
          
        }, error = function(e) {
          cat("ERROR processing rMATS transcript", txID, ":", e$message, "\n")
          return(NULL)
        })
        
      } else {
        cat("DEBUG: rMATS pipeline results not available for transcript:", txID, "\n")
      }
    } else {
      cat("DEBUG: rmats_pipeline_results not found in global environment\n")
    }
    
    # If we get here, transcript not found in known, novel, or rMATS data
    return(NULL)
  } else {
    # Known transcript - return the result we already got
    if (!is.null(known_result) && !is.null(known_result$genomic_ranges) && length(known_result$genomic_ranges) > 0) {
      return(known_result$genomic_ranges)
    } else {
      return(NULL)
    }
  }
} 