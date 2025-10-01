#===============================================================================
# VISUALIZATION FUNCTIONS
#===============================================================================

# Load required libraries for genomic ranges
if (!require(GenomicRanges, quietly = TRUE)) {
  stop("GenomicRanges package is required for visualization functions")
}

# Global hover text cleaning function to remove hex color codes
clean_hover_text <- function(text) {
  if (is.null(text) || length(text) == 0) return(text)
  
  # Check for color codes (silently)
  
  # Handle vectors properly - use vectorized is.na
  if (any(is.na(text))) {
    # Keep NA values as they are, only clean non-NA elements
    result <- text
    non_na_idx <- !is.na(text)
    if (any(non_na_idx)) {
      result[non_na_idx] <- gsub("#[A-Fa-f0-9]{6}\\b", "", text[non_na_idx])
      result[non_na_idx] <- gsub("\\s+", " ", result[non_na_idx])
      result[non_na_idx] <- trimws(result[non_na_idx])
    }
    return(result)
  } else {
    # All elements are non-NA, process normally
    cleaned <- gsub("#[A-Fa-f0-9]{6}\\b", "", text)
    cleaned <- gsub("\\s+", " ", cleaned)
    cleaned <- trimws(cleaned)
    return(cleaned)
  }
}

# Clean plotly object and enhance hover control to prevent color code display
clean_plotly_hover <- function(p) {
  # Extract and clean hover text from all traces
  for (i in seq_along(p$x$data)) {
    # Clean direct text properties
    if (!is.null(p$x$data[[i]]$text)) {
      p$x$data[[i]]$text <- clean_hover_text(p$x$data[[i]]$text)
    }
    if (!is.null(p$x$data[[i]]$hovertext)) {
      p$x$data[[i]]$hovertext <- clean_hover_text(p$x$data[[i]]$hovertext)
    }
    
    # Handle ggplotly-generated traces that may have fillcolor in hover
    if (!is.null(p$x$data[[i]]$fillcolor)) {
      # Check if fillcolor contains rgba values that might leak into hover
      if (is.character(p$x$data[[i]]$fillcolor) && 
          any(grepl("rgba\\(", p$x$data[[i]]$fillcolor))) {
        # Force specific hover behavior to prevent automatic fillcolor display
        if (is.null(p$x$data[[i]]$hovertemplate)) {
          if (!is.null(p$x$data[[i]]$text) || !is.null(p$x$data[[i]]$hovertext)) {
            p$x$data[[i]]$hoverinfo <- "text"
          } else {
            p$x$data[[i]]$hoverinfo <- "skip"
          }
        }
      }
    }
    
    # For traces without explicit hovertemplate, ensure no default hover shows color info
    if (is.null(p$x$data[[i]]$hovertemplate) && 
        (is.null(p$x$data[[i]]$hoverinfo) || p$x$data[[i]]$hoverinfo != "none")) {
      # Set hoverinfo to only show text if text is available, otherwise skip
      if (!is.null(p$x$data[[i]]$text) || !is.null(p$x$data[[i]]$hovertext)) {
        p$x$data[[i]]$hoverinfo <- "text"
      } else {
        p$x$data[[i]]$hoverinfo <- "skip"
      }
    }
    
    # Additional cleanup for ggplotly traces - remove any hover data containing color values
    if (!is.null(p$x$data[[i]]$customdata)) {
      # ggplotly sometimes stores additional data here that can leak into hover
      if (is.character(p$x$data[[i]]$customdata)) {
        p$x$data[[i]]$customdata <- clean_hover_text(p$x$data[[i]]$customdata)
      }
    }
  }
  
  # Add layout-level hover configuration to prevent color code display
  if (is.null(p$x$layout$hovermode) || p$x$layout$hovermode == "auto") {
    p <- p %>% plotly::layout(
      hovermode = "closest",
      hoverlabel = list(
        bgcolor = "rgba(255, 255, 255, 0.9)",
        bordercolor = "rgba(0, 0, 0, 0.2)",
        font = list(size = 12, color = "black")
      )
    )
  }
  
  return(p)
}

# Helper function to create optimized plotly plots with legend interactivity
create_optimized_plotly <- function(ggplot_obj, tooltip = "text", source = NULL) {
  # Use 'all' tooltip to preserve fill aesthetics for legend interactivity, then clean hover manually
  plotly_obj <- ggplotly(ggplot_obj, tooltip = "all", source = source) %>%
    layout(
      hoverlabel = list(bgcolor = "white", font = list(size = 10)),
      # Performance optimizations
      hovermode = "closest",
      showlegend = TRUE
    ) %>%
    config(
      displayModeBar = TRUE,
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("sendDataToCloud", "editInChartStudio", "lasso2d", "select2d"),
      # Performance settings
      doubleClick = "reset",
      scrollZoom = TRUE
    )
  
  # Clean hover text while preserving legend interactivity
  for (i in seq_along(plotly_obj$x$data)) {
    # Only clean hover text if there is custom text available
    if (!is.null(plotly_obj$x$data[[i]]$text)) {
      # Use custom text for hover, preserving fill for legend functionality
      plotly_obj$x$data[[i]]$hovertemplate <- paste0(
        clean_hover_text(plotly_obj$x$data[[i]]$text), 
        "<extra></extra>"
      )
    } else {
      # For traces without custom text, clean any automatic hover
      if (!is.null(plotly_obj$x$data[[i]]$hovertext)) {
        plotly_obj$x$data[[i]]$hovertext <- clean_hover_text(plotly_obj$x$data[[i]]$hovertext)
      }
    }
  }
  
  # Clean hover text but preserve legend functionality
  plotly_obj <- clean_plotly_hover(plotly_obj)
  
  # Register events if a source is provided (following pattern from as_analysis.R)
  if (!is.null(source)) {
    plotly_obj <- event_register(plotly_obj, "plotly_click")
  }
  
  return(plotly_obj)
}

# Function to create transcript structure plot data
create_transcript_plot_data <- function(exons_result, transcript_ids, gene_symbol, gene_id, as_events = NULL) {
  if (!exons_result$success) {
    return(list(
      success = FALSE,
      message = exons_result$message,
      plot = ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = exons_result$message, 
                size = 5, hjust = 0.5) +
        theme_void() +
        xlim(0, 1) + ylim(0, 1)
    ))
  }
  
  exons_by_transcript <- exons_result$exons
  cds_by_transcript <- exons_result$cds
  
  # Check if we have any exons to plot
  if (length(exons_by_transcript) == 0) {
    return(list(
      success = FALSE,
      message = "No exons found for any transcripts",
      plot = ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = "No exons found for any transcripts", 
                size = 5, hjust = 0.5) +
        theme_void() +
        xlim(0, 1) + ylim(0, 1)
    ))
  }
  
  # Collect all exons to calculate gene boundaries
  all_exons <- list()
  for (tx in names(exons_by_transcript)) {
    if (length(exons_by_transcript[[tx]]) > 0) {
      all_exons <- c(all_exons, list(exons_by_transcript[[tx]]))
    }
  }
  
  # Calculate gene boundaries with padding
  all_ranges <- do.call(c, all_exons)
  padding <- 5000  # Increased padding for better visualization
  gene_start <- min(start(all_ranges)) - padding
  gene_end <- max(end(all_ranges)) + padding
  chromosome <- as.character(seqnames(all_ranges)[1])
  
  # Determine strand (all exons should be on the same strand)
  strand_char <- as.character(strand(all_ranges)[1])
  if (strand_char == "+" || strand_char == "-") {
    strand_direction <- strand_char
  } else {
    strand_direction <- "unknown"
  }
  
  # Create strand direction display (plots are always shown 5' to 3')
  strand_display <- ""
  
  # Create data frames for plotting
  # One for transcript backbones
  transcript_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    stringsAsFactors = FALSE
  )
  
  # One for exons
  exon_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    exon_number = integer(),
    stringsAsFactors = FALSE
  )
  
  # One for CDS regions
  cds_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    cds_number = integer(),
    stringsAsFactors = FALSE
  )
  
  # Add new data frame for start codons
  start_codon_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    position = numeric(),
    strand = character(),
    stringsAsFactors = FALSE
  )
  
  # Fill data frames with transcript and exon data
  for (i in seq_along(transcript_ids)) {
    tx <- transcript_ids[i]
    
    # Add transcript to backbone data frame
    transcript_df <- rbind(transcript_df, data.frame(
      transcript = tx,
      y_position = i,
      stringsAsFactors = FALSE
    ))
    
    # Add exons if available
    if (!is.null(exons_by_transcript[[tx]]) && length(exons_by_transcript[[tx]]) > 0) {
      tx_exons <- exons_by_transcript[[tx]]
      for (j in seq_along(tx_exons)) {
        exon_df <- rbind(exon_df, data.frame(
          transcript = tx,
          y_position = i,
          start = start(tx_exons[j]),
          end = end(tx_exons[j]),
          exon_number = j,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add CDS regions if available
    if (!is.null(cds_by_transcript[[tx]]) && length(cds_by_transcript[[tx]]) > 0) {
      tx_cds <- cds_by_transcript[[tx]]
      for (j in seq_along(tx_cds)) {
        cds_df <- rbind(cds_df, data.frame(
          transcript = tx,
          y_position = i,
          start = start(tx_cds[j]),
          end = end(tx_cds[j]),
          cds_number = j,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add start codon if CDS is available
    if (!is.null(cds_by_transcript[[tx]]) && length(cds_by_transcript[[tx]]) > 0) {
      tx_cds <- cds_by_transcript[[tx]]
      
      # Determine strand
      tx_strand <- as.character(strand(tx_cds[1]))
      
      # Find start codon position based on strand
      if (tx_strand == "+") {
        # On positive strand, start codon is at beginning of first CDS
        tx_cds_sorted <- tx_cds[order(start(tx_cds))]
        start_pos <- start(tx_cds_sorted[1])
      } else {
        # On negative strand, start codon is at end of last CDS
        tx_cds_sorted <- tx_cds[order(end(tx_cds), decreasing = TRUE)]
        start_pos <- end(tx_cds_sorted[1])
      }
      
      # Add to start codon data frame
      start_codon_df <- rbind(start_codon_df, data.frame(
        transcript = tx,
        y_position = i,
        position = start_pos,
        strand = tx_strand,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Number of transcripts to determine top and bottom of plot
  num_transcripts <- nrow(transcript_df)
  max_y_pos <- max(transcript_df$y_position) + 1  # Add space for direction indicator
  
  # Create the plot
  p <- ggplot() +
    # Add transcript lines (introns)
    geom_segment(data = transcript_df, 
                aes(x = gene_start, xend = gene_end, 
                    y = y_position, yend = y_position),
                linewidth = 0.5, color = "grey70") +
    
    # Add exon blocks with green color
    geom_rect(data = exon_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.3, ymax = y_position + 0.3,
                 fill = "Transcript"),
             color = "black", alpha = 0.8) +
    
    # Add CDS overlay with yellow color
    geom_rect(data = cds_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.25, ymax = y_position + 0.25,
                 fill = "CDS"),
             color = "black", alpha = 0.9) +
    
    # Add transcript labels
    geom_text(data = transcript_df,
             aes(x = gene_start - 100, y = y_position, label = transcript),
             hjust = 1, size = 3.5) +
    
    # Add start codon markers
    geom_point(data = start_codon_df,
              aes(x = position, y = y_position),
              color = "green3", size = 3, shape = 24) +  # Triangle shape
    
    # Add start codon labels
    geom_text(data = start_codon_df,
             aes(x = position, y = y_position, 
                 label = "Start",
                 hjust = ifelse(strand == "+", -0.3, 1.3)),  # Adjust label position based on strand
             size = 3, color = "#8b0002", fontface = "bold")
  
  # Add 5' to 3' direction indicator at top
  # Create data frame for direction arrow
  direction_df <- data.frame(
    x = c(gene_start + padding, gene_end - padding),
    y = c(max_y_pos, max_y_pos),
    label = c("5'", "3'"),
    stringsAsFactors = FALSE
  )
  
  # Add 5' label
  p <- p + geom_text(
    data = direction_df[1,],
    aes(x = x, y = y, label = label),
    size = 4, fontface = "bold"
  )
  
  # Add 3' label
  p <- p + geom_text(
    data = direction_df[2,],
    aes(x = x, y = y, label = label),
    size = 4, fontface = "bold"
  )
  
  # Add direction arrow
  p <- p + geom_segment(
    data = data.frame(x = gene_start + padding + 50, xend = gene_end - padding - 50, 
                    y = max_y_pos, yend = max_y_pos),
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.3, "cm"), ends = "last", type = "closed"),
    color = "black", linewidth = 0.7
  )
  
  # Add AS event visualization
  as_df <- data.frame(
    start = numeric(),
    end = numeric(),
    y_position = numeric(),
    event_type = character(),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(as_events)) {
    gene_events <- as_events[as_events$geneID == gene_id,]
    if (nrow(gene_events) > 0) {
      for(i in seq_len(nrow(gene_events))) {
        event <- gene_events[i,]
        event_type <- event$AS_type
        
        # Get relevant ranges based on AS type
        ranges <- switch(event_type,
          "SE" = c(event$AS_range, event$AS_range1, event$AS_range2),
          "RI" = event$AS_range,
          "MX" = c(event$AS_range1, event$AS_range2),
          "A3" = event$AS_range,
          "A5" = event$AS_range,
          NULL
        )
        
        # Add valid ranges to dataframe
        for(rng in ranges) {
          if(length(rng) > 0 && !is.null(rng)) {
            as_df <- rbind(as_df, data.frame(
              start = start(rng),
              end = end(rng),
              y_position = max(transcript_df$y_position) + 1 + (i*0.2),
              event_type = event_type,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  # Unified fill values for consistent legends across all plots
  fill_values <- c(
    "Transcript" = "rgba(77, 175, 74, 0.8)",
    "CDS" = "rgba(255, 221, 0, 0.8)",
    "peptide" = "rgba(52, 152, 219, 0.9)",
    "junction spanning peptide" = "rgba(231, 76, 60, 0.9)",
    "SE" = "rgba(255, 0, 0, 0.8)",
    "RI" = "rgba(0, 255, 0, 0.8)",
    "MX" = "rgba(0, 0, 255, 0.8)",
    "A3" = "rgba(35, 94, 36, 0.8)",
    "A5" = "rgba(255, 0, 255, 0.8)"
  )

  # Create combined fill scale - show only relevant elements for transcript structure
  p <- p + 
    scale_fill_manual(
      values = fill_values,
      breaks = c("Transcript", "CDS", "peptide", "junction spanning peptide"),
      labels = c("Transcript", "CDS", "peptide", "junction spanning peptide"),
      name = ""
    )
  
  # Has CDS regions?
  has_cds <- nrow(cds_df) > 0
  
  # Complete the styling with manual fill scale for transcript and CDS
  p <- p + 
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(6, 6, 6, 6)
    ) +
    labs(
      x = paste0("Genomic Position (chromosome ", chromosome, ")"), 
      title = paste0("Transcript Structure - ", gene_symbol, " (", gene_id, ")"),
      subtitle = paste0("Chromosome ", chromosome, ": ", gene_start, " - ", gene_end)
    ) +
    coord_cartesian(xlim = c(gene_start, gene_end), 
                   ylim = c(0.5, max_y_pos + 0.5))  # Set y limits to show direction indicator
  
  # Always show consistent legend elements
  # (Peptide elements won't be visible in transcript plot but will be in legend for consistency)
  
  return(list(
    success = TRUE,
    plot = p,
    chromosome = chromosome,
    has_cds = has_cds
  ))
}

# Function to get peptides for a transcript with genomic positions
get_transcript_peptides_genomic <- function(txID, processed_data, protease = "trp") {
  # Check if processed_data exists and has required structure
  if (is.null(processed_data)) {
    return(NULL)
  }
  
  # Check for different data structures with better error handling
  tryCatch({
    if (!is.null(processed_data$original_peptides)) {
      # Standard structure
      original_peptides <- processed_data$original_peptides
      if (!is.null(original_peptides) && "txID" %in% names(original_peptides)) {
        tx_matches <- original_peptides$txID == txID
        if (length(tx_matches) > 0 && any(tx_matches)) {
          tx_data <- original_peptides[tx_matches, ]
        } else {
          tx_data <- NULL
        }
      } else {
        tx_data <- NULL
      }
    } else if (!is.null(processed_data$peptides) && !is.null(names(processed_data$peptides)) && txID %in% names(processed_data$peptides)) {
      # Alternative structure - direct isoform access
      isoform_data <- processed_data$peptides[[txID]]
      if (!is.null(isoform_data)) {
        # Try to create a compatible structure
        peptide_col <- paste0(protease, "Peps")
        mapped_ranges_col <- paste0(protease, "Peps_mapped_ranges")
        
        if (!is.null(names(isoform_data)) && peptide_col %in% names(isoform_data) && mapped_ranges_col %in% names(isoform_data)) {
          peptides <- isoform_data[[peptide_col]][[1]]
          mapped_ranges <- isoform_data[[mapped_ranges_col]][[1]]
          return(list(
            peptides = peptides,
            genomic_ranges = mapped_ranges
          ))
        }
      }
      return(NULL)
    } else {
      return(NULL)
    }
  }, error = function(e) {
    cat("Error in get_transcript_peptides_genomic:", e$message, "\n")
    return(NULL)
  })
  
  # Safely check nrow for standard structure
  if (is.null(tx_data) || length(tx_data) == 0) {
    return(NULL)
  }
  
  tryCatch({
    if (nrow(tx_data) == 0) {
      return(NULL)
    }
  }, error = function(e) {
    cat("Error checking nrow:", e$message, "\n")
    return(NULL)
  })
  
  # Extract peptides
  peptide_col <- paste0(protease, "Peps")
  peptides <- tx_data[[peptide_col]][[1]]
  
  # Extract mapped ranges
  mapped_ranges_col <- paste0(protease, "Peps_mapped_ranges")
  
  # Check if mapped ranges exist
  if (!mapped_ranges_col %in% names(tx_data) || is.null(tx_data[[mapped_ranges_col]][[1]])) {
    return(NULL)
  }
  
  mapped_ranges <- tx_data[[mapped_ranges_col]][[1]]
  return(list(
    peptides = peptides,
    genomic_ranges = mapped_ranges
  ))
}

# Function to create peptide visualization plot data
create_peptide_plot_data <- function(exons_result, transcript_ids, gene_symbol, gene_id, processed_data, protease, use_compression = FALSE) {
  
  
  if (!exons_result$success) {
    return(list(
      success = FALSE,
      message = exons_result$message,
      plot = ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = exons_result$message, 
                size = 5, hjust = 0.5) +
        theme_void() +
        xlim(0, 1) + ylim(0, 1)
    ))
  }
  
  exons_by_transcript <- exons_result$exons
  cds_by_transcript <- exons_result$cds
  
  # Check if we have any exons to plot
  if (length(exons_by_transcript) == 0) {
    return(list(
      success = FALSE,
      message = "No exons found for any transcripts",
      plot = ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = "No exons found for any transcripts", 
                size = 5, hjust = 0.5) +
        theme_void() +
        xlim(0, 1) + ylim(0, 1)
    ))
  }
  
  # Collect all exons to calculate gene boundaries
  all_exons <- list()
  for (tx in names(exons_by_transcript)) {
    if (length(exons_by_transcript[[tx]]) > 0) {
      all_exons <- c(all_exons, list(exons_by_transcript[[tx]]))
    }
  }
  
  # Calculate gene boundaries with padding
  all_ranges <- do.call(c, all_exons)
  padding <- 5000  # Increased padding for better visualization
  gene_start <- min(start(all_ranges)) - padding
  gene_end <- max(end(all_ranges)) + padding
  chromosome <- as.character(seqnames(all_ranges)[1])
  
  # Determine strand (all exons should be on the same strand)
  strand_char <- as.character(strand(all_ranges)[1])
  if (strand_char == "+" || strand_char == "-") {
    strand_direction <- strand_char
  } else {
    strand_direction <- "unknown"
  }
  
  # Create strand direction display (plots are always shown 5' to 3')
  strand_display <- ""
  
  # Create data frames for plotting
  # One for transcript backbones
  transcript_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    stringsAsFactors = FALSE
  )
  
  # One for exons
  exon_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    exon_number = integer(),
    stringsAsFactors = FALSE
  )
  
  # One for CDS regions
  cds_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    cds_number = integer(),
    stringsAsFactors = FALSE
  )
  
  # One for peptides with hover text
  peptide_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    peptide = character(),
    is_junction_spanning = logical(),
    hover_text = character(),
    stringsAsFactors = FALSE
  )
  
  # Fill data frames with transcript and exon data
  for (i in seq_along(transcript_ids)) {
    tx <- transcript_ids[i]
    
    # Add transcript to backbone data frame
    transcript_df <- rbind(transcript_df, data.frame(
      transcript = tx,
      y_position = i,
      stringsAsFactors = FALSE
    ))
    
    # Add exons if available
    if (!is.null(exons_by_transcript[[tx]]) && length(exons_by_transcript[[tx]]) > 0) {
      tx_exons <- exons_by_transcript[[tx]]
      for (j in seq_along(tx_exons)) {
        exon_df <- rbind(exon_df, data.frame(
          transcript = tx,
          y_position = i,
          start = start(tx_exons[j]),
          end = end(tx_exons[j]),
          exon_number = j,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add CDS regions if available
    if (!is.null(cds_by_transcript[[tx]]) && length(cds_by_transcript[[tx]]) > 0) {
      tx_cds <- cds_by_transcript[[tx]]
      for (j in seq_along(tx_cds)) {
        cds_df <- rbind(cds_df, data.frame(
          transcript = tx,
          y_position = i,
          start = start(tx_cds[j]),
          end = end(tx_cds[j]),
          cds_number = j,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add peptides if available
    peptide_data <- get_transcript_peptides_genomic(tx, processed_data, protease)
    if (!is.null(peptide_data) && length(peptide_data$genomic_ranges) > 0) {
      tx_peptides <- peptide_data$genomic_ranges
      
      # Count occurrences of each peptide sequence
      peptide_counts <- table(tx_peptides$peptide)
      junction_spanning <- names(peptide_counts[peptide_counts > 1])
      
      # Add all peptide ranges with their metadata and hover text
      peptide_df <- rbind(peptide_df, data.frame(
        transcript = tx,
        y_position = i,
        start = start(tx_peptides),
        end = end(tx_peptides),
        peptide = tx_peptides$peptide,
        is_junction_spanning = tx_peptides$peptide %in% junction_spanning,
        hover_text = clean_hover_text(paste0(
          "Peptide: ", tx_peptides$peptide,
          "<br>Position: ", start(tx_peptides), "-", end(tx_peptides),
          "<br>Length: ", width(tx_peptides), " bp",
          "<br>Type: ", ifelse(tx_peptides$peptide %in% junction_spanning, 
                              "Junction-spanning", "Regular")
        )),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Add hover text information to peptides
  if (nrow(peptide_df) == 0) {
    # Create empty hover_text column for consistency
    peptide_df$hover_text <- character(0)
  } else {
    peptide_df$hover_text <- paste0(
      "Peptide: ", peptide_df$peptide,
      "<br>Position: ", peptide_df$start, "-", peptide_df$end,
      "<br>Transcript: ", peptide_df$transcript,
      "<br>Type: ", ifelse(peptide_df$is_junction_spanning, "Junction-spanning", "Regular"),
      "<br>AS-affected: ", ifelse(peptide_df$is_junction_spanning, "Yes", "No")
    )
  }
  
  # Apply compression if requested
  intron_markers_df <- NULL
  intron_labels_df <- NULL
  x_axis_label <- paste0("Genomic Position (chromosome ", chromosome, ")")
  
  # Store original gene boundaries before compression
  original_gene_start <- gene_start
  original_gene_end <- gene_end
  
  if (use_compression) {
    # Source compression functions if not already loaded
    if (!exists("create_compression_map")) {
      source("R/coordinate_compression.R")
    }
    
    # Create compression map
    compression_map <- create_compression_map(exons_by_transcript)
    
    # Compress all coordinate data frames
    if (nrow(exon_df) > 0) {
      exon_df <- compress_coordinates_df(exon_df, compression_map)
    }
    if (nrow(cds_df) > 0) {
      cds_df <- compress_coordinates_df(cds_df, compression_map)
    }
    if (nrow(peptide_df) > 0) {
      peptide_df <- compress_coordinates_df(peptide_df, compression_map)
    }
    
    # Create intron markers for visualization between actual transcript exons
    # Note: We'll adjust Y positions per transcript later
    intron_markers_df <- create_intron_markers_for_transcripts(exons_by_transcript, compression_map, y_position = 0)
    
    # No intron labels needed - they clutter the visualization
    intron_labels_df <- NULL
    
    # Update gene boundaries for compressed view
    if (nrow(compression_map$coords) > 0) {
      gene_start <- min(compression_map$coords$compressed_start) - 100
      gene_end <- max(compression_map$coords$compressed_end) + 100
    }
    
    # Update axis label to indicate compression
    x_axis_label <- paste0("Genomic Position (compressed view) - chromosome ", chromosome, ")")
    
    # Update hover text to show original positions
    if (nrow(peptide_df) > 0) {
      peptide_df$hover_text <- paste0(
        peptide_df$hover_text,
        "<br><i>Display uses compressed scale</i>"
      )
    }
  }
  
  # Create the plot
  p <- ggplot()
  
  # Add intron lines - simple straight lines like in the overview visualization
  p <- p + geom_segment(data = transcript_df, 
                       aes(x = gene_start, xend = gene_end, 
                           y = y_position, yend = y_position),
                       linewidth = 0.5, color = "grey70")
  
  p <- p +
    # Add exon blocks with green color for transcript structure (matching transcript plot)
    geom_rect(data = exon_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.3, ymax = y_position + 0.3,
                 fill = "Transcript"),
             color = "black", alpha = 0.8) +
    
    # Add CDS overlay with yellow color using fill aesthetic for legend consistency
    geom_rect(data = cds_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.25, ymax = y_position + 0.25,
                 fill = "CDS"),
             color = "black", alpha = 0.9) +
    
    # Add peptides with hover information using fill aesthetic for legend consistency
    # Only add peptide layers if we have peptide data to avoid ggplotly conversion errors
    {if (nrow(peptide_df) > 0 && any(!peptide_df$is_junction_spanning)) {
      geom_rect(data = peptide_df[!peptide_df$is_junction_spanning, ],
               aes(xmin = start, xmax = end, 
                   ymin = y_position - 0.15, ymax = y_position + 0.15,
                   text = hover_text, fill = "peptide"),
               color = "black", alpha = 0.9)
    } else {
      geom_blank()
    }} +
    {if (nrow(peptide_df) > 0 && any(peptide_df$is_junction_spanning)) {
      geom_rect(data = peptide_df[peptide_df$is_junction_spanning, ],
               aes(xmin = start, xmax = end, 
                   ymin = y_position - 0.15, ymax = y_position + 0.15,
                   text = hover_text, fill = "junction spanning peptide"),
               color = "black", alpha = 0.9)
    } else {
      geom_blank()
    }} +
    
    # Add transcript labels
    geom_text(data = transcript_df,
             aes(x = gene_start - 100, y = y_position, label = transcript),
             hjust = 1, size = 3.5)
  
  # Add 5' to 3' direction indicator at top
  # Create data frame for direction arrow
  direction_df <- data.frame(
    x = c(gene_start + padding, gene_end - padding),
    y = c(max(transcript_df$y_position) + 1, max(transcript_df$y_position) + 1),
    label = c("5'", "3'"),
    stringsAsFactors = FALSE
  )
  
  # Add 5' label
  p <- p + geom_text(
    data = direction_df[1,],
    aes(x = x, y = y, label = label),
    size = 4, fontface = "bold"
  )
  
  # Add 3' label
  p <- p + geom_text(
    data = direction_df[2,],
    aes(x = x, y = y, label = label),
    size = 4, fontface = "bold"
  )
  
  # Add direction arrow
  p <- p + geom_segment(
    data = data.frame(x = gene_start + padding + 50, xend = gene_end - padding - 50, 
                    y = max(transcript_df$y_position) + 0.5, yend = max(transcript_df$y_position) + 0.5),
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.3, "cm"), ends = "last", type = "closed"),
    color = "black", linewidth = 0.7
  )
  
  # Has peptides?
  has_peptides <- nrow(peptide_df) > 0
  
  # Use the exact same unified fill scale as transcript structure plot
  fill_values <- c(
    "Transcript" = "rgba(77, 175, 74, 0.8)",
    "CDS" = "rgba(255, 221, 0, 0.8)",
    "peptide" = "rgba(52, 152, 219, 0.9)",
    "junction spanning peptide" = "rgba(231, 76, 60, 0.9)",
    "SE" = "rgba(255, 0, 0, 0.8)",
    "RI" = "rgba(0, 255, 0, 0.8)",
    "MX" = "rgba(0, 0, 255, 0.8)",
    "A3" = "rgba(35, 94, 36, 0.8)",
    "A5" = "rgba(255, 0, 255, 0.8)"
  )

  p <- p + 
    scale_fill_manual(
      values = fill_values,
      breaks = c("Transcript", "CDS", "peptide", "junction spanning peptide"),
      labels = c("Transcript", "CDS", "peptide", "junction spanning peptide"),
      name = ""
    )
  
  p <- p +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(6, 6, 6, 6)
    ) +
    labs(
      x = x_axis_label, 
      title = paste0("Peptide Visualization - ", gene_symbol, " (", gene_id, ")"),
      subtitle = paste0("Chromosome ", chromosome, ": ", 
                       ifelse(use_compression, "Compressed view", 
                              paste0(gene_start, " - ", gene_end)), 
                       " | Enzyme: ", protease)
    )
  
  # Apply coordinate system based on compression mode
  if (use_compression && exists("compression_map")) {
    # Use ACTUAL raw GTF coordinates at each exon position
    axis_breaks <- compression_map$coords$compressed_start  # Where each exon appears on plot
    axis_labels <- as.character(compression_map$coords$original_start)  # Raw GTF coordinates
    
    # Apply axis with actual GTF coordinates
    p <- p + 
      scale_x_continuous(
        breaks = axis_breaks,
        labels = axis_labels
      ) +
      coord_cartesian(ylim = c(0.5, max(transcript_df$y_position) + 0.5))
  } else {
    # Use original coordinate system for true scale
    p <- p + coord_cartesian(xlim = c(gene_start, gene_end), 
                            ylim = c(0.5, max(transcript_df$y_position) + 0.5))
  }
  
  # If no peptides, show a message
  if (!has_peptides) {
    # Calculate center position based on compression mode
    center_x <- (gene_start + gene_end) / 2
    
    p <- p + annotate("text", x = center_x, y = max(transcript_df$y_position)/2, 
                     label = paste0("No mapped peptides found for enzyme: ", protease),
                     size = 5, color = "red")
  }
  
  return(list(
    success = TRUE,
    plot = p,
    chromosome = chromosome,
    has_peptides = has_peptides
  ))
} 