#===============================================================================
# VISUALIZATION FUNCTIONS
#===============================================================================

# Helper function to create optimized plotly plots
create_optimized_plotly <- function(ggplot_obj, tooltip = "text", source = NULL) {
  plotly_obj <- ggplotly(ggplot_obj, tooltip = tooltip, source = source) %>%
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
  
  # Combine all fill values including AS events
  fill_values <- c(
    "Transcript" = "#4DAF4A",
    "CDS" = "#FFDD00",
    "SE" = "#FF0000",
    "RI" = "#00FF00",
    "MX" = "#0000FF",
    "A3" = "#235e24",
    "A5" = "#FF00FF"
  )

  # Create combined fill scale
  p <- p + 
    scale_fill_manual(
      values = fill_values,
      breaks = c("Transcript", "CDS", "SE", "RI", "MX", "A3", "A5"),
      labels = c("Transcript", "CDS", "Skipped Exon", "Retained Intron", 
                "Mutually Exclusive", "Alt 3'", "Alt 5'"),
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
  
  # If we don't have CDS, adjust the legend
  if (!has_cds) {
    p <- p + scale_fill_manual(values = c("Transcript" = "#4DAF4A"), name = "")
  }
  
  return(list(
    success = TRUE,
    plot = p,
    chromosome = chromosome,
    has_cds = has_cds
  ))
}

# Function to get peptides for a transcript with genomic positions
get_transcript_peptides_genomic <- function(txID, processed_data, protease = "trp") {
  # Find the row in peptides for this transcript
  tx_data <- processed_data$original_peptides[processed_data$original_peptides$txID == txID, ]
  
  if (nrow(tx_data) == 0) {
    return(NULL)
  }
  
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
create_peptide_plot_data <- function(exons_result, transcript_ids, gene_symbol, gene_id, processed_data, protease) {
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
        hover_text = paste0(
          "Peptide: ", tx_peptides$peptide,
          "<br>Position: ", start(tx_peptides), "-", end(tx_peptides),
          "<br>Length: ", width(tx_peptides), " bp",
          "<br>Type: ", ifelse(tx_peptides$peptide %in% junction_spanning, 
                              "Junction-spanning", "Regular")
        ),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Add hover text information to peptides
  peptide_df$hover_text <- paste0(
    "Peptide: ", peptide_df$peptide,
    "<br>Position: ", peptide_df$start, "-", peptide_df$end,
    "<br>Transcript: ", peptide_df$transcript,
    "<br>Type: ", ifelse(peptide_df$is_junction_spanning, "Junction-spanning", "Regular"),
    "<br>AS-affected: ", ifelse(peptide_df$is_junction_spanning, "Yes", "No")
  )
  
  # Create the plot
  p <- ggplot() +
    # Add transcript lines (introns)
    geom_segment(data = transcript_df, 
                aes(x = gene_start, xend = gene_end, 
                    y = y_position, yend = y_position),
                linewidth = 0.5, color = "grey70") +
    
    # Add exon blocks with light gray
    geom_rect(data = exon_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.3, ymax = y_position + 0.3),
             fill = "#CCCCCC", color = "black", alpha = 0.5) +
    
    # Add CDS overlay with yellow color
    geom_rect(data = cds_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.25, ymax = y_position + 0.25,
                 fill = "CDS"),
             color = "black", alpha = 0.7) +
    
    # Add peptides with hover information
    geom_rect(data = peptide_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.15, ymax = y_position + 0.15,
                 fill = factor(ifelse(is_junction_spanning, "Junction-spanning", "Regular")),
                 text = hover_text),
             color = "black", alpha = 0.9) +
    
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
  
  # Complete the styling with manual fill scale for CDS and peptides
  p <- p + 
    scale_fill_manual(
      values = c(
        "CDS" = "#F1C40F",         # Brighter yellow for CDS
        "Regular" = "#3498DB",      # Bright blue for regular peptides
        "Junction-spanning" = "#E74C3C"  # Bright red for junction-spanning
      ),
      name = "",
      labels = c(
        "CDS" = "CDS",
        "Regular" = paste0(protease, " peptides"),
        "Junction-spanning" = "Junction-spanning peptides"
      )
    ) +
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
      title = paste0("Peptide Visualization - ", gene_symbol, " (", gene_id, ")"),
      subtitle = paste0("Chromosome ", chromosome, ": ", gene_start, " - ", gene_end, 
                       " | Enzyme: ", protease)
    ) +
    coord_cartesian(xlim = c(gene_start, gene_end), 
                   ylim = c(0.5, max(transcript_df$y_position) + 0.5))
  
  # If no peptides, show a message
  if (!has_peptides) {
    p <- p + annotate("text", x = (gene_start + gene_end)/2, y = max(transcript_df$y_position)/2, 
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