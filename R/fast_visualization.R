#===============================================================================
# LIGHTNING-FAST VISUALIZATION FUNCTIONS
#===============================================================================

#' Create transcript plot data using pre-computed GTF cache (INSTANT)
#' 
create_fast_transcript_plot_data <- function(gene_id, gene_symbol, as_events = NULL) {
  # Load pre-computed GTF data (0.002 seconds vs 15+ seconds)
  gtf_data <- load_gtf_visualization_data(gene_id)
  
  if (!gtf_data$success) {
    return(list(
      success = FALSE,
      message = gtf_data$message,
      plot = ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = gtf_data$message, 
                size = 5, hjust = 0.5) +
        theme_void() +
        xlim(0, 1) + ylim(0, 1)
    ))
  }
  
  # Extract pre-computed data
  exons_by_transcript <- gtf_data$exons_by_transcript
  cds_by_transcript <- gtf_data$cds_by_transcript
  transcript_ids <- gtf_data$transcript_ids
  
  # Check if we have any transcripts
  if (length(transcript_ids) == 0) {
    return(list(
      success = FALSE,
      message = "No transcripts found for this gene",
      plot = ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = "No transcripts found for this gene", 
                size = 5, hjust = 0.5) +
        theme_void() +
        xlim(0, 1) + ylim(0, 1)
    ))
  }
  
  # Use pre-computed gene boundaries (no calculation needed)
  padding <- 5000
  gene_start <- gtf_data$gene_start - padding
  gene_end <- gtf_data$gene_end + padding
  chromosome <- gtf_data$chromosome
  strand_display <- ""
  
  # Create plot data frames efficiently
  transcript_df <- data.frame(
    transcript = transcript_ids,
    y_position = seq_along(transcript_ids),
    stringsAsFactors = FALSE
  )
  
  # Pre-compute exon data
  exon_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    exon_number = integer(),
    stringsAsFactors = FALSE
  )
  
  # Pre-compute CDS data
  cds_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    cds_number = integer(),
    stringsAsFactors = FALSE
  )
  
  # Process transcripts efficiently
  for (i in seq_along(transcript_ids)) {
    tx <- transcript_ids[i]
    
    # Add exons (already filtered by transcript)
    if (tx %in% names(exons_by_transcript) && length(exons_by_transcript[[tx]]) > 0) {
      tx_exons <- exons_by_transcript[[tx]]
      tx_exons <- tx_exons[order(start(tx_exons))]
      
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
    
    # Add CDS (already filtered by transcript)
    if (tx %in% names(cds_by_transcript) && length(cds_by_transcript[[tx]]) > 0) {
      tx_cds <- cds_by_transcript[[tx]]
      tx_cds <- tx_cds[order(start(tx_cds))]
      
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
  }
  
  # Create the plot (same as before but with pre-computed data)
  max_y_pos <- max(transcript_df$y_position) + 1
  
  p <- ggplot() +
    # Add transcript lines
    geom_segment(data = transcript_df, 
                aes(x = gene_start, xend = gene_end, 
                    y = y_position, yend = y_position),
                linewidth = 0.5, color = "grey70") +
    
    # Add exon blocks
    geom_rect(data = exon_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.3, ymax = y_position + 0.3,
                 fill = "Transcript"),
             color = "black", alpha = 0.8) +
    
    # Add CDS overlay
    geom_rect(data = cds_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.25, ymax = y_position + 0.25,
                 fill = "CDS"),
             color = "black", alpha = 0.9) +
    
    # Add transcript labels
    geom_text(data = transcript_df,
             aes(x = gene_start - 100, y = y_position, label = transcript),
             hjust = 1, size = 3.5)
  
  # Add direction indicator
  direction_df <- data.frame(
    x = c(gene_start + padding, gene_end - padding),
    y = c(max_y_pos, max_y_pos),
    label = c("5'", "3'"),
    stringsAsFactors = FALSE
  )
  
  p <- p + 
    geom_text(data = direction_df[1,], aes(x = x, y = y, label = label), size = 4, fontface = "bold") +
    geom_text(data = direction_df[2,], aes(x = x, y = y, label = label), size = 4, fontface = "bold") +
    geom_segment(data = data.frame(x = gene_start + padding + 50, xend = gene_end - padding - 50, 
                                  y = max_y_pos, yend = max_y_pos),
                aes(x = x, xend = xend, y = y, yend = yend),
                arrow = arrow(length = unit(0.3, "cm"), ends = "last", type = "closed"),
                color = "black", linewidth = 0.7)
  
  # Add fill scale for transcript structure - only show relevant elements
  has_cds <- nrow(cds_df) > 0
  fill_values <- c(
    "Transcript" = "rgba(77, 175, 74, 0.8)",
    "CDS" = "rgba(255, 221, 0, 0.8)"
  )
  
  p <- p + 
    scale_fill_manual(
      values = fill_values,
      breaks = c("Transcript", "CDS"),
      labels = c("Transcript", "CDS"),
      name = ""
    ) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
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
                   ylim = c(0.5, max_y_pos + 0.5))
  
  return(list(
    success = TRUE,
    plot = p,
    chromosome = chromosome,
    has_cds = has_cds
  ))
}

#' Create peptide plot data using pre-computed GTF cache (INSTANT)
#' 
create_fast_peptide_plot_data <- function(gene_id, gene_symbol, processed_data, protease, use_compression = FALSE) {
  # Load pre-computed GTF data (0.002 seconds vs 15+ seconds)
  gtf_data <- load_gtf_visualization_data(gene_id)
  
  if (!gtf_data$success) {
    return(list(
      success = FALSE,
      message = gtf_data$message,
      plot = ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = gtf_data$message, 
                size = 5, hjust = 0.5) +
        theme_void() +
        xlim(0, 1) + ylim(0, 1)
    ))
  }
  
  # Extract pre-computed data
  exons_by_transcript <- gtf_data$exons_by_transcript
  cds_by_transcript <- gtf_data$cds_by_transcript
  transcript_ids <- gtf_data$transcript_ids
  
  # Use pre-computed gene boundaries
  padding <- 5000
  gene_start <- gtf_data$gene_start - padding
  gene_end <- gtf_data$gene_end + padding
  chromosome <- gtf_data$chromosome
  strand_display <- ""
  
  # Create plot data frames (same structure as before)
  transcript_df <- data.frame(
    transcript = transcript_ids,
    y_position = seq_along(transcript_ids),
    stringsAsFactors = FALSE
  )
  
  # Pre-compute all data frames
  exon_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    exon_number = integer(),
    stringsAsFactors = FALSE
  )
  
  cds_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    cds_number = integer(),
    stringsAsFactors = FALSE
  )
  
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
  
  # Process efficiently (same logic as before but with pre-computed data)
  for (i in seq_along(transcript_ids)) {
    tx <- transcript_ids[i]
    
    # Add exons
    if (tx %in% names(exons_by_transcript) && length(exons_by_transcript[[tx]]) > 0) {
      tx_exons <- exons_by_transcript[[tx]]
      tx_exons <- tx_exons[order(start(tx_exons))]
      
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
    
    # Add CDS
    if (tx %in% names(cds_by_transcript) && length(cds_by_transcript[[tx]]) > 0) {
      tx_cds <- cds_by_transcript[[tx]]
      tx_cds <- tx_cds[order(start(tx_cds))]
      
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
    
    # Add peptides
    peptide_data <- get_transcript_peptides_genomic(tx, processed_data, protease)
    if (!is.null(peptide_data) && length(peptide_data$genomic_ranges) > 0) {
      tx_peptides <- peptide_data$genomic_ranges
      
      # Count peptide occurrences
      peptide_counts <- table(tx_peptides$peptide)
      junction_spanning <- names(peptide_counts[peptide_counts > 1])
      
      # Add all peptide ranges
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
          "<br>Transcript: ", tx,
          "<br>Type: ", ifelse(tx_peptides$peptide %in% junction_spanning, "Junction-spanning", "Regular")
        )),
        stringsAsFactors = FALSE
      ))
    }
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
    intron_markers_df <- create_intron_markers_for_transcripts(exons_by_transcript, compression_map, y_position = 0)
    # No intron labels needed - they clutter the visualization
    intron_labels_df <- NULL
    
    # Update gene boundaries for compressed view
    if (nrow(compression_map$coords) > 0) {
      gene_start <- min(compression_map$coords$compressed_start) - 100
      gene_end <- max(compression_map$coords$compressed_end) + 100
    }
    
    # Update axis label to indicate compression
    x_axis_label <- paste0("Genomic Position (compressed view) - chromosome ", chromosome, "")
    
    # Update hover text
    if (nrow(peptide_df) > 0) {
      peptide_df$hover_text <- paste0(
        peptide_df$hover_text,
        "<br><i>Display uses compressed scale</i>"
      )
    }
  }
  
  # Create the plot (same as before)
  p <- ggplot()
  
  # Add intron visualization based on compression mode
  if (use_compression && !is.null(intron_markers_df) && nrow(intron_markers_df) > 0) {
    # Create zigzag data by matching transcripts to their Y positions
    zigzag_data <- data.frame()
    
    for (j in seq_along(transcript_df$y_position)) {
      tx_id <- transcript_df$transcript[j]
      
      # Get zigzag markers for this specific transcript
      tx_markers <- intron_markers_df[intron_markers_df$transcript == tx_id, ]
      
      if (nrow(tx_markers) > 0) {
        # Adjust Y positions for this transcript's row
        tx_markers$y <- tx_markers$y + transcript_df$y_position[j]
        tx_markers$plot_transcript_id <- j  # Track which plot row this is
        zigzag_data <- rbind(zigzag_data, tx_markers)
      }
    }
    
    # Add all zigzag lines at once
    if (nrow(zigzag_data) > 0) {
      p <- p + geom_line(data = zigzag_data,
                        aes(x = x, y = y, group = group),
                        color = "grey50", linewidth = 0.5)
    }
    
    # Intron labels disabled - they clutter the visualization
    # if (!is.null(intron_labels_df) && nrow(intron_labels_df) > 0) {
    #   # Create label data for each transcript
    #   label_data <- data.frame()
    #   for (j in seq_along(transcript_df$y_position)) {
    #     transcript_labels <- intron_labels_df
    #     transcript_labels$y <- transcript_df$y_position[j] - 0.5
    #     label_data <- rbind(label_data, transcript_labels)
    #   }
    #   
    #   if (nrow(label_data) > 0) {
    #     p <- p + geom_text(data = label_data,
    #                       aes(x = x, y = y, label = label),
    #                       size = 2.5, color = "grey40", fontface = "italic")
    #   }
    # }
  } else {
    # Original intron lines for true scale
    p <- p + geom_segment(data = transcript_df, 
                         aes(x = gene_start, xend = gene_end, 
                             y = y_position, yend = y_position),
                         linewidth = 0.5, color = "grey70")
  }
  
  p <- p +
    
    geom_rect(data = exon_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.3, ymax = y_position + 0.3,
                 fill = "Transcript"),
             color = "black", alpha = 0.8) +
    
    geom_rect(data = cds_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.25, ymax = y_position + 0.25,
                 fill = "CDS"),
             color = "black", alpha = 0.9) +
    
    # Add peptides with fill aesthetic for legend consistency
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
    
    geom_text(data = transcript_df,
             aes(x = gene_start - 100, y = y_position, label = transcript),
             hjust = 1, size = 3.5)
  
  # Add direction indicator and styling (same as before)
  direction_df <- data.frame(
    x = c(gene_start + padding, gene_end - padding),
    y = c(max(transcript_df$y_position) + 1, max(transcript_df$y_position) + 1),
    label = c("5'", "3'"),
    stringsAsFactors = FALSE
  )
  
  p <- p + 
    geom_text(data = direction_df[1,], aes(x = x, y = y, label = label), size = 4, fontface = "bold") +
    geom_text(data = direction_df[2,], aes(x = x, y = y, label = label), size = 4, fontface = "bold") +
    geom_segment(data = data.frame(x = gene_start + padding + 50, xend = gene_end - padding - 50, 
                                  y = max(transcript_df$y_position) + 0.5, yend = max(transcript_df$y_position) + 0.5),
                aes(x = x, xend = xend, y = y, yend = yend),
                arrow = arrow(length = unit(0.3, "cm"), ends = "last", type = "closed"),
                color = "black", linewidth = 0.7) +
    
    # Use the exact same unified fill scale as transcript structure plot
    {
      fill_values <- c(
        "Transcript" = "rgba(77, 175, 74, 0.8)",
        "CDS" = "rgba(255, 221, 0, 0.8)",
        "peptide" = "rgba(52, 152, 219, 0.9)",
        "junction spanning peptide" = "rgba(231, 76, 60, 0.9)"
      )
      
      scale_fill_manual(
        values = fill_values,
        breaks = c("Transcript", "CDS", "peptide", "junction spanning peptide"),
        labels = c("Transcript", "CDS", "peptide", "junction spanning peptide"),
        name = ""
      )
    } +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
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
  
  has_peptides <- nrow(peptide_df) > 0
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