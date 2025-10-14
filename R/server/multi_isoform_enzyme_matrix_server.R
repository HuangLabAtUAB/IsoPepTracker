#===============================================================================
# MULTI-ISOFORM MULTI-ENZYME MATRIX ANALYSIS SERVER
# Combines patterns from multi_isoform_comparison_module.R and multi_enzyme_server.R
#===============================================================================

# Helper function for empty plotly messages (from multi_isoform_comparison_module.R)
empty_plotly_message <- function(message) {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = message, size = 5, hjust = 0.5) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "#f8f9fa", color = NA),
      plot.background = element_rect(fill = "#f8f9fa", color = NA)
    )
  
  ggplotly(p) %>%
    config(displayModeBar = FALSE) %>%
    layout(
      xaxis = list(visible = FALSE),
      yaxis = list(visible = FALSE)
    )
}

# Update matrix isoform choices when isoform analysis is loaded
observeEvent(input$load_all_isoforms, {
  req(selected_gene_transcripts())
  
  transcripts <- selected_gene_transcripts()
  if (length(transcripts) > 0) {
    transcript_choices <- setNames(transcripts, transcripts)
    updateSelectizeInput(session, "matrix_isoforms", 
                        choices = transcript_choices,
                        selected = character(0))  # Start with no selection
  }
})

# Main matrix analysis data processing reactive
matrix_analysis_data <- reactive({
  req(input$run_matrix_analysis > 0, input$matrix_isoforms, gene_data())
  
  selected_isoforms <- input$matrix_isoforms
  miscleavage_type <- input$matrix_miscleavage
  
  # Validate isoform selection
  if (length(selected_isoforms) < 2) {
    showNotification("Please select 2 or more isoforms for matrix analysis", type = "warning")
    return(NULL)
  }
  
  # Get selected enzymes following multi_enzyme_server.R pattern
  selected_enzymes <- c()
  if (input$matrix_trp) selected_enzymes <- c(selected_enzymes, "trp")
  if (input$matrix_chymo) selected_enzymes <- c(selected_enzymes, "chymo")
  if (input$matrix_aspn) selected_enzymes <- c(selected_enzymes, "aspn")
  if (input$matrix_lysc) selected_enzymes <- c(selected_enzymes, "lysc")
  if (input$matrix_lysn) selected_enzymes <- c(selected_enzymes, "lysn")
  if (input$matrix_gluc) selected_enzymes <- c(selected_enzymes, "gluc")
  
  # Validate enzyme selection
  if (length(selected_enzymes) < 2) {
    showNotification("Please select at least 2 enzymes for matrix analysis", type = "warning")
    return(NULL)
  }
  
  withProgress(message = 'Generating matrix analysis...', value = 0, {
    # Load gene data with correct miscleavage type if needed (from multi_enzyme_server.R)
    current_gene_data <- gene_data()
    if (current_gene_data$miscleavage_type != miscleavage_type) {
      incProgress(0.1, detail = 'Loading gene data with correct miscleavage settings...')
      success <- load_and_cache_gene_data(input$gene, miscleavage_type)
      if (!success) {
        showNotification("Failed to load gene data", type = "error")
        return(NULL)
      }
      current_gene_data <- gene_data()
    }
    
    incProgress(0.1, detail = 'Creating vis_data structure...')
    
    # Create standardized vis_data structure using core data module
    vis_data <- core_data_module$create_vis_data_structure()
    
    incProgress(0.2, detail = 'Extracting peptides for all isoform×enzyme combinations...')
    
    # Enzyme names mapping (from multi_enzyme_server.R)
    enzyme_names <- list(
      "trp" = "Trypsin",
      "chymo" = "Chymotrypsin", 
      "aspn" = "Asp-N",
      "lysc" = "Lys-C",
      "lysn" = "Lys-N",
      "gluc" = "Glu-C"
    )
    
    # Matrix data structure: isoforms × enzymes
    matrix_peptides_list <- list()
    
    # Extract peptides for all combinations using the same function as existing modules
    for (isoform in selected_isoforms) {
      for (enzyme in selected_enzymes) {
        combination_key <- paste(isoform, enzyme, sep = "_")
        
        # Use get_transcript_peptides_for_comparison like multi-isoform analysis
        peptides_gr <- get_transcript_peptides_for_comparison(isoform, vis_data, enzyme)
        
        if (!is.null(peptides_gr) && length(peptides_gr) > 0) {
          matrix_peptides_list[[combination_key]] <- data.frame(
            isoform = isoform,
            enzyme = enzyme,
            enzyme_name = enzyme_names[[enzyme]],
            combination = combination_key,
            start = start(peptides_gr),  # Genomic coordinates
            end = end(peptides_gr),      # Genomic coordinates
            peptide = peptides_gr$peptide,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    incProgress(0.2, detail = 'Calculating matrix statistics...')
    
    # Generate coverage matrix (isoforms × enzymes)
    coverage_matrix <- create_coverage_matrix(selected_isoforms, selected_enzymes, matrix_peptides_list, current_gene_data)
    
    incProgress(0.2, detail = 'Analyzing peptide specificity across both dimensions...')
    
    # Generate specificity analysis across both isoforms and enzymes
    specificity_analysis <- analyze_matrix_specificity(selected_isoforms, selected_enzymes, matrix_peptides_list)
    
    incProgress(0.2, detail = 'Finalizing visualization data...')
    
    # Combine all matrix data
    if (length(matrix_peptides_list) == 0) {
      showNotification("No peptide data found for selected combinations", type = "warning")
      return(list(success = FALSE, message = "No peptide data available"))
    }
    
    all_matrix_peptides <- do.call(rbind, matrix_peptides_list)
    
    # Add positional information for visualization
    all_matrix_peptides$isoform_y <- match(all_matrix_peptides$isoform, selected_isoforms)
    all_matrix_peptides$enzyme_x <- match(all_matrix_peptides$enzyme, selected_enzymes)
    
    # Add hover text following existing patterns
    all_matrix_peptides$hover_text <- paste0(
      "Peptide: ", all_matrix_peptides$peptide,
      "<br>Isoform: ", all_matrix_peptides$isoform,
      "<br>Enzyme: ", all_matrix_peptides$enzyme_name,
      "<br>Position: ", all_matrix_peptides$start, "-", all_matrix_peptides$end,
      "<br>Combination: ", all_matrix_peptides$combination
    )
    
    # Calculate gene boundaries from genomic coordinates (from multi_enzyme_server.R)
    if (nrow(all_matrix_peptides) > 0) {
      gene_start <- min(all_matrix_peptides$start) - 5000
      gene_end <- max(all_matrix_peptides$end) + 5000
    } else {
      gene_start <- 0
      gene_end <- 10000
    }
    
    incProgress(0.1, detail = 'Complete!')
    
    return(list(
      success = TRUE,
      selected_isoforms = selected_isoforms,
      selected_enzymes = selected_enzymes,
      enzyme_names = enzyme_names,
      matrix_peptides = all_matrix_peptides,
      coverage_matrix = coverage_matrix,
      specificity_analysis = specificity_analysis,
      gene_start = gene_start,
      gene_end = gene_end,
      miscleavage_type = miscleavage_type
    ))
  })
})

# Helper function: Create coverage matrix (percentage coverage per combination)
create_coverage_matrix <- function(isoforms, enzymes, peptides_list, gene_data) {
  coverage_matrix <- matrix(0, nrow = length(isoforms), ncol = length(enzymes),
                           dimnames = list(isoforms, enzymes))
  
  # Calculate coverage for each combination
  for (isoform in isoforms) {
    for (enzyme in enzymes) {
      combination_key <- paste(isoform, enzyme, sep = "_")
      
      if (combination_key %in% names(peptides_list)) {
        peptide_data <- peptides_list[[combination_key]]
        
        # Calculate proper protein coverage using amino acid positions (from multi_enzyme_server.R)
        isoform_data <- gene_data$peptides[gene_data$peptides$txID == isoform, ]
        
        if (nrow(isoform_data) > 0) {
          position_col <- paste0(enzyme, "Peps_positions")
          
          if (position_col %in% names(isoform_data) && !is.null(isoform_data[[position_col]][[1]])) {
            position_df <- isoform_data[[position_col]][[1]]
            
            if (nrow(position_df) > 0) {
              # Get unique amino acid positions covered by peptides
              covered_positions <- c()
              for (i in 1:nrow(position_df)) {
                aa_range <- position_df$aa_start[i]:position_df$aa_end[i]
                covered_positions <- c(covered_positions, aa_range)
              }
              covered_positions <- unique(covered_positions)
              
              # Calculate protein length and coverage percentage
              protein_length <- max(position_df$aa_end)
              coverage_percent <- round((length(covered_positions) / protein_length) * 100, 1)
              
              coverage_matrix[isoform, enzyme] <- coverage_percent
            }
          }
        }
      }
    }
  }
  
  return(as.data.frame(coverage_matrix))
}

# Helper function: Analyze peptide specificity across both dimensions
analyze_matrix_specificity <- function(isoforms, enzymes, peptides_list) {
  if (length(peptides_list) == 0) return(NULL)
  
  # Combine all peptide data
  all_peptides <- do.call(rbind, peptides_list)
  
  # Analyze peptide specificity
  peptide_analysis <- list()
  
  for (peptide_seq in unique(all_peptides$peptide)) {
    peptide_rows <- all_peptides[all_peptides$peptide == peptide_seq, ]
    
    # Count unique isoforms and enzymes for this peptide
    unique_isoforms <- unique(peptide_rows$isoform)
    unique_enzymes <- unique(peptide_rows$enzyme)
    
    # Classify peptide across both dimensions
    isoform_specificity <- if (length(unique_isoforms) == length(isoforms)) {
      "Universal"
    } else if (length(unique_isoforms) == 1) {
      "Unique"
    } else {
      "Shared"
    }
    
    enzyme_specificity <- if (length(unique_enzymes) == length(enzymes)) {
      "Universal"
    } else if (length(unique_enzymes) == 1) {
      "Unique"
    } else {
      "Shared"
    }
    
    # Get enzyme names for this peptide
    enzyme_names_for_peptide <- unique(peptide_rows$enzyme_name)
    enzyme_names_text <- paste(enzyme_names_for_peptide, collapse = ", ")
    
    peptide_analysis[[peptide_seq]] <- data.frame(
      peptide = peptide_seq,
      isoform_count = length(unique_isoforms),
      enzyme_name = enzyme_names_text,
      isoform_specificity = isoform_specificity,
      stringsAsFactors = FALSE
    )
  }
  
  return(do.call(rbind, peptide_analysis))
}

# Create genomic visualization (always overlay plot with exon/CDS structure)
output$matrix_analysis_plot <- renderPlotly({
  data <- matrix_analysis_data()
  
  if (is.null(data) || !data$success) {
    return(empty_plotly_message("No matrix data available. Please run analysis with valid selections."))
  }
  
  # Always use genomic overlay plot with exon/CDS structure
  tryCatch({
    # Add compression support - respect user preference (using isoform-specific control)
    use_compression_isoforms <- !is.null(input$isoform_intron_scale) && input$isoform_intron_scale == "compressed"
    
    create_overlay_plot(data, use_compression = use_compression_isoforms, intron_scale = input$isoform_intron_scale)
  }, error = function(e) {
    cat("Error in matrix_analysis_plot:", e$message, "\n")
    return(empty_plotly_message(paste("Error creating genomic plot:", e$message)))
  })
})


# Helper function: Create overlay plot with separate tracks for each isoform×enzyme combination
create_overlay_plot <- function(data, use_compression = FALSE, intron_scale = NULL) {
  all_peptides <- data$matrix_peptides
  selected_isoforms <- data$selected_isoforms
  selected_enzymes <- data$selected_enzymes
  gene_start <- data$gene_start
  gene_end <- data$gene_end
  gene_id <- input$gene
  
  # Load GTF data for exon and CDS boundaries (following multi_enzyme_server.R pattern)
  exons_result <- NULL
  if (dir.exists("data/gtf_cache")) {
    gtf_data <- load_gtf_visualization_data(gene_id)
    if (gtf_data$success) {
      exons_result <- list(
        success = TRUE,
        exons = gtf_data$exons_by_transcript,
        cds = gtf_data$cds_by_transcript
      )
    } else {
      # Fast GTF cache failed, use fallback
      gene_details <- load_gene_details(gene_id)
      exons_result <- load_transcript_exons(gene_details, selected_isoforms)
    }
  } else {
    # No GTF cache, use original method
    gene_details <- load_gene_details(gene_id)
    exons_result <- load_transcript_exons(gene_details, selected_isoforms)
  }
  
  if (!exons_result$success) {
    return(plotly::plot_ly() %>%
      plotly::add_annotations(
        x = 0.5, y = 0.5,
        text = "Unable to load transcript structure",
        showarrow = FALSE,
        font = list(size = 16)
      ))
  }
  
  # Add compression support
  compression_map_isoforms <- NULL
  
  # Create compression map if user selected compressed mode and exon data exists
  if (use_compression && exons_result$success && length(exons_result$exons) > 0) {
    compression_map_isoforms <- create_compression_map(exons_result$exons)
  }
  
  # Define colors for enzymes (consistent coloring across all isoforms)
  enzyme_colors <- list(
    "trp" = "rgba(31, 119, 180, 0.8)",     # Blue
    "chymo" = "rgba(255, 127, 14, 0.8)",   # Orange  
    "aspn" = "rgba(44, 160, 44, 0.8)",     # Green
    "lysc" = "rgba(214, 39, 40, 0.8)",     # Red
    "lysn" = "rgba(148, 103, 189, 0.8)",   # Purple
    "gluc" = "rgba(140, 86, 75, 0.8)"      # Brown
  )
  
  # Create combinations list - each isoform×enzyme gets its own track
  combinations <- expand.grid(
    isoform = selected_isoforms,
    enzyme = selected_enzymes,
    stringsAsFactors = FALSE
  )
  combinations$track_id <- paste(combinations$isoform, combinations$enzyme, sep = "_")
  combinations$y_position <- 1:nrow(combinations)
  
  # Create plotly object
  p <- plotly::plot_ly()
  
  # Calculate overall gene boundaries and get chromosome info
  if (nrow(all_peptides) > 0) {
    final_gene_start <- min(all_peptides$start) - 5000
    final_gene_end <- max(all_peptides$end) + 5000
  } else {
    final_gene_start <- gene_start
    final_gene_end <- gene_end
  }
  
  # Apply compression to gene boundaries and peptide coordinates
  if (use_compression && !is.null(compression_map_isoforms)) {
    plot_gene_start <- compression_map_isoforms$compress(final_gene_start)
    plot_gene_end <- compression_map_isoforms$compress(final_gene_end)
    
    # Apply compression to all peptide coordinates
    all_peptides$plot_start <- sapply(all_peptides$start, compression_map_isoforms$compress)
    all_peptides$plot_end <- sapply(all_peptides$end, compression_map_isoforms$compress)
  } else {
    plot_gene_start <- final_gene_start
    plot_gene_end <- final_gene_end
    all_peptides$plot_start <- all_peptides$start
    all_peptides$plot_end <- all_peptides$end
  }
  
  # Get chromosome information from first available isoform (fix duplication bug)
  chromosome_info <- "unknown"
  for (isoform in selected_isoforms) {
    transcript_exons <- exons_result$exons[[isoform]]
    if (!is.null(transcript_exons) && length(transcript_exons) > 0) {
      chromosome_info <- as.character(seqnames(transcript_exons)[1])
      break  # Use first valid chromosome info
    }
  }
  
  # Add tracks for each isoform×enzyme combination
  for (i in 1:nrow(combinations)) {
    isoform <- combinations$isoform[i]
    enzyme <- combinations$enzyme[i]
    y_pos <- combinations$y_position[i]
    track_label <- paste(isoform, "-", data$enzyme_names[[enzyme]])
    
    # Get transcript structure for this isoform
    transcript_exons <- exons_result$exons[[isoform]]
    transcript_cds <- exons_result$cds[[isoform]]
    
    # Add transcript line (backbone) using compressed coordinates
    p <- p %>% plotly::add_trace(
      type = "scatter",
      x = c(plot_gene_start, plot_gene_end),
      y = c(y_pos, y_pos),
      mode = "lines",
      line = list(color = "rgba(204, 204, 204, 0.8)", width = 2),
      showlegend = FALSE,
      hovertemplate = paste0("Track: ", track_label, " (", chromosome_info, ")<extra></extra>"),
      hoverinfo = "none"
    )
    
    # Add exon blocks for this isoform
    if (!is.null(transcript_exons) && length(transcript_exons) > 0) {
      for (j in seq_along(transcript_exons)) {
        exon_start <- start(transcript_exons[j])
        exon_end <- end(transcript_exons[j])
        
        # Apply compression to exon coordinates
        if (use_compression && !is.null(compression_map_isoforms)) {
          exon_start <- compression_map_isoforms$compress(exon_start)
          exon_end <- compression_map_isoforms$compress(exon_end)
        }
        
        p <- p %>% plotly::add_trace(
          type = "scatter",
          mode = "lines",
          x = c(exon_start, exon_end, exon_end, exon_start, exon_start),
          y = c(y_pos - 0.3, y_pos - 0.3, y_pos + 0.3, y_pos + 0.3, y_pos - 0.3),
          fill = "toself",
          fillcolor = "rgba(211, 211, 211, 0.3)",
          line = list(color = "rgba(128, 128, 128, 0.5)", width = 1),
          showlegend = FALSE,
          hovertemplate = paste0(track_label, " - Exon ", j, " (", exon_start, "-", exon_end, ")<extra></extra>"),
          hoverinfo = "none"
        )
      }
    }
    
    # Add CDS blocks for this isoform
    if (!is.null(transcript_cds) && length(transcript_cds) > 0) {
      for (j in seq_along(transcript_cds)) {
        cds_start <- start(transcript_cds[j])
        cds_end <- end(transcript_cds[j])
        
        # Apply compression to CDS coordinates
        if (use_compression && !is.null(compression_map_isoforms)) {
          cds_start <- compression_map_isoforms$compress(cds_start)
          cds_end <- compression_map_isoforms$compress(cds_end)
        }
        
        p <- p %>% plotly::add_trace(
          type = "scatter",
          mode = "lines",
          x = c(cds_start, cds_end, cds_end, cds_start, cds_start),
          y = c(y_pos - 0.25, y_pos - 0.25, y_pos + 0.25, y_pos + 0.25, y_pos - 0.25),
          fill = "toself",
          fillcolor = "rgba(241, 196, 15, 0.8)",  # Yellow for CDS
          line = list(color = "rgba(218, 165, 32, 1)", width = 1),
          showlegend = FALSE,
          hovertemplate = paste0(track_label, " - CDS ", j, " (", cds_start, "-", cds_end, ")<extra></extra>"),
          hoverinfo = "none"
        )
      }
    }
    
    # Add peptides for this specific isoform×enzyme combination
    combination_peptides <- all_peptides[all_peptides$isoform == isoform & all_peptides$enzyme == enzyme, ]
    
    if (nrow(combination_peptides) > 0) {
      enzyme_color <- enzyme_colors[[enzyme]]
      
      for (j in 1:nrow(combination_peptides)) {
        peptide_start <- combination_peptides$plot_start[j]
        peptide_end <- combination_peptides$plot_end[j]
        
        p <- p %>% plotly::add_trace(
          type = "scatter",
          mode = "lines",
          x = c(peptide_start, peptide_end, peptide_end, peptide_start, peptide_start),
          y = c(y_pos - 0.15, y_pos - 0.15, y_pos + 0.15, y_pos + 0.15, y_pos - 0.15),
          fill = "toself",
          fillcolor = enzyme_color,
          line = list(color = "black", width = 0.5),
          legendgroup = enzyme,
          name = data$enzyme_names[[enzyme]],
          showlegend = i == 1 && j == 1,  # Only show legend for first occurrence of each enzyme
          hovertemplate = paste0(clean_hover_text(combination_peptides$hover_text[j]), "<extra></extra>"),
          hoverinfo = "none"
        )
      }
    }
  }
  
  # Create track labels combining isoform and enzyme names
  track_labels <- paste(combinations$isoform, "-", 
                       unlist(data$enzyme_names[combinations$enzyme]))
  
  # Layout with separate tracks
  p <- p %>% plotly::layout(
    title = list(
      text = paste0("Multi-Isoform Multi-Enzyme Analysis: ", length(selected_isoforms), " Isoforms × ", length(selected_enzymes), " Enzymes = ", nrow(combinations), " Tracks (", gene_id, ", ", chromosome_info, ")"),
      font = list(size = 16)
    ),
    xaxis = if (use_compression && !is.null(compression_map_isoforms)) {
      # Use the same axis mapping pattern as other visualizations
      axis_breaks <- compression_map_isoforms$coords$compressed_start
      axis_labels <- as.character(compression_map_isoforms$coords$original_start)
      compressed_range <- range(c(compression_map_isoforms$coords$compressed_start, compression_map_isoforms$coords$compressed_end))
      
      list(
        title = paste0("Genomic Position (compressed view) - chromosome ", chromosome_info, ")"),
        tickvals = axis_breaks,
        ticktext = axis_labels,
        range = compressed_range,
        showgrid = TRUE,
        gridcolor = "rgba(128,128,128,0.2)",
        zeroline = FALSE
      )
    } else {
      list(
        title = paste0("Genomic Position (chromosome ", chromosome_info, ")"),
        range = c(final_gene_start, final_gene_end),
        showgrid = TRUE,
        gridcolor = "rgba(128,128,128,0.2)"
      )
    },
    yaxis = list(
      title = "Isoform × Enzyme Combinations",
      showticklabels = TRUE,
      showgrid = FALSE,
      range = c(0.5, nrow(combinations) + 0.5),
      tickvals = combinations$y_position,
      ticktext = track_labels,
      tickfont = list(size = 10),
      tickangle = 0
    ),
    hovermode = "closest",
    legend = list(
      orientation = "h",
      x = 0,
      y = -0.15,
      title = list(text = "Enzymes:")
    ),
    margin = list(l = 150, r = 50, t = 80, b = 100),
    plot_bgcolor = "white",
    paper_bgcolor = "white"
  )
  
  return(p %>% clean_plotly_hover())
}


# Render coverage matrix table - DEPRECATED (now integrated into enzyme_peptide_specificity_table)
# output$coverage_matrix_table <- DT::renderDataTable({
#   req(matrix_analysis_data())
#   
#   data <- matrix_analysis_data()
#   if (!data$success || is.null(data$coverage_matrix)) {
#     return(data.frame(Message = "No coverage matrix available"))
#   }
#   
#   coverage_df <- data$coverage_matrix
#   coverage_df$Isoform <- rownames(coverage_df)
#   coverage_df <- coverage_df[, c("Isoform", colnames(data$coverage_matrix))]
#   
#   # Add enzyme names to column headers
#   enzyme_names <- unlist(data$enzyme_names[colnames(data$coverage_matrix)])
#   colnames(coverage_df)[-1] <- enzyme_names
#   
#   # Format coverage values and add % sign
#   dt <- DT::datatable(
#     coverage_df,
#     options = list(
#       pageLength = 10,
#       scrollX = TRUE,
#       dom = 'Bfrtip',
#       buttons = c('copy', 'csv')
#     ),
#     class = 'cell-border stripe',
#     rownames = FALSE
#   )
#   
#   # Format each enzyme column to show % and apply color styling
#   for (col in enzyme_names) {
#     if (col %in% colnames(coverage_df)) {
#       dt <- dt %>%
#         DT::formatRound(col, digits = 1) %>%
#         DT::formatString(col, suffix = "%") %>%
#         DT::formatStyle(
#           col,
#           backgroundColor = DT::styleInterval(
#             cuts = c(0, 25, 50, 75),
#             values = c("#FFFFFF", "#E3F2FD", "#90CAF9", "#42A5F5", "#1E88E5")
#           )
#         )
#     }
#   }
#   
#   dt
# })

# Render specificity analysis table
output$specificity_matrix_table <- DT::renderDataTable({
  req(matrix_analysis_data())
  
  data <- matrix_analysis_data()
  if (!data$success || is.null(data$specificity_analysis)) {
    return(data.frame(Message = "No specificity analysis available"))
  }
  
  specificity_df <- data$specificity_analysis
  
  DT::datatable(
    specificity_df,
    options = list(
      pageLength = 10,
      scrollX = TRUE,
      dom = 'Bfrtip',
      buttons = c('copy', 'csv')
    ),
    class = 'cell-border stripe',
    rownames = FALSE
  ) %>%
  DT::formatStyle(
    "isoform_specificity",
    backgroundColor = DT::styleEqual(
      c("Universal", "Shared", "Unique"),
      c("#ccffcc", "#ffffcc", "#ffcccc")
    )
  )
})

# Isoform-specific enzyme performance matrix
output$enzyme_peptide_specificity_table <- DT::renderDataTable({
  req(matrix_analysis_data())
  
  data <- matrix_analysis_data()
  if (!data$success || is.null(data$matrix_peptides) || is.null(data$coverage_matrix)) {
    return(data.frame(Message = "No enzyme analysis available"))
  }
  
  # Get selected enzymes and isoforms
  selected_enzymes <- data$selected_enzymes
  selected_isoforms <- data$selected_isoforms
  enzyme_names <- c(
    "trp" = "Trypsin",
    "chymo" = "Chymotrypsin", 
    "aspn" = "Asp-N",
    "lysc" = "Lys-C",
    "lysn" = "Lys-N",
    "gluc" = "Glu-C"
  )
  
  # Create combined performance matrix with coverage, peptide count, and unique peptides
  performance_matrix <- data.frame(
    Isoform = selected_isoforms,
    stringsAsFactors = FALSE
  )
  
  # Get coverage data
  coverage_df <- data$coverage_matrix
  
  # Add columns for each enzyme: Coverage (%), Peptide Count, Unique Peptide Count
  for (enzyme in selected_enzymes) {
    enzyme_name <- enzyme_names[enzyme]
    enzyme_coverage <- c()
    enzyme_counts <- c()
    enzyme_unique_counts <- c()
    
    for (isoform in selected_isoforms) {
      # Get coverage percentage from coverage matrix
      isoform_idx <- which(rownames(coverage_df) == isoform)
      enzyme_idx <- which(colnames(coverage_df) == enzyme)
      coverage_val <- if (length(isoform_idx) > 0 && length(enzyme_idx) > 0) {
        coverage_df[isoform_idx, enzyme_idx]
      } else {
        0
      }
      enzyme_coverage <- c(enzyme_coverage, coverage_val)
      
      # Count total unique peptides for this isoform-enzyme combination
      isoform_enzyme_peptides <- data$matrix_peptides[
        data$matrix_peptides$isoform == isoform & 
        data$matrix_peptides$enzyme == enzyme, 
      ]
      total_peptides <- length(unique(isoform_enzyme_peptides$peptide))
      enzyme_counts <- c(enzyme_counts, total_peptides)
      
      # Count isoform-specific peptides for this enzyme
      unique_peptides <- unique(isoform_enzyme_peptides$peptide)
      unique_count <- 0
      
      for (peptide in unique_peptides) {
        # Check if this peptide appears in other isoforms for this same enzyme
        other_isoforms_same_enzyme <- unique(data$matrix_peptides$isoform[
          data$matrix_peptides$peptide == peptide & 
          data$matrix_peptides$enzyme == enzyme
        ])
        if (length(other_isoforms_same_enzyme) == 1 && other_isoforms_same_enzyme[1] == isoform) {
          unique_count <- unique_count + 1
        }
      }
      enzyme_unique_counts <- c(enzyme_unique_counts, unique_count)
    }
    
    # Add three columns for this enzyme: Coverage (%), Peptide Count, Unique Peptide Count
    performance_matrix[[paste0(enzyme_name, " Coverage (%)")]] <- enzyme_coverage
    performance_matrix[[paste0(enzyme_name, " Peptide Count")]] <- enzyme_counts
    performance_matrix[[paste0(enzyme_name, " Unique Peptide Count")]] <- enzyme_unique_counts
  }
  
  
  if (nrow(performance_matrix) == 0) {
    return(data.frame(Message = "No performance data available"))
  }
  
  # Create the datatable without color highlighting
  dt <- DT::datatable(
    performance_matrix,
    options = list(
      pageLength = 15,
      scrollX = TRUE,
      dom = 'Bfrtip',
      buttons = c('copy', 'csv'),
      columnDefs = list(
        list(className = 'dt-center', targets = '_all')
      )
    ),
    class = 'cell-border stripe',
    rownames = FALSE
  )
  
  # Format columns for each enzyme (only formatting, no colors)
  for (enzyme in selected_enzymes) {
    enzyme_name <- enzyme_names[enzyme]
    
    # Coverage column - format with % only (no color)
    coverage_col <- paste0(enzyme_name, " Coverage (%)")
    if (coverage_col %in% colnames(performance_matrix)) {
      dt <- dt %>% 
        DT::formatRound(coverage_col, digits = 1) %>%
        DT::formatString(coverage_col, suffix = "%")
    }
  }
  
  return(dt)
})

