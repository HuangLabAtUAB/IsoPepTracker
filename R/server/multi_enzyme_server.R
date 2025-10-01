#===============================================================================
# MULTI-ENZYME COVERAGE ANALYSIS
#===============================================================================

# Update enzyme coverage isoform choices when isoform analysis is loaded
observeEvent(input$load_all_isoforms, {
  req(selected_gene_transcripts())
  
  transcripts <- selected_gene_transcripts()
  if (length(transcripts) > 0) {
    transcript_choices <- setNames(transcripts, transcripts)
    updateSelectInput(session, "enzyme_coverage_isoform", 
                     choices = transcript_choices,
                     selected = transcripts[1])
  }
})

# Multi-enzyme coverage data processing
multi_enzyme_coverage_data <- reactive({
  req(input$run_enzyme_coverage > 0, input$enzyme_coverage_isoform, 
      input$enzyme_miscleavage_global, gene_data())
  
  selected_isoform <- input$enzyme_coverage_isoform
  miscleavage_type <- input$enzyme_miscleavage_global
  
  # Get selected enzymes
  selected_enzymes <- c()
  if (input$use_trp) selected_enzymes <- c(selected_enzymes, "trp")
  if (input$use_chymo) selected_enzymes <- c(selected_enzymes, "chymo")
  if (input$use_aspn) selected_enzymes <- c(selected_enzymes, "aspn")
  if (input$use_lysc) selected_enzymes <- c(selected_enzymes, "lysc")
  if (input$use_lysn) selected_enzymes <- c(selected_enzymes, "lysn")
  if (input$use_gluc) selected_enzymes <- c(selected_enzymes, "gluc")
  
  if (length(selected_enzymes) == 0) {
    showNotification("Please select at least one enzyme for comparison", type = "warning")
    return(NULL)
  }
  
  withProgress(message = 'Processing multi-enzyme coverage...', value = 0, {
    # Load gene data with correct miscleavage type if needed
    current_gene_data <- gene_data()
    if (current_gene_data$miscleavage_type != miscleavage_type) {
      incProgress(0.2, detail = 'Loading gene data with correct miscleavage settings...')
      success <- load_and_cache_gene_data(input$gene, miscleavage_type)
      if (!success) {
        showNotification("Failed to load gene data", type = "error")
        return(NULL)
      }
      current_gene_data <- gene_data()
    }
    
    incProgress(0.3, detail = 'Creating vis_data structure...')
    
    # Create standardized vis_data structure using core data module
    vis_data <- core_data_module$create_vis_data_structure()
    
    incProgress(0.5, detail = 'Extracting enzyme-specific peptides with genomic coordinates...')
    
    # Get enzyme-specific peptides with genomic coordinates (same as multi-isoform)
    enzyme_names <- list(
      "trp" = "Trypsin",
      "chymo" = "Chymotrypsin", 
      "aspn" = "Asp-N",
      "lysc" = "Lys-C",
      "lysn" = "Lys-N",
      "gluc" = "Glu-C"
    )
    
    # Extract peptides for each enzyme using same pattern as multi-isoform
    all_enzyme_peptides_list <- list()
    for (i in seq_along(selected_enzymes)) {
      enzyme_code <- selected_enzymes[i]
      enzyme_name <- enzyme_names[[enzyme_code]]
      
      # Use get_transcript_peptides_for_comparison like multi-isoform analysis
      enzyme_peptides_granges <- get_transcript_peptides_for_comparison(selected_isoform, vis_data, enzyme_code)
      
      if (!is.null(enzyme_peptides_granges) && length(enzyme_peptides_granges) > 0) {
        all_enzyme_peptides_list[[enzyme_code]] <- data.frame(
          transcript = selected_isoform,
          enzyme = enzyme_name,
          enzyme_code = enzyme_code,
          y_position = i,  # Assign y position for each enzyme track
          start = start(enzyme_peptides_granges),  # Genomic coordinates
          end = end(enzyme_peptides_granges),      # Genomic coordinates
          peptide = enzyme_peptides_granges$peptide,
          stringsAsFactors = FALSE
        )
      }
    }
    
    incProgress(0.7, detail = 'Combining enzyme data...')
    
    # Combine all enzyme data
    if (length(all_enzyme_peptides_list) == 0) {
      showNotification("No peptide data found for selected enzymes", type = "warning")
      return(list(success = FALSE, message = "No peptide data available"))
    }
    
    all_enzyme_peptides <- do.call(rbind, all_enzyme_peptides_list)
    
    # Create enzyme dataframe for track layout
    enzyme_df <- data.frame(
      enzyme = unlist(enzyme_names[selected_enzymes]),
      enzyme_code = selected_enzymes,
      y_position = seq_along(selected_enzymes),
      stringsAsFactors = FALSE
    )
    
    # Calculate gene boundaries from genomic coordinates
    if (nrow(all_enzyme_peptides) > 0) {
      gene_start <- min(all_enzyme_peptides$start) - 5000
      gene_end <- max(all_enzyme_peptides$end) + 5000
    } else {
      gene_start <- 0
      gene_end <- 10000
    }
    
    incProgress(0.8, detail = 'Calculating coverage statistics...')
    
    # Calculate coverage statistics per enzyme using PROPER amino acid positions
    coverage_stats <- list()
    isoform_data <- current_gene_data$peptides[current_gene_data$peptides$txID == selected_isoform, ]
    
    for (enzyme_code in selected_enzymes) {
      enzyme_data <- all_enzyme_peptides[all_enzyme_peptides$enzyme_code == enzyme_code, ]
      if (nrow(enzyme_data) > 0) {
        # Calculate basic coverage statistics
        total_peptides <- nrow(enzyme_data)
        
        # Calculate PROPER protein coverage using amino acid positions
        protein_coverage_percent <- 0
        
        # Get amino acid position data from original isoform data
        position_col <- paste0(enzyme_code, "Peps_positions")
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
            
            # Calculate protein length from maximum amino acid position
            protein_length <- max(position_df$aa_end)
            
            # Calculate proper coverage percentage (amino acids covered / total protein length)
            protein_coverage_percent <- round((length(covered_positions) / protein_length) * 100, 1)
          }
        }
        
        coverage_stats[[enzyme_code]] <- list(
          enzyme = enzyme_names[[enzyme_code]],
          total_peptides = total_peptides,
          protein_coverage_percent = protein_coverage_percent
        )
      }
    }
    
    incProgress(1.0, detail = 'Finalizing results...')
    
    return(list(
      isoform = selected_isoform,
      miscleavage_type = miscleavage_type,
      selected_enzymes = selected_enzymes,
      all_enzyme_peptides = all_enzyme_peptides,  # Structured data with genomic coordinates
      enzyme_df = enzyme_df,                      # Enzyme track layout data
      gene_start = gene_start,                    # Gene boundaries
      gene_end = gene_end,
      coverage_stats = coverage_stats,
      success = TRUE
    ))
  })
})

# Create multi-enzyme highlighted data (similar to multi-isoform analysis)
multi_enzyme_highlighted_data <- reactive({
  req(multi_enzyme_coverage_data())
  
  base_data <- multi_enzyme_coverage_data()
  if (is.null(base_data) || !base_data$success) return(NULL)
  
  all_enzyme_peptides <- base_data$all_enzyme_peptides
  selected_enzymes <- base_data$selected_enzymes
  
  # Add enzyme specificity categories (similar to isoform specificity)
  all_enzyme_peptides$enzyme_category <- ""
  
  for (i in 1:nrow(all_enzyme_peptides)) {
    peptide_seq <- all_enzyme_peptides$peptide[i]
    current_enzyme <- all_enzyme_peptides$enzyme_code[i]
    
    # Count other enzymes that have this peptide
    other_enzymes_with_peptide <- unique(all_enzyme_peptides$enzyme_code[
      all_enzyme_peptides$peptide == peptide_seq & all_enzyme_peptides$enzyme_code != current_enzyme
    ])
    
    # Classify enzyme specificity
    if (length(other_enzymes_with_peptide) == 0) {
      # Peptide unique to this enzyme
      all_enzyme_peptides$enzyme_category[i] <- "Unique"
    } else if (length(other_enzymes_with_peptide) == (length(selected_enzymes) - 1)) {
      # Peptide found in all enzymes
      all_enzyme_peptides$enzyme_category[i] <- "Universal"
    } else {
      # Peptide shared with some enzymes
      all_enzyme_peptides$enzyme_category[i] <- "Shared"
    }
  }
  
  # Add hover text
  all_enzyme_peptides$hover_text <- paste0(
    "Peptide: ", all_enzyme_peptides$peptide,
    "<br>Enzyme: ", all_enzyme_peptides$enzyme,
    "<br>Position: ", all_enzyme_peptides$start, "-", all_enzyme_peptides$end,
    "<br>Specificity: ", all_enzyme_peptides$enzyme_category
  )
  
  return(list(
    all_enzyme_peptides = all_enzyme_peptides,
    enzyme_df = base_data$enzyme_df,
    gene_start = base_data$gene_start,
    gene_end = base_data$gene_end,
    selected_enzymes = base_data$selected_enzymes,
    isoform = base_data$isoform
  ))
})

# Render enzyme coverage plot using TRUE genomic visualization like multi-isoform analysis
output$enzyme_coverage_plot <- renderPlotly({
  req(multi_enzyme_highlighted_data())
  
  data <- multi_enzyme_highlighted_data()
  if (is.null(data)) {
    p <- plotly::plot_ly() %>%
      plotly::add_annotations(
        x = 0.5, y = 0.5,
        text = "No peptide data available for selected enzymes",
        showarrow = FALSE,
        font = list(size = 16)
      ) %>%
      plotly::layout(
        xaxis = list(range = c(0, 1), showticklabels = FALSE),
        yaxis = list(range = c(0, 1), showticklabels = FALSE)
      )
    return(p)
  }
  
  selected_isoform <- data$isoform
  selected_enzymes <- data$selected_enzymes
  all_enzyme_peptides <- data$all_enzyme_peptides
  enzyme_df <- data$enzyme_df
  gene_start <- data$gene_start
  gene_end <- data$gene_end
  gene_id <- input$gene
  
  # Load GTF data for exon and CDS boundaries (lightning-fast cached GTF)
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
      exons_result <- load_transcript_exons(gene_details, c(selected_isoform))
    }
  } else {
    # No GTF cache, use original method
    gene_details <- load_gene_details(gene_id)
    exons_result <- load_transcript_exons(gene_details, c(selected_isoform))
  }
  
  if (!exons_result$success) {
    p <- plotly::plot_ly() %>%
      plotly::add_annotations(
        x = 0.5, y = 0.5,
        text = "Unable to load transcript structure",
        showarrow = FALSE,
        font = list(size = 16)
      )
    return(p)
  }
  
  # Add compression support - respect user preference (using isoform-specific control)
  use_compression_isoforms <- !is.null(input$isoform_intron_scale) && input$isoform_intron_scale == "compressed"
  compression_map_isoforms <- NULL
  
  # Create compression map if user selected compressed mode and exon data exists
  if (use_compression_isoforms && exons_result$success && length(exons_result$exons) > 0) {
    compression_map_isoforms <- create_compression_map(exons_result$exons)
  }
  
  # Get transcript structure
  transcript_exons <- exons_result$exons[[selected_isoform]]
  transcript_cds <- exons_result$cds[[selected_isoform]]
  
  # Check if transcript structure is valid
  if (is.null(transcript_exons) || length(transcript_exons) == 0) {
    p <- plotly::plot_ly() %>%
      plotly::add_annotations(
        x = 0.5, y = 0.5,
        text = "No transcript structure available for selected isoform",
        showarrow = FALSE,
        font = list(size = 16)
      )
    return(p)
  }
  
  # Calculate genomic boundaries with proper padding (same as multi-isoform)
  padding <- 5000
  gene_start <- min(start(transcript_exons)) - padding
  gene_end <- max(end(transcript_exons)) + padding
  chromosome <- as.character(seqnames(transcript_exons)[1])
  
  # Define colors for enzymes using rgba to prevent hex codes in hover
  enzyme_colors <- list(
    "trp" = "rgba(31, 119, 180, 0.8)",     # Blue
    "chymo" = "rgba(255, 127, 14, 0.8)",   # Orange  
    "aspn" = "rgba(44, 160, 44, 0.8)",     # Green
    "lysc" = "rgba(214, 39, 40, 0.8)",     # Red
    "lysn" = "rgba(148, 103, 189, 0.8)",   # Purple
    "gluc" = "rgba(140, 86, 75, 0.8)"      # Brown
  )
  
  
  # Create plotly object
  p <- plotly::plot_ly()
  
  # For each enzyme track, modify the base plot to show enzyme-specific peptides
  for (i in 1:nrow(enzyme_df)) {
    enzyme_name <- enzyme_df$enzyme[i]
    enzyme_code <- enzyme_df$enzyme_code[i]
    y_pos <- enzyme_df$y_position[i]
    
    # Add transcript line with compression support
    # Apply compression to gene boundaries
    plot_gene_start <- gene_start
    plot_gene_end <- gene_end
    if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
      plot_gene_start <- compression_map_isoforms$compress(gene_start)
      plot_gene_end <- compression_map_isoforms$compress(gene_end)
    }
    
    p <- p %>% plotly::add_trace(
      type = "scatter",
      x = c(plot_gene_start, plot_gene_end),
      y = c(y_pos, y_pos),
      mode = "lines",
      line = list(color = "rgba(204, 204, 204, 0.8)", width = 2),
      showlegend = FALSE,
      hovertemplate = paste0("Enzyme: ", enzyme_name, "<extra></extra>"),
      hoverinfo = "none"
    )
    
    # Add exon blocks (exactly like multi-isoform analysis)
    if (!is.null(transcript_exons) && length(transcript_exons) > 0) {
      for (j in seq_along(transcript_exons)) {
        exon_start <- start(transcript_exons[j])
        exon_end <- end(transcript_exons[j])
        
        # Apply compression to exon coordinates
        if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
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
          hovertemplate = paste0("Exon ", j, " (", exon_start, "-", exon_end, ")<extra></extra>"),
          hoverinfo = "none"
        )
      }
    }
    
    # Add CDS blocks (exactly like multi-isoform analysis)
    if (!is.null(transcript_cds) && length(transcript_cds) > 0) {
      for (j in seq_along(transcript_cds)) {
        cds_start <- start(transcript_cds[j])
        cds_end <- end(transcript_cds[j])
        
        # Apply compression to CDS coordinates
        if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
          cds_start <- compression_map_isoforms$compress(cds_start)
          cds_end <- compression_map_isoforms$compress(cds_end)
        }
        
        p <- p %>% plotly::add_trace(
          type = "scatter",
          mode = "lines",
          x = c(cds_start, cds_end, cds_end, cds_start, cds_start),
          y = c(y_pos - 0.25, y_pos - 0.25, y_pos + 0.25, y_pos + 0.25, y_pos - 0.25),
          fill = "toself",
          fillcolor = "rgba(241, 196, 15, 0.8)",
          line = list(color = "rgba(218, 165, 32, 1)", width = 1),
          showlegend = FALSE,
          hovertemplate = paste0("CDS ", j, " (", cds_start, "-", cds_end, ")<extra></extra>"),
          hoverinfo = "none"
        )
      }
    }
    
    # Add enzyme-specific peptides using structured data (SAME pattern as multi-isoform)
    enzyme_peptides <- all_enzyme_peptides[all_enzyme_peptides$enzyme_code == enzyme_code, ]
    
    if (nrow(enzyme_peptides) > 0) {
      enzyme_color <- enzyme_colors[[enzyme_code]]
      
      # Create peptide rectangles using EXACT same pattern as multi-isoform analysis
      for (k in 1:nrow(enzyme_peptides)) {
        peptide_start <- enzyme_peptides$start[k]
        peptide_end <- enzyme_peptides$end[k]
        
        # Apply compression to peptide coordinates
        if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
          peptide_start <- compression_map_isoforms$compress(peptide_start)
          peptide_end <- compression_map_isoforms$compress(peptide_end)
        }
        
        # Use exact same Plotly pattern as multi-isoform comparative analysis
        p <- p %>% plotly::add_trace(
          type = "scatter", 
          mode = "lines",
          x = c(peptide_start, peptide_end, peptide_end, peptide_start, peptide_start),
          y = c(y_pos - 0.15, y_pos - 0.15, y_pos + 0.15, y_pos + 0.15, y_pos - 0.15),
          fill = "toself",
          fillcolor = enzyme_color,
          line = list(color = "black", width = 0.5),
          marker = list(opacity = 0),
          legendgroup = enzyme_name,
          name = enzyme_name,  # Proper enzyme name for legend
          showlegend = k == 1,  # Only show legend for first peptide
          hovertemplate = paste0(clean_hover_text(enzyme_peptides$hover_text[k]), "<extra></extra>"),
          hoverinfo = "none"  # Disable default hover, use hovertemplate instead
        )
      }
    }
  }
  
  # Configure X-axis based on compression mode
  xaxis_config <- if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
    # Compressed view: real genomic coordinates at compressed positions
    axis_breaks <- compression_map_isoforms$coords$compressed_start
    axis_labels <- as.character(compression_map_isoforms$coords$original_start)
    compressed_range <- range(c(compression_map_isoforms$coords$compressed_start, compression_map_isoforms$coords$compressed_end))
    
    list(
      title = paste0("Genomic Position (compressed view) - chromosome ", chromosome),
      tickvals = axis_breaks,
      ticktext = axis_labels,
      range = compressed_range,
      showgrid = TRUE,
      gridcolor = "rgba(128,128,128,0.2)",
      zeroline = FALSE
    )
  } else {
    # Original view
    list(
      title = paste0("Genomic Position (chromosome ", chromosome, ")"),
      range = c(gene_start, gene_end),
      showgrid = TRUE,
      gridcolor = "rgba(128,128,128,0.2)",
      zeroline = FALSE
    )
  }
  
  # Add title and layout with compression-aware axis
  p <- p %>% plotly::layout(
    title = list(
      text = paste0("Multi-Enzyme Genomic Coverage: ", selected_isoform, " (Chr", chromosome, ")"),
      font = list(size = 16)
    ),
    xaxis = xaxis_config,
    yaxis = list(
      title = "Enzymes",
      showticklabels = TRUE,
      showgrid = FALSE,
      range = c(0.5, nrow(enzyme_df) + 0.5),
      tickvals = enzyme_df$y_position,
      ticktext = enzyme_df$enzyme,
      tickfont = list(size = 12)
    ),
    hovermode = "closest",
    legend = list(
      orientation = "h",
      x = 0,
      y = -0.1
    ),
    margin = list(l = 120, r = 50, t = 80, b = 80),
    plot_bgcolor = "white",
    paper_bgcolor = "white"
  )
  
  return(p %>% clean_plotly_hover())
})

# Render coverage statistics table
output$enzyme_coverage_stats_table <- DT::renderDataTable({
  req(multi_enzyme_coverage_data())
  
  data <- multi_enzyme_coverage_data()
  if (!data$success || length(data$coverage_stats) == 0) {
    return(data.frame(Message = "No coverage statistics available"))
  }
  
  # Convert to dataframe
  stats_df <- do.call(rbind, lapply(data$coverage_stats, function(x) {
    data.frame(
      Enzyme = x$enzyme,
      Total_Peptides = x$total_peptides,
      Protein_Coverage = paste0(x$protein_coverage_percent, "%"),
      stringsAsFactors = FALSE
    )
  }))
  
  DT::datatable(stats_df, 
                options = list(pageLength = 10, scrollX = TRUE),
                colnames = c("Enzyme", "Total Peptides", "Protein Coverage"))
})


# Download handler for enzyme coverage results
output$download_enzyme_coverage <- downloadHandler(
  filename = function() {
    miscleavage_suffix <- switch(
      input$enzyme_miscleavage_global,
      "no_miss_cleavage" = "no_miss",
      "upto_two_misscleavage" = "upto_2_miss"
    )
    paste0("enzyme_coverage_", input$enzyme_coverage_isoform, "_", 
          miscleavage_suffix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  },
  content = function(file) {
    req(multi_enzyme_coverage_data())
    
    data <- multi_enzyme_coverage_data()
    if (data$success && nrow(data$all_enzyme_peptides) > 0) {
      # Add metadata
      peptides_df <- data$all_enzyme_peptides
      peptides_df$gene <- input$gene
      peptides_df$isoform <- data$isoform
      peptides_df$miscleavage_type <- data$miscleavage_type
      peptides_df$analysis_timestamp <- Sys.time()
      
      write.csv(peptides_df, file, row.names = FALSE)
    }
  }
)