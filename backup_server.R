#===============================================================================
# SERVER LOGIC (from server.R)
#===============================================================================

# Define the NULL coalescing operator
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

server <- function(input, output, session) {
  
  # Clear file inputs on session start to prevent autoloading
  # Note: File inputs will be properly handled by robust conditional checks in UI
  
  # Source compression functions
  if (!exists("create_compression_map")) {
    source("R/coordinate_compression.R")
  }
  
  # ============================================================================
  # CARD-BASED UI NAVIGATION - NO SIDEBAR
  # ============================================================================
  
  # Navigation is handled by JavaScript onclick and current_view selectInput
  # No server-side navigation needed for pure card-based UI
  
  # ============================================================================
  # LOAD CORE DATA MANAGEMENT MODULE
  # ============================================================================
  source("R/modules/core_data_management_module.R", local = TRUE)
  core_data_module <- create_core_data_management_module(input, output, session)
  
  # Extract reactive values from core data management module
  processed_data <- core_data_module$processed_data
  gene_data <- core_data_module$gene_data
  gene_cache <- core_data_module$gene_cache
  selected_gene_transcripts <- core_data_module$selected_gene_transcripts
  selected_gene_as_events <- core_data_module$selected_gene_as_events
  selected_gene_peptides <- core_data_module$selected_gene_peptides
  load_and_cache_gene_data <- core_data_module$load_and_cache_gene_data
  
  # ============================================================================
  # LOAD NOVEL ISOFORM ANALYSIS MODULE
  # ============================================================================
  source("R/modules/novel_isoform_analysis_module.R", local = TRUE)
  novel_isoform_module <- create_novel_isoform_analysis_module(input, output, session)
  
  # Extract reactive values from novel isoform analysis module
  novel_pipeline_results <- novel_isoform_module$novel_pipeline_results
  novel_isoform_data <- novel_isoform_module$novel_isoform_data
  novel_analysis_results <- novel_isoform_module$novel_analysis_results
  novel_orf_results <- novel_isoform_module$novel_orf_results
  novel_gene_search_results <- novel_isoform_module$novel_gene_search_results
  novel_merged_data <- novel_isoform_module$novel_merged_data
  novel_all_isoforms_data <- novel_isoform_module$novel_all_isoforms_data
  novel_highlighted_isoform_data <- novel_isoform_module$novel_highlighted_isoform_data
  
  # ============================================================================
  # LOAD MULTI-ISOFORM COMPARISON MODULE
  # ============================================================================
  source("R/modules/multi_isoform_comparison_module.R", local = TRUE)
  multi_isoform_module <- create_multi_isoform_comparison_module(
    input, output, session, processed_data, gene_data, 
    selected_gene_transcripts, selected_gene_peptides, novel_pipeline_results
  )
  
  # Extract reactive values from multi-isoform comparison module
  multi_isoform_data <- multi_isoform_module$multi_isoform_data
  multi_isoform_highlighted_data <- multi_isoform_module$multi_isoform_highlighted_data
  comparison_isoforms <- multi_isoform_module$comparison_isoforms
  
  # Store visualization data
  peptide_visualization_data <- reactiveVal(NULL)
  exon_data <- reactiveVal(NULL)
  
  # Load common data needed for both visualizations
  observeEvent(input$update, {
    req(input$gene, input$miscleavage_type)
    
    withProgress(message = 'Updating visualizations...', {
      gene_id <- input$gene
      transcript_ids <- isolate(selected_gene_transcripts())
      gene_symbol <- isolate(core_data_module$current_gene_symbol())
      
      # Use cached visualization data
      message("🚀 Using cached visualization for gene: ", gene_id)
      incProgress(0.4, detail = "Creating transcript structure")
      plot_data <- create_fast_transcript_plot_data(gene_id, gene_symbol, selected_gene_as_events())
      
      # Create compatible exon_data for AS visualization using cached GTF data
      gtf_data <- load_gtf_visualization_data(gene_id)
      if (gtf_data$success) {
        # Create exon_data structure that's compatible with create_as_view
        exon_data(list(
          success = TRUE,
          exons = gtf_data$exons_by_transcript,
          cds = gtf_data$cds_by_transcript,
          message = "Using cached data"
        ))
      } else {
        exon_data(list(success = FALSE, message = gtf_data$message))
      }
      
      # Create standardized vis_data structure using core data module
      vis_data <- core_data_module$create_vis_data_structure()
      
      # Use cached peptide visualization
      incProgress(0.2, detail = "Creating peptide visualization")
      message("🚀 Using cached peptide visualization for gene: ", gene_id)
      # Get compression preference
      use_compression <- !is.null(input$intron_scale) && input$intron_scale == "compressed"
      peptide_data <- create_fast_peptide_plot_data(gene_id, gene_symbol, vis_data, input$protease, use_compression)
      peptide_visualization_data(peptide_data)
      
      incProgress(0.2, detail = "Creating AS visualization")
      as_data <- create_as_view(
        exon_data(),
        gene_id,
        gene_symbol,
        selected_gene_as_events()
      )
      
      # Update column names properly when creating the data
      if (!is.null(as_data$events)) {
        colnames(as_data$events) <- c("EventID", "Type", "Spliced-out Transcripts", "Spliced-in Transcripts")
      }
      
      as_view_data(as_data)
    })
  })
  
  # Update peptide visualization when protease, miscleavage type, or intron scale changes
  observeEvent(list(input$protease, input$miscleavage_type, input$intron_scale), {
    req(exon_data(), input$miscleavage_type)
    if (!exon_data()$success) return()
    
    gene_id <- input$gene
    transcript_ids <- isolate(selected_gene_transcripts())
    gene_symbol <- isolate(core_data_module$current_gene_symbol())
    
    # Create standardized vis_data structure using core data module
    vis_data <- core_data_module$create_vis_data_structure()
    
    # Create peptide visualization plot with new protease
    # Check if we have AS events (same logic as main update)
    gene_as_events <- selected_gene_as_events()
    
    # Get compression preference
    use_compression <- !is.null(input$intron_scale) && input$intron_scale == "compressed"
    
    if (nrow(gene_as_events) > 0) {
      # Use alternative splicing coordinates
      peptide_data <- tryCatch({
        create_peptide_plot_data(
          exon_data(), exon_data()$transcript_ids, gene_symbol, gene_id, vis_data, input$protease, use_compression
        )
      }, error = function(e) {
        if (grepl("unused argument.*use_compression", e$message)) {
          # Fallback: call without use_compression parameter
          cat("WARNING: create_peptide_plot_data doesn't support use_compression, using default\n")
          create_peptide_plot_data(
            exon_data(), exon_data()$transcript_ids, gene_symbol, gene_id, vis_data, input$protease
          )
        } else {
          stop(e)
        }
      })
    } else {
      # Use lightning-fast visualization if GTF cache exists (for genes without AS events)
      if (dir.exists("data/gtf_cache")) {
        peptide_data <- create_fast_peptide_plot_data(gene_id, gene_symbol, vis_data, input$protease, use_compression)
      } else {
        # Fallback to original method
        peptide_data <- tryCatch({
          create_peptide_plot_data(
            exon_data(), transcript_ids, gene_symbol, gene_id, vis_data, input$protease, use_compression
          )
        }, error = function(e) {
          if (grepl("unused argument.*use_compression", e$message)) {
            # Fallback: call without use_compression parameter
            cat("WARNING: create_peptide_plot_data doesn't support use_compression, using default\n")
            create_peptide_plot_data(
              exon_data(), transcript_ids, gene_symbol, gene_id, vis_data, input$protease
            )
          } else {
            stop(e)
          }
        })
      }
    }
    
    # Store result
    peptide_visualization_data(peptide_data)
  })
  
  # Dynamic title for peptide plot showing miscleavage type
  output$peptide_plot_title <- renderText({
    if (!is.null(input$miscleavage_type)) {
      miscleavage_label <- switch(input$miscleavage_type,
        "no_miss_cleavage" = "No Miscleavage",
        "upto_two_misscleavage" = "Up to 2 Miscleavages"
      )
      paste("Peptide Visualization -", miscleavage_label)
    } else {
      "Peptide Visualization"
    }
  })
  
  
  # Render peptide visualization plot
  output$peptide_plot <- renderPlotly({
    req(input$gene, input$protease)
    req(isolate(selected_gene_transcripts()))
    
    # Get plot data
    plot_data <- peptide_visualization_data()
    
    # If no data or error, show error message
    if (is.null(plot_data) || !plot_data$success) {
      message_text <- if (is.null(plot_data)) {
        "Click 'Update View' to load peptide visualization"
      } else {
        plot_data$message
      }
      
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                label = message_text, 
                size = 5, hjust = 0.5) +
        theme_void() +
        xlim(0, 1) + ylim(0, 1)
      
      return(ggplotly(p))
    }
    
    # Return plotly version for interactivity using optimized function
    create_optimized_plotly(plot_data$plot, tooltip = "text") %>%
      toWebGL()
  })
  
  # Store AS visualization data
  as_view_data <- reactiveVal(NULL)
  
  # Update AS view when gene changes
  observeEvent(input$update, {
    req(input$gene, input$miscleavage_type, exon_data())
    if (!exon_data()$success) return()
    
    gene_id <- input$gene
    gene_symbol <- core_data_module$current_gene_symbol()
    
    # Create AS visualization
    as_data <- create_as_view(
      exon_data(),
      gene_id,
      gene_symbol,
      selected_gene_as_events()
    )
    
    # Update column names properly when creating the data
    if (!is.null(as_data$events)) {
      colnames(as_data$events) <- c("EventID", "Type", "Spliced-out Transcripts", "Spliced-in Transcripts")
    }
    
    as_view_data(as_data)
  })
  
  # Render AS structure plot with properly registered click event
  output$as_structure_plot <- renderPlotly({
    req(as_view_data())
    
    # Get the plot
    p <- as_view_data()$plot
    
    # Return plot without click event registration
    return(p %>% clean_plotly_hover())
  })
  
  # Render AS events table
  output$as_events_table <- DT::renderDataTable({
    req(as_view_data())
    
    # Get the events data
    events_data <- as_view_data()$events
    
    # Create the datatable
    DT::datatable(
      events_data,
      options = list(
        pageLength = 5,
        scrollX = TRUE,
        dom = 'tip'  # Show table, info, and pagination only
      ),
      rownames = FALSE
    )
  })
  
  # Add reactive values for selected event
  selected_event <- reactiveVal(NULL)
  
  # Store whether an event is selected
  output$event_selected <- reactive({
    !is.null(selected_event())
  })
  outputOptions(output, "event_selected", suspendWhenHidden = FALSE)
  
  # Replace the problematic observer with this corrected version
  observeEvent(event_data("plotly_click", source = "as_plot"), {
    click_data <- event_data("plotly_click", source = "as_plot")
    
    # Only proceed if we have valid click data and events
    if (!is.null(click_data) && !is.null(as_view_data())) {
      event_df <- as_view_data()$events
      
      # Only proceed if we have events
      if (!is.null(event_df) && nrow(event_df) > 0) {
        # Get the clicked point number (plotly is 0-indexed)
        point_number <- click_data$pointNumber
        
        # Make sure we have a valid point number
        if (!is.null(point_number) && length(point_number) > 0) {
          point_number <- as.numeric(point_number)
          if (!is.na(point_number)) {
            # Calculate the actual event index (1-based)
            event_index <- point_number + 1
            
            # Verify the index is within bounds
            if (event_index > 0 && event_index <= nrow(event_df)) {
            # Get the event ID
            selected_event_id <- event_df$event_id[event_index]
            
            # Look up complete event data
            gene_id <- input$gene
            gene_events <- selected_gene_as_events()
            matching_events <- gene_events[gene_events$eventID == selected_event_id,]
            
            if (nrow(matching_events) > 0) {
              selected_event(matching_events[1,])
            }
            }
          }
        }
      }
    }
  }, ignoreNULL = FALSE)
  
  # Make the selection more reliable by using event_id instead of row numbers
  observeEvent(input$as_events_table_rows_selected, {
    if (is.null(as_view_data()) || length(input$as_events_table_rows_selected) == 0) return()
    
    # Get the event ID from the selected row
    selected_row <- input$as_events_table_rows_selected[1]
    events_df <- as_view_data()$events
    
    if (nrow(events_df) >= selected_row) {
      selected_event_id <- events_df$EventID[selected_row]
      
      # Look up complete event data
      gene_id <- input$gene
      gene_events <- selected_gene_as_events()
      matching_events <- gene_events[gene_events$eventID == selected_event_id,]
      
      if (nrow(matching_events) > 0) {
        # Store the selected event
        selected_event(matching_events[1,])
      }
    }
  })
  
  # Show information about the selected event
  output$selected_as_event_info <- renderUI({
    event_data <- selected_event()
    if (is.null(event_data)) {
      return(div(
        h5("Click on an event in the plot or table above to view peptide comparison"),
        style = "color: #888; font-style: italic; text-align: center; margin: 20px 0;"
      ))
    }
    
    # Get event details
    gene_id <- event_data$geneID
    event_id <- event_data$eventID
    
    # Get ALL rows for this event
    gene_events <- selected_gene_as_events()
    event_rows <- gene_events[gene_events$eventID == event_id,]
    
    # Extract spliced-out transcripts
    spliced_out_tx <- c()
    if ("RefTranscripts" %in% names(event_rows)) {
      spliced_out_tx <- unlist(strsplit(as.character(event_rows$RefTranscripts), ", "))
    } else if (!is.null(event_rows$refTx)) {
      spliced_out_tx <- as.character(event_rows$refTx)
    }
    
    # Extract spliced-in transcripts
    spliced_in_tx <- c()
    if ("AltTranscripts" %in% names(event_rows)) {
      spliced_in_tx <- unlist(strsplit(as.character(event_rows$AltTranscripts), ", "))
    } else if (!is.null(event_rows$asTx)) {
      spliced_in_tx <- as.character(event_rows$asTx)
    }
    
    # Remove NAs and ensure unique values
    spliced_out_tx <- unique(spliced_out_tx[!is.na(spliced_out_tx)])
    spliced_in_tx <- unique(spliced_in_tx[!is.na(spliced_in_tx)])
    
    tagList(
      h4(paste0("Selected Event: ", event_data$eventID, " (", event_data$AS_type, ")")),
      div(
        style = "display: flex; justify-content: space-between; margin-bottom: 15px;",
        div(
          h5("Spliced-out Transcripts:"),
          tags$ul(
            lapply(spliced_out_tx, function(tx) {
              tags$li(tx)
            })
          ),
          style = "flex: 1; margin-right: 20px;"
        ),
        div(
          h5("Spliced-in Transcripts:"),
          tags$ul(
            lapply(spliced_in_tx, function(tx) {
              tags$li(tx)
            })
          ),
          style = "flex: 1;"
        )
      )
    )
  })
  
  
  # Create peptide comparison data and visualization
  as_peptide_comparison_data <- reactive({
    # Requirements: selected event, protease, and miscleavage type
    req(selected_event(), input$protease, input$miscleavage_type)
    
    event_data <- selected_event()
    protease <- input$protease
    
    # Get reference and alternative transcripts
    spliced_out_tx <- unique(as.character(event_data$refTx))
    spliced_in_tx <- unique(as.character(event_data$asTx))
    
    # Create standardized vis_data structure using core data module
    vis_data <- core_data_module$create_vis_data_structure()
    
    # Get peptides for each transcript
    spliced_out_peptides_list <- lapply(spliced_out_tx, function(tx) {
      peptides <- get_transcript_peptides_for_comparison(tx, vis_data, protease)
      if (!is.null(peptides)) {
        return(data.frame(
          transcript = tx,
          chromosome = as.character(seqnames(peptides)),
          start = start(peptides),
          end = end(peptides),
          peptide = peptides$peptide,
          type = "Spliced-out",
          stringsAsFactors = FALSE
        ))
      }
      return(NULL)
    })
    
    spliced_in_peptides_list <- lapply(spliced_in_tx, function(tx) {
      peptides <- get_transcript_peptides_for_comparison(tx, vis_data, protease)
      if (!is.null(peptides)) {
        return(data.frame(
          transcript = tx,
          chromosome = as.character(seqnames(peptides)),
          start = start(peptides),
          end = end(peptides),
          peptide = peptides$peptide,
          type = "Spliced-in",
          stringsAsFactors = FALSE
        ))
      }
      return(NULL)
    })
    
    # Combine all peptides
    all_peptides_list <- c(spliced_out_peptides_list, spliced_in_peptides_list)
    all_peptides_list <- all_peptides_list[!sapply(all_peptides_list, is.null)]
    
    if (length(all_peptides_list) == 0) {
      return(NULL)
    }
    
    # Combine all peptides
    all_peptides <- data.table::rbindlist(all_peptides_list)
    
    # Extract all AS ranges (AS_range, AS_range1, AS_range2) and detect peptide overlaps
    # Create a named list to store all ranges
    all_ranges <- list()
    
    # Check and extract each range type
    range_cols <- c("AS_range", "AS_range1", "AS_range2")
    for (col in range_cols) {
      if (col %in% names(event_data) && !is.null(event_data[[col]]) && length(event_data[[col]][[1]]) > 0) {
        # Store the range with its name
        all_ranges[[col]] <- event_data[[col]][[1]]
      }
    }
    
    # Create a combined GRanges object for overlap detection
    combined_as_range <- NULL
    if (length(all_ranges) > 0) {
      # Unlist all ranges into a single GRanges object
      all_ranges_gr <- unlist(GRangesList(all_ranges))
      
      if (length(all_ranges_gr) > 0) {
        # Get the chromosomes of all ranges
        chromosomes <- as.character(unique(GenomicRanges::seqnames(all_ranges_gr)))
        
        # Process each chromosome separately to create combined ranges
        for (chr in chromosomes) {
          chr_ranges <- all_ranges_gr[GenomicRanges::seqnames(all_ranges_gr) == chr]
          
          if (length(chr_ranges) > 0) {
            # Get min start and max end to create a combined range
            min_start <- min(GenomicRanges::start(chr_ranges))
            max_end <- max(GenomicRanges::end(chr_ranges))
            
            # Create a combined range for this chromosome
            chr_combined <- GenomicRanges::GRanges(
              seqnames = chr,
              ranges = IRanges::IRanges(start = min_start, end = max_end)
            )
            
            # Add to the combined range
            if (is.null(combined_as_range)) {
              combined_as_range <- chr_combined
            } else {
              combined_as_range <- c(combined_as_range, chr_combined)
            }
          }
        }
      }
    }

    # Initialize event_overlap column
    all_peptides$event_overlap <- FALSE
    
    if (!is.null(combined_as_range) && length(combined_as_range) > 0) {
      # Create GRanges for all peptides at once
      peptide_ranges <- GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(combined_as_range)[1],
        ranges = IRanges::IRanges(
          start = all_peptides$start,
          end = all_peptides$end
        )
      )
      
      # Find overlaps between peptide ranges and the COMBINED AS region
      overlaps <- GenomicRanges::findOverlaps(peptide_ranges, combined_as_range)
      
      # Mark peptides that overlap with the combined AS region
      if (length(overlaps) > 0) {
        # Extract query hits (peptide indices) from overlaps
        overlapping_peptides <- unique(S4Vectors::queryHits(overlaps))
        all_peptides$event_overlap[overlapping_peptides] <- TRUE
      }
      
      # Also check the individual overlaps with each original AS region
      # This ensures we don't miss any overlaps due to chromosome differences
      all_ranges_gr <- unlist(GRangesList(all_ranges))
      if (length(all_ranges_gr) > 0) {
        for (i in 1:nrow(all_peptides)) {
          if (!all_peptides$event_overlap[i]) {
            pep_range <- GenomicRanges::GRanges(
              seqnames = GenomicRanges::seqnames(all_ranges_gr)[1],
              ranges = IRanges::IRanges(
                start = all_peptides$start[i],
                end = all_peptides$end[i]
              )
            )
            
            # Check against each AS range individually
            if (length(GenomicRanges::findOverlaps(pep_range, all_ranges_gr)) > 0) {
              all_peptides$event_overlap[i] <- TRUE
            }
          }
        }
      }
    }
    
    # Propagate AS flag to all instances of the same peptide sequence
    unique_peptides <- unique(all_peptides$peptide)
    for (pep in unique_peptides) {
      if (any(all_peptides$event_overlap[all_peptides$peptide == pep])) {
        all_peptides$event_overlap[all_peptides$peptide == pep] <- TRUE
      }
    }
    
    # Create set difference categories for AS-affected peptides
    # Initialize a new category column
    all_peptides$set_category <- "Not AS-affected"
    
    # Get all unique peptide sequences
    unique_peptide_seqs <- unique(all_peptides$peptide)
    
    # For each unique peptide sequence
    for (pep_seq in unique_peptide_seqs) {
      # Get all rows for this peptide
      pep_rows <- all_peptides[all_peptides$peptide == pep_seq, ]
      
      # Only process if it's AS-affected
      if (any(pep_rows$event_overlap)) {
        # Check if it appears in spliced-out transcripts
        has_spliced_out <- any(pep_rows$type == "Spliced-out")
        
        # Check if it appears in spliced-in transcripts
        has_spliced_in <- any(pep_rows$type == "Spliced-in")
        
        # Assign category based on presence in transcript types
        category <- if (has_spliced_out && has_spliced_in) {
          "Shared AS-affected"
        } else if (has_spliced_out) {
          "Spliced-out unique AS-affected"
        } else if (has_spliced_in) {
          "Spliced-in unique AS-affected"
        } else {
          "Not AS-affected" # Fallback, shouldn't happen
        }
        
        # Update all instances of this peptide
        all_peptides$set_category[all_peptides$peptide == pep_seq] <- category
      }
    }
    
    # Calculate gene boundaries for plot
    padding <- 5000
    gene_start <- max(1, min(all_peptides$start) - padding)
    gene_end <- max(all_peptides$end) + padding
    
    # Add hover text information to peptides
    all_peptides$hover_text <- paste0(
      "Peptide: ", all_peptides$peptide,
      "<br>Position: ", all_peptides$start, "-", all_peptides$end,
      "<br>Chromosome: ", all_peptides$chromosome,
      "<br>Transcript: ", all_peptides$transcript,
      "<br>Type: ", all_peptides$type,
      "<br>AS-affected: ", ifelse(all_peptides$event_overlap, "Yes", "No"),
      "<br>Category: ", all_peptides$set_category
    )
    
    # Create summary statistics for peptides
    peptide_summary <- data.frame(
      category = character(),
      count = numeric(),
      stringsAsFactors = FALSE
    )

    # Count based on set categories
    as_categories <- table(all_peptides$set_category)
    
    # Convert table to data frame for display
    for (cat in names(as_categories)) {
      peptide_summary <- rbind(peptide_summary, data.frame(
        category = cat,
        count = as_categories[cat],
        stringsAsFactors = FALSE
      ))
    }
    
    # Return the data
    return(list(
      peptides = all_peptides,
      gene_start = gene_start,
      gene_end = gene_end,
      event_region = all_ranges,  # Return the named list of ranges
      transcript_list = c(spliced_out_tx, spliced_in_tx),
      summary = peptide_summary
    ))
  })
  
  # Render peptide comparison visualization
  output$peptide_comparison_plot <- renderPlotly({
    # Get the peptides that are filtered for the selected pair
    peptides_to_plot <- filtered_peptides()

    # Get overall event context (gene boundaries, event regions)
    event_context_data <- as_peptide_comparison_data()
    
    # If no data or context, show a message
    if (is.null(peptides_to_plot) || is.null(event_context_data) || nrow(peptides_to_plot) == 0) {
      p <- plotly::plot_ly() %>%
        plotly::add_annotations(
          x = 0.5, y = 0.5,
          text = clean_hover_text(paste0("No mapped peptides found for selected pair using ", input$protease, " enzyme, or event data is missing.")),
          showarrow = FALSE,
          font = list(size = 16)
        ) %>%
        plotly::layout(
          xaxis = list(range = c(0, 1), showticklabels = FALSE),
          yaxis = list(range = c(0, 1), showticklabels = FALSE)
        )
      return(p)
    }
    
    # Determine the specific pair of transcripts being compared
    req(input$selected_pair)
    selected_pair_object <- transcript_pairs()[[input$selected_pair]]
    if (is.null(selected_pair_object)) {
      p <- plotly::plot_ly() %>%
        plotly::add_annotations(
          x = 0.5, y = 0.5, text = "Selected transcript pair not found.",
          showarrow = FALSE, font = list(size = 16)
        ) %>%
        plotly::layout(
          xaxis = list(range = c(0, 1), showticklabels = FALSE),
          yaxis = list(range = c(0, 1), showticklabels = FALSE)
        )
      return(p)
    }
    
    ref_tx_name <- selected_pair_object["ref"]
    alt_tx_name <- selected_pair_object["alt"]
    current_transcript_list <- c(ref_tx_name, alt_tx_name)
    
    # Extract other necessary data from event_context_data
    gene_start <- event_context_data$gene_start
    gene_end <- event_context_data$gene_end
    event_region <- event_context_data$event_region # This is from the overall event

    # Create transcript data frame for plot based on the selected pair
    transcript_df <- data.frame(
      transcript = current_transcript_list,
      y_position = seq_along(current_transcript_list), # Will be 1 and 2
      stringsAsFactors = FALSE
    )
    
    # Assign type based on the order in the pair (1st is "Spliced-out" context, 2nd is "Spliced-in")
    transcript_df$type <- ifelse(transcript_df$transcript == ref_tx_name, "Spliced-out", "Spliced-in")
    
    # Get event type for better labeling (from the overall selected event)
    event_type <- selected_event()$AS_type
    event_type_full <- switch(event_type,
      "SE" = "Skipped Exon",
      "RI" = "Retained Intron",
      "MX" = "Mutually Exclusive Exons",
      "A3" = "Alternative 3' Splice Site",
      "A5" = "Alternative 5' Splice Site",
      event_type
    )
    
    # Prepare exon data for the two selected transcripts using fast GTF cache
    gene_id <- selected_event()$geneID
    exons_by_transcript_pair <- list()
    
    # Initialize gene_details for compatibility
    gene_details <- NULL
    
    # Use lightning-fast GTF cache if available
    if (dir.exists("data/gtf_cache")) {
      gtf_data <- load_gtf_visualization_data(gene_id)
      if (gtf_data$success) {
        # Create exons_result compatible structure from fast GTF data
        exons_result <- list(
          success = TRUE,
          exons = gtf_data$exons_by_transcript,
          cds = gtf_data$cds_by_transcript
        )
        exons_by_transcript_pair <- exons_result$exons
        
        # Create gene_details compatible structure for downstream use
        gene_details <- list(
          success = TRUE,
          chromosome = gtf_data$chromosome,
          strand = gtf_data$strand
        )
      } else {
        # Fast GTF cache failed, use fallback
        gene_details <- load_gene_details(gene_id)
        if (gene_details$success) {
          exons_result <- load_transcript_exons(gene_details, current_transcript_list)
          if (exons_result$success) {
            exons_by_transcript_pair <- exons_result$exons
          }
        }
      }
    } else {
      # No GTF cache, use original method
      gene_details <- load_gene_details(gene_id)
      if (gene_details$success) {
        exons_result <- load_transcript_exons(gene_details, current_transcript_list)
        if (exons_result$success) {
          exons_by_transcript_pair <- exons_result$exons
        }
      }
    }
    
    # Add compression support for peptide comparison - respect user preference
    use_compression <- !is.null(input$intron_scale) && input$intron_scale == "compressed"
    compression_map <- NULL
    intron_markers <- NULL
    
    if (use_compression && length(exons_by_transcript_pair) > 0) {
      # Create compression map from the transcript pair exons
      compression_map <- create_compression_map(exons_by_transcript_pair)
      
      # Create intron markers for zigzag lines for each transcript at their specific Y positions
      intron_markers <- data.frame(x = numeric(), y = numeric(), group = character(), transcript = character(), stringsAsFactors = FALSE)
      for(tx_name in names(exons_by_transcript_pair)) {
        tx_y_pos <- ifelse(tx_name == current_transcript_list[1], 1, 2)  # Y position for each transcript
        tx_markers <- create_intron_markers_for_transcripts(list(tx_name = exons_by_transcript_pair[[tx_name]]), compression_map, y_position = tx_y_pos)
        if(nrow(tx_markers) > 0) {
          intron_markers <- rbind(intron_markers, tx_markers)
        }
      }
    }
    
    # Create exon dataframe for plotting (for the pair)
    exon_df_pair <- data.frame(
      transcript = character(),
      y_position = numeric(),
      start = numeric(),
      end = numeric(),
      stringsAsFactors = FALSE
    )
    
    if (length(exons_by_transcript_pair) > 0) {
      for (tx_name_in_pair in names(exons_by_transcript_pair)) {
        if (tx_name_in_pair %in% current_transcript_list && length(exons_by_transcript_pair[[tx_name_in_pair]]) > 0) {
          tx_exons <- exons_by_transcript_pair[[tx_name_in_pair]]
          y_pos <- transcript_df$y_position[transcript_df$transcript == tx_name_in_pair] # Get y_pos from transcript_df (1 or 2)
          
          for (j in seq_along(tx_exons)) {
            exon_df_pair <- rbind(exon_df_pair, data.frame(
              transcript = tx_name_in_pair,
              y_position = y_pos,
              start = start(tx_exons[j]),
              end = end(tx_exons[j]),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
    
    # Apply compression to exon coordinates if compression is enabled
    if (use_compression && !is.null(compression_map) && nrow(exon_df_pair) > 0) {
      exon_df_pair$compressed_start <- sapply(exon_df_pair$start, compression_map$compress)
      exon_df_pair$compressed_end <- sapply(exon_df_pair$end, compression_map$compress)
      
      # Use compressed coordinates for plotting
      exon_df_pair$plot_start <- exon_df_pair$compressed_start
      exon_df_pair$plot_end <- exon_df_pair$compressed_end
    } else {
      # Use original coordinates for plotting
      exon_df_pair$plot_start <- exon_df_pair$start
      exon_df_pair$plot_end <- exon_df_pair$end
    }

    # Create a new plotly object
    p <- plotly::plot_ly()
    
    # Define colors for event types - highly distinct colors
    event_colors <- list(
      "SE" = "#FF0000",  # Bright red
      "RI" = "#00FF00",  # Bright green  
      "MX" = "#0000FF",  # Bright blue
      "A3" = "#235e24",  # Bright yellow
      "A5" = "#FF00FF"   # Bright magenta
    )
    
    # Define colors for AS status
    as_status_colors <- list(
      "Spliced-out unique AS-affected" = "#9B59B6",
      "Spliced-in unique AS-affected" = "#F39C12",
      "Shared AS-affected" = "#2ECC71",
      "Not AS-affected" = "#CCCCCC"
    )
    
    # Add event regions (context from the overall event)
    has_valid_event_region <- FALSE
    if (!is.null(event_region)) {
      if (is(event_region, "GRanges") && length(event_region) > 0) {
        has_valid_event_region <- TRUE
      } else if (is.list(event_region) && length(event_region) > 0) {
        for (range_item in event_region) {
          if (is(range_item, "GRanges") && length(range_item) > 0) {
            has_valid_event_region <- TRUE
            break
          }
        }
      }
    }
    
    if (has_valid_event_region) {
      range_colors <- list("AS_range" = event_colors[[event_type]], "AS_range1" = "#D35400", "AS_range2" = "#16A085")
      range_labels <- list(
        "AS_range" = paste0(event_type, ": ", event_type_full),
        "AS_range1" = paste0(event_type, " Region 1"),
        "AS_range2" = paste0(event_type, " Region 2")
      )
      
      processed_event_region <- if (is(event_region, "GRanges")) list("AS_range" = event_region) else event_region
      
      for (range_name in names(processed_event_region)) {
        range_obj <- processed_event_region[[range_name]]
        if (!is(range_obj, "GRanges") || length(range_obj) == 0) next
        
        color <- range_colors[[range_name]] %||% event_colors[[event_type]] %||% "grey"
        label <- range_labels[[range_name]] %||% paste0(event_type, " ", gsub("AS_", "", range_name))
        
        # Apply compression to AS event coordinates if compression is enabled
        range_start <- start(range_obj)[1]
        range_end <- end(range_obj)[1]
        if (use_compression && !is.null(compression_map)) {
          range_start <- compression_map$compress(range_start)
          range_end <- compression_map$compress(range_end)
        }
        
        p <- p %>% plotly::add_trace(
          type = "scatter", mode = "lines",
          x = c(range_start, range_end, range_end, range_start, range_start),
          y = c(0.5, 0.5, length(current_transcript_list) + 0.5, length(current_transcript_list) + 0.5, 0.5), # Adjusted for 2 transcripts
          fill = "toself", fillcolor = paste0(color, "40"),
          line = list(color = color, dash = ifelse(range_name == "AS_range", "solid", ifelse(range_name == "AS_range1", "dash", "dot")), width = 1),
          marker = list(opacity = 0), showlegend = TRUE, name = label,
          hovertemplate = paste0("Event Type: ", event_type, " (", event_type_full, ")", "<br>Range: ", gsub("AS_", "", range_name), "<br>Position: ", start(range_obj)[1], "-", end(range_obj)[1], "<extra></extra>"),
          hoverinfo = "none"
        )
      }
    }
    
    # Add AS position indicator (green vertical line)
    if (has_valid_event_region && !is.null(event_region)) {
      # Get the AS position from the event region
      as_position <- NULL
      
      # Try to get position from the first valid range
      if (is(event_region, "GRanges") && length(event_region) > 0) {
        # Use midpoint of the AS range
        as_position <- (start(event_region)[1] + end(event_region)[1]) / 2
      } else if (is.list(event_region)) {
        # Check for AS_range in the list
        if (!is.null(event_region[["AS_range"]]) && is(event_region[["AS_range"]], "GRanges") && length(event_region[["AS_range"]]) > 0) {
          as_position <- (start(event_region[["AS_range"]])[1] + end(event_region[["AS_range"]])[1]) / 2
        } else {
          # Try to find any valid range
          for (range_item in event_region) {
            if (is(range_item, "GRanges") && length(range_item) > 0) {
              as_position <- (start(range_item)[1] + end(range_item)[1]) / 2
              break
            }
          }
        }
      }
      
      if (!is.null(as_position)) {
        # Apply compression if needed
        plot_as_position <- as_position
        if (use_compression && !is.null(compression_map)) {
          plot_as_position <- compression_map$compress(as_position)
        }
        
        # Add vertical line for AS position
        p <- p %>% plotly::add_trace(
          type = "scatter", mode = "lines",
          x = c(plot_as_position, plot_as_position),
          y = c(0.3, length(current_transcript_list) + 0.7),
          line = list(color = "#27AE60", width = 3, dash = "solid"),
          showlegend = TRUE,
          name = "AS Position",
          hovertemplate = paste0("AS Position<br>Genomic: ", round(as_position), "<extra></extra>"),
          hoverinfo = "none"
        )
      }
    }
    
    # Add transcript lines for the selected pair
    for (i in 1:nrow(transcript_df)) {
      tx_display_type <- transcript_df$type[i] # "Spliced-out" or "Spliced-in" for the pair
      tx_color <- ifelse(tx_display_type == "Spliced-out", "#2C3E50", "#E74C3C") # Colors for pair context
      
      p <- p %>% plotly::add_trace(
        type = "scatter", x = c(gene_start, gene_end), y = c(transcript_df$y_position[i], transcript_df$y_position[i]),
        mode = "lines", line = list(color = tx_color, width = 1),
        showlegend = TRUE, # Show legend for both transcripts in the pair
        name = paste0(tx_display_type, ": ", transcript_df$transcript[i]), # More descriptive legend
        legendgroup = paste0("tx_pair_type_", tx_display_type),
        hoverinfo = "text", text = clean_hover_text(paste0("Transcript: ", transcript_df$transcript[i], " (", tx_display_type, " of pair)"))
      )
    }
    
    # Add exons for the selected pair
    if (nrow(exon_df_pair) > 0) {
      for (i in 1:nrow(exon_df_pair)) {
        p <- p %>% plotly::add_trace(
          type = "scatter", mode = "lines",
          x = c(exon_df_pair$plot_start[i], exon_df_pair$plot_end[i], exon_df_pair$plot_end[i], exon_df_pair$plot_start[i], exon_df_pair$plot_start[i]),
          y = c(exon_df_pair$y_position[i] - 0.3, exon_df_pair$y_position[i] - 0.3, exon_df_pair$y_position[i] + 0.3, exon_df_pair$y_position[i] + 0.3, exon_df_pair$y_position[i] - 0.3),
          fill = "toself", fillcolor = "grey70", line = list(color = "black", width = 1),
          marker = list(opacity = 0), showlegend = FALSE,
          hovertemplate = paste0("Exon<br>Transcript: ", exon_df_pair$transcript[i], "<br>Position: ", exon_df_pair$start[i], "-", exon_df_pair$end[i], "<extra></extra>"),
          hoverinfo = "none"
        )
      }
    }
    
    # Add intron markers if compression is enabled - only for displayed transcripts
    if (use_compression && !is.null(intron_markers) && nrow(intron_markers) > 0) {
      # Only show zigzag markers for the transcripts actually being displayed
      displayed_transcripts <- transcript_df$transcript
      filtered_markers <- intron_markers[intron_markers$transcript %in% displayed_transcripts, ]
      
      if (nrow(filtered_markers) > 0) {
        # Group intron markers by group and create discrete zigzag segments
        unique_groups <- unique(filtered_markers$group)
        for (group in unique_groups) {
          group_data <- filtered_markers[filtered_markers$group == group, ]
          if (nrow(group_data) > 0) {
            # Create discrete zigzag segment (not continuous line)
            p <- p %>% plotly::add_trace(
              type = "scatter", mode = "lines",
              x = group_data$x, y = group_data$y,
              line = list(color = "grey50", width = 1),
              connectgaps = FALSE,  # Don't connect to other segments
              showlegend = FALSE,
              hoverinfo = "none"
            )
          }
        }
      }
    }
    
    # Add peptides (from peptides_to_plot)
    if (nrow(peptides_to_plot) > 0) {
      # Apply compression to peptide coordinates if compression is enabled
      if (use_compression && !is.null(compression_map)) {
        peptides_to_plot$compressed_start <- sapply(peptides_to_plot$start, compression_map$compress)
        peptides_to_plot$compressed_end <- sapply(peptides_to_plot$end, compression_map$compress)
        
        # Use compressed coordinates for plotting
        peptides_to_plot$plot_start <- peptides_to_plot$compressed_start
        peptides_to_plot$plot_end <- peptides_to_plot$compressed_end
      } else {
        # Use original coordinates for plotting
        peptides_to_plot$plot_start <- peptides_to_plot$start
        peptides_to_plot$plot_end <- peptides_to_plot$end
      }
      p <- p %>% plotly::add_trace(
        type = "scatter", mode = "lines", x = c(0), y = c(0),
        line = list(color = "rgba(0,0,0,0)", width = 0), showlegend = TRUE,
        name = "<b>AS Status:</b>", hoverinfo = "none", visible = TRUE
      )
      
      categories_in_plot <- unique(peptides_to_plot$set_category)
      for (category in c("Spliced-out unique AS-affected", "Spliced-in unique AS-affected", "Shared AS-affected", "Not AS-affected")) {
        if (category %in% categories_in_plot) {
          p <- p %>% plotly::add_trace(
            type = "scatter", mode = "markers", x = c(0), y = c(0),
            marker = list(color = as_status_colors[[category]], size = 10), showlegend = TRUE,
            name = switch(category,
                         "Spliced-out unique AS-affected" = "AS-affected only in Spliced-out Tx",
                         "Spliced-in unique AS-affected" = "AS-affected only in Spliced-in Tx",
                         "Shared AS-affected" = "AS-affected in both Txs",
                         "Not AS-affected" = "Not AS-affected"),
            legendgroup = paste0("category_", category), hoverinfo = "none", visible = TRUE
          )
        }
      }
      
      for (i in 1:nrow(peptides_to_plot)) {
        category <- peptides_to_plot$set_category[i]
        color <- as_status_colors[[category]] %||% "#CCCCCC"
        
        # Map peptide's transcript to its y_position in the current pair (1 or 2)
        tx_name_for_peptide <- peptides_to_plot$transcript[i]
        tx_pos_for_peptide <- transcript_df$y_position[transcript_df$transcript == tx_name_for_peptide]
        
        if (length(tx_pos_for_peptide) == 0 || is.na(tx_pos_for_peptide)) next # Should not happen if peptides_to_plot is correct

        p <- p %>% plotly::add_trace(
          type = "scatter", mode = "lines",
          x = c(peptides_to_plot$plot_start[i], peptides_to_plot$plot_end[i], peptides_to_plot$plot_end[i], peptides_to_plot$plot_start[i], peptides_to_plot$plot_start[i]),
          y = c(tx_pos_for_peptide - 0.2, tx_pos_for_peptide - 0.2, tx_pos_for_peptide + 0.2, tx_pos_for_peptide + 0.2, tx_pos_for_peptide - 0.2),
          fill = "toself", fillcolor = color, line = list(color = "black", width = 1),
          marker = list(opacity = 0), legendgroup = paste0("category_", category), showlegend = FALSE,
          hovertemplate = paste0(clean_hover_text(peptides_to_plot$hover_text[i]), "<extra></extra>"),
          hoverinfo = "none"
        )
      }
    }
    
    # Get chromosome and strand information
    chromosome <- if(gene_details$success) gene_details$chromosome else "Unknown"
    
    # Extract strand from event context (plots are always shown 5' to 3')
    strand_display <- ""
    
    # Set the layout with compression-aware axis
    xaxis_config <- if (use_compression && !is.null(compression_map)) {
      # Use compression: axis shows real coordinates positioned at compressed locations
      axis_breaks <- compression_map$coords$compressed_start
      axis_labels <- as.character(compression_map$coords$original_start)
      compressed_range <- range(c(compression_map$coords$compressed_start, compression_map$coords$compressed_end))
      
      list(
        title = paste0("Genomic Position (compressed view) - chromosome ", chromosome, ")"),
        tickvals = axis_breaks,
        ticktext = axis_labels,
        range = compressed_range,
        zeroline = FALSE
      )
    } else {
      # Use original coordinates
      list(
        title = paste0("Genomic Position (chromosome ", chromosome, ")"),
        range = c(gene_start, gene_end),
        zeroline = FALSE
      )
    }
    
    p <- p %>% plotly::layout(
      title = paste0('Peptide Comparison for Event: ', selected_event()$eventID, ' - Pair: ', input$selected_pair),
      xaxis = xaxis_config,
      yaxis = list(
        title = "",
        range = c(0.5, length(current_transcript_list) + 0.5), # Adjusted for 2 transcripts
        tickvals = transcript_df$y_position, # seq_along(current_transcript_list) which is 1, 2
        ticktext = current_transcript_list, # Names of the two transcripts in the pair
        zeroline = FALSE
      ),
      legend = list(orientation = "v", x = 1.02, y = 0.5, title = list(text = "Legend:"), itemsizing = "constant", tracegroupgap = 10, bgcolor = "rgba(255, 255, 255, 0.8)", bordercolor = "rgba(0, 0, 0, 0.2)", borderwidth = 1),
      margin = list(t = 50, b = 50, l = 100, r = 250) # Increased right margin for longer legend text
    ) %>%
    plotly::config(displayModeBar = TRUE, modeBarButtonsToRemove = c('lasso2d', 'select2d'))
    
    return(p %>% clean_plotly_hover())
  })
  
  # Render peptide comparison table
  output$peptide_comparison_table <- DT::renderDataTable({
    peptides <- filtered_peptides()
    
    if(nrow(peptides) == 0) {
      return(data.frame(
        Message = "No peptides found for this transcript pair"
      ))
    }
    
    # Group peptides to eliminate duplicates and merge transcript information
    display_df <- peptides %>%
      group_by(peptide) %>%
      summarise(
        chromosome = first(chromosome),
        all_positions = paste(paste(start, end, sep = "-"), collapse = ", "),
        transcript = paste(unique(transcript), collapse = ", "),
        type = paste(unique(type), collapse = ", "),
        is_junction_spanning = n() > 1,
        event_overlap = any(event_overlap),
        set_category = first(set_category),
        .groups = "drop"
      )
    
    # Rename columns for display
    colnames(display_df) <- c("Peptide", "Chromosome", "All Positions", "Transcript", "Splice Source", "Junction-spanning", "AS-affected", "Category")
    
    # Create custom category labels for display
    display_df$Category <- sapply(display_df$Category, function(cat) {
      switch(cat,
        "Spliced-out unique AS-affected" = "Only in Spliced-out",
        "Spliced-in unique AS-affected" = "Only in Spliced-in",
        "Shared AS-affected" = "In both transcripts",
        "Not AS-affected" = "Not AS-affected"
      )
    })
    
    DT::datatable(
      display_df,
      selection = 'single',
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      rownames = FALSE
    ) %>%
    DT::formatStyle(
      'Splice Source',
      backgroundColor = styleEqual(
        c("Spliced-out", "Spliced-in"),
        c("#3498DB", "#E74C3C")
      ),
      color = styleEqual(
        c("Spliced-out", "Spliced-in"),
        c("white", "white")
      )
    )
  })
  
  # Store selected peptide from table
  selected_peptide <- reactiveVal(NULL)
  
  # Observer for peptide table row selection
  observeEvent(input$peptide_comparison_table_rows_selected, {
    req(filtered_peptides())
    if (length(input$peptide_comparison_table_rows_selected) > 0) {
      selected_row <- input$peptide_comparison_table_rows_selected[1]
      peptides <- filtered_peptides()
      if (selected_row <= nrow(peptides)) {
        selected_peptide(peptides[selected_row, ])
      }
    } else {
      selected_peptide(NULL)
    }
  })
  
  # Render summary plot
  output$peptide_summary_plot <- renderPlot({
    req(as_peptide_comparison_data())
    
    # Get comparison data
    comparison_data <- as_peptide_comparison_data()
    
    # Get the summary data
    summary_data <- comparison_data$summary
    
    # Display only the AS-affected categories for clarity
    summary_data <- summary_data[grepl("AS-affected", summary_data$category), ]
    
    # Create better labels for display
    summary_data$display_category <- sapply(summary_data$category, function(cat) {
      switch(cat,
        "Spliced-out unique AS-affected" = "Only in Spliced-out",
        "Spliced-in unique AS-affected" = "Only in Spliced-in",
        "Shared AS-affected" = "In both transcripts",
        cat  # Default fallback
      )
    })
    
    # Create summary plot with set difference categories
    ggplot(summary_data, aes(x = display_category, y = count, fill = display_category)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(
        values = c(
          "Only in Spliced-out" = "#9B59B6",  # Purple
          "Only in Spliced-in" = "#F39C12",   # Orange
          "In both transcripts" = "#2ECC71"   # Green
        )
      ) +
      labs(
        title = "AS-affected Peptide Classification",
        x = "",
        y = "Count"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  })
  
  # Render summary table
  output$peptide_summary_table <- DT::renderDataTable({
    req(as_peptide_comparison_data())
    
    # Get comparison data
    comparison_data <- as_peptide_comparison_data()
    
    # Get the summary data
    summary_data <- comparison_data$summary
    
    # Create better labels for display
    display_df <- data.frame(
      Category = sapply(summary_data$category, function(cat) {
        switch(cat,
          "Spliced-out unique AS-affected" = "AS-affected only in Spliced-out",
          "Spliced-in unique AS-affected" = "AS-affected only in Spliced-in",
          "Shared AS-affected" = "AS-affected in both transcripts",
          "Not AS-affected" = "Not AS-affected",
          cat  # Default fallback
        )
      }),
      Count = summary_data$count,
      stringsAsFactors = FALSE
    )
    
    # Order the categories in a logical way
    display_order <- c(
      "AS-affected only in Spliced-out",
      "AS-affected only in Spliced-in",
      "AS-affected in both transcripts",
      "Not AS-affected"
    )
    
    # Match existing categories to order
    existing_cats <- intersect(display_order, display_df$Category)
    display_df <- display_df[match(existing_cats, display_df$Category), ]
    
    DT::datatable(
      display_df,
      options = list(
        dom = 't',
        ordering = FALSE
      ),
      rownames = FALSE
    ) %>%
    DT::formatStyle(
      'Category',
      backgroundColor = styleEqual(
        c(
          "AS-affected only in Spliced-out",
          "AS-affected only in Spliced-in",
          "AS-affected in both transcripts", 
          "Not AS-affected"
        ),
        c(
          "#9B59B6",  # Purple
          "#F39C12",  # Orange
          "#2ECC71",  # Green
          "#EEEEEE"   # Light gray
        )
      ),
      color = styleEqual(
        c(
          "AS-affected only in Spliced-out",
          "AS-affected only in Spliced-in",
          "AS-affected in both transcripts",
          "Not AS-affected"
        ),
        c("white", "white", "white", "black")
      )
    )
  })

  # Render peptide comparison summary
  output$peptide_comparison_summary <- renderText({
    peptides <- filtered_peptides()
    comparison_data <- as_peptide_comparison_data()
    
    if (is.null(peptides) || nrow(peptides) == 0 || is.null(comparison_data)) {
      return("No peptide data available for summary.")
    }
    
    # Get the selected pair information
    req(input$selected_pair)
    pair <- transcript_pairs()[[input$selected_pair]]
    if (is.null(pair)) {
      return("Selected transcript pair not found.")
    }
    
    ref_tx <- pair["ref"]
    alt_tx <- pair["alt"]
    
    # Calculate summary statistics based on UNIQUE peptide sequences
    # Group by peptide sequence to avoid double-counting junction-spanning peptides
    unique_peptides <- peptides %>%
      group_by(peptide) %>%
      summarise(
        is_as_affected = any(event_overlap),
        set_category = first(set_category),
        transcripts = paste(unique(transcript), collapse = ", "),
        .groups = "drop"
      )
    
    total_peptides <- nrow(unique_peptides)
    as_affected <- sum(unique_peptides$is_as_affected)
    not_as_affected <- total_peptides - as_affected
    
    # Count by category based on unique peptides
    ref_unique <- sum(unique_peptides$set_category == "Spliced-out unique AS-affected", na.rm = TRUE)
    alt_unique <- sum(unique_peptides$set_category == "Spliced-in unique AS-affected", na.rm = TRUE)
    shared <- sum(unique_peptides$set_category == "Shared AS-affected", na.rm = TRUE)
    
    # Create summary text
    summary_text <- paste0(
      "PEPTIDE COMPARISON SUMMARY\n",
      "=========================\n\n",
      "Selected Pair: ", input$selected_pair, "\n",
      "Reference Transcript: ", ref_tx, "\n", 
      "Alternative Transcript: ", alt_tx, "\n\n",
      "OVERALL STATISTICS:\n",
      "- Total peptides analyzed: ", total_peptides, "\n",
      "- AS-affected peptides: ", as_affected, " (", round(as_affected/total_peptides*100, 1), "%)\n",
      "- Not AS-affected: ", not_as_affected, " (", round(not_as_affected/total_peptides*100, 1), "%)\n\n",
      "AS-AFFECTED PEPTIDE BREAKDOWN:\n",
      "- Unique to Reference (", ref_tx, "): ", ref_unique, "\n",
      "- Unique to Alternative (", alt_tx, "): ", alt_unique, "\n",
      "- Shared between both: ", shared, "\n\n",
      "INTERPRETATION:\n",
      if (ref_unique > 0) paste0("- ", ref_unique, " peptides are lost in the alternative transcript\n") else "",
      if (alt_unique > 0) paste0("- ", alt_unique, " peptides are gained in the alternative transcript\n") else "",
      if (shared > 0) paste0("- ", shared, " peptides are present in both transcripts but affected by the AS event\n") else "",
      "\nNOTE:\n",
      "- Junction-spanning peptides are counted once per unique sequence\n",
      "- Multiple coordinate ranges for the same peptide are not double-counted\n",
      "\nProtease used: ", input$protease, "\n",
      "Miscleavage setting: ", input$miscleavage_type
    )
    
    return(summary_text)
  })

  # Download handler for peptide data
  output$download_peptide_data <- downloadHandler(
    filename = function() {
      miscleavage_suffix <- switch(input$miscleavage_type,
        "no_miss_cleavage" = "no_miss",
        "upto_two_misscleavage" = "upto_2_miss"
      )
      paste0("peptide_comparison_", selected_event()$eventID, "_", 
            miscleavage_suffix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      # Get peptide data
      comparison_data <- as_peptide_comparison_data()
      peptides <- comparison_data$peptides
      
      # Add miscleavage information to the data
      peptides$miscleavage_type <- input$miscleavage_type
      peptides$protease <- input$protease
      
      # Filter if needed
      if(!input$show_all_peptides) {
        peptides <- peptides[peptides$event_overlap, ]
      }
      
      # Write to CSV
      write.csv(peptides, file, row.names = FALSE)
    }
  )
  
  # Create all possible ref-alt transcript pairs for the selected event
  transcript_pairs <- reactive({
    req(selected_event())
    
    event_data <- selected_event()
    
    # Get the event ID to find all matching rows
    event_id <- event_data$eventID
    gene_id <- event_data$geneID
    
    # Get ALL rows for this event
    gene_events <- selected_gene_as_events()
    event_rows <- gene_events[gene_events$eventID == event_id,]
    
    # Extract spliced-out transcripts - try multiple methods to ensure we get all
    ref_tx <- c()
    
    # Method 1: Try RefTranscripts column
    if ("RefTranscripts" %in% names(event_rows)) {
      for (i in 1:nrow(event_rows)) {
        ref_tx_str <- as.character(event_rows$RefTranscripts[i])
        if (!is.na(ref_tx_str) && ref_tx_str != "") {
          ref_tx <- c(ref_tx, unlist(strsplit(ref_tx_str, ", ")))
        }
      }
    }
    
    # Method 2: Try refTx column directly
    if (length(ref_tx) == 0 || all(is.na(ref_tx))) {
      for (i in 1:nrow(event_rows)) {
        if (!is.null(event_rows$refTx[[i]]) && length(event_rows$refTx[[i]]) > 0) {
          ref_tx <- c(ref_tx, as.character(event_rows$refTx[[i]]))
        }
      }
    }
    
    # Extract spliced-in transcripts with same approach
    alt_tx <- c()
    
    # Method 1: Try AltTranscripts column
    if ("AltTranscripts" %in% names(event_rows)) {
      for (i in 1:nrow(event_rows)) {
        alt_tx_str <- as.character(event_rows$AltTranscripts[i])
        if (!is.na(alt_tx_str) && alt_tx_str != "") {
          alt_tx <- c(alt_tx, unlist(strsplit(alt_tx_str, ", ")))
        }
      }
    }
    
    # Method 2: Try asTx column directly
    if (length(alt_tx) == 0 || all(is.na(alt_tx))) {
      for (i in 1:nrow(event_rows)) {
        if (!is.null(event_rows$asTx[[i]]) && length(event_rows$asTx[[i]]) > 0) {
          alt_tx <- c(alt_tx, as.character(event_rows$asTx[[i]]))
        }
      }
    }
    
    # Ensure unique values and remove any NAs
    ref_tx <- unique(ref_tx[!is.na(ref_tx)])
    alt_tx <- unique(alt_tx[!is.na(alt_tx)])
    
    # Print debug info
    print("ALL Spliced-out transcripts:")
    print(ref_tx)
    print("ALL Spliced-in transcripts:")
    print(alt_tx)
    
    # Generate all possible spliced-out-spliced-in pairs
    pairs <- list()
    pair_labels <- character()
    
    for (r in ref_tx) {
      for (a in alt_tx) {
        pair_labels <- c(pair_labels, paste0("Spliced-out: ", r, " vs Spliced-in: ", a))
        pairs <- c(pairs, list(c(ref = r, alt = a)))
      }
    }
    
    # Print the generated pairs
    print("Generated ALL spliced-out-spliced-in pairs:")
    print(pair_labels)
    
    # Return named list for selectInput
    if (length(pairs) > 0) {
      names(pairs) <- pair_labels
      return(pairs)
    } else {
      return(NULL)
    }
  })

  # Render the comparison pair selector UI
  output$comparison_pair_selector <- renderUI({
    pairs <- transcript_pairs()
    if (is.null(pairs) || length(pairs) == 0) {
      return(div(
        style = "color: #999; font-style: italic;",
        "No transcript pairs available for comparison"
      ))
    }
    
    selectInput("selected_pair", 
               "Choose Transcript Pair to Compare:",
               choices = names(pairs),
               selected = names(pairs)[1],
               width = "100%")
  })

  # Update the pair selection dropdown when an event is selected
  observeEvent(selected_event(), {
    pairs <- transcript_pairs()
    if (!is.null(pairs) && length(pairs) > 0) {
      updateSelectInput(session, "selected_pair", 
                       choices = names(pairs),
                       selected = names(pairs)[1])
    }
  })

  # Filter peptides based on the selected pair
  filtered_peptides <- reactive({
    req(as_peptide_comparison_data(), input$selected_pair)
    
    withProgress(message = 'Loading peptide comparison...', {
      comparison_data <- as_peptide_comparison_data()
      all_peptides <- comparison_data$peptides
      
      incProgress(0.4, detail = "Filtering peptides")
      pair <- transcript_pairs()[[input$selected_pair]]
      if (is.null(pair)) return(all_peptides)
      
      filtered <- all_peptides[all_peptides$transcript %in% c(pair["ref"], pair["alt"]), ]
      
      incProgress(0.6, detail = "Applying AS filter")
      if (!input$show_all_peptides) {
        filtered <- filtered[filtered$event_overlap, ]
      }
      
      return(filtered)
    })
  })

  # Add this to update the selected pair display
  output$selected_pair_display <- renderUI({
    req(input$selected_pair)
    
    # Get the selected pair name
    pair_name <- input$selected_pair
    
    div(
      strong("Currently comparing:"),
      p(pair_name)
    )
  })

  output$compare_button <- renderUI({
    req(input$selected_pair)
    actionButton("compare_peptides", "Compare Peptides",
                icon = icon("chart-bar"),
                style = "margin: 15px 0; color: #fff; background-color: #3c8dbc; border-color: #367fa9;")
  })










  # Reactive data for all-events peptide comparison
  as_all_peptide_comparison_data <- reactive({
    req(input$load_all_events > 0, input$gene, input$protease, input$miscleavage_type)
    withProgress(message = 'Loading all-event peptide comparison...', value = 0, {
      # Step 1: Gathering AS events & transcripts (0-20%)
      incProgress(0.1, detail = 'Gathering AS events...')
      gene_id <- input$gene
      protease <- input$protease
      
      # Check gene_events structure first
      gene_events <- selected_gene_as_events()
      
      # Add debug messages to check structure
      cat("Gene events structure:", paste(names(gene_events), collapse=", "), "\n")
      cat("Number of AS events found:", nrow(gene_events), "\n")
      
      incProgress(0.1, detail = 'Collecting transcripts...')
      # Defensively extract transcripts from gene_events
      all_transcripts <- c()
      if (!is.null(gene_events$refTx)) {
        ref_tx <- unlist(gene_events$refTx)
        ref_tx <- ref_tx[!is.na(ref_tx)]
        all_transcripts <- c(all_transcripts, ref_tx)
      }
      if (!is.null(gene_events$asTx)) {
        as_tx <- unlist(gene_events$asTx)
        as_tx <- as_tx[!is.na(as_tx)]
        all_transcripts <- c(all_transcripts, as_tx)
      }
      
      # Make sure we have unique transcripts
      transcripts <- unique(all_transcripts)
      cat("Number of transcripts found:", length(transcripts), "\n")

      # Step 2: Loading exons (20-40%)
      incProgress(0.1, detail = 'Retrieving gene details...')
      # Use lightning-fast GTF cache if available
      exons_by_tx <- list()
      
      incProgress(0.1, detail = 'Loading transcript exons (lightning-fast)...')
      if (dir.exists("data/gtf_cache")) {
        gtf_data <- load_gtf_visualization_data(gene_id)
        if (gtf_data$success) {
          exons_by_tx <- gtf_data$exons_by_transcript
        } else {
          # Fast GTF cache failed, use fallback
          gene_details <- load_gene_details(gene_id)
          if (gene_details$success) {
            ex <- load_transcript_exons(gene_details, transcripts)
            if (ex$success) exons_by_tx <- ex$exons
          }
        }
      } else {
        # No GTF cache, use original method
        gene_details <- load_gene_details(gene_id)
        if (gene_details$success) {
          ex <- load_transcript_exons(gene_details, transcripts)
          if (ex$success) exons_by_tx <- ex$exons
        }
      }
      cat("Number of transcripts with exons:", length(exons_by_tx), "\n")
      
      # Add compression support for all events comparison - respect user preference
      use_compression_all <- !is.null(input$intron_scale) && input$intron_scale == "compressed"
      compression_map_all <- NULL
      intron_markers_all <- NULL
      
      if (use_compression_all && length(exons_by_tx) > 0) {
        # Create compression map from all transcript exons
        compression_map_all <- create_compression_map(exons_by_tx)
        
        # Create intron markers for zigzag lines - defer this until after transcript_df is created
        intron_markers_all <- NULL
      }

      # Step 3: Fetching peptides (40-60%)
      # Create a modified processed_data structure for visualization functions
      vis_data <- list(
        genes = processed_data()$genes,
        gene_symbols = processed_data()$gene_symbols,
        gene_lookup = processed_data()$gene_lookup,
        proteases = processed_data()$proteases,
        original_peptides = selected_gene_peptides()
      )
      
      incProgress(0.1, detail = 'Fetching transcript peptides...')
      peptide_list <- lapply(seq_along(transcripts), function(i) {
        tx <- transcripts[i]
        gr <- get_transcript_peptides_for_comparison(tx, vis_data, protease)
        if (is.null(gr)) return(NULL)
        data.frame(transcript = tx, y = i, start = start(gr), end = end(gr), peptide = gr$peptide, stringsAsFactors = FALSE)
      })
      
      incProgress(0.1, detail = 'Combining peptide data...')
      # Filter out NULL results before rbind
      peptide_list <- peptide_list[!sapply(peptide_list, is.null)]
      peptide_df <- if (length(peptide_list) > 0) data.table::rbindlist(peptide_list) else data.frame()
      cat("Number of peptides found:", nrow(peptide_df), "\n")
      
      # Apply compression to all events peptide coordinates  
      if (use_compression_all && !is.null(compression_map_all) && nrow(peptide_df) > 0) {
        peptide_df$compressed_start <- sapply(peptide_df$start, compression_map_all$compress)
        peptide_df$compressed_end <- sapply(peptide_df$end, compression_map_all$compress)
        
        # Use compressed coordinates for plotting
        peptide_df$plot_start <- peptide_df$compressed_start
        peptide_df$plot_end <- peptide_df$compressed_end
      } else if (nrow(peptide_df) > 0) {
        # Use original coordinates for plotting
        peptide_df$plot_start <- peptide_df$start
        peptide_df$plot_end <- peptide_df$end
      }

      # Step 4: Determining AS overlaps (60-80%)
      incProgress(0.1, detail = 'Collecting AS event regions...')
      # Handle the AS ranges more carefully
      event_ranges <- tryCatch({
        # Check if the columns exist
        range_cols <- intersect(c('AS_range','AS_range1','AS_range2'), names(gene_events))
        cat("Available range columns:", paste(range_cols, collapse=", "), "\n")
        
        if (length(range_cols) == 0) {
          # No range columns found
          cat("No range columns found in gene_events\n")
          GRanges()
        } else {
          # Extract ranges directly without using data.table syntax
          all_ranges <- list()
          
          # Process each event
          for (i in 1:nrow(gene_events)) {
            event_id <- gene_events$eventID[i]
            event_type <- gene_events$AS_type[i]
            
            # Process each range column for this event
            event_range_list <- list()
            
            for (col in range_cols) {
              if (!is.null(gene_events[[col]][i]) && length(gene_events[[col]][[i]]) > 0) {
                # Store the range with its associated event type and column
                current_range <- gene_events[[col]][[i]]
                
                # Add metadata
                mcols(current_range)$event_type <- event_type
                mcols(current_range)$range_type <- col
                mcols(current_range)$event_id <- event_id
                mcols(current_range)$for_visualization <- TRUE
                
                # Add to the list of ranges
                all_ranges <- c(all_ranges, list(current_range))
                
                # Also add to event-specific list for overlap calculation
                event_range_list <- c(event_range_list, list(current_range))
              }
            }
            
            # If we found any ranges for this event, create a combined range for overlap calculation
            if (length(event_range_list) > 0) {
              combined_ranges <- unlist(GRangesList(event_range_list))
              
              # Only proceed if we have valid ranges
              if (length(combined_ranges) > 0) {
                min_start <- min(start(combined_ranges))
                max_end <- max(end(combined_ranges))
                
                # Create a combined range that spans all coordinates for this event
                combined_range <- GRanges(
                  seqnames = seqnames(combined_ranges)[1],
                  ranges = IRanges(start = min_start, end = max_end)
                )
                
                # Add metadata
                mcols(combined_range)$event_type <- event_type
                mcols(combined_range)$event_id <- event_id
                mcols(combined_range)$for_visualization <- FALSE
                mcols(combined_range)$for_overlap <- TRUE
                
                # Add to the list of ranges
                all_ranges <- c(all_ranges, list(combined_range))
              }
            }
          }
          
          # Convert to GRanges if we have ranges
          if (length(all_ranges) > 0) {
            unlist(GRangesList(all_ranges))
          } else {
            GRanges()
          }
        }
      }, error = function(e) {
        cat("Error creating event ranges:", e$message, "\n")
        GRanges()
      })
      cat("Number of event ranges found:", length(event_ranges), "\n")
      
      # Initialize event_overlap column
      if (nrow(peptide_df) > 0) {
      peptide_df$event_overlap <- FALSE
      }
      
      incProgress(0.1, detail = 'Finding AS-affected peptides...')
      if (nrow(peptide_df) > 0 && length(event_ranges) > 0) {
        # Add safer handling of peptide and event ranges
        tryCatch({
          # Make sure we have a valid seqname before creating GRanges
          if (length(seqnames(event_ranges)) > 0) {
            event_seq <- as.character(seqnames(event_ranges)[1])
            
            # Create peptide ranges with proper error handling
            peptide_ranges <- GRanges(
              seqnames = event_seq,
              ranges = IRanges(start = peptide_df$start, end = peptide_df$end)
            )
            
            # Find overlap ranges - these are the combined ranges we created for overlap detection
            overlap_ranges <- tryCatch({
              # Filter to get only the ranges marked for overlap
              if ("for_overlap" %in% names(mcols(event_ranges))) {
                event_ranges[!is.na(mcols(event_ranges)$for_overlap) & mcols(event_ranges)$for_overlap]
              } else {
                # If no metadata, use all ranges
                event_ranges
              }
            }, error = function(e) {
              cat("Error filtering overlap ranges:", e$message, "\n")
              event_ranges
            })
            
            # Find overlaps with the combined ranges
            if (length(overlap_ranges) > 0) {
              hits <- findOverlaps(peptide_ranges, overlap_ranges)
              if (length(hits) > 0) {
                peptide_df$event_overlap[queryHits(hits)] <- TRUE
              }
            }
            
            cat("Number of overlapping peptides found:", sum(peptide_df$event_overlap), "\n")
          }
        }, error = function(e) {
          # Log the error but continue gracefully
          cat("Error finding peptide-AS overlaps:", e$message, "\n")
        })
      }

      # Step 5: Building plot (80-100%)
      incProgress(0.1, detail = 'Preparing visualization data...')
      # gene bounds - handle empty exon case
      gene_start <- 1
      gene_end <- 1000
      
      if (length(exons_by_tx) > 0) {
      sup_ex <- unlist(GRangesList(exons_by_tx))
        if (length(sup_ex) > 0) {
          gene_start <- min(start(sup_ex))
          gene_end <- max(end(sup_ex))
        }
      }
      
      # build data frames
      transcript_df <- data.frame(transcript=transcripts, y=seq_along(transcripts), stringsAsFactors=FALSE)
      
      # Create intron markers for zigzag lines with correct Y positions for each transcript
      if (use_compression_all && !is.null(compression_map_all) && length(exons_by_tx) > 0) {
        intron_markers_all <- data.frame(x = numeric(), y = numeric(), group = character(), transcript = character(), stringsAsFactors = FALSE)
        for(i in seq_along(transcripts)) {
          tx_name <- transcripts[i]
          tx_y_pos <- transcript_df$y[i]  # Y position for this transcript
          if(tx_name %in% names(exons_by_tx)) {
            tx_markers <- create_intron_markers_for_transcripts(list(tx_name = exons_by_tx[[tx_name]]), compression_map_all, y_position = tx_y_pos)
            if(nrow(tx_markers) > 0) {
              intron_markers_all <- rbind(intron_markers_all, tx_markers)
            }
          }
        }
      }
      
      # Build exon data frame safely
      exon_df <- data.frame(
        transcript = character(),
        y = numeric(),
        start = numeric(),
        end = numeric(),
        stringsAsFactors = FALSE
      )
      
      for (i in seq_along(transcripts)) {
        tx <- transcripts[i]
        if (tx %in% names(exons_by_tx) && length(exons_by_tx[[tx]]) > 0) {
          tx_exons <- exons_by_tx[[tx]]
          for (j in seq_along(tx_exons)) {
            exon_df <- rbind(exon_df, data.frame(
              transcript = tx,
              y = i,
              start = start(tx_exons[j]),
              end = end(tx_exons[j]),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
      cat("Number of exons in visualization:", nrow(exon_df), "\n")
      
      # Apply compression to all events exon coordinates
      if (use_compression_all && !is.null(compression_map_all) && nrow(exon_df) > 0) {
        exon_df$compressed_start <- sapply(exon_df$start, compression_map_all$compress)
        exon_df$compressed_end <- sapply(exon_df$end, compression_map_all$compress)
        
        # Use compressed coordinates for plotting
        exon_df$plot_start <- exon_df$compressed_start
        exon_df$plot_end <- exon_df$compressed_end
      } else {
        # Use original coordinates for plotting
        exon_df$plot_start <- exon_df$start
        exon_df$plot_end <- exon_df$end
      }
      
      # Build event region data frame for visualization
      event_df <- data.frame(
        start = numeric(),
        end = numeric(),
        event_type = character(),
        range_type = character(),
        event_id = character(),
        stringsAsFactors = FALSE
      )
      
      if (length(event_ranges) > 0) {
        # Get only the ranges marked for visualization
        vis_ranges <- tryCatch({
          if ("for_visualization" %in% names(mcols(event_ranges))) {
            # Create a safe logical vector with no NAs
            is_vis <- rep(FALSE, length(event_ranges))
            vis_flag <- mcols(event_ranges)$for_visualization
            
            # Mark ranges explicitly for visualization
            for (i in seq_along(vis_flag)) {
              if (!is.na(vis_flag[i]) && vis_flag[i] == TRUE) {
                is_vis[i] <- TRUE
              }
            }
            
            # Return visualization ranges
            event_ranges[is_vis]
          } else {
            # If no metadata, use all ranges
            event_ranges
          }
        }, error = function(e) {
          cat("Error filtering visualization ranges:", e$message, "\n")
          event_ranges
        })
        
        if (length(vis_ranges) > 0) {
          event_df <- data.frame(
            start = start(vis_ranges),
            end = end(vis_ranges),
            event_type = mcols(vis_ranges)$event_type,
            range_type = mcols(vis_ranges)$range_type,
            event_id = mcols(vis_ranges)$event_id,
            stringsAsFactors = FALSE
          )
        }
      }
      cat("Number of event regions in visualization:", nrow(event_df), "\n")
      
      # Apply compression to all events AS event region coordinates
      if (use_compression_all && !is.null(compression_map_all) && nrow(event_df) > 0) {
        event_df$compressed_start <- sapply(event_df$start, compression_map_all$compress)
        event_df$compressed_end <- sapply(event_df$end, compression_map_all$compress)
        
        # Use compressed coordinates for plotting
        event_df$plot_start <- event_df$compressed_start
        event_df$plot_end <- event_df$compressed_end
      } else if (nrow(event_df) > 0) {
        # Use original coordinates for plotting
        event_df$plot_start <- event_df$start
        event_df$plot_end <- event_df$end
      }

      incProgress(0.1, detail = 'Creating plot...')
      
      # Create direct plotly plot instead of using ggplotly
      if (nrow(peptide_df) == 0 || length(transcripts) == 0) {
        # Create empty plot with message
        p <- plotly::plot_ly() %>%
          plotly::add_annotations(
            x = 0.5, y = 0.5,
            text = "No peptide or transcript data available",
            showarrow = FALSE,
            font = list(size = 16)
          ) %>%
          plotly::layout(
            xaxis = list(range = c(0, 1), showticklabels = FALSE),
            yaxis = list(range = c(0, 1), showticklabels = FALSE)
          )
      } else {
        # Create a new plotly object
        p <- plotly::plot_ly()
        
        # Define colors for event types - highly distinct colors
        event_colors <- list(
          "SE" = "#FF0000",  # Bright red
          "RI" = "#00FF00",  # Bright green  
          "MX" = "#0000FF",  # Bright blue
          "A3" = "#235e24",  # Bright yellow
          "A5" = "#FF00FF"   # Bright magenta
        )

        # Define colors for AS status - NEUTRAL COLORS ONLY
        as_status_colors <- list(
          "TRUE" = "#FF0000",   # Red for AS-affected
          "FALSE" = "#0000FF"   # Blue for non-AS-affected
        )
        
        # Define dash types for range types
        dash_types <- list(
          "AS_range" = "solid",
          "AS_range1" = "dash",
          "AS_range2" = "dot"
        )
        
        # Create a list to collect all shapes for the layout
        shapes <- list()
        
        # Add event regions if available
        if (nrow(event_df) > 0) {
          # Initialize tracking for labeled events
          labeled_events <- c()
          
          # Group event regions by event type for legend control
          event_types <- unique(event_df$event_type)
          
          # Add a header trace for the "Event Types" section
          p <- p %>% plotly::add_trace(
            type = "scatter",
            mode = "lines",
            x = c(0), y = c(0),
            line = list(color = "rgba(0,0,0,0)", width = 0),
            showlegend = TRUE,
            name = "<b>Event Types:</b>",
            hoverinfo = "none",
            visible = TRUE
          )
          
          # Process each event type separately
          for (event_type in event_types) {
            # Get events of this type
            type_events <- event_df[event_df$event_type == event_type,]
            color <- if (!is.null(event_colors[[event_type]])) event_colors[[event_type]] else "gray"
            
            # Add a trace for this event type
            p <- p %>% plotly::add_trace(
              type = "scatter",
              mode = "lines",
              x = c(0), y = c(0),
              line = list(color = color, width = 4),
              showlegend = TRUE,
              name = event_type,
              legendgroup = paste0("event_type_", event_type),
              hoverinfo = "none",
              visible = TRUE
            )
            
            # Add rectangular regions for each event of this type
            for (i in 1:nrow(type_events)) {
              range_type <- type_events$range_type[i]
              event_id <- type_events$event_id[i]
              
              # Get the dash type based on range type
              dash <- if (!is.null(dash_types[[range_type]])) dash_types[[range_type]] else "solid"
              
              # Add rectangle as a filled area trace
              p <- p %>% plotly::add_trace(
                type = "scatter",
                mode = "lines",
                x = c(type_events$plot_start[i], type_events$plot_end[i], type_events$plot_end[i], type_events$plot_start[i], type_events$plot_start[i]),
                y = c(0.5, 0.5, length(transcripts) + 0.5, length(transcripts) + 0.5, 0.5),
                fill = "toself",
                fillcolor = paste0(color, "40"), # 40 is hex for 25% opacity
                line = list(color = color, dash = dash, width = 1),
                marker = list(opacity = 0),
                legendgroup = paste0("event_type_", event_type),
                showlegend = FALSE,
                hoverinfo = "text",
                text = clean_hover_text(paste0(
                  "Event ID: ", event_id,
                  "<br>Event Type: ", event_type, 
                  "<br>Range: ", gsub("AS_", "", range_type),
                  "<br>Position: ", type_events$start[i], "-", type_events$end[i]
                ))
              )
              
              # Add event ID label as annotation (one per unique event ID)
              if (input$show_event_labels) {
                # Only add label for the first occurrence of each unique event ID
                if (!event_id %in% labeled_events) {
                  p <- p %>% plotly::add_trace(
                    type = "scatter",
                    x = c(mean(c(type_events$plot_start[i], type_events$plot_end[i]))),
                    y = c(length(transcripts) + 0.8),
                    mode = "text",
                    text = paste0(event_id, " (", event_type, ")"),
                    textfont = list(color = color, size = 10, family = "Arial, sans-serif"),
                    legendgroup = paste0("event_type_", event_type),
                    showlegend = FALSE,
                    hoverinfo = "none"
                  )
                  labeled_events <- c(labeled_events, event_id)
                }
              }
            }
          }
          
          # No need to add invisible hover points anymore since the rectangle traces handle this
        }
        
        # Add transcript lines
        for (i in 1:nrow(transcript_df)) {
          p <- p %>% plotly::add_trace(
            type = "scatter",
            x = c(gene_start, gene_end),
            y = c(transcript_df$y[i], transcript_df$y[i]),
            mode = "lines",
            line = list(color = "black", width = 1),
            showlegend = FALSE,
            hoverinfo = "text",
            text = clean_hover_text(paste0("Transcript: ", transcript_df$transcript[i]))
          )
        }
        
        # Add exons if available
        if (nrow(exon_df) > 0) {
          for (i in 1:nrow(exon_df)) {
            # Add exon as filled trace
            p <- p %>% plotly::add_trace(
              type = "scatter",
              mode = "lines",
              x = c(exon_df$plot_start[i], exon_df$plot_end[i], exon_df$plot_end[i], exon_df$plot_start[i], exon_df$plot_start[i]),
              y = c(exon_df$y[i] - 0.3, exon_df$y[i] - 0.3, exon_df$y[i] + 0.3, exon_df$y[i] + 0.3, exon_df$y[i] - 0.3),
              fill = "toself",
              fillcolor = "grey70",
              line = list(color = "black", width = 1),
              marker = list(opacity = 0),
              showlegend = FALSE,
              hoverinfo = "text",
              text = clean_hover_text(paste0(
                "Exon<br>Transcript: ", exon_df$transcript[i],
                "<br>Position: ", exon_df$start[i], "-", exon_df$end[i]
              ))
            )
          }
        }
        
        # Add peptides if available
        if (nrow(peptide_df) > 0) {
          # Define colors for AS status
          as_status_colors <- list(
            "TRUE" = "#FF0000",   # Red for AS-affected
            "FALSE" = "#0000FF"   # Blue for non-AS-affected
          )
          
          # Create placeholder traces for AS status legend
          p <- p %>% 
            # Add a header trace for the "AS Status" section
            plotly::add_trace(
              type = "scatter",
              mode = "lines",
              x = c(0), y = c(0),
              line = list(color = "rgba(0,0,0,0)", width = 0),
              showlegend = TRUE,
              name = "<b>AS Status:</b>",
              hoverinfo = "none",
              visible = TRUE
            )
          
          # Add AS-affected trace (with red color)
          p <- p %>% plotly::add_trace(
            type = "scatter",
            mode = "markers",
            x = c(0), y = c(0),
            marker = list(color = as_status_colors[["TRUE"]], size = 10),
            showlegend = TRUE,
            name = "AS-affected",
            legendgroup = "as_status_affected",
            hoverinfo = "none",
            visible = TRUE
          )
          
          # Add non-AS-affected trace (with blue color)
          p <- p %>% plotly::add_trace(
            type = "scatter",
            mode = "markers",
            x = c(0), y = c(0),
            marker = list(color = as_status_colors[["FALSE"]], size = 10),
            showlegend = TRUE,
            name = "Not AS-affected",
            legendgroup = "as_status_not_affected",
            hoverinfo = "none",
            visible = TRUE
          )
          
          # Add peptides with AS status colors
          for (i in 1:nrow(peptide_df)) {
            # Get AS status
            is_as <- peptide_df$event_overlap[i]
            color <- if (is_as) as_status_colors[["TRUE"]] else as_status_colors[["FALSE"]]
            legend_group <- if (is_as) "as_status_affected" else "as_status_not_affected"
            
            # Add the peptide as trace instead of shape
            p <- p %>% plotly::add_trace(
              type = "scatter",
              mode = "lines",
              x = c(peptide_df$plot_start[i], peptide_df$plot_end[i], peptide_df$plot_end[i], peptide_df$plot_start[i], peptide_df$plot_start[i]),
              y = c(peptide_df$y[i] - 0.2, peptide_df$y[i] - 0.2, peptide_df$y[i] + 0.2, peptide_df$y[i] + 0.2, peptide_df$y[i] - 0.2),
              fill = "toself",
              fillcolor = color,
              line = list(color = "black", width = 1),
              marker = list(opacity = 0),
              legendgroup = legend_group,
              showlegend = FALSE,
              hoverinfo = "text",
              text = clean_hover_text(paste0(
                "Peptide: ", peptide_df$peptide[i],
                "<br>Position: ", peptide_df$start[i], "-", peptide_df$end[i],
                "<br>AS-affected: ", if (is_as) "Yes" else "No"
              ))
            )
          }
        }
        
        # Add intron markers if compression is enabled - only for displayed transcripts
        if (use_compression_all && !is.null(intron_markers_all) && nrow(intron_markers_all) > 0) {
          # Only show zigzag markers for the transcripts actually being displayed
          displayed_transcripts <- transcript_df$transcript
          filtered_markers_all <- intron_markers_all[intron_markers_all$transcript %in% displayed_transcripts, ]
          
          if (nrow(filtered_markers_all) > 0) {
            # Group intron markers by group and create discrete zigzag segments
            unique_groups <- unique(filtered_markers_all$group)
            for (group in unique_groups) {
              group_data <- filtered_markers_all[filtered_markers_all$group == group, ]
              if (nrow(group_data) > 0) {
                # Create discrete zigzag segment (not continuous line)
                p <- p %>% plotly::add_trace(
                  type = "scatter", mode = "lines",
                  x = group_data$x, y = group_data$y,
                  line = list(color = "grey50", width = 1),
                  connectgaps = FALSE,  # Don't connect to other segments
                  showlegend = FALSE,
                  hoverinfo = "none"
                )
              }
            }
          }
        }
        
        # Get chromosome and strand information (lightning-fast)
        chromosome <- "Unknown"
        if (dir.exists("data/gtf_cache")) {
          gtf_data <- load_gtf_visualization_data(gene_id)
          if (gtf_data$success) {
            chromosome <- gtf_data$chromosome
          } else {
            # Fast GTF cache failed, use fallback
            gene_details <- load_gene_details(gene_id)
            chromosome <- if(gene_details$success) gene_details$chromosome else "Unknown"
          }
        } else {
          # No GTF cache, use original method
          gene_details <- load_gene_details(gene_id)
          chromosome <- if(gene_details$success) gene_details$chromosome else "Unknown"
        }
        
        # Extract strand from gene events (plots are always shown 5' to 3')
        strand_display <- ""
        
        # Set the layout with compression-aware axis
        xaxis_config_all <- if (use_compression_all && !is.null(compression_map_all)) {
          # Use compression: axis shows real coordinates positioned at compressed locations
          axis_breaks <- compression_map_all$coords$compressed_start
          axis_labels <- as.character(compression_map_all$coords$original_start)
          compressed_range <- range(c(compression_map_all$coords$compressed_start, compression_map_all$coords$compressed_end))
          
          list(
            title = paste0("Genomic Position (compressed view) - chromosome ", chromosome, ")"),
            tickvals = axis_breaks,
            ticktext = axis_labels,
            range = compressed_range,
            zeroline = FALSE
          )
        } else {
          # Use original coordinates
          list(
            title = paste0("Genomic Position (chromosome ", chromosome, ")"),
            range = c(gene_start, gene_end),
            zeroline = FALSE
          )
        }
        
        p <- p %>% plotly::layout(
          title = paste0('All-event Comparison - ', gene_id),
          xaxis = xaxis_config_all,
          yaxis = list(
            title = "",
            range = c(0.5, length(transcripts) + 2.5),  # Extended range to accommodate annotations
            tickvals = seq_along(transcripts),
            ticktext = transcripts,
            zeroline = FALSE
          ),
          legend = list(
            orientation = "v",
            x = 1.02,
            y = 0.5,
            title = list(text = "Filter visualization:"),
            itemsizing = "constant",
            tracegroupgap = 10,
            bgcolor = "rgba(255, 255, 255, 0.8)",
            bordercolor = "rgba(0, 0, 0, 0.2)",
            borderwidth = 1
          ),
          margin = list(t = 50, b = 50, l = 100, r = 50)  # Increased margins
        ) %>%
        # Configure legend click behavior to make it more intuitive
        plotly::config(
          displayModeBar = TRUE,
          modeBarButtonsToRemove = c('lasso2d', 'select2d')
        )
      }
      
      # Complete progress and return the plot
      setProgress(1, detail = 'Completed')
      return(p)
    })
  })
  
  # Extract table data for all events peptide comparison
  as_all_events_table_data <- reactive({
    req(input$load_all_events > 0, input$gene, input$protease, input$miscleavage_type)
    
    withProgress(message = 'Preparing table data...', value = 0, {
      # Get the same base data as the plot
      gene_id <- input$gene
      protease <- input$protease
      
      # Get AS events and transcripts
      gene_events <- selected_gene_as_events()
      if (is.null(gene_events) || nrow(gene_events) == 0) {
        return(data.frame(Message = "No AS events found for this gene"))
      }
      
      # Get transcripts involved in AS events
      all_transcripts <- c()
      if (!is.null(gene_events$refTx)) {
        ref_tx <- unlist(gene_events$refTx)
        ref_tx <- ref_tx[!is.na(ref_tx)]
        all_transcripts <- c(all_transcripts, ref_tx)
      }
      if (!is.null(gene_events$asTx)) {
        as_tx <- unlist(gene_events$asTx)
        as_tx <- as_tx[!is.na(as_tx)]
        all_transcripts <- c(all_transcripts, as_tx)
      }
      
      transcripts <- unique(all_transcripts)
      if (length(transcripts) == 0) {
        return(data.frame(Message = "No transcripts found for AS events"))
      }
      
      incProgress(0.3, detail = 'Getting peptides...')
      
      # Create vis_data structure
      vis_data <- list(
        genes = processed_data()$genes,
        gene_symbols = processed_data()$gene_symbols,
        gene_lookup = processed_data()$gene_lookup,
        proteases = processed_data()$proteases,
        original_peptides = selected_gene_peptides()
      )
      
      # Get peptides for all transcripts
      peptide_list <- lapply(seq_along(transcripts), function(i) {
        tx <- transcripts[i]
        gr <- get_transcript_peptides_for_comparison(tx, vis_data, protease)
        if (is.null(gr)) return(NULL)
        data.frame(
          transcript = tx, 
          start = start(gr), 
          end = end(gr), 
          peptide = gr$peptide, 
          stringsAsFactors = FALSE
        )
      })
      
      # Filter out NULL results and combine
      peptide_list <- peptide_list[!sapply(peptide_list, is.null)]
      if (length(peptide_list) == 0) {
        return(data.frame(Message = "No peptides found for transcripts"))
      }
      
      peptide_df <- data.table::rbindlist(peptide_list)
      peptide_df$event_overlap <- FALSE
      
      incProgress(0.4, detail = 'Determining AS overlaps...')
      
      # Get AS event ranges for overlap detection
      event_ranges <- list()
      range_cols <- intersect(c('AS_range','AS_range1','AS_range2'), names(gene_events))
      
      for (i in 1:nrow(gene_events)) {
        event_id <- gene_events$eventID[i]
        event_type <- gene_events$AS_type[i]
        
        for (col in range_cols) {
          if (!is.null(gene_events[[col]][i]) && length(gene_events[[col]][[i]]) > 0) {
            current_range <- gene_events[[col]][[i]]
            mcols(current_range)$event_type <- event_type
            mcols(current_range)$event_id <- event_id
            event_ranges <- c(event_ranges, list(current_range))
          }
        }
      }
      
      # Find overlaps if we have both peptides and event ranges
      if (nrow(peptide_df) > 0 && length(event_ranges) > 0) {
        event_ranges_gr <- unlist(GRangesList(event_ranges))
        
        # Add chromosome information to peptide_df
        chromosome <- as.character(seqnames(event_ranges_gr)[1])
        peptide_df$chromosome <- chromosome
        
        # Create peptide ranges
        peptide_ranges <- GRanges(
          seqnames = seqnames(event_ranges_gr)[1],  # Use same chromosome
          ranges = IRanges(start = peptide_df$start, end = peptide_df$end)
        )
        
        # Find overlaps
        hits <- findOverlaps(peptide_ranges, event_ranges_gr)
        if (length(hits) > 0) {
          peptide_df$event_overlap[queryHits(hits)] <- TRUE
          
          # Add event information to overlapping peptides
          peptide_df$event_id <- ""
          peptide_df$event_type <- ""
          
          for (hit_idx in seq_along(hits)) {
            peptide_idx <- queryHits(hits)[hit_idx]
            event_idx <- subjectHits(hits)[hit_idx]
            
            event_id <- mcols(event_ranges_gr)$event_id[event_idx]
            event_type <- mcols(event_ranges_gr)$event_type[event_idx]
            
            if (peptide_df$event_id[peptide_idx] == "") {
              peptide_df$event_id[peptide_idx] <- event_id
              peptide_df$event_type[peptide_idx] <- event_type
            } else {
              # Append additional events
              if (!grepl(event_id, peptide_df$event_id[peptide_idx])) {
                peptide_df$event_id[peptide_idx] <- paste(peptide_df$event_id[peptide_idx], event_id, sep = ", ")
                peptide_df$event_type[peptide_idx] <- paste(peptide_df$event_type[peptide_idx], event_type, sep = ", ")
              }
            }
          }
        } else {
          peptide_df$event_id <- ""
          peptide_df$event_type <- ""
        }
      } else {
        # No event ranges available - try to get chromosome from gene_events directly
        chromosome <- NULL
        
        # Try to extract chromosome from any available AS range in gene_events
        range_cols <- intersect(c('AS_range','AS_range1','AS_range2'), names(gene_events))
        for (col in range_cols) {
          for (i in 1:nrow(gene_events)) {
            if (!is.null(gene_events[[col]][i]) && length(gene_events[[col]][[i]]) > 0) {
              chromosome <- as.character(seqnames(gene_events[[col]][[i]])[1])
              break
            }
          }
          if (!is.null(chromosome)) break
        }
        
        if (!is.null(chromosome)) {
          peptide_df$chromosome <- chromosome
        } else {
          # If we still can't determine chromosome, return an error
          return(data.frame(Message = "Cannot determine chromosome information for peptides"))
        }
        
        peptide_df$event_id <- ""
        peptide_df$event_type <- ""
      }
      
      # Ensure chromosome column exists even if no overlaps were processed
      if (!"chromosome" %in% names(peptide_df)) {
        # Try one more time to get chromosome from gene_events
        chromosome <- NULL
        range_cols <- intersect(c('AS_range','AS_range1','AS_range2'), names(gene_events))
        for (col in range_cols) {
          for (i in 1:nrow(gene_events)) {
            if (!is.null(gene_events[[col]][i]) && length(gene_events[[col]][[i]]) > 0) {
              chromosome <- as.character(seqnames(gene_events[[col]][[i]])[1])
              break
            }
          }
          if (!is.null(chromosome)) break
        }
        
        if (!is.null(chromosome)) {
          peptide_df$chromosome <- chromosome
        } else {
          return(data.frame(Message = "Cannot determine chromosome information for peptides"))
        }
      }
      
      incProgress(0.2, detail = 'Formatting table...')
      
      # Group peptides like existing peptide comparison table
      display_df <- peptide_df %>%
        group_by(peptide) %>%
        summarise(
          all_positions = paste(unique(paste0(first(chromosome), ":", start, "-", end)), collapse = ", "),
          is_junction_spanning = n() > 1,
          transcripts = paste(unique(transcript), collapse = ", "),
          event_ids = paste(unique(event_id[event_id != ""]), collapse = ", "),
          event_types = paste(unique(event_type[event_type != ""]), collapse = ", "),
          as_affected = any(event_overlap),
          .groups = "drop"
        )
      
      # Clean up empty event data
      display_df$event_ids[display_df$event_ids == ""] <- "-"
      display_df$event_types[display_df$event_types == ""] <- "-"
      
      incProgress(0.1, detail = 'Complete!')
      
      return(display_df)
    })
  })
  
  # Render the all-events peptide comparison plot
  output$as_all_peptide_comparison <- renderPlotly({
    req(as_all_peptide_comparison_data())
    as_all_peptide_comparison_data()
  })
  
  # Render the all-events peptide table
  output$as_all_events_table <- DT::renderDataTable({
    data <- as_all_events_table_data()
    
    if ("Message" %in% names(data)) {
      return(data)
    }
    
    # Format column names for display
    colnames(data) <- c("Peptide", "All Positions", "Junction-spanning", 
                       "Transcript(s)", "Event ID(s)", "Event Type(s)", "AS-affected")
    
    DT::datatable(
      data,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel'),
        order = list(list(6, 'desc'), list(2, 'desc'))  # Sort by AS-affected, then junction-spanning
      ),
      rownames = FALSE
    ) %>%
    # Color coding for junction-spanning peptides
    DT::formatStyle(
      'Junction-spanning',
      backgroundColor = styleEqual(c(TRUE, FALSE), c("#FFE6E6", "#E6F3FF"))
    ) %>%
    # Color coding for AS-affected peptides
    DT::formatStyle(
      'AS-affected', 
      backgroundColor = styleEqual(c(TRUE, FALSE), c("#FFEBE6", "#F0F8F0"))
    ) %>%
    # Color coding for event types
    DT::formatStyle(
      'Event Type(s)',
      backgroundColor = styleEqual(
        c("SE", "RI", "MX", "A3", "A5", "-"),
        c("#FFE6E6", "#E6FFE6", "#E6E6FF", "#FFFFE6", "#FFE6FF", "#F8F9FA")
      )
    )
  })
  
  # Download handler for all events table
  output$download_all_events_table <- downloadHandler(
    filename = function() {
      gene_symbol <- processed_data()$gene_lookup[input$gene] %||% input$gene
      paste0("all_events_peptides_", gene_symbol, "_", input$protease, "_", 
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      data <- as_all_events_table_data()
      if (!"Message" %in% names(data)) {
        # Add metadata columns
        data$gene_id <- input$gene
        data$gene_symbol <- processed_data()$gene_lookup[input$gene] %||% input$gene
        data$enzyme <- input$protease
        data$miscleavage_type <- input$miscleavage_type
        data$analysis_timestamp <- Sys.time()
        
        write.csv(data, file, row.names = FALSE)
      }
    }
  )
  
  #===============================================================================
  # ISOFORM ANALYSIS TAB
  #===============================================================================
  
  # Update isoform selector when gene OR miscleavage type changes
  observeEvent(list(input$gene, input$miscleavage_type), {
    if (!is.null(input$gene) && !is.null(input$miscleavage_type)) {
      # Get transcripts for the selected gene (this will use the gene-by-gene loading with correct miscleavage)
      transcripts <- isolate(selected_gene_transcripts())
      cat("Debug: Found transcripts for gene", input$gene, "with miscleavage", input$miscleavage_type, ":", length(transcripts), "\n")
      cat("Debug: Transcript IDs:", paste(transcripts, collapse=", "), "\n")
      
      if (!is.null(transcripts) && length(transcripts) > 0) {
        # Create nice labels for the dropdown
        transcript_choices <- setNames(transcripts, transcripts)
        updateSelectInput(session, "highlight_isoform", 
                         choices = transcript_choices,
                         selected = transcripts[1])
      } else {
        # Clear choices if no transcripts found
        updateSelectInput(session, "highlight_isoform", 
                         choices = character(0),
                         selected = NULL)
      }
    }
  })
  
  # Reactive data for all isoforms analysis (like all events peptide comparison)
  all_isoforms_data <- reactive({
    req(input$load_all_isoforms > 0, input$gene, input$protease, input$miscleavage_type)
    
    # The gene data should already be loaded by the observeEvent - don't call load_and_cache_gene_data here!
    req(gene_data()) # Now, we can safely require the data
    
    withProgress(message = 'Loading all isoforms analysis...', value = 0, {
      gene_id <- input$gene
      protease <- input$protease
      miscleavage_type <- input$miscleavage_type
      
      # Step 1: Get all transcripts for this gene (isolated to prevent reactive loops)
      incProgress(0.2, detail = 'Getting all gene transcripts...')
      all_transcripts <- isolate(selected_gene_transcripts())
      
      if (is.null(all_transcripts) || length(all_transcripts) == 0) {
        showNotification(paste("No transcripts found for gene", gene_id, "with miscleavage type", miscleavage_type), type = "warning")
        return(NULL)
      }
      
      cat("All transcripts for isoform analysis:", paste(all_transcripts, collapse=", "), "\n")
      cat("Using miscleavage type:", miscleavage_type, "\n")
      
      # Step 2: Get peptides for each transcript (using gene-by-gene loaded data)
      incProgress(0.3, detail = 'Loading peptides for all transcripts...')
      # Use the gene-by-gene loaded peptides which already have the correct miscleavage type (isolated to prevent reactive loops)
      gene_peptides <- isolate(selected_gene_peptides())
      if (is.null(gene_peptides) || nrow(gene_peptides) == 0) {
        showNotification(paste("No peptides found for gene", gene_id, "with miscleavage type", miscleavage_type), type = "warning")
        return(NULL)
      }
      
      # Create standardized vis_data structure using core data module  
      vis_data <- core_data_module$create_vis_data_structure()
      
      # Get peptides for all transcripts
      all_peptides_list <- list()
      for (i in seq_along(all_transcripts)) {
        tx <- all_transcripts[i]
        tx_peptides <- get_transcript_peptides_for_comparison(tx, vis_data, protease)
        if (!is.null(tx_peptides) && length(tx_peptides) > 0) {
          all_peptides_list[[tx]] <- data.frame(
            transcript = tx,
            y_position = i,  # Assign y position here
            start = start(tx_peptides),
            end = end(tx_peptides),
            peptide = tx_peptides$peptide,
            stringsAsFactors = FALSE
          )
        }
      }
      
      if (length(all_peptides_list) == 0) {
        return(NULL)
      }
      
      # Step 3: Combine all peptides
      incProgress(0.2, detail = 'Combining peptide data...')
      all_peptides_df <- data.table::rbindlist(all_peptides_list)
      
      # Calculate gene boundaries
      gene_start <- min(all_peptides_df$start) - 1000
      gene_end <- max(all_peptides_df$end) + 1000
      
      # Create transcript position mapping
      transcript_df <- data.frame(
        transcript = all_transcripts,
        y_position = seq_along(all_transcripts),
        stringsAsFactors = FALSE
      )
      
          # Add hover text for all peptides
    all_peptides_df$hover_text <- paste0(
      "Peptide: ", all_peptides_df$peptide,
      "<br>Position: ", all_peptides_df$start, "-", all_peptides_df$end,
      "<br>Transcript: ", all_peptides_df$transcript,
      "<br>Miscleavage: ", miscleavage_type
    )
      
      incProgress(0.2, detail = 'Completed')
      
      return(list(
        all_peptides = all_peptides_df,
        transcript_df = transcript_df,
        gene_start = gene_start,
        gene_end = gene_end,
        all_transcripts = all_transcripts
      ))
    })
  })
  
  # Reactive data for highlighting specific isoform with chromosome info and merged peptides
  highlighted_isoform_data <- reactive({
    req(all_isoforms_data(), input$highlight_isoform, input$miscleavage_type)
    
    base_data <- all_isoforms_data()
    highlight_isoform <- input$highlight_isoform
    miscleavage_type <- input$miscleavage_type
    
    if (is.null(base_data)) return(NULL)
    
    all_peptides_df <- base_data$all_peptides
    all_transcripts <- base_data$all_transcripts
    
    # Calculate specificity for highlighted isoform peptides
    highlight_peptides <- all_peptides_df[all_peptides_df$transcript == highlight_isoform, ]
    
    if (nrow(highlight_peptides) == 0) {
      return(base_data)
    }
    
    # For each peptide in highlighted isoform, calculate specificity
    highlight_peptides$other_transcript_count <- 0
    highlight_peptides$other_transcripts <- ""
    highlight_peptides$specificity_category <- ""
    highlight_peptides$chromosome <- ""
    highlight_peptides$genomic_positions <- ""
    
    # Get chromosome info from genomic ranges using the RDS data structure
    gene_id <- input$gene
    protease <- input$protease
    
    tryCatch({
      # Get gene data for the highlighted isoform
      gene_data <- selected_gene_peptides()  # This contains the current gene data
      if (!is.null(gene_data) && highlight_isoform %in% gene_data$txID) {
        tx_row <- which(gene_data$txID == highlight_isoform)[1]
        
        # Get the mapped ranges column for this protease
        mapped_ranges_col <- paste0(protease, "Peps_mapped_ranges")
        
        if (mapped_ranges_col %in% names(gene_data) && !is.null(gene_data[[mapped_ranges_col]][[tx_row]])) {
          genomic_ranges <- gene_data[[mapped_ranges_col]][[tx_row]]
          
          if (length(genomic_ranges) > 0) {
            # Extract chromosome using seqnames
            chromosomes <- unique(as.character(seqnames(genomic_ranges)))
            if (length(chromosomes) > 0) {
              chromosome_info <- chromosomes[1]
              
              # Map peptides to their genomic locations
              for (i in 1:nrow(highlight_peptides)) {
                peptide_seq <- highlight_peptides$peptide[i]
                
                # Find matching ranges for this peptide
                matching_ranges_idx <- which(genomic_ranges$peptide == peptide_seq)
                
                if (length(matching_ranges_idx) > 0) {
                  # Get all genomic positions for this peptide
                  peptide_ranges <- genomic_ranges[matching_ranges_idx]
                  genomic_positions <- paste(
                    paste0(start(peptide_ranges), "-", end(peptide_ranges)),
                    collapse = "; "
                  )
                  
                  highlight_peptides$chromosome[i] <- chromosome_info
                  highlight_peptides$genomic_positions[i] <- genomic_positions
                }
              }
            }
          }
        }
      }
    }, error = function(e) {
      cat("Warning: Could not extract chromosome info:", e$message, "\n")
    })
    
    for (i in 1:nrow(highlight_peptides)) {
      peptide_seq <- highlight_peptides$peptide[i]
      # Count transcripts (excluding highlight) that have this peptide
      other_tx_with_peptide <- unique(all_peptides_df$transcript[
        all_peptides_df$peptide == peptide_seq & all_peptides_df$transcript != highlight_isoform
      ])
      
      highlight_peptides$other_transcript_count[i] <- length(other_tx_with_peptide)
      highlight_peptides$other_transcripts[i] <- paste(other_tx_with_peptide, collapse = ", ")
      
      # Classify specificity into three categories: Unique, Shared, Universal
      total_isoforms <- length(all_transcripts)
      other_count <- length(other_tx_with_peptide)
      
      if (other_count == 0) {
        # Peptide found only in this isoform
        highlight_peptides$specificity_category[i] <- "Unique"
      } else if (other_count == (total_isoforms - 1)) {
        # Peptide found in ALL isoforms
        highlight_peptides$specificity_category[i] <- "Universal"
      } else {
        # Peptide found in some but not all isoforms
        highlight_peptides$specificity_category[i] <- "Shared"
      }
    }
    
    # Merge duplicate peptides (same peptide sequence appears multiple times due to multiple genomic locations)
    merged_peptides <- highlight_peptides %>%
      dplyr::group_by(peptide, specificity_category, other_transcript_count, other_transcripts, chromosome) %>%
      dplyr::summarise(
        genomic_positions = paste(unique(genomic_positions[genomic_positions != ""]), collapse = "; "),
        start_pos = min(start, na.rm = TRUE),
        end_pos = max(end, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Create final merged table with proper column names
    highlight_peptides_merged <- data.frame(
      peptide = merged_peptides$peptide,
      chromosome = merged_peptides$chromosome,
      genomic_positions = ifelse(merged_peptides$genomic_positions == "", "N/A", merged_peptides$genomic_positions),
      start = merged_peptides$start_pos,
      end = merged_peptides$end_pos,
      specificity_category = merged_peptides$specificity_category,
      other_transcript_count = merged_peptides$other_transcript_count,
      other_transcripts = merged_peptides$other_transcripts,
      stringsAsFactors = FALSE
    )
    
    # Add specificity info to all peptides for visualization
    all_peptides_df$is_highlighted <- all_peptides_df$transcript == highlight_isoform
    all_peptides_df$specificity_category <- ""
    all_peptides_df$specificity_category[all_peptides_df$is_highlighted] <- 
      highlight_peptides$specificity_category[match(all_peptides_df$peptide[all_peptides_df$is_highlighted], highlight_peptides$peptide)]
    
    # Update hover text
    all_peptides_df$hover_text <- paste0(
      "Peptide: ", all_peptides_df$peptide,
      "<br>Position: ", all_peptides_df$start, "-", all_peptides_df$end,
      "<br>Transcript: ", all_peptides_df$transcript,
      "<br>Miscleavage: ", miscleavage_type,
      ifelse(all_peptides_df$is_highlighted, 
             paste0("<br>Specificity: ", all_peptides_df$specificity_category), 
             "")
    )
    
    # Create summary statistics for the three categories using merged data
    peptide_categories_summary <- as.data.frame(table(highlight_peptides_merged$specificity_category))
    names(peptide_categories_summary) <- c("Category", "Count")
    peptide_categories_summary$Percentage <- round(peptide_categories_summary$Count / nrow(highlight_peptides_merged) * 100, 1)
    
    # Create isoform coverage analysis
    isoform_coverage <- data.frame(
      Isoform = all_transcripts,
      Total_Peptides = sapply(all_transcripts, function(tx) {
        sum(all_peptides_df$transcript == tx)
      }),
      Unique_Peptides = sapply(all_transcripts, function(tx) {
        tx_peptides <- all_peptides_df[all_peptides_df$transcript == tx, "peptide"]
        sum(sapply(tx_peptides, function(pep) {
          sum(all_peptides_df$peptide == pep) == 1
        }))
      }),
      Shared_Peptides = sapply(all_transcripts, function(tx) {
        tx_peptides <- all_peptides_df[all_peptides_df$transcript == tx, "peptide"]
        sum(sapply(tx_peptides, function(pep) {
          count <- sum(all_peptides_df$peptide == pep)
          count > 1 && count < length(all_transcripts)
        }))
      }),
      Universal_Peptides = sapply(all_transcripts, function(tx) {
        tx_peptides <- all_peptides_df[all_peptides_df$transcript == tx, "peptide"]
        sum(sapply(tx_peptides, function(pep) {
          sum(all_peptides_df$peptide == pep) == length(all_transcripts)
        }))
      }),
      stringsAsFactors = FALSE
    )
    
    return(list(
      all_peptides = all_peptides_df,
      highlight_peptides = highlight_peptides_merged,  # Use merged data
      transcript_df = base_data$transcript_df,
      gene_start = base_data$gene_start,
      gene_end = base_data$gene_end,
      highlight_isoform = highlight_isoform,
      peptide_categories_summary = peptide_categories_summary,
      isoform_coverage = isoform_coverage
    ))
  })
  
  # Render all isoforms plot with GTF-based exon/CDS boundaries
  output$all_isoforms_plot <- renderPlotly({
    data <- highlighted_isoform_data()
    
    if (is.null(data)) {
      p <- plotly::plot_ly() %>%
        plotly::add_annotations(
          x = 0.5, y = 0.5,
          text = "No peptide data available. Please ensure gene data is loaded.",
          showarrow = FALSE,
          font = list(size = 16)
        ) %>%
        plotly::layout(
          xaxis = list(range = c(0, 1), showticklabels = FALSE),
          yaxis = list(range = c(0, 1), showticklabels = FALSE)
        )
      return(p)
    }
    
    # Extract data
    all_peptides <- data$all_peptides
    transcript_df <- data$transcript_df
    gene_start <- data$gene_start
    gene_end <- data$gene_end
    highlight_isoform <- data$highlight_isoform
    all_transcripts <- transcript_df$transcript
    
    # Load GTF data for exon and CDS boundaries (lightning-fast cached GTF)
    exons_result <- NULL
    if (dir.exists("data/gtf_cache")) {
      gtf_data <- load_gtf_visualization_data(input$gene)
      if (gtf_data$success) {
        exons_result <- list(
          success = TRUE,
          exons = gtf_data$exons_by_transcript,
          cds = gtf_data$cds_by_transcript
        )
      } else {
        # Fast GTF cache failed, use fallback
        gene_details <- load_gene_details(input$gene)
        exons_result <- load_transcript_exons(gene_details, all_transcripts)
      }
    } else {
      # No GTF cache, use original method
      gene_details <- load_gene_details(input$gene)
      exons_result <- load_transcript_exons(gene_details, all_transcripts)
    }
    
    # Add compression support - respect user preference
    use_compression_isoforms <- !is.null(input$intron_scale) && input$intron_scale == "compressed"
    compression_map_isoforms <- NULL
    
    # Create compression map if user selected compressed mode and exon data exists
    if (use_compression_isoforms && exons_result$success && length(exons_result$exons) > 0) {
      compression_map_isoforms <- create_compression_map(exons_result$exons)
    }
    
    # Create plotly object
    p <- plotly::plot_ly()
    
    # Define colors for specificity categories using rgba to prevent hex codes in hover
    specificity_colors <- list(
      "Unique" = "rgba(255, 0, 0, 0.8)",      # Red - peptides unique to this isoform
      "Shared" = "rgba(243, 156, 18, 0.8)",   # Orange - peptides shared with some isoforms
      "Universal" = "rgba(46, 204, 113, 0.8)" # Green - peptides found in all isoforms
    )
    
    # Add transcript lines (ALL transcripts)
    for (i in 1:nrow(transcript_df)) {
      line_color <- if (transcript_df$transcript[i] == highlight_isoform) "#000000" else "#CCCCCC"
      line_width <- if (transcript_df$transcript[i] == highlight_isoform) 3 else 1
      
      # Apply compression to transcript line coordinates
      plot_gene_start <- gene_start
      plot_gene_end <- gene_end
      if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
        plot_gene_start <- compression_map_isoforms$compress(gene_start)
        plot_gene_end <- compression_map_isoforms$compress(gene_end)
      }
      
      p <- p %>% plotly::add_trace(
        type = "scatter",
        x = c(plot_gene_start, plot_gene_end),
        y = c(transcript_df$y_position[i], transcript_df$y_position[i]),
        mode = "lines",
        line = list(color = line_color, width = line_width),
        showlegend = FALSE,
        hoverinfo = "text",
        text = clean_hover_text(paste0("Transcript: ", transcript_df$transcript[i]))
      )
    }
    
    # Add exon and CDS boundaries if GTF data is available
    if (exons_result$success && length(exons_result$exons) > 0) {
      exons_by_transcript <- exons_result$exons
      cds_by_transcript <- exons_result$cds
      
      # Add exon blocks for each transcript
      for (i in 1:nrow(transcript_df)) {
        tx <- transcript_df$transcript[i]
        y_pos <- transcript_df$y_position[i]
        
                 # Add exons if available (make them larger to contain peptides)
         if (!is.null(exons_by_transcript[[tx]]) && length(exons_by_transcript[[tx]]) > 0) {
           tx_exons <- exons_by_transcript[[tx]]
           for (j in seq_along(tx_exons)) {
             # Apply compression to exon coordinates
             exon_start <- start(tx_exons[j])
             exon_end <- end(tx_exons[j])
             if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
               exon_start <- compression_map_isoforms$compress(exon_start)
               exon_end <- compression_map_isoforms$compress(exon_end)
             }
             
             p <- p %>% plotly::add_trace(
               type = "scatter", mode = "lines",
               x = c(exon_start, exon_end, exon_end, exon_start, exon_start),
               y = c(y_pos - 0.35, y_pos - 0.35, y_pos + 0.35, y_pos + 0.35, y_pos - 0.35),
               fill = "toself",
               fillcolor = "rgba(211, 211, 211, 0.3)",  # Transparent grey for exons
               line = list(color = "rgba(128, 128, 128, 0.5)", width = 1),  # Transparent grey border
              legendgroup = "structure",
              showlegend = FALSE,
              hovertemplate = paste0("Exon ", j, " (", start(tx_exons[j]), "-", end(tx_exons[j]), ")<br>Transcript: ", tx, "<extra></extra>"),
              hoverinfo = "none"
            )
          }
        }
        
                 # Add CDS overlay if available (slightly smaller than exons)
         if (!is.null(cds_by_transcript[[tx]]) && length(cds_by_transcript[[tx]]) > 0) {
           tx_cds <- cds_by_transcript[[tx]]
           for (j in seq_along(tx_cds)) {
             # Apply compression to CDS coordinates
             cds_start <- start(tx_cds[j])
             cds_end <- end(tx_cds[j])
             if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
               cds_start <- compression_map_isoforms$compress(cds_start)
               cds_end <- compression_map_isoforms$compress(cds_end)
             }
             
             p <- p %>% plotly::add_trace(
               type = "scatter", mode = "lines",
               x = c(cds_start, cds_end, cds_end, cds_start, cds_start),
               y = c(y_pos - 0.25, y_pos - 0.25, y_pos + 0.25, y_pos + 0.25, y_pos - 0.25),
               fill = "toself",
               fillcolor = "#F1C40F",  # Gold/yellow for CDS (consistent with visualization tab)
               line = list(color = "#DAA520", width = 1),  # Goldenrod border
              legendgroup = "structure",
              showlegend = FALSE,
              hovertemplate = paste0("CDS ", j, " (", start(tx_cds[j]), "-", end(tx_cds[j]), ")<br>Transcript: ", tx, "<extra></extra>"),
              hoverinfo = "none"
            )
          }
        }
      }
    }
    
    # Add legend headers
    p <- p %>% plotly::add_trace(
      type = "scatter", mode = "lines",
      x = c(0), y = c(0),
      line = list(color = "rgba(0,0,0,0)", width = 0),
      showlegend = TRUE,
      name = paste0("<b>Highlighted: ", highlight_isoform, "</b>"),
      hoverinfo = "none",
      visible = TRUE
    )
    
    # Add legend entries for each category with counts and percentages
    peptide_summary <- data$peptide_categories_summary
    for (category in names(specificity_colors)) {
      # Get count and percentage for this category
      category_row <- peptide_summary[peptide_summary$Category == category, ]
      if (nrow(category_row) > 0) {
        count <- category_row$Count
        percentage <- category_row$Percentage
        legend_name <- paste0(category, ": ", count, " (", percentage, "%)")
      } else {
        # Handle case where category has no peptides
        legend_name <- paste0(category, ": 0 (0%)")
      }
      
      p <- p %>% plotly::add_trace(
        type = "scatter", mode = "markers",
        x = c(0), y = c(0),
        marker = list(color = specificity_colors[[category]], size = 10),
        showlegend = TRUE,
        name = legend_name,
        legendgroup = paste0("category_", category),
        hoverinfo = "none",
        visible = TRUE
      )
    }
    
         # Add legend entries for gene structure
     p <- p %>% plotly::add_trace(
       type = "scatter", mode = "markers",
       x = c(0), y = c(0),
       marker = list(color = "rgba(211, 211, 211, 0.3)", size = 10, symbol = "square"),
       showlegend = TRUE,
       name = "Exons",
       legendgroup = "structure",
       hoverinfo = "none",
       visible = TRUE
     ) %>% plotly::add_trace(
       type = "scatter", mode = "markers",
       x = c(0), y = c(0),
       marker = list(color = "#F1C40F", size = 10, symbol = "square"),
       showlegend = TRUE,
       name = "CDS",
       legendgroup = "structure",
       hoverinfo = "none",
       visible = TRUE
     )
    
    # Add peptides for ALL transcripts
    if (nrow(all_peptides) > 0) {
      for (i in 1:nrow(all_peptides)) {
        # Determine color based on highlighting
        if (all_peptides$is_highlighted[i] && all_peptides$specificity_category[i] != "") {
          # Highlighted isoform peptides with specificity colors
          color <- specificity_colors[[all_peptides$specificity_category[i]]]
          legend_group <- paste0("category_", all_peptides$specificity_category[i])
        } else {
          # Non-highlighted peptides in light gray
          color <- "#DDDDDD"
          legend_group <- "other"
        }
        
        # Add peptide as filled rectangle (smaller to fit inside exons)
        # Apply compression to peptide coordinates
        peptide_start <- all_peptides$start[i]
        peptide_end <- all_peptides$end[i]
        if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
          peptide_start <- compression_map_isoforms$compress(peptide_start)
          peptide_end <- compression_map_isoforms$compress(peptide_end)
        }
        
        p <- p %>% plotly::add_trace(
          type = "scatter", mode = "lines",
          x = c(peptide_start, peptide_end, peptide_end, peptide_start, peptide_start),
          y = c(all_peptides$y_position[i] - 0.15, all_peptides$y_position[i] - 0.15, all_peptides$y_position[i] + 0.15, all_peptides$y_position[i] + 0.15, all_peptides$y_position[i] - 0.15),
          fill = "toself",
          fillcolor = color,
          line = list(color = "black", width = 0.5),
          marker = list(opacity = 0),
          legendgroup = legend_group,
          showlegend = FALSE,
          hovertemplate = paste0(clean_hover_text(all_peptides$hover_text[i]), "<extra></extra>"),
          hoverinfo = "none"
        )
      }
    }
    
    
    # Get chromosome and strand information from the base data
    chromosome <- if(exists("base_data") && !is.null(base_data$chromosome)) {
      base_data$chromosome
    } else {
      # Try to get from fast GTF cache
      if (dir.exists("data/gtf_cache")) {
        gtf_data <- load_gtf_visualization_data(input$gene)
        if (gtf_data$success) {
          gtf_data$chromosome
        } else {
          # Fast GTF cache failed, use fallback
          gene_details_local <- load_gene_details(input$gene)
          if(gene_details_local$success) gene_details_local$chromosome else "Unknown"
        }
      } else {
        # No GTF cache, use original method
        gene_details_local <- load_gene_details(input$gene)
        if(gene_details_local$success) gene_details_local$chromosome else "Unknown"
      }
    }
    
    # Extract strand from peptide data if available (lightning-fast)
    strand_display <- ""  # Default to standard display
    if(nrow(all_peptides) > 0) {
      # Get first transcript's peptides to check strand (assuming all same gene)
      gene_id <- input$gene
      
      # Use fast GTF cache if available
      if (dir.exists("data/gtf_cache")) {
        gtf_data <- load_gtf_visualization_data(gene_id)
        if (gtf_data$success) {
          strand_display <- ""  # Already computed in GTF data
        } else {
          # Fast GTF cache failed, use fallback
          gene_details_local <- load_gene_details(gene_id)
          if(gene_details_local$success) {
            first_tx <- all_transcripts[1]
            exons_result <- load_transcript_exons(gene_details_local, c(first_tx))
            if(exons_result$success && length(exons_result$exons) > 0) {
              first_exons <- exons_result$exons[[1]]
              if(length(first_exons) > 0) {
                strand_char <- as.character(strand(first_exons)[1])
                strand_display <- ""
              }
            }
          }
        }
      } else {
        # No GTF cache, use original method
        gene_details_local <- load_gene_details(gene_id)
        if(gene_details_local$success) {
          first_tx <- all_transcripts[1]
          exons_result <- load_transcript_exons(gene_details_local, c(first_tx))
          if(exons_result$success && length(exons_result$exons) > 0) {
            first_exons <- exons_result$exons[[1]]
            if(length(first_exons) > 0) {
              strand_char <- as.character(strand(first_exons)[1])
              strand_display <- ""
            }
          }
        }
      }
    }
    
    # Add triangle marker for selected isoform-specific peptide
    selected_isoform_pep <- selected_isoform_peptide()
    if (!is.null(selected_isoform_pep) && nrow(all_peptides) > 0) {
      # Get the highlighted isoform name
      highlight_isoform <- data$highlight_isoform
      
      # Find the y_position of the highlighted isoform in transcript_df
      highlight_y_position <- transcript_df$y_position[transcript_df$transcript == highlight_isoform]
      
      if (length(highlight_y_position) > 0) {
        # Find the selected peptide on the highlighted isoform
        matching_peptide <- all_peptides[all_peptides$peptide == selected_isoform_pep$peptide & 
                                        all_peptides$transcript == highlight_isoform, ]
        
        if (nrow(matching_peptide) > 0) {
          # Use the peptide coordinates but force y_position BELOW the highlighted isoform
          # Apply compression to triangle marker coordinates
          triangle_x <- (matching_peptide$start[1] + matching_peptide$end[1]) / 2
          if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
            triangle_x <- compression_map_isoforms$compress(triangle_x)
          }
          triangle_y <- highlight_y_position[1] - 0.35
          
          # Add triangle marker
          p <- p %>% plotly::add_trace(
            type = "scatter", mode = "markers",
            x = triangle_x, y = triangle_y,
            marker = list(symbol = "triangle-up", size = 15, color = "#FF0000", 
                         line = list(color = "black", width = 1)),
            showlegend = TRUE, name = "Selected Peptide",
            hovertemplate = paste0("Selected: ", selected_isoform_pep$peptide, "<br>On Isoform: ", highlight_isoform, "<extra></extra>"),
            hoverinfo = "none"
          )
        }
      }
    }
    
    # Set layout
    miscleavage_label <- switch(input$miscleavage_type,
      "no_miss_cleavage" = "No Miscleavage",
      "upto_two_misscleavage" = "Up to 2 Miscleavages"
    )
    # Configure X-axis based on compression mode
    xaxis_config <- if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
      # Compressed view: let plotly auto-generate cleaner axis labels
      compressed_range <- range(c(compression_map_isoforms$coords$compressed_start, compression_map_isoforms$coords$compressed_end))
      
      list(
        title = paste0("Genomic Position (compressed view) - chromosome ", chromosome, ")"),
        range = compressed_range,
        zeroline = FALSE
      )
    } else {
      # Original view
      list(
        title = paste0("Genomic Position (chromosome ", chromosome, ")"),
        range = c(gene_start, gene_end + 2000),
        zeroline = FALSE
      )
    }
    
    p <- p %>% plotly::layout(
      title = paste0('All Isoforms Analysis (', miscleavage_label, ') - Highlighted: ', highlight_isoform),
      xaxis = xaxis_config,
      yaxis = list(
        title = "",
        range = c(0.3, nrow(transcript_df) + 0.7),
        tickvals = transcript_df$y_position,
        ticktext = transcript_df$transcript,
        zeroline = FALSE
      ),
      legend = list(
        orientation = "v", x = 1.02, y = 0.5,
        title = list(text = "Legend:"),
        itemsizing = "constant",
        tracegroupgap = 10,
        bgcolor = "rgba(255, 255, 255, 0.8)",
        bordercolor = "rgba(0, 0, 0, 0.2)",
        borderwidth = 1
      ),
      margin = list(t = 50, b = 50, l = 100, r = 200)
    ) %>%
    plotly::config(displayModeBar = TRUE, modeBarButtonsToRemove = c('lasso2d', 'select2d'))
    
    return(p %>% clean_plotly_hover())
  })
  
  # Render highlighted isoform peptides table with chromosome info and merged peptides
  output$highlighted_isoform_table <- DT::renderDataTable({
    data <- highlighted_isoform_data()
    
    if (is.null(data)) {
      return(data.frame(Message = "No data available"))
    }
    
    highlight_peptides <- data$highlight_peptides
    
    # Create display table with chromosome and genomic positions
    if ("chromosome" %in% names(highlight_peptides) && "genomic_positions" %in% names(highlight_peptides)) {
      display_df <- highlight_peptides[, c("peptide", "chromosome", "genomic_positions", "specificity_category", "other_transcript_count", "other_transcripts")]
      colnames(display_df) <- c("Peptide", "Chromosome", "Genomic Positions", "Specificity", "Other Isoforms Count", "Other Isoforms")
    } else {
      # Fallback if chromosome info not available
    display_df <- highlight_peptides[, c("peptide", "start", "end", "specificity_category", "other_transcript_count", "other_transcripts")]
    colnames(display_df) <- c("Peptide", "Start", "End", "Specificity", "Other Isoforms Count", "Other Isoforms")
    }
    
    DT::datatable(
      display_df,
      selection = 'single',
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(which(names(display_df) == "Specificity") - 1, 'asc'))  # Sort by specificity (0-indexed)
      ),
      rownames = FALSE
    ) %>%
    DT::formatStyle(
      'Specificity',
      backgroundColor = styleEqual(
        c("Unique", "Shared", "Universal"),
        c("#FF0000", "#F39C12", "#2ECC71")
      ),
      color = styleEqual(
        c("Unique", "Shared", "Universal"),
        c("white", "white", "white")
      )
    )
  })
  
  # Store selected peptide from isoform-specific table
  selected_isoform_peptide <- reactiveVal(NULL)
  
  # Observer for isoform-specific table row selection
  observeEvent(input$highlighted_isoform_table_rows_selected, {
    req(highlighted_isoform_data())
    if (length(input$highlighted_isoform_table_rows_selected) > 0) {
      selected_row <- input$highlighted_isoform_table_rows_selected[1]
      data <- highlighted_isoform_data()
      highlight_peptides <- data$highlight_peptides
      if (selected_row <= nrow(highlight_peptides)) {
        selected_isoform_peptide(highlight_peptides[selected_row, ])
      }
    } else {
      selected_isoform_peptide(NULL)
    }
  })
  
  
  # Render isoform coverage analysis table
  output$isoform_coverage_table <- DT::renderDataTable({
    data <- highlighted_isoform_data()
    
    if (is.null(data)) {
      return(data.frame(Message = "No coverage data available"))
    }
    
    coverage_df <- data$isoform_coverage
    
    DT::datatable(
      coverage_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(1, 'desc'))  # Sort by total peptides descending
      ),
      rownames = FALSE
    ) %>%
    DT::formatStyle(
      'Unique_Peptides',
      backgroundColor = styleInterval(
        c(1, 5),
        c("#F8F9FA", "#FFE4E1", "#FF6B6B")
      )
    ) %>%
    DT::formatStyle(
      'Universal_Peptides',
      backgroundColor = styleInterval(
        c(1, 10),
        c("#F8F9FA", "#E8F5E8", "#4CAF50")
      )
    )
  })
  

  
  # Download handler for isoform analysis data
  output$download_isoform_analysis <- downloadHandler(
    filename = function() {
      miscleavage_suffix <- switch(input$miscleavage_type,
        "no_miss_cleavage" = "no_miss",
        "upto_two_misscleavage" = "upto_2_miss"
      )
      paste0("isoform_analysis_", input$highlight_isoform, "_", 
            miscleavage_suffix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      data <- highlighted_isoform_data()
      if (!is.null(data)) {
        highlight_peptides <- data$highlight_peptides
        
        # Add metadata
        highlight_peptides$gene <- input$gene
        highlight_peptides$highlight_isoform <- input$highlight_isoform
        highlight_peptides$protease <- input$protease
        highlight_peptides$miscleavage_type <- input$miscleavage_type
        
        write.csv(highlight_peptides, file, row.names = FALSE)
      }
    }
  )
  
  #===============================================================================
  # PEPTIDE SEARCH TAB
  #===============================================================================
  
  # Load peptide search server logic
  source("R/server/peptide_search_server.R", local = TRUE)
  
  # Load multi-enzyme coverage server logic
  source("R/server/multi_enzyme_server.R", local = TRUE)
  
  # Load multi-isoform multi-enzyme matrix analysis server logic
  source("R/server/multi_isoform_enzyme_matrix_server.R", local = TRUE)
  
  # Load rMATS peptide analysis server logic
  source("R/server/rmats_peptide_server.R", local = TRUE)
  
  # Alternative splicing analysis available through rMATS interface
  
  #===============================================================================
  # NOVEL ISOFORM DISCOVERY SERVER LOGIC
  #===============================================================================
  # DISABLED - Functionality moved to novel_isoform_analysis_module.R to avoid conflicts
  # source("R/server/novel_isoform_server.R", local = TRUE)
  
  # ============================================================================
  # NOVEL MULTI-ISOFORM COMPARATIVE ANALYSIS
  # ============================================================================
  source("R/server/novel_multi_isoform_server.R", local = TRUE)
  
  #===============================================================================
  # NEW AS PEPTIDE MAPPING SERVER LOGIC
  #===============================================================================
  
  # Source new AS peptide mapping functions
  
  
  # Reactive values for new AS peptide mapping
  new_as_mapping_results <- reactiveVal(NULL)
  new_as_mapping_processing <- reactiveVal(FALSE)
  
  # Run new AS peptide mapping
  # Note: Alternative splicing peptide mapping available through other analysis modes

  #===============================================================================
  # SPLADDER ANALYSIS MODULE (Reusing 100% of rMATS logic)
  #===============================================================================
  
  # Source SplAdder GFF3 converter
  source("R/spladder_gff3_converter.R", local = TRUE)
  
  # Reactive values for SplAdder analysis (copying rMATS structure)
  spladder_pipeline_results <- reactiveVal(NULL)
  spladder_isoform_data <- reactiveVal(NULL)
  spladder_merged_data <- reactiveVal(NULL)
  spladder_multi_isoform_data <- reactiveVal(NULL)
  
  # SplAdder pipeline state (copying rMATS structure)
  spladder_pipeline_state <- reactiveValues(
    spladder_data = NULL,
    event_type = NULL,
    analysis_results = NULL,
    pipeline_completed = FALSE
  )
  
  # File upload validation
  output$spladder_file_uploaded <- reactive({
    # Only return TRUE when file has actual content AND event type selected
    if (is.null(input$spladder_file) || is.null(input$spladder_file$datapath) || 
        !file.exists(input$spladder_file$datapath) ||
        is.null(input$spladder_event_type) || input$spladder_event_type == "") {
      return(FALSE)
    }
    return(TRUE)
  })
  outputOptions(output, "spladder_file_uploaded", suspendWhenHidden = FALSE)
  
  # Event type display
  output$spladder_event_type_display2 <- renderText({
    switch(input$spladder_event_type,
           "exon_skip" = "SE (Skipped Exon)",
           "alt_3prime" = "A3SS (Alternative 3' Splice Site)",
           "alt_5prime" = "A5SS (Alternative 5' Splice Site)",
           "mutex_exons" = "MXE (Mutually Exclusive Exons)",
           "intron_retention" = "RI (Retained Intron)",
           "")
  })
  
  output$spladder_expected_format <- renderText({
    "SplAdder GFF3 files contain hierarchical gene → mRNA → exon structures that will be converted to rMATS-compatible coordinate format for analysis."
  })
  
  # Step 1: Parse SplAdder events (THE ONLY DIFFERENT STEP)
  observeEvent(input$spladder_parse_events, {
    req(input$spladder_file)
    req(input$spladder_event_type != "")
    
    tryCatch({
      # Parse SplAdder GFF3 file and convert to rMATS format
      cat("🔄 Converting SplAdder GFF3 to rMATS format...\n")
      gff3_file <- input$spladder_file$datapath
      rmats_format_data <- parse_spladder_gff3(gff3_file)
      
      if (length(rmats_format_data) == 0) {
        stop("No events found in the uploaded file")
      }
      
      # Extract the appropriate event type based on what was found
      if (input$spladder_event_type == "exon_skip" && "SE" %in% names(rmats_format_data)) {
        selected_events <- rmats_format_data[["SE"]]
        event_type <- "SE"
      } else if ("A3SS" %in% names(rmats_format_data)) {
        selected_events <- rmats_format_data[["A3SS"]]
        event_type <- "A3SS"
      } else if ("A5SS" %in% names(rmats_format_data)) {
        selected_events <- rmats_format_data[["A5SS"]]
        event_type <- "A5SS"
      } else if ("MXE" %in% names(rmats_format_data)) {
        selected_events <- rmats_format_data[["MXE"]]
        event_type <- "MXE"
      } else if ("RI" %in% names(rmats_format_data)) {
        selected_events <- rmats_format_data[["RI"]]
        event_type <- "RI"
      } else {
        stop("No compatible event types found in the uploaded file")
      }
      
      # Store in pipeline state (now in rMATS format!)
      spladder_pipeline_state$spladder_data <- selected_events
      spladder_pipeline_state$event_type <- event_type
      
      cat("✅ Successfully converted", nrow(selected_events), "SplAdder", event_type, "events to rMATS format\n")
      
    }, error = function(e) {
      cat("❌ Error parsing SplAdder GFF3 file:", conditionMessage(e), "\n")
      showNotification(paste("Error parsing file:", conditionMessage(e)), type = "error")
    })
  })
  
  # Step 1 completion indicator
  output$spladder_step1_completed <- reactive({
    return(!is.null(spladder_pipeline_state$spladder_data) && 
           !is.null(spladder_pipeline_state$analysis_results) &&
           spladder_pipeline_state$pipeline_completed == TRUE)
  })
  outputOptions(output, "spladder_step1_completed", suspendWhenHidden = FALSE)
  
  # Events table (identical to rMATS)
  output$spladder_events_table <- DT::renderDataTable({
    req(spladder_pipeline_state$spladder_data)
    
    display_cols <- c("ID", "GeneID", "geneSymbol", "chr", "strand")
    available_cols <- intersect(display_cols, colnames(spladder_pipeline_state$spladder_data))
    
    # FIX: Ensure all data is atomic (not lists) before rendering
    spladder_data <- spladder_pipeline_state$spladder_data
    
    # Convert any list columns to atomic values
    for (col in available_cols) {
      if (is.list(spladder_data[[col]])) {
        cat("DEBUG: Converting list column", col, "to atomic\n")
        spladder_data[[col]] <- sapply(spladder_data[[col]], function(x) {
          if (length(x) == 1) return(as.character(x))
          else return(paste(x, collapse = ","))
        })
      }
    }
    
    DT::datatable(
      spladder_data[, available_cols, drop = FALSE],
      selection = "single",
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'rtip'
      ),
      rownames = FALSE
    )
  })
  
  # Step 1 results summary (identical to rMATS)
  output$spladder_step1_results <- renderText({
    req(spladder_pipeline_state$spladder_data)
    
    data <- spladder_pipeline_state$spladder_data
    paste("📊 Parsed", nrow(data), "SplAdder events converted to rMATS format")
  })
  
  # Selected event information (identical to rMATS)
  output$spladder_selected_event_info <- renderText({
    req(input$spladder_events_table_rows_selected)
    req(spladder_pipeline_state$spladder_data)
    
    selected_row <- input$spladder_events_table_rows_selected
    selected_event <- spladder_pipeline_state$spladder_data[selected_row, ]
    
    # FIX: Ensure selected_event has atomic values
    selected_event[] <- lapply(selected_event, function(x) {
      if (is.list(x) && length(x) == 1) return(x[[1]])
      return(x)
    })
    
    paste("🎯 Selected Event:", selected_event$ID,
          "\n🧬 Gene:", selected_event$GeneID, "Symbol:", selected_event$geneSymbol,
          "\n📍 Location:", selected_event$chr, selected_event$strand)
  })
  
  # COMPREHENSIVE ANALYSIS - Uses EXACT same rMATS logic!
  observeEvent(input$run_spladder_comprehensive_analysis, {
    req(input$spladder_events_table_rows_selected, spladder_pipeline_state$spladder_data)
    
    selected_row <- input$spladder_events_table_rows_selected
    if (length(selected_row) == 0) return()
    
    rmats_data <- spladder_pipeline_state$spladder_data
    selected_event <- rmats_data[selected_row, ]
    event_type <- spladder_pipeline_state$event_type
    
    # FIX: Ensure selected_event has atomic values before comprehensive analysis
    selected_event[] <- lapply(selected_event, function(x) {
      if (is.list(x) && length(x) == 1) return(x[[1]])
      return(x)
    })
    
    cat("DEBUG: About to run comprehensive analysis on event with data types:\n")
    str(selected_event)
    
    withProgress(message = 'Running comprehensive SplAdder analysis...', value = 0, {
      tryCatch({
        # Use EXACT same comprehensive analysis as rMATS!
        analysis_results <- analyze_rmats_event_comprehensive(selected_event, event_type)
        
        if (analysis_results$success) {
          # Store pipeline results for visualization system
          spladder_pipeline_results(analysis_results$pipeline_results)
          
          # Store analysis results in pipeline state
          spladder_pipeline_state$analysis_results <- analysis_results
          spladder_pipeline_state$pipeline_completed <- TRUE
          
          # Load SplAdder isoform data from analysis results (same as rMATS)
          cat("Loading 24-column peptide file from analysis results:", analysis_results$pipeline_results$dataframe_file, "\n")
          if (!is.null(analysis_results$pipeline_results$dataframe_file) && 
              file.exists(analysis_results$pipeline_results$dataframe_file)) {
            spladder_data_24col <- readRDS(analysis_results$pipeline_results$dataframe_file)
          } else {
            showNotification("No SplAdder data file found in analysis results", type = "error")
            return()
          }
          
          if (!is.null(spladder_data_24col)) {
            spladder_isoform_data(spladder_data_24col)
            
            cat("✅ Loaded data contains gene IDs:", paste(unique(spladder_data_24col$geneID), collapse = ", "), "\n")
            
            # Use SAME gene merging logic as rMATS
            gene_id <- extract_gene_id_from_rmats_data(spladder_data_24col)
            
            if (!is.null(gene_id)) {
              cat("🔄 Loading existing gene transcripts for", gene_id, "to merge with SplAdder isoforms...\n")
              
              tryCatch({
                merged_data <- load_and_merge_gene_data(
                  gene_id = gene_id,
                  novel_data = spladder_data_24col,  # SplAdder inclusion/exclusion data
                  miscleavage_type = input$spladder_miscleavage_type,
                  rds_dir = "data/genes"
                )
                
                spladder_merged_data(merged_data)
                # CRITICAL FIX: Also store in rmats_merged_data for visualization system
                rmats_merged_data(merged_data)
                
                cat("✅ Successfully merged SplAdder data with existing gene data\n")
                
                # Update dropdown choices with ALL transcripts (SplAdder + existing)
                transcript_choices <- unique(merged_data$txID)
                spladder_transcripts <- transcript_choices[grepl("\\.(inclusion|exclusion)$", transcript_choices)]
                existing_transcripts <- transcript_choices[!grepl("\\.(inclusion|exclusion)$", transcript_choices)]
                
                updateSelectInput(session, "spladder_highlight_isoform", 
                                 choices = setNames(transcript_choices, transcript_choices),
                                 selected = spladder_transcripts[1])
                
                updateSelectizeInput(session, "spladder_compare_isoforms",
                                    choices = setNames(transcript_choices, transcript_choices),
                                    selected = transcript_choices)
                
              }, error = function(e) {
                cat("Warning: Could not merge with existing gene data:", e$message, "\n")
                cat("Falling back to SplAdder data only\n")
                
                spladder_merged_data(spladder_data_24col)
                # CRITICAL FIX: Also store in rmats_merged_data for visualization system
                rmats_merged_data(spladder_data_24col)
                
                transcript_choices <- unique(spladder_data_24col$txID)
                updateSelectInput(session, "spladder_highlight_isoform", 
                                 choices = setNames(transcript_choices, transcript_choices),
                                 selected = transcript_choices[1])
                
                updateSelectizeInput(session, "spladder_compare_isoforms",
                                    choices = setNames(transcript_choices, transcript_choices),
                                    selected = transcript_choices)
              })
            } else {
              cat("Warning: Could not extract gene ID from SplAdder data\n")
              spladder_merged_data(spladder_data_24col)
              # CRITICAL FIX: Also store in rmats_merged_data for visualization system
              rmats_merged_data(spladder_data_24col)
              
              transcript_choices <- unique(spladder_data_24col$txID)
              updateSelectInput(session, "spladder_highlight_isoform", 
                               choices = setNames(transcript_choices, transcript_choices),
                               selected = transcript_choices[1])
              
              updateSelectizeInput(session, "spladder_compare_isoforms",
                                  choices = setNames(transcript_choices, transcript_choices),
                                  selected = transcript_choices)
            }
          }
          
          showNotification("SplAdder comprehensive analysis completed successfully!", type = "message")
          
        } else {
          showNotification(paste("SplAdder analysis failed:", analysis_results$error), type = "error")
        }
        
      }, error = function(e) {
        cat("❌ Error in SplAdder comprehensive analysis:", conditionMessage(e), "\n")
        showNotification(paste("Analysis failed:", conditionMessage(e)), type = "error")
      })
    })
  })
  
  # Comprehensive analysis completion indicator
  output$spladder_comprehensive_completed <- reactive({
    return(spladder_pipeline_state$pipeline_completed)
  })
  outputOptions(output, "spladder_comprehensive_completed", suspendWhenHidden = FALSE)
  
  # Comprehensive analysis status
  output$spladder_comprehensive_status <- renderText({
    if (!spladder_pipeline_state$pipeline_completed) {
      "⏳ Click to start comprehensive analysis..."
    } else {
      results <- spladder_pipeline_state$analysis_results
      paste("✅ Analysis completed successfully!",
            paste("📊 Generated", nrow(results$peptide_data), "protein isoforms"),
            paste("📁 Results saved in:", results$pipeline_results$output_dir),
            sep = "\n")
    }
  })
  
  #===============================================================================
  # SPLADDER GENE LOADING AND MERGING LOGIC (IDENTICAL TO rMATS)
  #===============================================================================

  # SplAdder Load Gene Data observer - enables comparative analysis
  observeEvent(input$load_spladder_gene, {
    req(spladder_pipeline_state$pipeline_completed, spladder_pipeline_state$analysis_results$success)
    
    withProgress(message = 'Loading SplAdder gene data for comparative analysis...', value = 0, {
      tryCatch({
        incProgress(0.2, detail = "Extracting SplAdder data...")
        
        # Get results from comprehensive analysis
        analysis_results <- spladder_pipeline_state$analysis_results
        gene_id <- analysis_results$pipeline_results$event_info$gene_id
        gene_symbol <- analysis_results$pipeline_results$event_info$gene_symbol
        
        cat("Loading SplAdder gene data for:", gene_symbol, "(", gene_id, ")\n")
        
        # Extract SplAdder isoform data from temp files
        incProgress(0.5, detail = "Processing SplAdder isoforms...")
        pipeline_results <- analysis_results$pipeline_results
        
        # Load the RDS data generated by the comprehensive analysis
        if (file.exists(pipeline_results$dataframe_file)) {
          spladder_rds_data <- readRDS(pipeline_results$dataframe_file)
          
          # Store in spladder_isoform_data for access by other reactives
          spladder_isoform_data(spladder_rds_data)
          
          cat("✅ Loaded", nrow(spladder_rds_data), "SplAdder isoform entries from comprehensive analysis\n")
          
          # Merge with existing gene data if available (same approach as rMATS)
          incProgress(0.8, detail = "Merging with existing gene data...")
          
          tryCatch({
            merged_data <- load_and_merge_gene_data_with_spladder_isoforms(
              novel_data = spladder_rds_data,
              miscleavage_type = input$spladder_miscleavage_type,
              gene_id = gene_id,
              gene_symbol = gene_symbol
            )
            
            spladder_merged_data(merged_data)
            
            # Update transcript selection choices with ALL transcripts (SplAdder + existing)
            transcript_choices <- unique(merged_data$txID)
            spladder_transcripts <- transcript_choices[grepl("\\.(inclusion|exclusion)$", transcript_choices)]
            
            # Update UI dropdowns
            updateSelectInput(session, "spladder_highlight_isoform", 
                             choices = transcript_choices,
                             selected = spladder_transcripts[1])
            
            updateSelectizeInput(session, "spladder_compare_isoforms",
                                choices = transcript_choices,
                                selected = spladder_transcripts[1:min(3, length(spladder_transcripts))])
            
            cat("✅ Successfully merged SplAdder data with existing gene data\n")
            cat("📊 Total transcripts available:", length(transcript_choices), "\n")
            cat("🧬 SplAdder transcripts:", length(spladder_transcripts), "\n")
            
          }, error = function(merge_error) {
            cat("Warning: Could not merge with existing gene data:", merge_error$message, "\n")
            cat("Using SplAdder data only\n")
            
            # Fallback: use SplAdder data only
            spladder_merged_data(spladder_rds_data)
            
            transcript_choices <- unique(spladder_rds_data$txID)
            updateSelectInput(session, "spladder_highlight_isoform", 
                             choices = transcript_choices,
                             selected = transcript_choices[1])
            
            updateSelectizeInput(session, "spladder_compare_isoforms",
                                choices = transcript_choices,
                                selected = transcript_choices[1:min(3, length(transcript_choices))])
          })
          
        } else {
          stop("SplAdder dataframe file not found: ", pipeline_results$dataframe_file)
        }
        
        incProgress(1)
        showNotification("SplAdder gene data loaded successfully! You can now run comparative analysis.", type = "message")
        
      }, error = function(e) {
        cat("❌ Error loading SplAdder gene data:", e$message, "\n")
        showNotification(paste("Failed to load gene data:", e$message), type = "error")
      })
    })
  })

  #===============================================================================
  # SPLADDER MULTI-ISOFORM ANALYSIS LOGIC (COPYING rMATS PATTERN)
  #===============================================================================

  # SplAdder Multi-Isoform data loader (copied from rMATS system)
  spladder_multi_isoform_data <- reactive({
    cat("🔍 DEBUGGING SPLADDER MULTI-ISOFORM REQUIREMENTS:\n")
    cat("  run_spladder_comparative_analysis:", input$run_spladder_comparative_analysis, "\n")
    cat("  spladder_protease:", input$spladder_protease, "\n")
    cat("  spladder_miscleavage_type:", input$spladder_miscleavage_type, "\n")
    cat("  spladder_merged_data() exists:", !is.null(spladder_merged_data()), "\n")
    cat("  spladder_compare_isoforms:", paste(input$spladder_compare_isoforms, collapse = ", "), "\n")
    cat("  spladder_compare_isoforms length:", length(input$spladder_compare_isoforms), "\n")
    cat("  Length >= 2:", length(input$spladder_compare_isoforms) >= 2, "\n")
    cat("  Length <= 8:", length(input$spladder_compare_isoforms) <= 8, "\n")
    
    req(input$run_spladder_comparative_analysis > 0, input$spladder_protease, input$spladder_miscleavage_type, 
        spladder_merged_data(), input$spladder_compare_isoforms, 
        length(input$spladder_compare_isoforms) >= 2, length(input$spladder_compare_isoforms) <= 8)
    
    withProgress(message = 'Loading SplAdder multi-isoform analysis...', value = 0, {
      protease <- input$spladder_protease
      miscleavage_type <- input$spladder_miscleavage_type
      selected_transcripts <- input$spladder_compare_isoforms
      merged_data <- spladder_merged_data()
      
      cat("Selected transcripts for SplAdder multi-isoform analysis:", paste(selected_transcripts, collapse=", "), "\n")
      cat("Using enzyme:", protease, "and miscleavage type:", miscleavage_type, "\n")
      
      incProgress(0.2, detail = 'Loading peptides for selected transcripts...')
      
      # Debug: Check what data we actually have
      cat("DEBUG: spladder_merged_data() structure:\n")
      cat("DEBUG: Number of rows:", nrow(merged_data), "\n")
      cat("DEBUG: Column names:", paste(names(merged_data), collapse=", "), "\n")
      cat("DEBUG: First few column names containing 'mapped':", paste(names(merged_data)[grepl("mapped", names(merged_data))], collapse=", "), "\n")
      
      # Create vis_data structure for compatibility with existing functions
      vis_data <- list(
        genes = c(),
        gene_symbols = c(),
        gene_lookup = c(),
        proteases = c("trp", "chymo", "aspn", "lysc", "lysn", "gluc"),
        original_peptides = merged_data
      )
      
      # Get peptides for SELECTED transcripts only using existing function
      all_peptides_list <- list()
      for (i in seq_along(selected_transcripts)) {
        tx <- selected_transcripts[i]
        tx_peptides <- get_transcript_peptides_for_comparison(tx, vis_data, protease)
        if (!is.null(tx_peptides) && length(tx_peptides) > 0) {
          all_peptides_list[[tx]] <- data.frame(
            transcript = tx,
            y_position = i,
            start = start(tx_peptides),
            end = end(tx_peptides),
            chromosome = as.character(seqnames(tx_peptides)),
            peptide = tx_peptides$peptide,
            stringsAsFactors = FALSE
          )
        }
      }
      
      if (length(all_peptides_list) == 0) {
        return(NULL)
      }
      
      # Combine all peptides
      all_peptides_df <- do.call(rbind, all_peptides_list)
      rownames(all_peptides_df) <- NULL
      
      # Apply 6-60 AA length filter consistently
      all_peptides_df$aa_length <- nchar(as.character(all_peptides_df$peptide))
      all_peptides_df <- all_peptides_df[all_peptides_df$aa_length >= 6 & all_peptides_df$aa_length <= 60, ]
      
      if (nrow(all_peptides_df) == 0) {
        showNotification("No peptides in 6-60 AA range found for selected transcripts", type = "warning")
        return(NULL)
      }
      
      incProgress(0.8, detail = 'Calculating peptide specificity...')
      
      # Calculate specificity for ALL peptides across selected transcripts (EXACT COPY from rMATS system)
      all_peptides_df$specificity_category <- ""
      
      for (i in 1:nrow(all_peptides_df)) {
        peptide_seq <- all_peptides_df$peptide[i]
        current_transcript <- all_peptides_df$transcript[i]
        
        # Count transcripts (excluding current) that have this peptide
        other_tx_with_peptide <- unique(all_peptides_df$transcript[
          all_peptides_df$peptide == peptide_seq & all_peptides_df$transcript != current_transcript
        ])
        
        total_selected_isoforms <- length(selected_transcripts)
        other_count <- length(other_tx_with_peptide)
        
        if (other_count == 0) {
          all_peptides_df$specificity_category[i] <- "Unique"
        } else if (other_count == (total_selected_isoforms - 1)) {
          all_peptides_df$specificity_category[i] <- "Universal"
        } else {
          all_peptides_df$specificity_category[i] <- "Shared"
        }
      }
      
      # Set up transcript positions and boundaries
      transcript_df <- data.frame(
        transcript = selected_transcripts,
        y_position = seq_along(selected_transcripts),
        stringsAsFactors = FALSE
      )
      
      # Calculate gene boundaries from all peptides
      gene_start <- min(all_peptides_df$start)
      gene_end <- max(all_peptides_df$end)
      
      incProgress(1)
      
      return(list(
        all_peptides = all_peptides_df,
        transcript_df = transcript_df,
        gene_start = gene_start,
        gene_end = gene_end,
        all_transcripts = selected_transcripts
      ))
    })
  })

  # SplAdder Multi-Isoform highlighted data (copied from rMATS system)
  spladder_multi_isoform_highlighted_data <- reactive({
    req(spladder_multi_isoform_data(), input$spladder_highlight_isoform)
    
    base_data <- spladder_multi_isoform_data()
    highlight_isoform <- input$spladder_highlight_isoform
    
    if (is.null(base_data) || is.null(highlight_isoform) || highlight_isoform == "") {
      return(base_data)
    }
    
    # Return the data with specificity already calculated
    all_peptides_df <- base_data$all_peptides
    
    # Update hover text with specificity
    all_peptides_df$hover_text <- paste0(
      "Peptide: ", all_peptides_df$peptide,
      "<br>Position: ", all_peptides_df$start, "-", all_peptides_df$end,
      "<br>Transcript: ", all_peptides_df$transcript,
      "<br>Specificity: ", all_peptides_df$specificity_category,
      "<br>Enzyme: ", input$spladder_protease,
      "<br>Miscleavage: ", input$spladder_miscleavage_type
    )
    
    return(list(
      all_peptides = all_peptides_df,
      transcript_df = base_data$transcript_df,
      gene_start = base_data$gene_start,
      gene_end = base_data$gene_end,
      all_transcripts = base_data$all_transcripts
    ))
  })

  # Render SplAdder comparative plot (EXACT COPY of rMATS plot logic)
  output$spladder_comparative_plot <- renderPlotly({
    # Pure SplAdder renderer - only uses SplAdder data
    req(spladder_multi_isoform_highlighted_data())
    
    data <- spladder_multi_isoform_highlighted_data()
    pipeline_results <- spladder_pipeline_results()
    cat("Using SplAdder data for visualization\n")
    
    tryCatch({
      
      if (is.null(data) || is.null(data$all_peptides) || nrow(data$all_peptides) == 0) {
        p <- plotly::plot_ly() %>%
          plotly::add_annotations(
            x = 0.5, y = 0.5,
            text = "No comparative data available. Please select 2-8 isoforms and run analysis.",
            showarrow = FALSE,
            font = list(size = 16)
          ) %>%
          plotly::layout(
            xaxis = list(range = c(0, 1), showticklabels = FALSE),
            yaxis = list(range = c(0, 1), showticklabels = FALSE)
          )
        return(p)
      }
      
      all_peptides_df <- data$all_peptides
      transcript_df <- data$transcript_df
      gene_start <- data$gene_start
      gene_end <- data$gene_end
      all_transcripts <- data$all_transcripts
      
      # Separate SplAdder transcripts from reference transcripts
      spladder_transcripts <- all_transcripts[grepl("\\.inclusion$|\\.exclusion$", all_transcripts)]
      reference_transcripts <- all_transcripts[!grepl("\\.inclusion$|\\.exclusion$", all_transcripts)]
      
      cat("SplAdder transcripts:", paste(spladder_transcripts, collapse = ", "), "\n")
      cat("Reference transcripts:", paste(reference_transcripts, collapse = ", "), "\n")
      
      # Initialize combined exon/CDS data
      exons_by_transcript <- list()
      cds_by_transcript <- list()
      
      # Load SplAdder GTF data for SplAdder transcripts (if available)
      spladder_results <- spladder_pipeline_results()
      cat("=== SPLADDER GTF LOADING DEBUG ===\n")
      cat("spladder_results is null:", is.null(spladder_results), "\n")
      if (!is.null(spladder_results)) {
        cat("spladder_results$gtf_file:", spladder_results$gtf_file, "\n")
        cat("GTF file exists:", file.exists(spladder_results$gtf_file), "\n")
      }
      cat("spladder_transcripts length:", length(spladder_transcripts), "\n")
      cat("spladder_transcripts:", paste(spladder_transcripts, collapse = ", "), "\n")
      
      if (!is.null(spladder_results) && !is.null(spladder_results$gtf_file) && 
          file.exists(spladder_results$gtf_file) && length(spladder_transcripts) > 0) {
        
        spladder_gtf_file <- spladder_results$gtf_file
        cat("🔍 DEBUGGING GTF LOADING:\n")
        cat("GTF file path:", spladder_gtf_file, "\n")
        cat("Looking for transcripts:", paste(spladder_transcripts, collapse = ", "), "\n")
        
        # CHECK GTF FILE CONTENTS FIRST
        cat("📄 Checking GTF file contents...\n")
        tryCatch({
          gtf_lines <- readLines(spladder_gtf_file, n = 50)  # Read first 50 lines
          transcript_lines <- gtf_lines[grepl("transcript_id", gtf_lines)]
          cat("Found", length(transcript_lines), "lines with transcript_id in first 50 lines\n")
          if (length(transcript_lines) > 0) {
            cat("Sample transcript lines:\n")
            cat(head(transcript_lines, 3), sep = "\n")
          }
          
          # Look for our specific transcript names
          for (tx in spladder_transcripts) {
            matching_lines <- gtf_lines[grepl(tx, gtf_lines, fixed = TRUE)]
            cat("Lines containing '", tx, "':", length(matching_lines), "\n")
            if (length(matching_lines) > 0) {
              cat("Sample:", head(matching_lines, 1), "\n")
            }
          }
        }, error = function(e) {
          cat("Error reading GTF file:", e$message, "\n")
        })
        
        cat("📋 Now calling load_rmats_transcript_exons...\n")
        spladder_exons_result <- load_rmats_transcript_exons(spladder_gtf_file, spladder_transcripts)
        
        cat("🔍 RESULT ANALYSIS:\n")
        cat("Success:", spladder_exons_result$success, "\n")
        cat("Message:", spladder_exons_result$message, "\n")
        cat("Exons found:", length(spladder_exons_result$exons), "\n")
        cat("CDS found:", length(spladder_exons_result$cds), "\n")
        if (length(spladder_exons_result$exons) > 0) {
          cat("Exon transcript names:", paste(names(spladder_exons_result$exons), collapse = ", "), "\n")
        }
        if (length(spladder_exons_result$cds) > 0) {
          cat("CDS transcript names:", paste(names(spladder_exons_result$cds), collapse = ", "), "\n")
        }
        
        if (spladder_exons_result$success && length(spladder_exons_result$exons) > 0) {
          # Add SplAdder exons and CDS to combined data
          exons_by_transcript <- c(exons_by_transcript, spladder_exons_result$exons)
          cds_by_transcript <- c(cds_by_transcript, spladder_exons_result$cds)
          cat("✅ Loaded", length(spladder_exons_result$exons), "SplAdder transcript exon sets\n")
        } else {
          cat("❌ Failed to load SplAdder exons:", spladder_exons_result$message, "\n")
        }
      } else {
        cat("❌ SPLADDER GTF LOADING CONDITIONS NOT MET:\n")
        if (length(spladder_transcripts) > 0) {
          cat("SplAdder GTF file not available for transcripts:", paste(spladder_transcripts, collapse = ", "), "\n")
        } else {
          cat("No SplAdder transcripts to load\n")
        }
      }
      cat("=== END SPLADDER GTF DEBUG ===\n")
      
      # Load reference GTF data for reference transcripts (if available)
      if (length(reference_transcripts) > 0) {
        cat("Loading reference GTF for transcripts:", paste(reference_transcripts, collapse = ", "), "\n")
        
        # Use cached GTF system if available
        if (dir.exists("data/gtf_cache")) {
          # Extract gene ID from SplAdder data to find GTF cache file
          gene_id <- extract_gene_id_from_rmats_data(spladder_merged_data())
          if (!is.null(gene_id)) {
            
            # Load GTF cache using version-agnostic lookup
            cache_file <- find_gtf_cache_file(gene_id)
            
            if (!is.null(cache_file)) {
              tryCatch({
                gene_details <- readRDS(cache_file)
                
                if (!is.null(gene_details) && 
                    !is.null(gene_details$exons_by_transcript) && 
                    !is.null(gene_details$cds_by_transcript)) {
                  
                  cat("DEBUG: GTF cache loaded successfully for gene", gene_id, "\n")
                  cat("DEBUG: Available transcript IDs from cache:", paste(gene_details$transcript_ids, collapse = ", "), "\n")
                  cat("DEBUG: Available CDS transcript names:", paste(names(gene_details$cds_by_transcript), collapse = ", "), "\n")
                  cat("DEBUG: Available exon transcript names:", paste(names(gene_details$exons_by_transcript), collapse = ", "), "\n")
                  cat("DEBUG: Looking for reference transcripts:", paste(reference_transcripts, collapse = ", "), "\n")
                  
                  # Handle GenomicRanges objects properly
                  # Filter exons for reference transcripts (if names exist)
                  ref_exons <- list()
                  if (!is.null(names(gene_details$exons_by_transcript))) {
                    available_exon_tx <- names(gene_details$exons_by_transcript)
                    matching_exon_tx <- available_exon_tx[available_exon_tx %in% reference_transcripts]
                    if (length(matching_exon_tx) > 0) {
                      ref_exons <- gene_details$exons_by_transcript[matching_exon_tx]
                      # Convert to list format for compatibility
                      ref_exons <- as.list(ref_exons)
                    }
                  } else {
                    cat("WARNING: Exon data has no transcript names - may need different approach\n")
                  }
                  
                  # Filter CDS for reference transcripts
                  ref_cds <- list()
                  if (!is.null(names(gene_details$cds_by_transcript))) {
                    available_cds_tx <- names(gene_details$cds_by_transcript)
                    matching_cds_tx <- available_cds_tx[available_cds_tx %in% reference_transcripts]
                    if (length(matching_cds_tx) > 0) {
                      ref_cds <- gene_details$cds_by_transcript[matching_cds_tx]
                      # Convert to list format for compatibility  
                      ref_cds <- as.list(ref_cds)
                    }
                  }
                  
                  # Add reference exons and CDS to combined data
                  exons_by_transcript <- c(exons_by_transcript, ref_exons)
                  cds_by_transcript <- c(cds_by_transcript, ref_cds)
                  
                  cat("✅ Loaded", length(ref_exons), "reference transcript exon sets from GTF cache\n")
                  cat("✅ Loaded", length(ref_cds), "reference transcript CDS sets from GTF cache\n")
                  
                } else {
                  cat("WARNING: GTF cache structure is invalid for gene", gene_id, "\n")
                }
                
              }, error = function(e) {
                cat("ERROR loading GTF cache:", e$message, "\n")
              })
            } else {
              cat("WARNING: GTF cache file not found:", cache_file, "\n")
            }
          } else {
            cat("WARNING: Could not extract gene ID for GTF cache loading\n")
          }
        } else {
          cat("WARNING: GTF cache directory not found\n")
        }
      }
      
      # Debug: show final combined exon/CDS data
      cat("=== GTF LOADING DEBUG ===\n")
      cat("All transcripts in visualization:", paste(all_transcripts, collapse = ", "), "\n")
      cat("SplAdder transcripts:", paste(spladder_transcripts, collapse = ", "), "\n")
      cat("Reference transcripts:", paste(reference_transcripts, collapse = ", "), "\n")
      cat("Total exon sets loaded:", length(exons_by_transcript), "\n")
      cat("Total CDS sets loaded:", length(cds_by_transcript), "\n")
      cat("Exon transcript IDs:", paste(names(exons_by_transcript), collapse = ", "), "\n")
      cat("CDS transcript IDs:", paste(names(cds_by_transcript), collapse = ", "), "\n")
      cat("========================\n")
      
      # ============================================================================
      # COORDINATE COMPRESSION (DEFAULT FOR SplAdder VISUALIZATIONS)
      # ============================================================================
      
      compression_map <- NULL
      use_compression <- TRUE  # Always use compression for SplAdder visualizations
      
      if (use_compression && length(exons_by_transcript) > 0) {
        # Create compression map from exon data
        cat("Creating coordinate compression map for SplAdder visualization...\n")
        compression_map <- create_compression_map(exons_by_transcript)
        
        # Apply compression to gene boundaries
        gene_start_compressed <- compression_map$compress(gene_start)
        gene_end_compressed <- compression_map$compress(gene_end)
        
        cat("Original gene range:", gene_start, "-", gene_end, "\n")
        cat("Compressed gene range:", gene_start_compressed, "-", gene_end_compressed, "\n")
        
        # Create intron markers for visualization (triangles/lines between exons)
        intron_markers_df <- create_intron_markers_for_transcripts(exons_by_transcript, compression_map, y_position = 0)
        cat("Created", nrow(intron_markers_df), "intron marker points\n")
      } else {
        # No compression - use original coordinates
        gene_start_compressed <- gene_start
        gene_end_compressed <- gene_end
        intron_markers_df <- data.frame()  # No intron markers for uncompressed view
      }
      
      # Extract chromosome information (similar to rMATS approach)
      chromosome <- "Unknown"
      # Try to get chromosome from cached GTF data if available
      if (dir.exists("data/gtf_cache")) {
        gene_id <- extract_gene_id_from_rmats_data(spladder_merged_data())
        if (!is.null(gene_id)) {
          cache_file <- find_gtf_cache_file(gene_id)
          if (!is.null(cache_file)) {
            tryCatch({
              gene_details <- readRDS(cache_file)
              if (!is.null(gene_details$chromosome)) {
                chromosome <- gene_details$chromosome
                cat("✅ Chromosome information:", chromosome, "\n")
              }
            }, error = function(e) {
              cat("Warning: Could not load chromosome info:", e$message, "\n")
            })
          }
        }
      }
      
      # Create the plot (same as rMATS system)
      p <- plotly::plot_ly()
      
      # Color mapping for specificity
      specificity_colors <- c("Unique" = "#FF0000", "Shared" = "#F39C12", "Universal" = "#2ECC71")
      
      # Add clean legend entries once (no clutter)
      for (category in names(specificity_colors)) {
        p <- p %>% plotly::add_trace(
          type = "scatter", mode = "markers",
          x = c(0), y = c(0),
          marker = list(color = specificity_colors[[category]], size = 10),
          showlegend = TRUE,
          name = category,
          legendgroup = paste0("category_", category),
          hoverinfo = "none",
          visible = "legendonly"
        )
      }
      
      # Add transcript lines for all transcripts - UPDATED FOR COMPRESSION
      for (i in 1:nrow(transcript_df)) {
        p <- p %>% plotly::add_trace(
          type = "scatter",
          x = c(gene_start_compressed, gene_end_compressed),
          y = c(transcript_df$y_position[i], transcript_df$y_position[i]),
          mode = "lines",
          line = list(color = "#666666", width = 2),
          showlegend = FALSE,
          hoverinfo = "text",
          text = paste0("Transcript: ", transcript_df$transcript[i])
        )
      }
      
      # Add intron markers (triangle/zigzag lines between exons) - SplAdder FEATURE
      if (nrow(intron_markers_df) > 0) {
        # Group intron markers by transcript and add them
        for (tx_name in unique(intron_markers_df$transcript)) {
          tx_markers <- intron_markers_df[intron_markers_df$transcript == tx_name, ]
          if (nrow(tx_markers) > 0) {
            # Find the y position for this transcript
            tx_y_pos <- transcript_df$y_position[transcript_df$transcript == tx_name]
            if (length(tx_y_pos) > 0) {
              # Add intron markers for this transcript
              for (group_name in unique(tx_markers$group)) {
                group_data <- tx_markers[tx_markers$group == group_name, ]
                if (nrow(group_data) > 0) {
                  p <- p %>% plotly::add_trace(
                    type = "scatter",
                    x = group_data$x,
                    y = rep(tx_y_pos, nrow(group_data)),
                    mode = "lines",
                    line = list(color = "grey50", width = 0.5),
                    showlegend = FALSE,
                    hoverinfo = "text",
                    text = paste0("Intron marker - ", tx_name)
                  )
                }
              }
            }
          }
        }
      }
      
      # Add exon and CDS boundaries if GTF data is available
      if (length(exons_by_transcript) > 0) {
        # Add exon blocks for each transcript (same as rMATS system)
        for (i in 1:nrow(transcript_df)) {
          tx <- transcript_df$transcript[i]
          y_pos <- transcript_df$y_position[i]
          
          # Add exons if available (transparent grey) - UPDATED FOR COMPRESSION
          if (!is.null(exons_by_transcript[[tx]]) && length(exons_by_transcript[[tx]]) > 0) {
            tx_exons <- exons_by_transcript[[tx]]
            for (j in seq_along(tx_exons)) {
              # Apply compression to exon coordinates
              exon_start <- if (!is.null(compression_map)) compression_map$compress(start(tx_exons[j])) else start(tx_exons[j])
              exon_end <- if (!is.null(compression_map)) compression_map$compress(end(tx_exons[j])) else end(tx_exons[j])
              
              p <- p %>% plotly::add_trace(
                type = "scatter", mode = "lines",
                x = c(exon_start, exon_end, exon_end, exon_start, exon_start),
                y = c(y_pos - 0.35, y_pos - 0.35, y_pos + 0.35, y_pos + 0.35, y_pos - 0.35),
                fill = "toself",
                fillcolor = "rgba(211, 211, 211, 0.3)",  # Light grey for exons
                line = list(color = "#D3D3D3", width = 1),
                legendgroup = "structure",
                showlegend = FALSE,
                hoverinfo = "text",
                text = paste0("Exon ", j, " (", start(tx_exons[j]), "-", end(tx_exons[j]), ")<br>Transcript: ", tx)
              )
            }
          }
          
          # Add CDS overlay if available (yellow) - UPDATED FOR COMPRESSION
          if (!is.null(cds_by_transcript[[tx]]) && length(cds_by_transcript[[tx]]) > 0) {
            tx_cds <- cds_by_transcript[[tx]]
            for (j in seq_along(tx_cds)) {
              # Apply compression to CDS coordinates
              cds_start <- if (!is.null(compression_map)) compression_map$compress(start(tx_cds[j])) else start(tx_cds[j])
              cds_end <- if (!is.null(compression_map)) compression_map$compress(end(tx_cds[j])) else end(tx_cds[j])
              
              p <- p %>% plotly::add_trace(
                type = "scatter", mode = "lines",
                x = c(cds_start, cds_end, cds_end, cds_start, cds_start),
                y = c(y_pos - 0.25, y_pos - 0.25, y_pos + 0.25, y_pos + 0.25, y_pos - 0.25),
                fill = "toself",
                fillcolor = "#F1C40F",  # Gold/yellow for CDS
                line = list(color = "#DAA520", width = 1),
                legendgroup = "structure",
                showlegend = FALSE,
                hoverinfo = "text",
                text = paste0("CDS ", j, " (", start(tx_cds[j]), "-", end(tx_cds[j]), ")<br>Transcript: ", tx)
              )
            }
          }
        }
      }
      
      # Add peptide rectangles with no legend entries (clean)
      for (i in 1:nrow(all_peptides_df)) {
        peptide_row <- all_peptides_df[i, ]
        
        # Get the correct color for this peptide's specificity
        peptide_category <- as.character(peptide_row$specificity_category)
        peptide_color <- specificity_colors[[peptide_category]]
        
        # Debug output
        cat("Peptide", i, "category:", peptide_category, "color:", peptide_color, "\n")
        
        # Fallback to default color if category not found
        if (is.null(peptide_color) || is.na(peptide_color)) {
          peptide_color <- "#CCCCCC"  # Grey fallback
          cat("WARNING: Unknown specificity category:", peptide_category, "\n")
        }
        
        # Apply compression to peptide coordinates
        peptide_start <- if (!is.null(compression_map)) compression_map$compress(peptide_row$start) else peptide_row$start
        peptide_end <- if (!is.null(compression_map)) compression_map$compress(peptide_row$end) else peptide_row$end
        
        p <- p %>% plotly::add_trace(
          type = "scatter",
          mode = "lines",
          x = c(peptide_start, peptide_end, peptide_end, peptide_start, peptide_start),
          y = c(peptide_row$y_position - 0.15, peptide_row$y_position - 0.15, 
                peptide_row$y_position + 0.15, peptide_row$y_position + 0.15, peptide_row$y_position - 0.15),
          fill = "toself",
          fillcolor = peptide_color,
          line = list(color = "black", width = 0.5),
          marker = list(opacity = 0),
          text = peptide_row$hover_text,
          hoverinfo = "text",
          legendgroup = paste0("category_", peptide_category),
          showlegend = FALSE
        )
      }
      
      # Add structure legend entries if exon/CDS data is available
      if (length(exons_by_transcript) > 0) {
        # Add Exons legend entry
        p <- p %>% plotly::add_trace(
          type = "scatter", mode = "markers",
          x = c(0), y = c(0),
          marker = list(color = "rgba(211, 211, 211, 0.8)", size = 12, symbol = "square"),
          showlegend = TRUE,
          name = "Exons",
          legendgroup = "structure_exons",
          hoverinfo = "none",
          visible = "legendonly"
        )
        
        # Add CDS legend entry if any CDS data exists
        if (length(cds_by_transcript) > 0) {
          p <- p %>% plotly::add_trace(
            type = "scatter", mode = "markers",
            x = c(0), y = c(0),
            marker = list(color = "#F1C40F", size = 12, symbol = "square"),
            showlegend = TRUE,
            name = "CDS",
            legendgroup = "structure_cds",
            hoverinfo = "none",
            visible = "legendonly"
          )
        }
      }
      
      # Configure layout
      p <- p %>% plotly::layout(
        title = list(
          text = paste("SplAdder Multi-Isoform Comparative Analysis -", spladder_results$event_info$gene_id, "- Compressed View - Chromosome", chromosome),
          font = list(size = 16)
        ),
        xaxis = list(
          title = paste0("Genomic Position (Compressed Introns) - chromosome ", chromosome),
          range = c(gene_start_compressed, gene_end_compressed),
          showgrid = TRUE,
          zeroline = FALSE
        ),
        yaxis = list(
          title = "Transcripts",
          tickvals = transcript_df$y_position,
          ticktext = transcript_df$transcript,
          range = c(0.5, max(transcript_df$y_position) + 0.5),
          showgrid = TRUE,
          zeroline = FALSE
        ),
        showlegend = TRUE,
        legend = list(
          orientation = "h",
          x = 0,
          y = 1.1,
          bgcolor = 'rgba(255, 255, 255, 0.8)'
        ),
        margin = list(t = 100, b = 80, l = 150, r = 80),
        hovermode = "closest"
      )
      
      return(p)
      
    }, error = function(e) {
      cat("Error in spladder_comparative_plot:", e$message, "\n")
      print(traceback())
      
      # Provide helpful error message to user
      error_msg <- if (grepl("object.*not found", e$message)) {
        "Data not available. Please ensure comprehensive analysis completed successfully."
      } else if (grepl("subscription", e$message)) {
        "Data format issue. Please re-run comprehensive analysis."
      } else {
        paste("Visualization error:", e$message, "\nCheck console for details.")
      }
      
      p <- plotly::plot_ly() %>%
        plotly::add_annotations(
          x = 0.5, y = 0.5,
          text = error_msg,
          showarrow = FALSE,
          font = list(size = 14, color = "red")
        ) %>%
        plotly::layout(
          xaxis = list(range = c(0, 1), showticklabels = FALSE),
          yaxis = list(range = c(0, 1), showticklabels = FALSE)
        )
      return(p)
    })
  })

  # SplAdder comparative analysis observer
  observeEvent(input$run_spladder_comparative_analysis, {
    req(input$spladder_compare_isoforms, length(input$spladder_compare_isoforms) >= 2, 
        length(input$spladder_compare_isoforms) <= 8, spladder_merged_data())
    
    withProgress(message = 'Running SplAdder comparative analysis...', value = 0, {
      selected_isoforms <- input$spladder_compare_isoforms
      
      cat("Running comparative analysis for SplAdder isoforms:", paste(selected_isoforms, collapse = ", "), "\n")
      
      incProgress(0.5, detail = "Processing selected isoforms...")
      
      # The reactive system (spladder_multi_isoform_data) will automatically update
      # when this button is clicked, triggering the plot generation
      
      incProgress(1, detail = "Complete!")
      
      showNotification(
        paste("✅ Comparative analysis completed for", length(selected_isoforms), "isoforms"), 
        type = "message"
      )
    })
  })

  # SplAdder highlighted isoform table - uses EXACT same logic as rMATS
  output$spladder_highlighted_isoform_table <- DT::renderDataTable({
    req(spladder_multi_isoform_highlighted_data(), input$spladder_highlight_isoform)
    
    data <- spladder_multi_isoform_highlighted_data()
    highlight_isoform <- input$spladder_highlight_isoform
    
    if (is.null(data) || is.null(data$all_peptides) || nrow(data$all_peptides) == 0) {
      return(data.frame(Message = "No peptide data available"))
    }
    
    # Filter for highlighted isoform
    highlight_peptides <- data$all_peptides[data$all_peptides$transcript == highlight_isoform, ]
    
    if (nrow(highlight_peptides) == 0) {
      return(data.frame(Message = paste("No peptides found for", highlight_isoform)))
    }
    
    # Display peptides
    display_peptides <- highlight_peptides[, c("peptide", "start", "end", "aa_length")]
    names(display_peptides) <- c("Peptide", "Start", "End", "Length (AA)")
    
    DT::datatable(
      display_peptides,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      rownames = FALSE,
      caption = paste("Peptides for", highlight_isoform)
    )
  }, server = TRUE)

}  # End of server function
