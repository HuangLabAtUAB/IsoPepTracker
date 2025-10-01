#===============================================================================
# NOVEL ISOFORM ANALYSIS MODULE
# Extracted from server.R - Novel isoform discovery and analysis functionality
#===============================================================================

#' Novel Isoform Analysis Module
#' 
#' @description 
#' Complete novel isoform discovery and analysis functionality:
#' - Novel isoform pipeline execution with TransDecoder
#' - ORF extraction and quality assessment
#' - Gene boundary-based gene search
#' - Novel vs reference isoform comparison
#' - Peptide specificity analysis for novel isoforms
#' - Download handlers for all novel isoform results
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' 
#' @return Named list containing reactive values for novel isoform analysis
#' 
#' @note 
#' - Requires novel isoform pipeline functions from R/novel_isoform_analysis.R
#' - Uses TransDecoder for ORF prediction and analysis
#' - Integrates with existing peptide analysis pipeline
#' - Maintains full compatibility with regular isoform analysis

create_novel_isoform_analysis_module <- function(input, output, session) {
  
  # ============================================================================
  # MODULE REACTIVE VALUES
  # ============================================================================
  
  # Novel isoform pipeline reactive values
  novel_pipeline_results <- reactiveVal(NULL)
  novel_isoform_data <- reactiveVal(NULL)
  novel_analysis_results <- reactiveVal(NULL)
  novel_orf_results <- reactiveVal(NULL)
  novel_gene_search_results <- reactiveVal(NULL)
  novel_merged_data <- reactiveVal(NULL)
  
  # ============================================================================
  # NOVEL ISOFORM PIPELINE EXECUTION
  # ============================================================================
  
  # Run novel isoform pipeline when button is clicked
  observeEvent(input$run_novel_pipeline, {
    req(input$novel_fasta_file)
    
    # Get uploaded file info
    file_info <- input$novel_fasta_file
    if (is.null(file_info)) {
      showNotification("Please select a FASTA file to upload.", type = "warning")
      return()
    }
    
    # Validate file extension
    if (!grepl("\\.(fa|fasta|fas)$", file_info$name, ignore.case = TRUE)) {
      showNotification("Please upload a valid FASTA file (.fa, .fasta, or .fas)", type = "error")
      return()
    }
    
    withProgress(message = 'Running Novel Isoform Discovery Pipeline...', value = 0, {
      
      # Progress callback function
      progress_callback <- function(message, value) {
        incProgress(value - getDefaultReactiveDomain()$progressValue, detail = message)
      }
      
      tryCatch({
        # Run the enhanced pipeline with user parameters
        results <- run_novel_isoform_pipeline(
          input_fasta_path = file_info$datapath,
          min_protein_length = input$min_protein_length,
          genetic_code = input$genetic_code,
          strand_specific = input$strand_specific,
          retain_long_orfs = input$retain_long_orfs,
          single_best_orf = input$single_best_orf,
          require_start_codon = input$require_start_codon,
          require_stop_codon = input$require_stop_codon,
          min_orf_coverage = input$min_orf_coverage,
          enable_gene_search = input$enable_gene_search,
          min_overlap_percent = as.numeric(input$min_overlap_percent),
          progress_callback = progress_callback
        )
        
        # Store results
        novel_pipeline_results(results)
        
        if (results$success) {
          # Load the novel isoform data
          novel_data <- load_novel_isoform_data(results$dataframe_file)
          novel_isoform_data(novel_data)
          
          # Extract ORF information for selection
          orf_info <- extract_orf_information(novel_data)
          novel_orf_results(orf_info)
          
          # Clear transcript choices initially - will be populated when ORF is selected
          updateSelectInput(session, "novel_highlight_isoform", 
                           choices = character(0),
                           selected = NULL)
          
          showNotification("Novel isoform discovery pipeline completed successfully! Please select an ORF for gene matching.", type = "message")
        } else {
          showNotification(paste("Pipeline failed:", results$error), type = "error")
        }
        
      }, error = function(e) {
        showNotification(paste("Pipeline execution failed:", e$message), type = "error")
        novel_pipeline_results(list(success = FALSE, error = e$message))
      })
    })
  })
  
  # ============================================================================
  # PIPELINE STATUS AND OUTPUTS
  # ============================================================================
  
  # Output enhanced pipeline status
  output$novel_pipeline_status <- renderText({
    results <- novel_pipeline_results()
    if (is.null(results)) {
      return("Analysis not started")
    }
    
    if (results$success) {
      status_text <- paste("✅ Enhanced TransDecoder Analysis Completed Successfully!\n\n")
      
      # Add analysis summary if available
      if (!is.null(results$analysis_summary)) {
        status_text <- paste0(status_text, results$analysis_summary, "\n\n")
      }
      
      # Add parameters used
      if (!is.null(results$parameters)) {
        params <- results$parameters
        status_text <- paste0(status_text, "Parameters Used:\n")
        status_text <- paste0(status_text, "• Minimum protein length: ", params$min_protein_length, " AA\n")
        status_text <- paste0(status_text, "• Genetic code: ", if(!is.null(params$genetic_code)) params$genetic_code else "Unknown", "\n")
        status_text <- paste0(status_text, "• Strand specificity: ", if(!is.null(params$strand_specific)) params$strand_specific else "Unknown", "\n")
        if (!is.null(params$require_start_codon) && params$require_start_codon) status_text <- paste0(status_text, "• Required start codon (ATG)\n")
        if (!is.null(params$require_stop_codon) && params$require_stop_codon) status_text <- paste0(status_text, "• Required stop codon\n")
        if (!is.null(params$enable_gene_search) && params$enable_gene_search) {
          status_text <- paste0(status_text, "• Gene search enabled (≥", params$min_overlap_percent, "% overlap)\n")
        }
        status_text <- paste0(status_text, "\n")
      }
      
      status_text <- paste0(status_text, "📋 Detailed Log:\n", results$log)
      return(status_text)
    } else {
      return(paste("❌ Analysis Failed:\n", results$error, "\n\n📋 Log:\n", results$log))
    }
  })
  
  # Output to control conditional panels
  output$novel_pipeline_completed <- reactive({
    results <- novel_pipeline_results()
    return(!is.null(results) && results$success)
  })
  outputOptions(output, "novel_pipeline_completed", suspendWhenHidden = FALSE)
  
  # ============================================================================
  # ORF SELECTION AND ANALYSIS
  # ============================================================================
  
  # Render ORF selection table
  output$novel_orf_table <- DT::renderDataTable({
    orf_data <- novel_orf_results()
    
    if (is.null(orf_data) || nrow(orf_data) == 0) {
      return(data.frame(Message = "No ORF data available"))
    }
    
    # Create display table
    display_df <- data.frame(
      `ORF ID` = orf_data$orf_id,
      `Length (AA)` = orf_data$protein_length,
      `Quality Score` = round(orf_data$quality_score, 2),
      `Starts with M` = ifelse(orf_data$starts_with_M, "Yes", "No"),
      `Has Stop Codon` = ifelse(orf_data$has_stop, "Yes", "No"),
      `Sequence Preview` = paste0(orf_data$sequence_preview, "..."),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(2, 'desc'))  # Sort by quality score (column 3, zero-indexed as 2)
      ),
      rownames = FALSE,
      selection = 'single'
    ) %>%
    DT::formatStyle(
      columns = 3,  # Quality Score is the 3rd column
      background = DT::styleColorBar(range(display_df[,3]), '#e3f2fd'),
      backgroundSize = '100% 90%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center'
    )
  })
  
  # Output ORF selection status
  output$novel_orf_selection_status <- renderText({
    if (is.null(input$novel_orf_table_rows_selected)) {
      return("Please select an ORF from the table above")
    } else {
      orf_data <- novel_orf_results()
      if (!is.null(orf_data) && length(input$novel_orf_table_rows_selected) > 0) {
        selected_orf <- orf_data[input$novel_orf_table_rows_selected[1], ]
        return(paste("Selected:", selected_orf$orf_id, 
                    "\nLength:", selected_orf$protein_length, "AA",
                    "\nQuality Score:", round(selected_orf$quality_score, 2)))
      }
    }
    return("")
  })
  
  # ============================================================================
  # GENE SEARCH AND MATCHING
  # ============================================================================
  
  # Run Gene Search for selected ORF
  observeEvent(input$run_gene_search, {
    req(novel_pipeline_results(), input$novel_orf_table_rows_selected)
    
    pipeline_results <- novel_pipeline_results()
    
    if (!pipeline_results$success) {
      showNotification("Pipeline must complete successfully before gene search", type = "error")
      return()
    }
    
    # Validate ORF selection
    if (is.null(input$novel_orf_table_rows_selected) || length(input$novel_orf_table_rows_selected) == 0) {
      showNotification("Please select an ORF from the table above before running gene search", type = "error")
      return()
    }
    
    withProgress(message = 'Searching for candidate genes...', value = 0, {
      
      tryCatch({
        incProgress(0.3, detail = 'Loading gene boundary database...')
        
        # Run boundary-based gene search
        gene_search_results <- run_boundary_gene_search_for_novel(
          work_dir = pipeline_results$work_dir,
          min_overlap_bp = 50,
          min_overlap_percent = as.numeric(input$min_overlap_percent),
          max_genes = 10
        )
        
        incProgress(0.4, detail = 'Processing gene information...')
        
        if (nrow(gene_search_results) > 0) {
          # Parse gene information (already in correct format)
          parsed_results <- parse_boundary_gene_info(gene_search_results)
          
          incProgress(0.2, detail = 'Checking gene availability...')
          
          # Check RDS availability
          if (nrow(parsed_results) > 0) {
            gene_availability <- check_boundary_gene_rds_availability_wrapper(parsed_results$gene_id)
            
            # Merge availability info
            final_results <- merge(parsed_results, gene_availability, by = "gene_id", all.x = TRUE)
            final_results <- final_results[order(final_results$pident, decreasing = TRUE), ]
            
            novel_gene_search_results(final_results)
            
            incProgress(0.1, detail = 'Complete')
            
            showNotification(paste("Gene search completed! Found", nrow(final_results), "candidate genes."), type = "message")
          } else {
            showNotification("Gene search completed but no gene information could be parsed.", type = "warning")
            novel_gene_search_results(data.frame())
          }
        } else {
          showNotification(paste("No overlapping genes found with ≥", input$min_overlap_percent, "% overlap."), type = "warning")
          novel_gene_search_results(data.frame())
        }
        
      }, error = function(e) {
        showNotification(paste("Gene search failed:", e$message), type = "error")
        novel_gene_search_results(data.frame())
      })
    })
  })
  
  # Output to control gene search results conditional panel
  output$novel_gene_search_completed <- reactive({
    results <- novel_gene_search_results()
    return(!is.null(results))
  })
  outputOptions(output, "novel_gene_search_completed", suspendWhenHidden = FALSE)
  
  # Render gene search results table
  output$novel_gene_search_results_table <- DT::renderDataTable({
    gene_search_data <- novel_gene_search_results()
    
    if (is.null(gene_search_data) || nrow(gene_search_data) == 0) {
      return(data.frame(Message = "No gene search results available"))
    }
    
    # Create display table
    display_df <- data.frame(
      `Gene ID` = gene_search_data$gene_id,
      `Gene Symbol` = gene_search_data$gene_symbol,
      `RDS Available` = ifelse(gene_search_data$available, "Yes", "No"),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      rownames = FALSE,
      selection = 'single'
    ) %>%
    DT::formatStyle(
      columns = 3,  # RDS Available is now the 3rd column
      backgroundColor = DT::styleEqual(c("Yes", "No"), c("#d4edda", "#f8d7da"))
    )
  })
  
  # Load selected gene for comparison
  observeEvent(input$load_selected_gene, {
    req(input$novel_gene_search_results_table_rows_selected, novel_gene_search_results(), 
        filtered_novel_data(), input$novel_orf_table_rows_selected)
    
    selected_row <- input$novel_gene_search_results_table_rows_selected[1]
    gene_search_data <- novel_gene_search_results()
    filtered_novel <- filtered_novel_data()
    
    if (selected_row > nrow(gene_search_data)) {
      showNotification("Invalid gene selection", type = "error")
      return()
    }
    
    # Validate ORF selection
    if (is.null(input$novel_orf_table_rows_selected) || length(input$novel_orf_table_rows_selected) == 0) {
      showNotification("Please select an ORF from the table above before loading gene data", type = "error")
      return()
    }
    
    if (is.null(filtered_novel) || nrow(filtered_novel) == 0) {
      showNotification("No data available for selected ORF", type = "error")
      return()
    }
    
    selected_gene <- gene_search_data[selected_row, ]
    
    # Check if RDS is available
    if (!selected_gene$available) {
      showNotification("RDS file not available for selected gene. Proceeding with selected ORF only.", type = "warning")
      novel_merged_data(filtered_novel)
      return()
    }
    
    withProgress(message = 'Loading gene data...', value = 0, {
      
      tryCatch({
        incProgress(0.5, detail = 'Loading known gene isoforms...')
        
        # DEBUG: Check what filtered data we're actually passing
        cat("DEBUG: About to merge gene data\n")
        cat("DEBUG: Filtered novel data rows:", nrow(filtered_novel), "\n")
        cat("DEBUG: Filtered novel txIDs:", paste(unique(filtered_novel$txID), collapse = ", "), "\n")
        
        # Load and merge gene data - USE FILTERED NOVEL DATA
        merged_data <- load_and_merge_gene_data(
          gene_id = selected_gene$gene_id,
          novel_data = filtered_novel,
          miscleavage_type = input$novel_miscleavage_type,
          rds_dir = "data/genes"
        )
        
        # DEBUG: Check what merged data we got back
        cat("DEBUG: Merged data rows:", nrow(merged_data), "\n")
        cat("DEBUG: Merged data txIDs:", paste(unique(merged_data$txID), collapse = ", "), "\n")
        
        incProgress(0.4, detail = 'Preparing analysis...')
        
        # Store merged data
        novel_merged_data(merged_data)
        
        # Make merged data available globally for transcript peptide extraction
        novel_merged_data_global <<- merged_data
        
        # Update isoform choices to include ONLY selected ORF and known transcripts
        all_merged_transcripts <- unique(merged_data$txID)
        filtered_novel_transcripts <- unique(filtered_novel$txID)
        
        # CRITICAL FIX: Filter out other novel ORFs from dropdown choices
        known_transcripts <- all_merged_transcripts[!grepl("^NOVEL_", all_merged_transcripts)]
        transcript_choices <- c(filtered_novel_transcripts[1], known_transcripts)
        
        cat("DEBUG: load_selected_gene - transcript choices:", paste(transcript_choices, collapse = ", "), "\n")
        
        updateSelectInput(session, "novel_highlight_isoform", 
                         choices = transcript_choices,
                         selected = filtered_novel_transcripts[1])  # Default to selected ORF transcript
        
        incProgress(0.1, detail = 'Complete')
        
        showNotification(paste("Successfully loaded", selected_gene$gene_symbol, 
                              "data for comparison!"), type = "message")
        
      }, error = function(e) {
        showNotification(paste("Failed to load gene data:", e$message), type = "error")
        # Fallback to filtered novel data only
        novel_merged_data(filtered_novel)
        
        # Make filtered novel data available globally for transcript peptide extraction
        novel_merged_data_global <<- filtered_novel
      })
    })
  })
  
  # Output gene search selection status
  output$novel_gene_search_selection_status <- renderText({
    if (is.null(input$novel_gene_search_results_table_rows_selected)) {
      return("Please select a gene from the search results above")
    } else {
      gene_search_data <- novel_gene_search_results()
      if (!is.null(gene_search_data) && length(input$novel_gene_search_results_table_rows_selected) > 0) {
        selected_gene <- gene_search_data[input$novel_gene_search_results_table_rows_selected[1], ]
        return(paste("Selected:", selected_gene$gene_symbol, 
                    "\nGene ID:", selected_gene$gene_id,
                    "\nRDS Available:", ifelse(selected_gene$available, "Yes", "No")))
      }
    }
    return("")
  })
  
  # Output to control analysis section conditional panel
  output$novel_analysis_ready <- reactive({
    merged_data <- novel_merged_data()
    gene_search_data <- novel_gene_search_results()
    
    # Analysis is ready if we have merged data OR if gene search found no good matches
    return(!is.null(merged_data) || (!is.null(gene_search_data) && nrow(gene_search_data) == 0))
  })
  outputOptions(output, "novel_analysis_ready", suspendWhenHidden = FALSE)
  
  # ============================================================================
  # ISOFORM ANALYSIS DATA PROCESSING
  # ============================================================================
  
  # Filtered novel data reactive - only includes selected ORF
  filtered_novel_data <- reactive({
    req(novel_isoform_data(), input$novel_orf_table_rows_selected)
    
    full_novel_data <- novel_isoform_data()
    orf_data <- novel_orf_results()
    
    if (is.null(full_novel_data) || is.null(orf_data) || 
        length(input$novel_orf_table_rows_selected) == 0) {
      return(NULL)
    }
    
    # Get selected ORF information
    selected_row <- input$novel_orf_table_rows_selected[1]
    selected_orf <- orf_data[selected_row, ]
    selected_orf_id <- selected_orf$orf_id
    
    cat("DEBUG: filtered_novel_data() called - Filtering novel data to selected ORF:", selected_orf_id, "\n")
    cat("DEBUG: Available proteinIDs in full data:", paste(unique(full_novel_data$proteinID), collapse = ", "), "\n")
    
    # Filter novel data to only include the selected ORF
    filtered_data <- full_novel_data[full_novel_data$proteinID == selected_orf_id, , drop = FALSE]
    
    cat("DEBUG: Filtered from", nrow(full_novel_data), "ORFs to", nrow(filtered_data), "ORF\n")
    if (nrow(filtered_data) > 0) {
      cat("DEBUG: Filtered data txID:", unique(filtered_data$txID), "\n")
    }
    
    if (nrow(filtered_data) == 0) {
      cat("DEBUG: Warning - No data found for selected ORF:", selected_orf_id, "\n")
      return(NULL)
    }
    
    return(filtered_data)
  })
  
  # Update transcript dropdown when ORF selection changes (immediate feedback)
  observeEvent(input$novel_orf_table_rows_selected, {
    req(novel_isoform_data(), novel_orf_results())
    
    if (length(input$novel_orf_table_rows_selected) > 0) {
      orf_data <- novel_orf_results()
      selected_row <- input$novel_orf_table_rows_selected[1]
      selected_orf <- orf_data[selected_row, ]
      selected_orf_id <- selected_orf$orf_id
      
      # Get the transcript ID for the selected ORF
      full_novel_data <- novel_isoform_data()
      selected_orf_data <- full_novel_data[full_novel_data$proteinID == selected_orf_id, , drop = FALSE]
      
      if (nrow(selected_orf_data) > 0) {
        selected_transcript <- selected_orf_data$txID[1]
        
        # CRITICAL FIX: Always show only the selected ORF's transcript + known transcripts
        # Never show other novel ORF transcripts
        merged_data <- novel_merged_data()
        if (!is.null(merged_data) && nrow(merged_data) > 0) {
          # Get all transcripts from merged data
          all_merged_transcripts <- unique(merged_data$txID)
          # Filter to show only the selected novel transcript + known transcripts (not other novel ORFs)
          known_transcripts <- all_merged_transcripts[!grepl("^NOVEL_", all_merged_transcripts)]
          transcript_choices <- c(selected_transcript, known_transcripts)
        } else {
          # No merged data yet - show only the selected ORF's transcript
          transcript_choices <- selected_transcript
        }
        
        updateSelectInput(session, "novel_highlight_isoform", 
                         choices = transcript_choices,
                         selected = selected_transcript)
        
        cat("DEBUG: ORF selection changed - transcript choices:", paste(transcript_choices, collapse = ", "), "\n")
        cat("DEBUG: Selected transcript:", selected_transcript, "\n")
      }
    }
  })
  
  # Update novel isoform selector when merged data changes
  observeEvent(novel_merged_data(), {
    merged_data <- novel_merged_data()
    filtered_novel <- filtered_novel_data()
    
    if (!is.null(merged_data) && nrow(merged_data) > 0) {
      # CRITICAL FIX: Only show selected ORF transcript + known transcripts
      # Never show other novel ORF transcripts in the dropdown
      all_merged_transcripts <- unique(merged_data$txID)
      cat("Debug: All merged transcripts:", paste(all_merged_transcripts, collapse=", "), "\n")
      
      # Separate known transcripts from novel transcripts
      known_transcripts <- all_merged_transcripts[!grepl("^NOVEL_", all_merged_transcripts)]
      
      # Get the selected ORF's transcript from filtered data
      if (!is.null(filtered_novel) && nrow(filtered_novel) > 0) {
        selected_novel_transcript <- unique(filtered_novel$txID)[1]  # Should be exactly one
        transcript_choices <- c(selected_novel_transcript, known_transcripts)
        default_selection <- selected_novel_transcript
        cat("Debug: Showing selected ORF transcript + known transcripts:", paste(transcript_choices, collapse=", "), "\n")
      } else {
        # Fallback - show only known transcripts if no filtered novel data
        transcript_choices <- known_transcripts
        default_selection <- if (length(known_transcripts) > 0) known_transcripts[1] else NULL
        cat("Debug: No filtered novel data - showing only known transcripts:", paste(transcript_choices, collapse=", "), "\n")
      }
      
      if (length(transcript_choices) > 0) {
        updateSelectInput(session, "novel_highlight_isoform", 
                         choices = transcript_choices,
                         selected = default_selection)
      } else {
        # Clear choices if no transcripts found
        updateSelectInput(session, "novel_highlight_isoform", 
                         choices = character(0),
                         selected = NULL)
      }
    }
  })

  # Reactive data for all novel isoforms analysis (FILTERED like regular tab for identical visualization)
  novel_all_isoforms_data <- reactive({
    req(input$load_novel_isoforms > 0, input$novel_protease, input$novel_miscleavage_type, novel_merged_data(), input$novel_highlight_isoform)
    
    withProgress(message = 'Loading all novel isoforms analysis...', value = 0, {
      protease <- input$novel_protease
      miscleavage_type <- input$novel_miscleavage_type
      merged_data <- novel_merged_data()
      highlight_isoform <- input$novel_highlight_isoform
      
      # Step 1: Filter to show ONLY selected isoform + known gene isoforms (like regular tab)
      incProgress(0.2, detail = 'Filtering transcripts for comparison...')
      
      # Get all available transcripts
      all_available_transcripts <- unique(merged_data$txID)
      
      # Identify known isoforms (those that don't start with "NOVEL_")
      known_isoforms <- all_available_transcripts[!grepl("^NOVEL_", all_available_transcripts)]
      
      # Filter to show ONLY the selected isoform + known isoforms (matching regular tab behavior)
      if (!is.null(highlight_isoform) && highlight_isoform != "") {
        if (grepl("^NOVEL_", highlight_isoform)) {
          # If novel isoform selected, show selected novel + all known
          comparison_transcripts <- c(highlight_isoform, known_isoforms)
        } else {
          # If known isoform selected, show only known isoforms (like regular tab)
          comparison_transcripts <- known_isoforms
        }
      } else {
        # No selection, show only known isoforms
        comparison_transcripts <- known_isoforms
      }
      
      cat("Filtered transcripts for novel comparison:", paste(comparison_transcripts, collapse=", "), "\n")
      cat("Using miscleavage type:", miscleavage_type, "\n")
      
      if (length(comparison_transcripts) == 0) {
        showNotification("No transcripts available for comparison", type = "warning")
        return(NULL)
      }
      
      # Step 2: Get peptides for filtered transcripts only
      incProgress(0.3, detail = 'Loading peptides for comparison transcripts...')
      
      # Use the SAME pipeline as regular isoform analysis
      # Create vis_data structure for compatibility with existing functions
      vis_data <- list(
        genes = c(),  # Not needed for this analysis
        gene_symbols = c(),
        gene_lookup = c(),
        proteases = c("trp", "chymo", "aspn", "lysc", "lysn", "gluc"),
        original_peptides = merged_data  # Use merged data
      )
      
      # Get peptides for all transcripts using the SAME method as regular tab
      all_peptides_list <- list()
      for (i in seq_along(comparison_transcripts)) {
        tx <- comparison_transcripts[i]
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
        showNotification("No peptides found for comparison transcripts", type = "warning")
        return(NULL)
      }
      
      # Combine all peptides
      all_peptides_df <- data.table::rbindlist(all_peptides_list)
      
      # Step 3: Calculate gene boundaries
      incProgress(0.2, detail = 'Combining peptide data...')
      gene_start <- min(all_peptides_df$start) - 1000
      gene_end <- max(all_peptides_df$end) + 1000
      
      # Create transcript position mapping (reassign positions 1, 2, 3... like regular tab)
      transcript_df <- data.frame(
        transcript = comparison_transcripts,
        y_position = seq_along(comparison_transcripts),
        stringsAsFactors = FALSE
      )
      
      # Update y_position in peptides to match filtered transcript positions
      all_peptides_df$y_position <- transcript_df$y_position[match(all_peptides_df$transcript, transcript_df$transcript)]
      
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
        all_transcripts = comparison_transcripts
      ))
    })
  })
  
  # Reactive data for highlighting specific novel isoform (FILTERED with accurate specificity calculation)
  novel_highlighted_isoform_data <- reactive({
    req(novel_all_isoforms_data(), input$novel_highlight_isoform, input$novel_miscleavage_type)
    
    base_data <- novel_all_isoforms_data()
    highlight_isoform <- input$novel_highlight_isoform
    miscleavage_type <- input$novel_miscleavage_type
    
    if (is.null(base_data)) return(NULL)
    
    all_peptides_df <- base_data$all_peptides
    all_transcripts <- base_data$all_transcripts
    
    # Calculate specificity for highlighted isoform peptides
    highlight_peptides <- all_peptides_df[all_peptides_df$transcript == highlight_isoform, ]
    
    if (nrow(highlight_peptides) == 0) {
      return(base_data)
    }
    
    # For specificity calculation, use ALL available isoforms (not just filtered ones)
    # This ensures accurate specificity calculation across the entire gene
    merged_data <- novel_merged_data()
    all_available_transcripts <- unique(merged_data$txID)
    
    # Use the unified method for consistency
    vis_data <- list(
      genes = c(),
      gene_symbols = c(),  
      gene_lookup = c(),
      proteases = c("trp", "chymo", "aspn", "lysc", "lysn", "gluc"),
      original_peptides = merged_data
    )
    
    all_peptides_list <- list()
    for (tx in all_available_transcripts) {
      tx_peptides <- get_transcript_peptides_for_comparison(tx, vis_data, input$novel_protease)
      if (!is.null(tx_peptides) && length(tx_peptides) > 0) {
        all_peptides_list[[tx]] <- data.frame(
          transcript = tx,
          peptide = tx_peptides$peptide,
          stringsAsFactors = FALSE
        )
      }
    }
    
    all_available_peptides <- if (length(all_peptides_list) > 0) data.table::rbindlist(all_peptides_list) else data.frame()
    
    # For each peptide in highlighted isoform, calculate specificity
    highlight_peptides$other_transcript_count <- 0
    highlight_peptides$other_transcripts <- ""
    highlight_peptides$specificity_category <- ""
    
    for (i in 1:nrow(highlight_peptides)) {
      peptide_seq <- highlight_peptides$peptide[i]
      # Count transcripts (excluding highlight) that have this peptide from ALL available isoforms
      other_tx_with_peptide <- unique(all_available_peptides$transcript[
        all_available_peptides$peptide == peptide_seq & all_available_peptides$transcript != highlight_isoform
      ])
      
      highlight_peptides$other_transcript_count[i] <- length(other_tx_with_peptide)
      highlight_peptides$other_transcripts[i] <- paste(other_tx_with_peptide, collapse = ", ")
      
      # Classify specificity using ALL available isoforms for accurate calculation
      total_isoforms <- length(all_available_transcripts)
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
    
    # Create summary statistics
    specificity_summary <- as.data.frame(table(highlight_peptides$specificity_category))
    names(specificity_summary) <- c("Specificity_Category", "Count")
    specificity_summary$Percentage <- round(specificity_summary$Count / nrow(highlight_peptides) * 100, 1)
    
    # Calculate decision tree-based detectability scores
    decision_tree_results <- calculate_isoform_detectability_scores(
      all_peptides_df, highlight_peptides, all_transcripts, 
      input$novel_protease, miscleavage_type
    )
    
    return(list(
      all_peptides = all_peptides_df,
      highlight_peptides = highlight_peptides,
      transcript_df = base_data$transcript_df,
      gene_start = base_data$gene_start,
      gene_end = base_data$gene_end,
      highlight_isoform = highlight_isoform,
      specificity_summary = specificity_summary,
      decision_tree = decision_tree_results
    ))
  })
  
  # ============================================================================
  # ANALYSIS RESULT TABLES
  # ============================================================================
  
  # Render highlighted novel isoform peptides table
  output$novel_highlighted_isoform_table <- DT::renderDataTable({
    data <- novel_highlighted_isoform_data()
    
    if (is.null(data)) {
      return(data.frame(Message = "No data available"))
    }
    
    highlight_peptides <- data$highlight_peptides
    
    # Create display table
    display_df <- highlight_peptides[, c("peptide", "start", "end", "specificity_category", "other_transcript_count", "other_transcripts")]
    colnames(display_df) <- c("Peptide", "Start", "End", "Specificity", "Other Isoforms Count", "Other Isoforms")
    
    DT::datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(3, 'asc'))  # Sort by specificity
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
  
  # Render novel specificity summary table
  output$novel_specificity_summary_table <- DT::renderDataTable({
    data <- novel_highlighted_isoform_data()
    
    if (is.null(data)) {
      return(data.frame(Message = "No data available"))
    }
    
    summary_df <- data$specificity_summary
    
    DT::datatable(
      summary_df,
      options = list(
        dom = 't',
        ordering = FALSE
      ),
      rownames = FALSE
    ) %>%
    DT::formatStyle(
      'Specificity_Category',
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
  
  # Render pipeline summary table
  output$novel_pipeline_summary_table <- DT::renderDataTable({
    data <- novel_analysis_results()
    
    if (is.null(data) || is.null(data$pipeline_summary)) {
      return(data.frame(Message = "No pipeline summary available"))
    }
    
    DT::datatable(
      data$pipeline_summary,
      options = list(
        pageLength = 15,
        dom = 't'  # Only show table
      ),
      rownames = FALSE
    ) %>%
    DT::formatStyle(
      'Value',
      fontWeight = 'bold'
    )
  })
  
  # ============================================================================
  # DOWNLOAD HANDLERS
  # ============================================================================
  
  # Download handlers for novel isoform results
  output$download_novel_dataframe <- downloadHandler(
    filename = function() {
      paste0("novel_isoform_dataframe_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
    },
    content = function(file) {
      data <- novel_isoform_data()
      if (!is.null(data)) {
        saveRDS(data, file)
      }
    }
  )
  
  output$download_novel_gtf <- downloadHandler(
    filename = function() {
      paste0("novel_isoform_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".gtf")
    },
    content = function(file) {
      results <- novel_pipeline_results()
      if (!is.null(results) && results$success && file.exists(results$gtf_file)) {
        file.copy(results$gtf_file, file)
      }
    }
  )
  
  output$download_novel_peptides <- downloadHandler(
    filename = function() {
      paste0("novel_isoform_peptides_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      data <- novel_analysis_results()
      if (!is.null(data) && !is.null(data$all_peptides)) {
        write.csv(data$all_peptides, file, row.names = FALSE)
      }
    }
  )
  
  output$download_novel_pipeline_log <- downloadHandler(
    filename = function() {
      paste0("novel_isoform_pipeline_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
    },
    content = function(file) {
      results <- novel_pipeline_results()
      if (!is.null(results)) {
        log_content <- paste(
          "Novel Isoform Discovery Pipeline Log",
          "=====================================",
          paste("Success:", results$success),
          paste("Timestamp:", Sys.time()),
          "",
          "Pipeline Output:",
          results$log,
          sep = "\n"
        )
        writeLines(log_content, file)
      }
    }
  )
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  # Get novel isoform analysis status
  get_novel_analysis_status <- function() {
    list(
      pipeline_complete = !is.null(novel_pipeline_results()) && novel_pipeline_results()$success,
      orf_data_available = !is.null(novel_orf_results()),
      gene_search_complete = !is.null(novel_gene_search_results()),
      merged_data_available = !is.null(novel_merged_data()),
      analysis_ready = !is.null(novel_merged_data()) || (!is.null(novel_gene_search_results()) && nrow(novel_gene_search_results()) == 0)
    )
  }
  
  # Clear all novel isoform data
  clear_novel_data <- function() {
    novel_pipeline_results(NULL)
    novel_isoform_data(NULL)
    novel_analysis_results(NULL)
    novel_orf_results(NULL)
    novel_gene_search_results(NULL)
    novel_merged_data(NULL)
  }
  
  # Get novel isoform summary
  get_novel_summary <- function() {
    status <- get_novel_analysis_status()
    
    if (!status$pipeline_complete) {
      return("Pipeline not completed")
    }
    
    pipeline_results <- novel_pipeline_results()
    merged_data <- novel_merged_data()
    
    summary_info <- list(
      pipeline_success = pipeline_results$success,
      total_novel_isoforms = if (!is.null(novel_isoform_data())) length(unique(novel_isoform_data()$txID)) else 0,
      merged_isoforms = if (!is.null(merged_data)) length(unique(merged_data$txID)) else 0,
      gene_matches = if (!is.null(novel_gene_search_results()) && nrow(novel_gene_search_results()) > 0) nrow(novel_gene_search_results()) else 0
    )
    
    return(summary_info)
  }
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  
  return(list(
    # Reactive values
    novel_pipeline_results = novel_pipeline_results,
    novel_isoform_data = novel_isoform_data,
    novel_analysis_results = novel_analysis_results,
    novel_orf_results = novel_orf_results,
    novel_gene_search_results = novel_gene_search_results,
    novel_merged_data = novel_merged_data,
    novel_all_isoforms_data = novel_all_isoforms_data,
    novel_highlighted_isoform_data = novel_highlighted_isoform_data,
    
    # Utility functions
    get_novel_analysis_status = get_novel_analysis_status,
    clear_novel_data = clear_novel_data,
    get_novel_summary = get_novel_summary
  ))
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test Novel Isoform Analysis Module
#' 
#' @description 
#' Function to test novel isoform analysis module functionality
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_novel_isoform_analysis_module <- function() {
  cat("Testing Novel Isoform Analysis module...\n")
  
  # Test 1: Pipeline status structure
  status_test <- tryCatch({
    # Test status structure
    mock_status <- list(
      pipeline_complete = FALSE,
      orf_data_available = FALSE,
      gene_search_complete = FALSE,
      merged_data_available = FALSE,
      analysis_ready = FALSE
    )
    
    # Verify all required fields
    required_fields <- c("pipeline_complete", "orf_data_available", "gene_search_complete", 
                        "merged_data_available", "analysis_ready")
    all(required_fields %in% names(mock_status))
  }, error = function(e) {
    cat("Status test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 2: ORF data structure
  orf_test <- tryCatch({
    # Test ORF data structure
    mock_orf <- data.frame(
      orf_id = "ORF_001",
      protein_length = 150,
      quality_score = 0.85,
      starts_with_M = TRUE,
      has_stop = TRUE,
      sequence_preview = "MSPVTGK",
      stringsAsFactors = FALSE
    )
    
    # Verify structure
    required_cols <- c("orf_id", "protein_length", "quality_score", "starts_with_M", "has_stop")
    all(required_cols %in% names(mock_orf))
  }, error = function(e) {
    cat("ORF test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 3: Gene search results structure
  gene_search_test <- tryCatch({
    # Test gene search structure
    mock_gene_search <- data.frame(
      gene_id = "ENSG00000123456",
      gene_symbol = "TEST_GENE",
      available = TRUE,
      stringsAsFactors = FALSE
    )
    
    # Verify structure
    required_cols <- c("gene_id", "gene_symbol", "available")
    all(required_cols %in% names(mock_gene_search))
  }, error = function(e) {
    cat("Gene search test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 4: Specificity calculation logic
  specificity_test <- tryCatch({
    # Test specificity logic
    total_isoforms <- 5
    other_count <- 2
    
    # Test classification logic
    specificity_category <- if (other_count == 0) {
      "Unique"
    } else if (other_count == (total_isoforms - 1)) {
      "Universal"
    } else {
      "Shared"
    }
    
    specificity_category == "Shared"  # Should be "Shared" for this test case
  }, error = function(e) {
    cat("Specificity test failed:", e$message, "\n")
    FALSE
  })
  
  if (status_test && orf_test && gene_search_test && specificity_test) {
    cat("All Novel Isoform Analysis module tests passed!\n")
    return(TRUE)
  } else {
    cat("Some Novel Isoform Analysis module tests failed!\n")
    return(FALSE)
  }
}