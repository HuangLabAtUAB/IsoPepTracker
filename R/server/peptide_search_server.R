  #===============================================================================
  # PEPTIDE SEARCH TAB
  #===============================================================================
  
  # Helper function to find GTF cache file with version-agnostic matching
  find_gtf_cache_file <- function(gene_id) {
    if (is.null(gene_id) || gene_id == "") {
      return(NULL)
    }
    
    cat("DEBUG: Looking for GTF cache file for gene ID:", gene_id, "\n")
    
    # Handle both versioned and non-versioned gene IDs
    base_gene_id <- gsub("\\.\\d+$", "", gene_id)  # Remove version if present
    
    # Try exact match first (most efficient)
    exact_file <- file.path("data/gtf_cache", paste0(gene_id, ".rds"))
    if (file.exists(exact_file)) {
      cat("DEBUG: Found exact GTF cache match:", exact_file, "\n")
      return(exact_file)
    }
    
    # Try version-agnostic pattern matching
    pattern <- paste0("^", base_gene_id, "(\\.\\d+)?\\.rds$")
    cache_files <- list.files("data/gtf_cache", pattern = pattern, full.names = TRUE)
    
    if (length(cache_files) > 0) {
      cat("DEBUG: Found GTF cache match via pattern:", cache_files[1], "\n")
      return(cache_files[1])
    }
    
    cat("WARNING: No GTF cache file found for gene ID:", gene_id, "\n")
    return(NULL)
  }
  
  # Reactive values for peptide search results
  peptide_search_results <- reactiveVal(NULL)
  
  # Store peptide for auto-selection after BLASTP navigation
  blast_peptide_for_selection <- reactiveVal(NULL)
  
  # Check if search databases are available
  search_databases_available <- reactive({
    search_dir <- "data/search"
    if (!dir.exists(search_dir)) return(FALSE)
    
    # Check if at least some search database files exist - simplified approach
    files <- list.files(search_dir, pattern = "peptide_search_", full.names = TRUE)
    rds_files <- files[grepl("\\.rds$", files)]
    return(length(rds_files) > 0)
  })
  
  # Run peptide search when button is clicked
  observeEvent(input$run_peptide_search, {
    req(input$peptide_search_query)
    
    # Clear any existing results to force fresh search with updated parameters
    peptide_search_results(NULL)
    
    # Validate BLAST database
    blast_db_check <- validate_blast_database()
    if (!blast_db_check$valid) {
      showNotification(
        paste("BLAST database not found:", blast_db_check$message),
        type = "error",
        duration = 10
      )
      return()
    }
    
    # Validate peptide query
    peptide_query <- trimws(input$peptide_search_query)
    if (peptide_query == "") {
      showNotification("Please enter a peptide sequence to search.", type = "warning")
      return()
    }
    
    withProgress(message = 'Searching peptides with BLASTP...', {
      tryCatch({
        # Get BLAST parameters from input (higher default E-value for short peptide sensitivity)
        evalue_threshold <- if (!is.null(input$peptide_search_evalue)) input$peptide_search_evalue else 20000
        identity_threshold <- if (!is.null(input$peptide_search_identity)) input$peptide_search_identity else 70
        max_targets <- if (!is.null(input$peptide_search_max_targets)) input$peptide_search_max_targets else 500
        
        # Progress callback function
        progress_callback <- function(message, value) {
          incProgress(value * 0.8, detail = message)
        }
        
        # Perform BLASTP search
        blast_results <- run_blastp_peptide_search(
          peptide_query = peptide_query,
          evalue = evalue_threshold,
          max_target_seqs = max_targets,
          identity_threshold = identity_threshold,
          progress_callback = progress_callback
        )
        
        incProgress(0.1, detail = "Processing BLAST results...")
        
        if (blast_results$success && nrow(blast_results$results) > 0) {
          # Convert BLAST results to format compatible with existing UI
          results <- blast_results$results
          
          # Cross-reference with peptide databases
          incProgress(0.05, detail = "Cross-referencing with peptide databases...")
          enhanced_results <- cross_reference_blast_with_peptide_databases(results, peptide_query)
          
          # Debug: Check what columns we have before storing results
          cat("DEBUG: Columns in enhanced_results before storing:", paste(names(enhanced_results), collapse = ", "), "\n")
          cat("DEBUG: Query coverage present before storing:", "query_coverage" %in% names(enhanced_results), "\n")
          
          # Add default values for compatibility with navigation
          enhanced_results$protease_used <- "blastp"
          enhanced_results$miscleavage_type_used <- "none"
          enhanced_results$peptide <- peptide_query  # Use the searched peptide sequence
          enhanced_results$txID <- enhanced_results$transcript_id
          enhanced_results$geneID <- enhanced_results$gene_id
          enhanced_results$geneSymbol <- enhanced_results$gene_symbol
          
          # Store enhanced results
          peptide_search_results(enhanced_results)
          
          # Debug: Check what was actually stored
          stored_results <- peptide_search_results()
          cat("DEBUG: Columns in stored results:", paste(names(stored_results), collapse = ", "), "\n")
          cat("DEBUG: Query coverage present in stored results:", "query_coverage" %in% names(stored_results), "\n")
          
          incProgress(0.1, detail = "Complete")
          
          # Show notification with results count
          showNotification(
            paste("BLASTP found", nrow(enhanced_results), "matches across", 
                  length(unique(enhanced_results$gene_id)), "genes and", 
                  length(unique(enhanced_results$transcript_id)), "transcript isoforms.",
                  "Best hit:", round(max(enhanced_results$identity_percent), 1), "% identity"),
            type = "message",
            duration = 8
          )
        } else {
          peptide_search_results(data.frame())
          
          if (!blast_results$success) {
            # Actual system error
            error_msg <- blast_results$error
            showNotification(paste("BLASTP search failed:", error_msg), type = "error")
          } else {
            # No matches found (successful search with 0 results)
            no_match_msg <- if (!is.null(blast_results$message)) {
              blast_results$message
            } else {
              paste("No matches found for peptide", blast_results$query, "in the protein database.")
            }
            
            # Add helpful suggestions
            suggestions <- "Try: • Using a longer peptide sequence • Increasing E-value threshold • Checking peptide sequence for typos"
            full_msg <- paste(no_match_msg, suggestions, sep = "\n\n")
            
            showNotification(
              full_msg,
              type = "message",
              duration = 10
            )
          }
        }
        
      }, error = function(e) {
        showNotification(paste("BLASTP search failed:", e$message), type = "error")
        peptide_search_results(NULL)
      })
    })
  })
  
  # Output for indicating if search results are available
  output$peptide_search_results_available <- reactive({
    results <- peptide_search_results()
    return(!is.null(results) && nrow(results) > 0)
  })
  outputOptions(output, "peptide_search_results_available", suspendWhenHidden = FALSE)
  
  # Render search results summary
  output$peptide_search_summary <- renderText({
    results <- peptide_search_results()
    if (is.null(results) || nrow(results) == 0) {
      return("No results")
    }
    
    # Enhanced summary with BLAST statistics
    if ("identity_percent" %in% names(results)) {
      best_identity <- max(results$identity_percent, na.rm = TRUE)
      best_evalue <- min(results$evalue, na.rm = TRUE)
      paste0("Found ", nrow(results), " BLASTP hits in ", 
             length(unique(results$geneID)), " genes and ", 
             length(unique(results$txID)), " transcript isoforms. ",
             "Best hit: ", round(best_identity, 1), "% identity (E=", 
             format(best_evalue, scientific = TRUE, digits = 2), ")")
    } else {
      paste0("Found ", nrow(results), " peptide matches in ", 
             length(unique(results$geneID)), " genes and ", 
             length(unique(results$txID)), " transcript isoforms")
    }
  })
  
  # Render search results table
  output$peptide_search_results_table <- DT::renderDataTable({
    results <- peptide_search_results()
    
    if (is.null(results) || nrow(results) == 0) {
      return(data.frame(Message = "No search results available"))
    }
    
    # Create display table with BLAST statistics if available
    if ("identity_percent" %in% names(results)) {
      # BLAST results format with enzyme availability
      base_cols <- c("peptide", "geneID", "geneSymbol", "txID", "identity_percent", "evalue", "bit_score")
      base_colnames <- c("Query Peptide", "Gene ID", "Gene Symbol", "Transcript ID", 
                        "Identity %", "E-value", "Bit Score")
      
      # Add exact match column if available
      if ("exact_match_found" %in% names(results)) {
        base_cols <- c(base_cols, "exact_match_found")
        base_colnames <- c(base_colnames, "Exact Match")
      }
      
      # Add all 12 enzyme availability columns if they exist
      enzyme_cols <- c("trp_no_miss", "trp_upto2miss", "chymo_no_miss", "chymo_upto2miss", 
                       "aspn_no_miss", "aspn_upto2miss", "lysc_no_miss", "lysc_upto2miss",
                       "lysn_no_miss", "lysn_upto2miss", "gluc_no_miss", "gluc_upto2miss")
      available_enzyme_cols <- enzyme_cols[enzyme_cols %in% names(results)]
      
      display_cols <- c(base_cols, available_enzyme_cols)
      display_df <- results[, display_cols]
      
      # Add enzyme column names
      enzyme_colnames <- c()
      if ("trp_no_miss" %in% available_enzyme_cols) enzyme_colnames <- c(enzyme_colnames, "Trypsin (No Miss)")
      if ("trp_upto2miss" %in% available_enzyme_cols) enzyme_colnames <- c(enzyme_colnames, "Trypsin (2 Miss)")
      if ("chymo_no_miss" %in% available_enzyme_cols) enzyme_colnames <- c(enzyme_colnames, "Chymo (No Miss)")
      if ("chymo_upto2miss" %in% available_enzyme_cols) enzyme_colnames <- c(enzyme_colnames, "Chymo (2 Miss)")
      if ("aspn_no_miss" %in% available_enzyme_cols) enzyme_colnames <- c(enzyme_colnames, "AspN (No Miss)")
      if ("aspn_upto2miss" %in% available_enzyme_cols) enzyme_colnames <- c(enzyme_colnames, "AspN (2 Miss)")
      if ("lysc_no_miss" %in% available_enzyme_cols) enzyme_colnames <- c(enzyme_colnames, "LysC (No Miss)")
      if ("lysc_upto2miss" %in% available_enzyme_cols) enzyme_colnames <- c(enzyme_colnames, "LysC (2 Miss)")
      if ("lysn_no_miss" %in% available_enzyme_cols) enzyme_colnames <- c(enzyme_colnames, "LysN (No Miss)")
      if ("lysn_upto2miss" %in% available_enzyme_cols) enzyme_colnames <- c(enzyme_colnames, "LysN (2 Miss)")
      if ("gluc_no_miss" %in% available_enzyme_cols) enzyme_colnames <- c(enzyme_colnames, "GluC (No Miss)")
      if ("gluc_upto2miss" %in% available_enzyme_cols) enzyme_colnames <- c(enzyme_colnames, "GluC (2 Miss)")
      
      colnames(display_df) <- c(base_colnames, enzyme_colnames)
      
      # Sort by identity percentage (descending) and E-value (ascending)
      display_df <- display_df[order(-display_df$`Identity %`, display_df$`E-value`), ]
      
    } else {
      # Fallback for non-BLAST results
      display_df <- results[, c("peptide", "geneID", "geneSymbol", "txID", "protease_used", "miscleavage_type_used")]
      colnames(display_df) <- c("Peptide", "Gene ID", "Gene Symbol", "Transcript ID", "Protease", "Miscleavage Type")
    }
    
    dt <- DT::datatable(
      display_df,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        order = list(list(4, 'desc'))  # Sort by Identity % (descending)
      ),
      rownames = FALSE,
      selection = 'single'
    )
    
    # Add formatting based on whether we have BLAST columns
    if ("Identity %" %in% colnames(display_df)) {
      dt <- dt %>%
        DT::formatStyle('Query Peptide', backgroundColor = '#e8f4fd', fontWeight = 'bold') %>%
        DT::formatRound('Identity %', 1) %>%
        DT::formatSignif('E-value', digits = 2) %>%
        DT::formatRound('Bit Score', 1)
      
      # Add formatting for enzyme availability columns
      enzyme_display_cols <- enzyme_colnames
      for (col in enzyme_display_cols) {
        if (col %in% colnames(display_df)) {
          dt <- dt %>%
            DT::formatStyle(col, 
                           backgroundColor = DT::styleEqual("✅ Found", "#d4edda"),
                           color = DT::styleEqual(c("✅ Found", "❌ Not Found"), 
                                                c("#155724", "#721c24")),
                           fontWeight = 'bold')
        }
      }
    } else {
      dt <- dt %>%
        DT::formatStyle('Peptide', backgroundColor = '#e8f4fd', fontWeight = 'bold')
    }
    
    return(dt)
  })
  
  # Download handler for peptide search results
  output$download_peptide_search <- downloadHandler(
    filename = function() {
      query_clean <- gsub("[^A-Za-z0-9]", "_", input$peptide_search_query)
      paste0("peptide_search_", query_clean, "_", input$peptide_search_protease, "_", 
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      results <- peptide_search_results()
      if (!is.null(results)) {
        write.csv(results, file, row.names = FALSE)
      }
    }
  )
  
  # Navigation button from search results - ONLY to isoform analysis
  
  observeEvent(input$goto_isoform_analysis, {
    cat("DEBUG: Navigation button clicked\n")
    
    # Check if row is selected
    if (is.null(input$peptide_search_results_table_rows_selected) || length(input$peptide_search_results_table_rows_selected) == 0) {
      showNotification("Please select a row from the search results table first", type = "warning")
      return()
    }
    
    selected_row <- input$peptide_search_results_table_rows_selected[1]
    results <- peptide_search_results()
    
    cat("DEBUG: Selected row:", selected_row, "\n")
    cat("DEBUG: Results available:", !is.null(results), "\n")
    if (!is.null(results)) cat("DEBUG: Results rows:", nrow(results), "\n")
    
    if (!is.null(results) && selected_row <= nrow(results)) {
      selected_gene <- results$geneID[selected_row]
      selected_transcript <- results$txID[selected_row]
      blast_peptide <- results$peptide[selected_row]
      
      # Find the best enzyme/miscleavage combination where peptide exists
      enzymes <- c("trp", "chymo", "aspn", "lysc", "lysn", "gluc")
      miscleavages <- c("no_miss", "upto2miss")
      
      found_enzyme <- "trp"  # Default fallback
      found_miscleavage <- "no_miss_cleavage"  # Default fallback
      
      # Search for first available combination
      enzyme_found <- FALSE
      for (enzyme in enzymes) {
        for (miscleavage in miscleavages) {
          col_name <- paste0(enzyme, "_", miscleavage)
          if (col_name %in% names(results) && results[selected_row, col_name] == "✅ Found") {
            found_enzyme <- enzyme
            found_miscleavage <- ifelse(miscleavage == "no_miss", "no_miss_cleavage", "upto_two_misscleavage")
            enzyme_found <- TRUE
            break
          }
        }
        if (enzyme_found) break  # Found a match, stop searching
      }
      
      # Log enzyme selection result
      if (enzyme_found) {
        cat("DEBUG: Found peptide in", found_enzyme, found_miscleavage, "\n")
      } else {
        cat("DEBUG: Peptide not found in any enzyme database, using defaults\n")
      }
      
      cat("DEBUG: Attempting navigation with:\n")
      cat("  Gene:", selected_gene, "\n")
      cat("  Transcript:", selected_transcript, "\n")
      cat("  Enzyme:", found_enzyme, "\n")
      cat("  Miscleavage:", found_miscleavage, "\n")
      
      # CORE NAVIGATION: Update inputs in the correct order to trigger data loading
      cat("DEBUG: Updating inputs to trigger proper data loading\n")
      
      # 1. First update enzyme and miscleavage (these need to be set before gene selection)
      updateSelectInput(session, "protease", selected = found_enzyme)
      updateSelectInput(session, "miscleavage_type", selected = found_miscleavage)
      
      # 2. Add the gene to dropdown choices and select it
      gene_symbol <- results$geneSymbol[selected_row]
      gene_choice_label <- paste0(gene_symbol, " (", selected_gene, ")")
      gene_choices <- setNames(selected_gene, gene_choice_label)
      
      updateSelectizeInput(session, "gene", 
                          choices = gene_choices,
                          selected = selected_gene)
      
      cat("DEBUG: Updated inputs - Gene:", selected_gene, "(", gene_symbol, ") Enzyme:", found_enzyme, "Miscleavage:", found_miscleavage, "\n")
      
      # CORE NAVIGATION: Switch to canonical analysis page, then isoform analysis tab
      updateTabItems(session, "tabs", "canonical_analysis")
      updateTabsetPanel(session, "canonical_tabs", "isoform_analysis")
      
      # Delay transcript selection to ensure gene data loads completely
      shinyjs::delay(3000, {
        cat("DEBUG: Attempting transcript selection after gene data should be loaded\n")
        updateSelectInput(session, "highlight_isoform", selected = selected_transcript)
        cat("DEBUG: Transcript selection attempted:", selected_transcript, "\n")
        
        # Also trigger peptide auto-selection after transcript is set
        shinyjs::delay(1000, {
          blast_peptide_for_selection(blast_peptide)
          cat("DEBUG: Peptide auto-selection triggered for:", blast_peptide, "\n")
        })
      })
      
      cat("DEBUG: Navigation commands sent\n")
      
      showNotification(paste("Navigated to", results$geneSymbol[selected_row], "isoform analysis with", 
                           toupper(found_enzyme), "enzyme.", "Look for peptide:", blast_peptide), 
                           type = "message", duration = 10)
    }
  })
  
  #===============================================================================
  # BLAST PERFECT MATCH VISUALIZATION (SAFE ADDITION)
  #===============================================================================
  
  # Source the genomic mapping module
  source("R/blast_genomic_mapper.R", local = TRUE)
  
  # Helper function for empty plotly messages
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
  
  #===============================================================================
  # MULTI-GENE BLAST VISUALIZATION - CLEANED UP
  #===============================================================================
  
  # Create combined isoform-centric BLAST plot with sequential gene faceting
  create_combined_isoform_blast_plot <- function(transcript_data_list, gene_transcript_mapping, 
                                                query_peptide, enzyme, use_compression = TRUE) {
    tryCatch({
      cat("DEBUG: Creating combined plot for", length(transcript_data_list), "transcripts\n")
      
      # Prepare combined data frame for ggplot
      plot_data_list <- list()
      
      # Process each gene sequentially (NO global y_offset - faceting will handle gene separation)
      unique_genes <- names(gene_transcript_mapping)
      
      for (gene_id in unique_genes) {
        gene_info <- gene_transcript_mapping[[gene_id]]
        gene_transcripts <- gene_info$all_transcripts
        
        cat("DEBUG: Processing gene", gene_id, "with", length(gene_transcripts), "transcripts\n")
        
        # Process each transcript in this gene (reset y-position for each gene)
        for (j in seq_along(gene_transcripts)) {
          transcript_id <- gene_transcripts[j]
          transcript_key <- paste0(gene_id, "_", transcript_id)
          
          if (!transcript_key %in% names(transcript_data_list)) {
            cat("WARNING: Transcript data not found for", transcript_key, "\n")
            next
          }
          
          transcript_data <- transcript_data_list[[transcript_key]]
          current_y <- j  # Reset y-position for each gene (1, 2, 3, etc.)
          
          cat("DEBUG: Transcript", transcript_id, "assigned y-position:", current_y, "in gene", gene_id, "\n")
          
          # Add exon/CDS structure data
          if (!is.null(transcript_data$exons) && length(transcript_data$exons) > 0) {
            # Create exon rectangles
            exon_data <- data.frame(
              gene_id = gene_id,
              gene_symbol = transcript_data$gene_symbol,
              transcript_id = transcript_id,
              transcript_label = transcript_id,  # For y-axis labeling
              element_type = "exon",
              start = start(transcript_data$exons),
              end = end(transcript_data$exons),
              y_position = current_y,
              y_min = current_y - 0.15,
              y_max = current_y + 0.15,
              has_blast_match = transcript_data$has_blast_match,
              stringsAsFactors = FALSE
            )
            plot_data_list[[paste0(transcript_key, "_exons")]] <- exon_data
          }
          
          # Add CDS data if available
          if (!is.null(transcript_data$cds) && length(transcript_data$cds) > 0) {
            cds_data <- data.frame(
              gene_id = gene_id,
              gene_symbol = transcript_data$gene_symbol,
              transcript_id = transcript_id,
              transcript_label = transcript_id,
              element_type = "cds",
              start = start(transcript_data$cds),
              end = end(transcript_data$cds),
              y_position = current_y,
              y_min = current_y - 0.1,
              y_max = current_y + 0.1,
              has_blast_match = transcript_data$has_blast_match,
              stringsAsFactors = FALSE
            )
            plot_data_list[[paste0(transcript_key, "_cds")]] <- cds_data
          }
          
          # Add digestible peptides (isoform-centric view)
          if (!is.null(transcript_data$gene_peptides) && !is.null(transcript_data$gene_peptides$peptides)) {
            gene_peptides_data <- transcript_data$gene_peptides$peptides
            
            # Get peptides for this transcript and enzyme
            tx_rows <- which(gene_peptides_data$txID == transcript_id)
            if (length(tx_rows) > 0) {
              # Get enzyme-specific mapped ranges
              enzyme_mapped_ranges_col <- paste0(enzyme, "Peps_mapped_ranges")
              
              if (enzyme_mapped_ranges_col %in% names(gene_peptides_data)) {
                mapped_ranges_list <- gene_peptides_data[[enzyme_mapped_ranges_col]]
                if (!is.null(mapped_ranges_list) && length(mapped_ranges_list) >= tx_rows[1]) {
                  genomic_ranges <- mapped_ranges_list[[tx_rows[1]]]
                  
                  if (!is.null(genomic_ranges) && length(genomic_ranges) > 0) {
                    # Create peptide data
                    peptide_data <- data.frame(
                      gene_id = gene_id,
                      gene_symbol = transcript_data$gene_symbol,
                      transcript_id = transcript_id,
                      transcript_label = transcript_id,
                      element_type = "peptide",
                      start = start(genomic_ranges),
                      end = end(genomic_ranges),
                      y_position = current_y,
                      y_min = current_y - 0.05,
                      y_max = current_y + 0.05,
                      has_blast_match = transcript_data$has_blast_match,
                      stringsAsFactors = FALSE
                    )
                    plot_data_list[[paste0(transcript_key, "_peptides")]] <- peptide_data
                  }
                }
              }
            }
          }
          
          # Add BLAST peptide overlay for matching transcripts
          if (transcript_data$has_blast_match) {
            # Map BLAST peptide to this transcript
            transcript_structure <- list(
              success = TRUE,
              exons = gene_info$exons_by_transcript,
              cds = gene_info$cds_by_transcript
            )
            
            peptide_mapping <- map_blast_peptide_to_transcript(
              blast_peptide = query_peptide,
              transcript_id = transcript_id,
              gene_id = gene_id,
              transcript_structure = transcript_structure
            )
            
            if (!is.null(peptide_mapping) && peptide_mapping$success && 
                length(peptide_mapping$genomic_ranges) > 0) {
              
              blast_ranges <- peptide_mapping$genomic_ranges
              blast_data <- data.frame(
                gene_id = gene_id,
                gene_symbol = transcript_data$gene_symbol,
                transcript_id = transcript_id,
                transcript_label = transcript_id,
                element_type = "blast_overlay",
                start = start(blast_ranges),
                end = end(blast_ranges),
                y_position = current_y,
                y_min = current_y - 0.08,
                y_max = current_y + 0.08,
                has_blast_match = TRUE,
                stringsAsFactors = FALSE
              )
              plot_data_list[[paste0(transcript_key, "_blast")]] <- blast_data
            }
          }
        }
        
        # Gene separation handled by faceting, no need for y_offset
      }
      
      # Combine all plot data
      if (length(plot_data_list) == 0) {
        return(NULL)
      }
      
      combined_data <- do.call(rbind, plot_data_list)
      combined_data$gene_symbol <- factor(combined_data$gene_symbol, levels = unique(combined_data$gene_symbol))
      
      cat("DEBUG: Combined data has", nrow(combined_data), "elements\n")
      cat("DEBUG: Unique genes:", paste(unique(combined_data$gene_symbol), collapse = ", "), "\n")
      cat("DEBUG: Y-position range:", min(combined_data$y_position), "to", max(combined_data$y_position), "\n")
      
      # Create transcript labels for y-axis (using base R to avoid dplyr dependency)
      transcript_labels_data <- unique(combined_data[, c("gene_symbol", "transcript_id", "y_position")])
      transcript_labels_data <- transcript_labels_data[order(transcript_labels_data$gene_symbol, transcript_labels_data$y_position), ]
      
      cat("DEBUG: Created labels for", nrow(transcript_labels_data), "transcripts\n")
      
      # Create the ggplot with proper faceting and transcript labels
      p <- ggplot(combined_data) +
        theme_minimal() +
        theme(
          strip.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 11),
          axis.text.y = element_text(size = 9),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing = unit(0.5, "lines")  # Space between gene facets
        )
      
      # Add different elements with isoform-centric styling (matching your example)
      exon_data <- combined_data[combined_data$element_type == "exon", ]
      if (nrow(exon_data) > 0) {
        p <- p + geom_rect(
          data = exon_data,
          aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max),
          fill = "#E8F4FD", color = "#B0C4DE", alpha = 0.7, size = 0.3
        )
      }
      
      cds_data <- combined_data[combined_data$element_type == "cds", ]
      if (nrow(cds_data) > 0) {
        p <- p + geom_rect(
          data = cds_data,
          aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max),
          fill = "#FFE4B5", color = "#DAA520", alpha = 0.9, size = 0.3
        )
      }
      
      peptide_data <- combined_data[combined_data$element_type == "peptide", ]
      if (nrow(peptide_data) > 0) {
        p <- p + geom_rect(
          data = peptide_data,
          aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max),
          fill = "orange", color = "darkorange", alpha = 0.8, size = 0.2
        )
      }
      
      blast_data <- combined_data[combined_data$element_type == "blast_overlay", ]
      if (nrow(blast_data) > 0) {
        p <- p + geom_rect(
          data = blast_data,
          aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max),
          fill = "#FF6B6B", color = "#D63031", alpha = 0.9, size = 0.6
        )
      }
      
      # Add transcript ID labels on y-axis 
      y_breaks <- transcript_labels_data$y_position
      y_labels <- transcript_labels_data$transcript_id
      
      p <- p + scale_y_continuous(
        breaks = y_breaks,
        labels = y_labels,
        name = "Transcripts"
      )
      
      # Add gene faceting (sequential, one below another)
      p <- p + facet_grid(gene_symbol ~ ., 
                          scales = "free", 
                          space = "free_y",
                          labeller = labeller(gene_symbol = function(x) paste0(x, " Gene")))
      
      # Add labels and title
      p <- p + labs(
        title = paste0("Multi-Gene BLAST Visualization (100% Identity Matches)"),
        subtitle = paste0("Query: ", substr(query_peptide, 1, 50), 
                         ifelse(nchar(query_peptide) > 50, "...", ""), 
                         " | Enzyme: ", toupper(enzyme)),
        x = "Genomic Position (bp)",
        y = "Transcripts"
      )
      
      # Add a simple legend by adding invisible points
      p <- p + 
        geom_point(aes(x = Inf, y = Inf, color = "Exons"), alpha = 0) +
        geom_point(aes(x = Inf, y = Inf, color = "CDS"), alpha = 0) +
        geom_point(aes(x = Inf, y = Inf, color = "Peptides"), alpha = 0) +
        geom_point(aes(x = Inf, y = Inf, color = "BLAST Match"), alpha = 0) +
        scale_color_manual(
          name = "Elements",
          values = c("Exons" = "#B0C4DE", "CDS" = "#DAA520", 
                    "Peptides" = "darkorange", "BLAST Match" = "#D63031"),
          guide = guide_legend(override.aes = list(alpha = 1, size = 3))
        )
      
      # Convert to plotly with improved settings
      plotly_obj <- ggplotly(p, tooltip = c("text")) %>%
        config(
          displayModeBar = TRUE,
          displaylogo = FALSE,
          modeBarButtonsToRemove = c("pan2d", "select2d", "lasso2d", "autoScale2d")
        ) %>%
        layout(
          title = list(
            text = paste0("Multi-Gene BLAST Visualization<br>Query: ", 
                         substr(query_peptide, 1, 60), 
                         ifelse(nchar(query_peptide) > 60, "...", "")),
            font = list(size = 14),
            x = 0.05,
            xanchor = 'left'
          ),
          margin = list(l = 150, r = 50, t = 120, b = 80),  # More space for transcript labels
          showlegend = TRUE,
          legend = list(
            x = 1.02,
            y = 1,
            xanchor = 'left',
            bgcolor = 'rgba(255,255,255,0.8)',
            bordercolor = 'rgba(0,0,0,0.2)',
            borderwidth = 1
          )
        )
      
      return(plotly_obj)
      
    }, error = function(e) {
      cat("ERROR in create_combined_isoform_blast_plot:", e$message, "\n")
      return(NULL)
    })
  }
  
  # Cached GTF visualization using gene-first approach for BLAST transcripts
  create_cached_gtf_visualization <- function(match_info, blast_peptide) {
    tryCatch({
      cat("DEBUG: Creating cached GTF visualization for", match_info$transcript_id, "\n")
      
      withProgress(message = "Loading transcript from cached GTF...", value = 0, {
        incProgress(0.2, detail = "Loading gene details...")
        
        # Step 1: Use gene-first approach - load cached GTF directly with version-agnostic lookup
        cache_file <- find_gtf_cache_file(match_info$gene_id)
        if (is.null(cache_file)) {
          return(empty_plotly_message(paste("Gene", match_info$gene_id, "not found in cached GTF system")))
        }
        
        gene_details <- readRDS(cache_file)
        if (is.null(gene_details)) {
          return(empty_plotly_message(paste("Failed to load gene details for", match_info$gene_id)))
        }
        
        cat("DEBUG: Successfully loaded gene details for", match_info$gene_id, "\n")
        
        incProgress(0.3, detail = "Loading transcript structure...")
        
        # Step 2: Extract transcript structure from cached GTF data
        if (!match_info$transcript_id %in% gene_details$transcript_ids) {
          return(empty_plotly_message(paste("Transcript", match_info$transcript_id, "not found in gene", match_info$gene_id)))
        }
        
        transcript_structure <- list(
          success = TRUE,
          exons = gene_details$exons_by_transcript,
          cds = gene_details$cds_by_transcript,
          gene_id = match_info$gene_id,
          transcript_ids = match_info$transcript_id
        )
        
        cat("DEBUG: Successfully loaded transcript structure for", match_info$transcript_id, "\n")
        
        # Step 3: Extract exons and CDS from cached structure
        exons_by_transcript <- transcript_structure$exons
        cds_by_transcript <- transcript_structure$cds
        
        # Get specific transcript data
        exons <- exons_by_transcript[[match_info$transcript_id]]
        cds <- cds_by_transcript[[match_info$transcript_id]]
        
        if (is.null(exons) || length(exons) == 0) {
          return(empty_plotly_message("No exons found for transcript in cached data"))
        }
        
        cat("DEBUG: Found", length(exons), "exons and", ifelse(is.null(cds), 0, length(cds)), "CDS segments from cached GTF\n")
        
        incProgress(0.2, detail = "Mapping BLAST peptide...")
        
        # Step 4: Map BLAST peptide to genomic coordinates using existing system
        peptide_mapping <- map_blast_peptide_to_transcript(
          blast_peptide = blast_peptide,
          transcript_id = match_info$transcript_id,
          gene_id = match_info$gene_id,
          transcript_structure = transcript_structure
        )
        
        incProgress(0.1, detail = "Creating visualization...")
        
        # Step 5: Create visualization in isoform analysis style with compression
        gene_start <- min(start(exons)) - 1000
        gene_end <- max(end(exons)) + 1000
        
        # Extract chromosome information
        chromosome <- as.character(seqnames(exons)[1])
        cat("DEBUG: Extracted chromosome:", chromosome, "\n")
        
        # Create compression map for better visualization (default enabled)
        use_compression <- TRUE
        compression_map <- NULL
        
        if (use_compression && length(exons) > 0) {
          # Source compression functions if not already loaded
          if (!exists("create_compression_map")) {
            source("R/coordinate_compression.R")
          }
          
          # Create compression map from exons
          exons_by_transcript <- list()
          exons_by_transcript[[match_info$transcript_id]] <- exons
          compression_map <- create_compression_map(exons_by_transcript)
          
          cat("DEBUG: Created compression map for", length(exons), "exons\n")
        }
        
        # Create visualization data matching isoform analysis style
        transcript_y <- 1
        
        # Create plot data structures
        exon_starts <- start(exons)
        exon_ends <- end(exons)
        
        # Apply compression to exon coordinates if enabled
        if (use_compression && !is.null(compression_map)) {
          exon_starts_compressed <- sapply(exon_starts, compression_map$compress)
          exon_ends_compressed <- sapply(exon_ends, compression_map$compress)
        } else {
          exon_starts_compressed <- exon_starts
          exon_ends_compressed <- exon_ends
        }
        
        exon_plot_data <- data.frame(
          start = exon_starts_compressed,
          end = exon_ends_compressed,
          original_start = exon_starts,
          original_end = exon_ends,
          y_min = transcript_y - 0.15,
          y_max = transcript_y + 0.15,
          type = "exon",
          transcript = match_info$transcript_id,
          hover_text = paste0("Exon | ", exon_starts, "-", exon_ends, " | Transcript: ", match_info$transcript_id),
          stringsAsFactors = FALSE
        )
        
        # CDS plot data
        cds_plot_data <- data.frame()
        if (length(cds) > 0) {
          cds_starts <- start(cds)
          cds_ends <- end(cds)
          
          # Apply compression to CDS coordinates if enabled
          if (use_compression && !is.null(compression_map)) {
            cds_starts_compressed <- sapply(cds_starts, compression_map$compress)
            cds_ends_compressed <- sapply(cds_ends, compression_map$compress)
          } else {
            cds_starts_compressed <- cds_starts
            cds_ends_compressed <- cds_ends
          }
          
          cds_plot_data <- data.frame(
            start = cds_starts_compressed,
            end = cds_ends_compressed,
            original_start = cds_starts,
            original_end = cds_ends,
            y_min = transcript_y - 0.1,
            y_max = transcript_y + 0.1,
            type = "CDS",
            transcript = match_info$transcript_id,
            hover_text = paste0("CDS | ", cds_starts, "-", cds_ends, " | Transcript: ", match_info$transcript_id),
            stringsAsFactors = FALSE
          )
        }
        
        # BLAST peptide plot data
        peptide_plot_data <- data.frame()
        if (!is.null(peptide_mapping) && peptide_mapping$success && length(peptide_mapping$genomic_ranges) > 0) {
          blast_ranges <- peptide_mapping$genomic_ranges
          peptide_starts <- start(blast_ranges)
          peptide_ends <- end(blast_ranges)
          
          # Apply compression to peptide coordinates if enabled
          if (use_compression && !is.null(compression_map)) {
            peptide_starts_compressed <- sapply(peptide_starts, compression_map$compress)
            peptide_ends_compressed <- sapply(peptide_ends, compression_map$compress)
          } else {
            peptide_starts_compressed <- peptide_starts
            peptide_ends_compressed <- peptide_ends
          }
          
          peptide_plot_data <- data.frame(
            start = peptide_starts_compressed,
            end = peptide_ends_compressed,
            original_start = peptide_starts,
            original_end = peptide_ends,
            y_min = transcript_y - 0.05,
            y_max = transcript_y + 0.05,
            type = "blast_peptide",
            transcript = match_info$transcript_id,
            hover_text = paste0("BLAST Peptide: ", blast_peptide, " | ", peptide_starts, "-", peptide_ends),
            stringsAsFactors = FALSE
          )
          cat("DEBUG: Mapped BLAST peptide to", nrow(peptide_plot_data), "genomic segments\n")
        } else {
          cat("DEBUG: BLAST peptide mapping failed, using approximate location\n")
          # Fallback: place peptide in middle CDS
          if (length(cds) > 0) {
            mid_cds <- cds[ceiling(length(cds)/2)]
            fallback_start <- start(mid_cds)
            fallback_end <- min(end(mid_cds), start(mid_cds) + 27)  # 9 AA * 3 bp ≈ 27 bp
            
            # Apply compression to fallback coordinates if enabled
            if (use_compression && !is.null(compression_map)) {
              fallback_start_compressed <- compression_map$compress(fallback_start)
              fallback_end_compressed <- compression_map$compress(fallback_end)
            } else {
              fallback_start_compressed <- fallback_start
              fallback_end_compressed <- fallback_end
            }
            
            peptide_plot_data <- data.frame(
              start = fallback_start_compressed,
              end = fallback_end_compressed,
              original_start = fallback_start,
              original_end = fallback_end,
              y_min = transcript_y - 0.05,
              y_max = transcript_y + 0.05,
              type = "blast_peptide_approx",
              transcript = match_info$transcript_id,
              hover_text = paste0("BLAST Peptide (approx): ", blast_peptide),
              stringsAsFactors = FALSE
            )
          }
        }
        
        # Step 6: Create the plot with compression-aware coordinates
        if (use_compression && !is.null(compression_map)) {
          # Use compressed coordinates for plot limits
          compressed_gene_start <- compression_map$compress(gene_start)
          compressed_gene_end <- compression_map$compress(gene_end)
          plot_start <- compressed_gene_start - (compressed_gene_end - compressed_gene_start) * 0.15
          
          plot_xlim_start <- plot_start
          plot_xlim_end <- compressed_gene_end
          backbone_start <- compressed_gene_start
          backbone_end <- compressed_gene_end
        } else {
          # Use original coordinates for plot limits
          plot_start <- gene_start - (gene_end - gene_start) * 0.15
          plot_xlim_start <- plot_start
          plot_xlim_end <- gene_end
          backbone_start <- gene_start
          backbone_end <- gene_end
        }
        
        p <- ggplot() +
          theme_minimal() +
          xlim(plot_xlim_start, plot_xlim_end) +
          ylim(0.5, 1.5)
        
        # Add transcript backbone line (only for the gene region, not the extended area)
        p <- p + geom_segment(
          aes(x = backbone_start, xend = backbone_end, y = transcript_y, yend = transcript_y),
          color = "gray50", size = 0.8
        )
        
        # Add exons (light colored rectangles like in isoform analysis)
        p <- p + geom_rect(
          data = exon_plot_data,
          aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max, text = hover_text),
          fill = "#E8F4FD", color = "#4A90E2", alpha = 0.8, size = 0.3
        )
        
        # Add CDS regions (darker rectangles like in isoform analysis)
        if (nrow(cds_plot_data) > 0) {
          p <- p + geom_rect(
            data = cds_plot_data,
            aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max, text = hover_text),
            fill = "#FFE4B5", color = "#DAA520", alpha = 0.9, size = 0.3
          )
        }
        
        # Add BLAST peptide overlay (red like selected peptide in isoform analysis)
        if (nrow(peptide_plot_data) > 0) {
          p <- p + geom_rect(
            data = peptide_plot_data,
            aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max, text = hover_text),
            fill = "#FF6B6B", color = "#D63031", alpha = 0.9, size = 0.5
          )
        }
        
        # Add transcript label on the left (like isoform analysis)
        transcript_label_x <- plot_start + (gene_start - plot_start) * 0.8  # Position in the extended left area
        p <- p + annotate(
          "text",
          x = transcript_label_x,
          y = transcript_y,
          label = match_info$transcript_id,
          hjust = 1, size = 4, color = "#2c3e50", fontface = "bold"
        )
        
        # Add BLAST peptide label
        if (nrow(peptide_plot_data) > 0) {
          p <- p + annotate(
            "text",
            x = mean(c(min(peptide_plot_data$start), max(peptide_plot_data$end))),
            y = transcript_y + 0.3,
            label = blast_peptide,
            hjust = 0.5, size = 3, fontface = "bold", color = "#D63031"
          )
        }
        
        # Style the plot like isoform analysis with chromosome information
        title_text <- paste0("BLAST Perfect Match: ", match_info$gene_symbol, " (chromosome ", chromosome, ")")
        x_axis_title <- if (use_compression && !is.null(compression_map)) {
          paste0("Genomic Position (compressed view) - chromosome ", chromosome)
        } else {
          paste0("Genomic Position (chromosome ", chromosome, ")")
        }
        
        p <- p + 
          labs(
            title = title_text,
            subtitle = paste0("Transcript: ", match_info$transcript_id, " | 100% Identity Match"),
            x = x_axis_title,
            y = ""
          ) +
          theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 12),
            axis.title = element_text(size = 11),
            panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
          )
        
        # Add compression-aware axis labels if compression is enabled
        if (use_compression && !is.null(compression_map)) {
          # Set custom axis breaks and labels to show original coordinates
          axis_breaks <- compression_map$coords$compressed_start
          axis_labels <- as.character(compression_map$coords$original_start)
          
          p <- p + scale_x_continuous(
            name = x_axis_title,
            breaks = axis_breaks,
            labels = axis_labels
          )
        }
        
        # Convert to plotly with clean hover and proper margins for transcript ID
        plotly_obj <- ggplotly(p, tooltip = "text") %>%
          clean_plotly_hover() %>%
          config(
            displayModeBar = TRUE,
            displaylogo = FALSE,
            modeBarButtonsToRemove = c("pan2d", "select2d", "lasso2d")
          ) %>%
          layout(
            title = list(
              text = paste0("BLAST Perfect Match: ", match_info$gene_symbol, " (", match_info$transcript_id, ")"),
              font = list(size = 16)
            ),
            margin = list(l = 150, r = 50, t = 100, b = 50),  # Extra left margin for transcript ID
            showlegend = FALSE
          )
        
        return(plotly_obj)
      })
      
    }, error = function(e) {
      cat("ERROR in create_cached_gtf_visualization:", e$message, "\n")
      return(empty_plotly_message(paste("Error creating cached GTF visualization:", e$message)))
    })
  }
  
  # State management for multi-gene visualization display
  multi_gene_blast_viz_state <- reactiveVal(FALSE)
  multi_gene_blast_data <- reactiveVal(NULL)
  
  # Detect 100% perfect matches availability for multi-gene visualization
  output$blast_perfect_matches_available <- reactive({
    tryCatch({
      results <- peptide_search_results()
      
      if (is.null(results) || nrow(results) == 0) {
        multi_gene_blast_viz_state(FALSE)
        return(FALSE)
      }
      
      # Check for any 100% identity matches
      perfect_matches <- results[!is.na(results$identity_percent) & results$identity_percent == 100, ]
      has_perfect_matches <- nrow(perfect_matches) > 0
      
      if (!has_perfect_matches) {
        multi_gene_blast_viz_state(FALSE)
      }
      
      return(has_perfect_matches)
    }, error = function(e) {
      cat("DEBUG: blast_perfect_matches_available error:", e$message, "\n")
      multi_gene_blast_viz_state(FALSE)
      return(FALSE)
    })
  })
  outputOptions(output, "blast_perfect_matches_available", suspendWhenHidden = FALSE)
  
  # Gene tabs UI for multi-gene BLAST visualization
  output$gene_tabs_ui <- renderUI({
    tryCatch({
      # Use stored data instead of recreating to ensure consistency
      stored_data <- multi_gene_blast_data()
      
      if (is.null(stored_data) || !stored_data$success) {
        return(div("No BLAST visualization data available"))
      }
      
      # Get data from stored results
      perfect_matches <- stored_data$perfect_matches
      unique_genes <- data.frame(
        gene_id = stored_data$unique_genes,
        gene_symbol = stored_data$gene_symbols,
        stringsAsFactors = FALSE
      )
      
      if (nrow(perfect_matches) == 0 || nrow(unique_genes) == 0) {
        return(div("No valid gene information available"))
      }
      
      cat("DEBUG: Creating tabs for", nrow(unique_genes), "genes\n")
      
      # Create tab panels for each gene
      tab_panels <- lapply(1:nrow(unique_genes), function(i) {
        gene_info <- unique_genes[i, ]
        gene_id <- gene_info$gene_id
        gene_symbol <- gene_info$gene_symbol
        
        # Create unique plot output ID for this gene
        plot_output_id <- paste0("gene_blast_plot_", gsub("[^A-Za-z0-9]", "_", gene_id))
        
        tabPanel(
          title = paste0(gene_symbol, " (", gene_id, ")"),
          value = paste0("tab_", gene_id),
          fluidRow(
            column(12,
              div(
                style = "margin-top: 15px;",
                plotlyOutput(plot_output_id, height = "500px")
              )
            )
          )
        )
      })
      
      # Create the tabsetPanel
      do.call(tabsetPanel, c(
        list(id = "gene_blast_tabs", type = "tabs"),
        tab_panels
      ))
      
    }, error = function(e) {
      cat("DEBUG: gene_tabs_ui error:", e$message, "\n")
      return(div(paste("Error creating gene tabs:", e$message)))
    })
  })
  
  # Dynamic server logic for individual gene plots - REACTIVE VERSION
  observeEvent(list(input$blast_viz_enzyme, input$blast_viz_miscleavage, input$blast_viz_intron_scale, multi_gene_blast_data()), {
    tryCatch({
      # Use stored data from multi_gene_blast_data instead of recreating
      stored_data <- multi_gene_blast_data()
      
      if (is.null(stored_data) || !stored_data$success) {
        return()
      }
      
      # Get data from stored results (ensures consistency)
      perfect_matches <- stored_data$perfect_matches
      unique_genes <- data.frame(
        gene_id = stored_data$unique_genes,
        gene_symbol = stored_data$gene_symbols,
        stringsAsFactors = FALSE
      )
      
      if (nrow(perfect_matches) == 0 || nrow(unique_genes) == 0) {
        return()
      }
      
      cat("DEBUG: Recreating renderPlotly outputs due to input changes\n")
      
      # Create plot outputs for each gene - using local() to capture variables by value
      for (i in 1:nrow(unique_genes)) {
        local({
          gene_info <- unique_genes[i, ]
          gene_id <- gene_info$gene_id
          gene_symbol <- gene_info$gene_symbol
          
          # Create plot output ID
          plot_output_id <- paste0("gene_blast_plot_", gsub("[^A-Za-z0-9]", "_", gene_id))
          
          cat("DEBUG: Creating renderPlotly for gene:", gene_id, "with output ID:", plot_output_id, "\n")
          
          # Create the renderPlotly for this gene - now recreated when inputs change
          output[[plot_output_id]] <- renderPlotly({
            req(input$blast_viz_enzyme, input$blast_viz_miscleavage)
            
            create_gene_blast_visualization(gene_id, gene_symbol, perfect_matches)
          })
        })
      }
      
    }, error = function(e) {
      cat("DEBUG: Dynamic gene plot creation error:", e$message, "\n")
    })
  }, ignoreNULL = FALSE)
  
  # Control multi-gene visualization visibility 
  output$multi_gene_blast_viz_available <- reactive({
    return(multi_gene_blast_viz_state())
  })
  outputOptions(output, "multi_gene_blast_viz_available", suspendWhenHidden = FALSE)
  
  # Hide multi-gene visualization when hide button is clicked
  observeEvent(input$hide_multi_gene_blast_viz, {
    cat("DEBUG: Hide multi-gene visualization button clicked\n")
    multi_gene_blast_viz_state(FALSE)
  })
  
  # Summary text for perfect matches
  output$blast_perfect_matches_summary <- renderText({
    tryCatch({
      results <- peptide_search_results()
      
      if (is.null(results) || nrow(results) == 0) {
        return("No BLAST results available")
      }
      
      # Filter to 100% identity matches
      perfect_matches <- results[!is.na(results$identity_percent) & results$identity_percent == 100, ]
      
      if (nrow(perfect_matches) == 0) {
        return("No 100% identity matches found")
      }
      
      # Count unique genes and transcripts
      unique_genes <- length(unique(perfect_matches$gene_id))
      unique_transcripts <- length(unique(perfect_matches$transcript_id))
      
      return(paste0("Found ", nrow(perfect_matches), " perfect matches across ", 
                   unique_genes, " genes and ", unique_transcripts, " transcript isoforms"))
      
    }, error = function(e) {
      cat("DEBUG: blast_perfect_matches_summary error:", e$message, "\n")
      return("Error displaying match summary")
    })
  })
  
  # Reset visualization state when new search is performed
  observeEvent(input$run_peptide_search, {
    cat("DEBUG: New search - resetting visualization state\n")
    multi_gene_blast_viz_state(FALSE)
    multi_gene_blast_data(NULL)
  })
  
  # Generate multi-gene BLAST visualization
  observeEvent(input$generate_multi_gene_blast_viz, {
    req(input$blast_viz_enzyme, input$blast_viz_miscleavage, input$peptide_search_query)
    
    # Source compression functions if not already loaded
    if (!exists("create_compression_map")) {
      source("R/coordinate_compression.R")
    }
    
    results <- peptide_search_results()
    if (is.null(results) || nrow(results) == 0) {
      showNotification("No BLAST results available", type = "warning")
      return()
    }
    
    # Filter to 100% identity matches only
    perfect_matches <- results[!is.na(results$identity_percent) & results$identity_percent == 100, ]
    if (nrow(perfect_matches) == 0) {
      showNotification("No 100% identity matches found", type = "warning")
      return()
    }
    
    unique_genes <- unique(perfect_matches$gene_id)
    query_peptide <- input$peptide_search_query
    
    withProgress(message = "Creating multi-gene BLAST visualization...", value = 0, {
      incProgress(0.1, detail = paste("Processing", length(unique_genes), "genes..."))
      
      tryCatch({
        # Create gene tabs directly using the working tab system
        incProgress(0.8, detail = "Creating gene tabs...")
        
        # Prepare gene data for tabs
        gene_symbols <- sapply(unique_genes, function(gene_id) {
          gtf_data <- load_gtf_visualization_data(gene_id)
          if (gtf_data$success) gtf_data$gene_symbol else gene_id
        })
        
        # Store data for tabs to use
        multi_gene_blast_data(list(
          success = TRUE,
          unique_genes = unique_genes,
          gene_symbols = gene_symbols,
          perfect_matches = perfect_matches
        ))
        multi_gene_blast_viz_state(TRUE)
        
        showNotification(
          paste("Successfully created gene tabs for", length(unique_genes), "genes"),
          type = "message",
          duration = 5
        )
        
      }, error = function(e) {
        cat("ERROR in generate_multi_gene_blast_viz:", e$message, "\n")
        multi_gene_blast_data(NULL)
        multi_gene_blast_viz_state(FALSE)
        showNotification(paste("Visualization error:", e$message), type = "error")
      })
    })
  })
  
  # Multi-gene BLAST plot output
  output$multi_gene_blast_plot <- renderPlotly({
    tryCatch({
      req(multi_gene_blast_viz_state() == TRUE)
      
      viz_data <- multi_gene_blast_data()
      if (is.null(viz_data) || !viz_data$success) {
        return(empty_plotly_message("No visualization data available"))
      }
      
      cat("DEBUG: Rendering multi-gene BLAST plot\n")
      
      return(viz_data$plot)
      
    }, error = function(e) {
      cat("ERROR in multi_gene_blast_plot:", e$message, "\n")
      return(empty_plotly_message("Error rendering multi-gene plot"))
    })
  })
  
  # Show selected match info - SAFE: Independent display function
  output$blast_selected_info <- renderText({
    tryCatch({
      req(input$peptide_search_results_table_rows_selected)
      results <- peptide_search_results()
      selected_row <- input$peptide_search_results_table_rows_selected[1]
      
      if (!is.null(results) && selected_row <= nrow(results)) {
        match_info <- results[selected_row, ]
        return(paste0(
          "Gene: ", match_info$gene_symbol, " (", match_info$gene_id, ") | ",
          "Transcript: ", match_info$transcript_id, " | ",
          "Peptide: ", input$peptide_search_query, " | ",
          "Identity: ", match_info$identity_percent, "% | ",
          "E-value: ", format(match_info$evalue, scientific = TRUE, digits = 2)
        ))
      }
      return("")
    }, error = function(e) {
      cat("DEBUG: blast_selected_info error:", e$message, "\n")
      return("Error displaying match info")
    })
  })
  
  # Main transcript visualization - SAFE: Isolated from existing functionality
  output$blast_transcript_plot <- renderPlotly({
    tryCatch({
      req(input$peptide_search_results_table_rows_selected)
      req(blast_visualization_state() == TRUE)  # Only render when visualization is requested
      
      # Get selected match details
      results <- peptide_search_results()
      selected_row <- input$peptide_search_results_table_rows_selected[1]
      
      if (is.null(results) || selected_row > nrow(results)) {
        return(empty_plotly_message("Invalid selection"))
      }
      
      match_info <- results[selected_row, ]
      
      # Only proceed for 100% matches
      if (is.na(match_info$identity_percent) || match_info$identity_percent != 100) {
        return(empty_plotly_message("Select a 100% perfect match to see transcript visualization"))
      }
      
      cat("DEBUG: Creating visualization for 100% match:", match_info$transcript_id, "\n")
      
      withProgress(message = "Loading transcript visualization...", value = 0, {
        incProgress(0.2, detail = "Loading transcript structure...")
        
        # Step 1: Load transcript structure using existing system
        transcript_structure <- tryCatch({
          load_transcript_exons(
            c(match_info$transcript_id), 
            match_info$gene_id
          )
        }, error = function(e) {
          cat("DEBUG: load_transcript_exons error:", e$message, "\n")
          return(list(success = FALSE, message = paste("Transcript loading error:", e$message)))
        })
        
        # Handle different return types from load_transcript_exons
        if (is.null(transcript_structure)) {
          return(empty_plotly_message("Error: transcript structure is NULL"))
        }
        
        if (is.atomic(transcript_structure)) {
          return(empty_plotly_message("Error: transcript structure is atomic vector - transcript not found"))
        }
        
        if (!is.list(transcript_structure) || is.null(transcript_structure$success) || !transcript_structure$success) {
          # Transcript not found in internal database - use cached GTF via gene-first approach
          cat("DEBUG: Transcript not in internal database, using cached GTF via gene-first approach\n")
          
          return(create_cached_gtf_visualization(
            match_info = match_info,
            blast_peptide = input$peptide_search_query
          ))
        }
        
        incProgress(0.3, detail = "Mapping BLAST peptide to genome...")
        
        # Step 2: Map BLAST peptide to genomic coordinates
        peptide_mapping <- map_blast_peptide_to_transcript(
          blast_peptide = input$peptide_search_query,
          transcript_id = match_info$transcript_id,
          gene_id = match_info$gene_id,
          transcript_structure = transcript_structure
        )
        
        if (is.null(peptide_mapping) || !peptide_mapping$success) {
          return(empty_plotly_message("Error mapping peptide to genomic coordinates"))
        }
        
        incProgress(0.3, detail = "Creating visualization...")
        
        # Step 3: Create visualization using existing system
        # Create a minimal processed_data structure for compatibility
        blast_processed_data <- list(
          genes = match_info$gene_id,
          gene_symbols = match_info$gene_symbol,
          gene_lookup = setNames(match_info$gene_symbol, match_info$gene_id)
        )
        
        # Create plot data using existing visualization function
        plot_result <- create_peptide_plot_data(
          exons_result = transcript_structure,
          transcript_ids = match_info$transcript_id,
          gene_symbol = match_info$gene_symbol,
          gene_id = match_info$gene_id,
          processed_data = blast_processed_data,
          protease = "blastp"  # Special marker for BLAST results
        )
        
        if (!plot_result$success) {
          return(empty_plotly_message(paste("Visualization error:", plot_result$message)))
        }
        
        incProgress(0.2, detail = "Rendering plot...")
        
        # Step 4: Customize plot for BLAST visualization
        p <- plot_result$plot
        
        # Add BLAST peptide overlay to the plot
        if (!is.null(peptide_mapping$genomic_ranges) && length(peptide_mapping$genomic_ranges) > 0) {
          blast_ranges <- peptide_mapping$genomic_ranges
          
          # Create BLAST peptide overlay data
          blast_peptide_df <- data.frame(
            start = start(blast_ranges),
            end = end(blast_ranges),
            y_min = 0.8,  # Position above other elements
            y_max = 1.2,
            peptide = blast_ranges$peptide,
            stringsAsFactors = FALSE
          )
          
          # Add BLAST peptide as highlighted overlay
          p <- p + 
            geom_rect(
              data = blast_peptide_df,
              aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max),
              fill = "#ff6b6b",  # Red for BLAST match
              color = "#d63031",
              alpha = 0.8,
              size = 0.5
            ) +
            annotate(
              "text", 
              x = mean(c(min(blast_peptide_df$start), max(blast_peptide_df$end))),
              y = 1.4,
              label = paste("BLAST Match:", input$peptide_search_query),
              hjust = 0.5,
              size = 3.5,
              color = "#d63031",
              fontface = "bold"
            )
        }
        
        # Update plot title and labels for BLAST context
        p <- p + 
          labs(
            title = paste0("BLAST Perfect Match Visualization"),
            subtitle = paste0("Gene: ", match_info$gene_symbol, " | Transcript: ", match_info$transcript_id, " | 100% Identity"),
            caption = "Red overlay shows BLAST peptide mapping"
          )
        
        # Convert to plotly
        plotly_obj <- ggplotly(p, tooltip = "text") %>%
          config(
            displayModeBar = TRUE,
            displaylogo = FALSE,
            modeBarButtonsToRemove = c("pan2d", "select2d", "lasso2d")
          ) %>%
          layout(
            title = list(
              text = paste0("BLAST Perfect Match: ", match_info$gene_symbol, " (", match_info$transcript_id, ")"),
              font = list(size = 16)
            ),
            margin = list(l = 50, r = 50, t = 100, b = 50)
          )
        
        return(plotly_obj)
      })
      
    }, error = function(e) {
      cat("DEBUG: blast_transcript_plot error:", e$message, "\n")
      cat("DEBUG: Error details:", toString(e), "\n")
      return(empty_plotly_message("Visualization error - check console for details"))
    })
  })
  
  # Download handler for transcript plot - SAFE: Independent functionality
  output$download_blast_plot <- downloadHandler(
    filename = function() {
      paste0("blast_transcript_", input$peptide_search_query, "_", 
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      # TODO: Implement actual plot download
      # For now, create placeholder
      cat("Download functionality will be implemented with visualization\n")
    }
  )
  
  #===============================================================================
  # GENE-SPECIFIC BLAST VISUALIZATION FUNCTION
  #===============================================================================
  
  create_gene_blast_visualization <- function(gene_id, gene_symbol, perfect_matches) {
    tryCatch({
      cat("DEBUG: Creating BLAST visualization for gene:", gene_id, "\n")
      
      # Step 1: Get selected enzyme and miscleavage from UI controls
      selected_enzyme <- if(is.null(input$blast_viz_enzyme)) "trp" else input$blast_viz_enzyme
      selected_miscleavage <- if(is.null(input$blast_viz_miscleavage)) "no_miss_cleavage" else input$blast_viz_miscleavage
      
      # Get intron display setting
      use_compression <- !is.null(input$blast_viz_intron_scale) && input$blast_viz_intron_scale == "compressed"
      cat("DEBUG: Intron display mode:", ifelse(use_compression, "compressed", "true_scale"), "\n")
      
      cat("DEBUG: Using enzyme:", selected_enzyme, "miscleavage:", selected_miscleavage, "\n")
      
      # Step 2: Load pre-computed GTF data for gene structure
      gtf_data <- load_gtf_visualization_data(gene_id)
      
      if (!gtf_data$success) {
        return(empty_plotly_message(gtf_data$message))
      }
      
      # Step 3: Extract GTF structure data
      exons_by_transcript <- gtf_data$exons_by_transcript
      cds_by_transcript <- gtf_data$cds_by_transcript
      transcript_ids <- gtf_data$transcript_ids
      
      if (length(transcript_ids) == 0) {
        return(empty_plotly_message("No transcripts found for this gene"))
      }
      
      # Step 4: Load gene peptide data from RDS file
      gene_file <- paste0("data/genes/", selected_miscleavage, "/", gene_id, ".rds")
      gene_peptides <- NULL
      
      if (file.exists(gene_file)) {
        cat("DEBUG: Loading gene peptides from:", gene_file, "\n")
        gene_peptides <- readRDS(gene_file)
        cat("DEBUG: Loaded gene data with", nrow(gene_peptides), "rows\n")
      } else {
        cat("DEBUG: Gene peptide file not found:", gene_file, "\n")
      }
      
      # Step 5: Create compression map if needed
      compression_map <- NULL
      if (use_compression && length(exons_by_transcript) > 0) {
        cat("DEBUG: Creating compression map for exons\n")
        compression_map <- create_compression_map(exons_by_transcript)
      }
      
      # Step 6: Use gene boundaries with proper scaling
      padding <- 1000  # Reduced padding for better scaling
      gene_start <- gtf_data$gene_start - padding
      gene_end <- gtf_data$gene_end + padding
      chromosome <- gtf_data$chromosome
      
      # Apply compression to gene boundaries if needed
      plot_gene_start <- gene_start
      plot_gene_end <- gene_end
      if (use_compression && !is.null(compression_map)) {
        plot_gene_start <- compression_map$compress(gene_start)
        plot_gene_end <- compression_map$compress(gene_end)
        cat("DEBUG: Compressed gene range from", gene_start, "-", gene_end, "to", plot_gene_start, "-", plot_gene_end, "\n")
      }
      
      # Step 7: Create plot data frames
      transcript_df <- data.frame(
        transcript = transcript_ids,
        y_position = seq_along(transcript_ids),
        stringsAsFactors = FALSE
      )
      
      # Initialize data frames
      exon_df <- data.frame(
        transcript = character(), y_position = numeric(), start = numeric(), 
        end = numeric(), exon_number = integer(), stringsAsFactors = FALSE
      )
      
      cds_df <- data.frame(
        transcript = character(), y_position = numeric(), start = numeric(),
        end = numeric(), cds_number = integer(), stringsAsFactors = FALSE
      )
      
      peptide_df <- data.frame(
        transcript = character(), y_position = numeric(), start = numeric(),
        end = numeric(), peptide = character(), hover_text = character(), 
        stringsAsFactors = FALSE
      )
      
      # Step 7: Process transcripts and load peptide data
      for (i in seq_along(transcript_ids)) {
        tx <- transcript_ids[i]
        
        # Add exons
        if (tx %in% names(exons_by_transcript) && length(exons_by_transcript[[tx]]) > 0) {
          tx_exons <- exons_by_transcript[[tx]]
          tx_exons <- tx_exons[order(start(tx_exons))]
          
          for (j in seq_along(tx_exons)) {
            exon_start <- start(tx_exons[j])
            exon_end <- end(tx_exons[j])
            
            # Apply compression if enabled
            if (use_compression && !is.null(compression_map)) {
              exon_start <- compression_map$compress(exon_start)
              exon_end <- compression_map$compress(exon_end)
            }
            
            exon_df <- rbind(exon_df, data.frame(
              transcript = tx, y_position = i, start = exon_start,
              end = exon_end, exon_number = j, stringsAsFactors = FALSE
            ))
          }
        }
        
        # Add CDS
        if (tx %in% names(cds_by_transcript) && length(cds_by_transcript[[tx]]) > 0) {
          tx_cds <- cds_by_transcript[[tx]]
          tx_cds <- tx_cds[order(start(tx_cds))]
          
          for (j in seq_along(tx_cds)) {
            cds_start <- start(tx_cds[j])
            cds_end <- end(tx_cds[j])
            
            # Apply compression if enabled
            if (use_compression && !is.null(compression_map)) {
              cds_start <- compression_map$compress(cds_start)
              cds_end <- compression_map$compress(cds_end)
            }
            
            cds_df <- rbind(cds_df, data.frame(
              transcript = tx, y_position = i, start = cds_start,
              end = cds_end, cds_number = j, stringsAsFactors = FALSE
            ))
          }
        }
        
        # Add peptides from gene data
        if (!is.null(gene_peptides)) {
          # Find which row in gene_peptides corresponds to this transcript
          tx_row <- which(gene_peptides$txID == tx)
          
          if (length(tx_row) > 0) {
            mapped_ranges_col <- paste0(selected_enzyme, "Peps_mapped_ranges")
            
            if (mapped_ranges_col %in% names(gene_peptides)) {
              mapped_ranges_list <- gene_peptides[[mapped_ranges_col]]
              
              if (!is.null(mapped_ranges_list) && length(mapped_ranges_list) >= tx_row[1]) {
                # Get the mapped ranges for THIS specific transcript
                mapped_ranges <- mapped_ranges_list[[tx_row[1]]]
                
                if (class(mapped_ranges)[1] == "GRanges" && length(mapped_ranges) > 0) {
                  # Get peptide sequences from metadata
                  peptide_seqs <- mcols(mapped_ranges)$peptide
                  
                  for (k in 1:length(mapped_ranges)) {
                    range_k <- mapped_ranges[k]
                    peptide_seq <- if (!is.null(peptide_seqs)) peptide_seqs[k] else paste("Peptide", k)
                    
                    pep_start <- start(range_k)
                    pep_end <- end(range_k)
                    original_start <- pep_start  # Keep original for hover text
                    original_end <- pep_end
                    
                    # Apply compression if enabled
                    if (use_compression && !is.null(compression_map)) {
                      pep_start <- compression_map$compress(pep_start)
                      pep_end <- compression_map$compress(pep_end)
                    }
                    
                    peptide_df <- rbind(peptide_df, data.frame(
                      transcript = tx,
                      y_position = i,
                      start = pep_start,
                      end = pep_end,
                      peptide = peptide_seq,
                      hover_text = clean_hover_text(paste0(
                        "Peptide: ", peptide_seq,
                        "<br>Position: ", original_start, "-", original_end,
                        "<br>Transcript: ", tx
                      )),
                      stringsAsFactors = FALSE
                    ))
                  }
                }
              }
            }
          }
        }
      }
      
      # Step 8: Create the base plot
      p <- ggplot() +
        # Add transcript lines
        geom_segment(data = transcript_df, 
                    aes(x = gene_start, xend = gene_end, y = y_position, yend = y_position),
                    linewidth = 0.5, color = "grey70") +
        
        # Add exon blocks
        geom_rect(data = exon_df,
                 aes(xmin = start, xmax = end, ymin = y_position - 0.3, ymax = y_position + 0.3,
                     fill = "Exons"), color = "black") +
        
        # Add CDS overlay
        geom_rect(data = cds_df,
                 aes(xmin = start, xmax = end, ymin = y_position - 0.25, ymax = y_position + 0.25,
                     fill = "CDS"), color = "black") +
        
        # Add peptide overlays (only if peptide data is available)
        {if (nrow(peptide_df) > 0) {
          geom_rect(data = peptide_df,
                   aes(xmin = start, xmax = end, ymin = y_position - 0.15, ymax = y_position + 0.15,
                       fill = "Peptides", text = hover_text), color = "black", linewidth = 0.2)
        } else {
          NULL
        }} +
        
        # Add transcript labels
        geom_text(data = transcript_df,
                 aes(x = plot_gene_start - 100, y = y_position, label = transcript),
                 hjust = 1, size = 3.5)
      
      # Step 9: Collect BLAST peptide overlay data using R_backup working pattern
      gene_perfect_matches <- perfect_matches[perfect_matches$gene_id == gene_id, ]
      matching_transcripts <- unique(gene_perfect_matches$transcript_id)
      
      # Initialize BLAST overlay data frame to collect all data first
      blast_overlay_df <- data.frame(
        start = numeric(), end = numeric(), y_min = numeric(), y_max = numeric(),
        transcript = character(), hover_text = character(), stringsAsFactors = FALSE
      )
      
      if (length(matching_transcripts) > 0) {
        cat("DEBUG: Collecting BLAST overlay data for", length(matching_transcripts), "transcripts with 100% matches\n")
        
        # Load gene details using R_backup pattern
        cache_file <- find_gtf_cache_file(gene_id)
        if (!is.null(cache_file)) {
          gene_details <- readRDS(cache_file)
          
          if (!is.null(gene_details)) {
            # Create transcript_structure using R_backup pattern
            transcript_structure <- list(
              success = TRUE,
              exons = gene_details$exons_by_transcript,
              cds = gene_details$cds_by_transcript
            )
            
            cat("DEBUG: Successfully created transcript structure for", gene_id, "\n")
            
            # Collect BLAST peptide data for each transcript with 100% match
            for (transcript_id in matching_transcripts) {
              cat("DEBUG: Mapping BLAST peptide to transcript", transcript_id, "\n")
              
              peptide_mapping <- map_blast_peptide_to_transcript(
                blast_peptide = input$peptide_search_query,
                transcript_id = transcript_id,
                gene_id = gene_id,
                transcript_structure = transcript_structure
              )
              
              if (!is.null(peptide_mapping) && peptide_mapping$success) {
                blast_ranges <- peptide_mapping$genomic_ranges
                transcript_y <- which(transcript_ids == transcript_id)
                
                if (length(transcript_y) > 0 && length(blast_ranges) > 0) {
                  cat("DEBUG: Collecting BLAST overlay data for transcript", transcript_id, "at y-position", transcript_y, "\n")
                  
                  # Get coordinates and apply compression like R_backup
                  blast_starts <- start(blast_ranges)
                  blast_ends <- end(blast_ranges)
                  
                  if (use_compression && !is.null(compression_map)) {
                    blast_starts <- sapply(blast_starts, compression_map$compress)
                    blast_ends <- sapply(blast_ends, compression_map$compress)
                  }
                  
                  # Add to BLAST overlay data frame (thicker than exons)
                  transcript_blast_data <- data.frame(
                    start = blast_starts,
                    end = blast_ends,
                    y_min = transcript_y - 0.35,  # Thicker than exons (±0.3)
                    y_max = transcript_y + 0.35,
                    transcript = transcript_id,
                    hover_text = paste0("BLAST: ", input$peptide_search_query, " | ", blast_starts, "-", blast_ends),
                    stringsAsFactors = FALSE
                  )
                  
                  blast_overlay_df <- rbind(blast_overlay_df, transcript_blast_data)
                }
              }
            }
          }
        }
      }
      
      # Add single BLAST overlay layer if we have data (R_backup pattern)
      if (nrow(blast_overlay_df) > 0) {
        cat("DEBUG: Adding BLAST overlay layer with", nrow(blast_overlay_df), "segments\n")
        p <- p + geom_rect(
          data = blast_overlay_df,
          aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max, fill = "BLAST Match", text = hover_text),
          color = "#d63031",
          linewidth = 0.5
        )
      }
      
      # Step 10: Styling and legend
      has_cds <- nrow(cds_df) > 0
      has_peptides <- nrow(peptide_df) > 0
      has_blast <- length(matching_transcripts) > 0
      
      fill_values <- c(
        "Exons" = "rgba(77, 175, 74, 0.8)",
        "CDS" = "rgba(255, 221, 0, 0.8)", 
        "Peptides" = "rgba(255, 140, 0, 0.8)",  # Deeper orange color
        "BLAST Match" = "rgba(255, 107, 107, 0.9)"
      )
      
      fill_breaks <- c("Exons")
      fill_labels <- c("Exons")
      
      if (has_cds) {
        fill_breaks <- c(fill_breaks, "CDS")
        fill_labels <- c(fill_labels, "CDS")
      }
      if (has_peptides) {
        fill_breaks <- c(fill_breaks, "Peptides")
        fill_labels <- c(fill_labels, "Peptides")
      }
      if (has_blast) {
        fill_breaks <- c(fill_breaks, "BLAST Match")
        fill_labels <- c(fill_labels, "BLAST Match")
      }
      
      p <- p + 
        scale_fill_manual(values = fill_values, breaks = fill_breaks, labels = fill_labels, name = "")
      
      # Add proper axis ticks for compression
      if (use_compression && !is.null(compression_map)) {
        axis_breaks <- compression_map$coords$compressed_start
        axis_labels <- as.character(compression_map$coords$original_start)
        
        p <- p + scale_x_continuous(
          name = paste0("Genomic Position (chromosome ", chromosome, " - compressed view)"),
          breaks = axis_breaks,
          labels = axis_labels
        )
      }
      
      p <- p +
        theme_minimal() +
        theme(
          axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(), legend.position = "bottom", legend.box = "horizontal",
          legend.margin = margin(6, 6, 6, 6)
        ) +
        labs(
          x = if (!use_compression || is.null(compression_map)) paste0("Genomic Position (chromosome ", chromosome, ")") else NULL,
          title = paste0("Gene: ", gene_symbol, " (", gene_id, ")"),
          subtitle = paste0("Chromosome ", chromosome, " | Enzyme: ", selected_enzyme, 
                           " | Miscleavage: ", selected_miscleavage,
                           if(has_blast) " | Red = BLAST matches" else "")
        ) +
        coord_cartesian(xlim = c(plot_gene_start, plot_gene_end), ylim = c(0.5, length(transcript_ids) + 0.5))
      
      # Step 11: Convert to plotly using optimized function for legend interactivity
      plotly_obj <- create_optimized_plotly(p) %>%
        layout(
          title = list(text = paste0("Gene: ", gene_symbol, " (", chromosome, ")"), font = list(size = 14)),
          margin = list(l = 120, r = 50, t = 80, b = 50), yaxis = list(title = "Transcripts")
        )
      
      return(plotly_obj)
      
    }, error = function(e) {
      cat("DEBUG: create_gene_blast_visualization error:", e$message, "\n")
      return(empty_plotly_message(paste("Error creating visualization:", e$message)))
    })
  }
  
