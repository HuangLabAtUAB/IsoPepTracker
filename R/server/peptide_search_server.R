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
        # Get BLAST parameters from input
        evalue_threshold <- if (!is.null(input$peptide_search_evalue)) input$peptide_search_evalue else 10
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
          
          # Add default values for compatibility with navigation
          enhanced_results$protease_used <- "blastp"
          enhanced_results$miscleavage_type_used <- "none"
          enhanced_results$peptide <- peptide_query  # Use the searched peptide sequence
          enhanced_results$txID <- enhanced_results$transcript_id
          enhanced_results$geneID <- enhanced_results$gene_id
          enhanced_results$geneSymbol <- enhanced_results$gene_symbol
          
          # Store enhanced results
          peptide_search_results(enhanced_results)
          
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
      
      # Add all 12 enzyme availability columns if they exist
      enzyme_cols <- c("trp_no_miss", "trp_upto2miss", "chymo_no_miss", "chymo_upto2miss", 
                       "aspn_no_miss", "aspn_upto2miss", "lysc_no_miss", "lysc_upto2miss",
                       "lysn_no_miss", "lysn_upto2miss", "gluc_no_miss", "gluc_upto2miss")
      available_enzyme_cols <- enzyme_cols[enzyme_cols %in% names(results)]
      
      display_cols <- c(base_cols, available_enzyme_cols)
      display_df <- results[, display_cols]
      
      base_colnames <- c("Query Peptide", "Gene ID", "Gene Symbol", "Transcript ID", 
                        "Identity %", "E-value", "Bit Score")
      
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
      
      # Sort by bit score (descending) and E-value (ascending)
      display_df <- display_df[order(-display_df$`Bit Score`, display_df$`E-value`), ]
      
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
        order = list(list(2, 'asc'))  # Sort by Gene Symbol
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
  
  # State management for visualization display
  blast_visualization_state <- reactiveVal(FALSE)
  
  # Detect 100% match selection - SAFE: No changes to existing functionality
  output$blast_perfect_match_selected <- reactive({
    tryCatch({
      selected_row <- input$peptide_search_results_table_rows_selected
      results <- peptide_search_results()
      
      if (is.null(selected_row) || is.null(results) || length(selected_row) == 0) {
        # Reset visualization state when no row selected
        blast_visualization_state(FALSE)
        return(FALSE)
      }
      
      if (selected_row[1] > nrow(results)) {
        blast_visualization_state(FALSE)
        return(FALSE)
      }
      
      identity <- results$identity_percent[selected_row[1]]
      is_perfect_match <- !is.na(identity) && identity == 100
      
      # Reset visualization state when different row is selected
      if (!is_perfect_match) {
        blast_visualization_state(FALSE)
      }
      
      return(is_perfect_match)
    }, error = function(e) {
      cat("DEBUG: blast_perfect_match_selected error:", e$message, "\n")
      blast_visualization_state(FALSE)
      return(FALSE)  # Fail safely - existing functionality unaffected
    })
  })
  outputOptions(output, "blast_perfect_match_selected", suspendWhenHidden = FALSE)
  
  # Control visualization visibility based on button clicks
  output$blast_visualization_visible <- reactive({
    return(blast_visualization_state())
  })
  outputOptions(output, "blast_visualization_visible", suspendWhenHidden = FALSE)
  
  # Show visualization when button is clicked
  observeEvent(input$show_blast_visualization, {
    cat("DEBUG: Show visualization button clicked\n")
    blast_visualization_state(TRUE)
  })
  
  # Hide visualization when hide button is clicked
  observeEvent(input$hide_blast_visualization, {
    cat("DEBUG: Hide visualization button clicked\n")
    blast_visualization_state(FALSE)
  })
  
  # Reset visualization state when new search is performed
  observeEvent(input$run_peptide_search, {
    cat("DEBUG: New search - resetting visualization state\n")
    blast_visualization_state(FALSE)
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
  
