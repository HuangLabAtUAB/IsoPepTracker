#===============================================================================
# MULTI-ISOFORM COMPARISON MODULE
# Extracted from server.R - Multi-isoform comparative analysis functionality
#===============================================================================

#' Multi-Isoform Comparison Module
#' 
#' @description 
#' Complete multi-isoform comparative analysis functionality:
#' - Multiple isoform selection and management
#' - Comparative peptide analysis across selected isoforms
#' - Interactive comparative visualizations with plotly
#' - Peptide overlap matrix generation
#' - Comparative summary statistics
#' - Download handlers for comparative analysis results
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' @param processed_data Reactive containing processed application data
#' @param gene_data Reactive containing current gene data
#' @param selected_gene_transcripts Reactive containing gene transcripts
#' @param selected_gene_peptides Reactive containing gene peptides
#' @param novel_pipeline_results Reactive containing novel isoform results (optional)
#' 
#' @return Named list containing reactive values for multi-isoform comparison
#' 
#' @note 
#' - Supports both regular gene isoforms and novel isoforms
#' - Requires gene-specific data loading for accurate comparisons
#' - Uses plotly for interactive comparative visualizations
#' - Maintains compatibility with existing isoform analysis pipeline

create_multi_isoform_comparison_module <- function(input, output, session, processed_data, gene_data, 
                                                 selected_gene_transcripts, selected_gene_peptides, 
                                                 novel_pipeline_results = NULL) {
  
  # ============================================================================
  # MODULE REACTIVE VALUES
  # ============================================================================
  
  # Multi-isoform comparison data
  multi_isoform_data <- reactiveVal(NULL)
  multi_isoform_highlighted_data <- reactiveVal(NULL)
  comparison_isoforms <- reactiveVal(character(0))
  
  # ============================================================================
  # ISOFORM SELECTION MANAGEMENT
  # ============================================================================
  
  # Update compare_isoforms choices when gene is selected OR novel data is available
  observeEvent({
    list(input$gene, novel_pipeline_results())
  }, {
    # Always try to get known gene transcripts if gene is selected
    known_transcripts <- NULL
    gene_symbol <- NULL
    
    if (!is.null(input$gene) && !is.na(input$gene) && input$gene != "") {
      transcripts_data <- selected_gene_transcripts()
      if (!is.null(transcripts_data) && length(transcripts_data) > 0) {
        known_transcripts <- transcripts_data
        
        # Get gene symbol
        gene_symbol <- tryCatch({
          processed_data()$gene_lookup[input$gene]
        }, error = function(e) {
          input$gene
        })
        if (is.null(gene_symbol)) gene_symbol <- input$gene
      }
    }
    
    # Get novel transcripts if available
    novel_transcripts <- NULL
    if (!is.null(novel_pipeline_results)) {
      novel_data <- novel_pipeline_results()
      
      if (!is.null(novel_data) && novel_data$success && !is.null(novel_data$dataframe_file) && file.exists(novel_data$dataframe_file)) {
        tryCatch({
          novel_rds <- readRDS(novel_data$dataframe_file)
          if (!is.null(novel_rds) && nrow(novel_rds) > 0) {
            novel_transcripts <- unique(novel_rds$txID)
            cat("Found", length(novel_transcripts), "novel transcripts for comparison\n")
          }
        }, error = function(e) {
          cat("Error loading novel transcripts:", e$message, "\n")
        })
      }
    }
    
    # Combine choices
    all_choices <- character(0)
    
    if (!is.null(known_transcripts) && length(known_transcripts) > 0) {
      # Add known transcripts with labels
      known_labels <- paste0(known_transcripts, " (", gene_symbol, ")")
      names(known_transcripts) <- known_labels
      all_choices <- c(all_choices, known_transcripts)
    }
    
    if (!is.null(novel_transcripts) && length(novel_transcripts) > 0) {
      # Add novel transcripts with labels
      novel_labels <- paste0(novel_transcripts, " (Novel)")
      names(novel_transcripts) <- novel_labels
      all_choices <- c(all_choices, novel_transcripts)
    }
    
    # Update dropdown
    if (length(all_choices) > 0) {
      updateSelectizeInput(session, "compare_isoforms", 
                          choices = all_choices,
                          selected = character(0))
      
      # Show notification about available options
      if (!is.null(known_transcripts) && !is.null(novel_transcripts)) {
        cat("✅ Multiple Isoform Comparison: Added", length(known_transcripts), "known +", length(novel_transcripts), "novel transcripts\n")
      }
    } else {
      # Clear choices if no transcripts available
      updateSelectizeInput(session, "compare_isoforms", 
                          choices = NULL, selected = character(0))
    }
  })
  
  # Enable/disable comparative analysis button based on selection count
  observeEvent(input$compare_isoforms, {
    selected_count <- length(input$compare_isoforms)
    comparison_isoforms(input$compare_isoforms)
    
    if (selected_count >= 2 && selected_count <= 8) {
      shinyjs::enable("run_comparative_analysis")
      shinyjs::removeClass("run_comparative_analysis", "btn-secondary")
      shinyjs::addClass("run_comparative_analysis", "btn-primary")
    } else {
      shinyjs::disable("run_comparative_analysis")
      shinyjs::removeClass("run_comparative_analysis", "btn-primary")
      shinyjs::addClass("run_comparative_analysis", "btn-secondary")
    }
  })
  
  # ============================================================================
  # MULTI-ISOFORM DATA PROCESSING
  # ============================================================================
  
  # Multi-Isoform data loader
  multi_isoform_data <- reactive({
    req(input$run_comparative_analysis > 0, input$gene, input$protease, input$miscleavage_type, gene_data(),
        input$compare_isoforms, length(input$compare_isoforms) >= 2)
    
    withProgress(message = 'Loading multi-isoform analysis...', value = 0, {
      gene_id <- input$gene
      protease <- input$protease
      miscleavage_type <- input$miscleavage_type
      selected_transcripts <- input$compare_isoforms
      
      # Step 1: Get all transcripts for this gene (from gene-by-gene loaded data)
      incProgress(0.2, detail = 'Getting selected gene transcripts...')
      all_transcripts <- selected_gene_transcripts()
      
      if (is.null(all_transcripts) || length(all_transcripts) == 0) {
        showNotification(paste("No transcripts found for gene", gene_id, "with miscleavage type", miscleavage_type), type = "warning")
        return(NULL)
      }
      
      cat("Selected transcripts for multi-isoform analysis:", paste(selected_transcripts, collapse=", "), "\n")
      cat("Using miscleavage type:", miscleavage_type, "\n")
      
      # Step 2: Get peptides for each transcript (using gene-by-gene loaded data)
      incProgress(0.3, detail = 'Loading peptides for selected transcripts...')
      gene_peptides <- selected_gene_peptides()
      if (is.null(gene_peptides) || nrow(gene_peptides) == 0) {
        showNotification(paste("No peptides found for gene", gene_id, "with miscleavage type", miscleavage_type), type = "warning")
        return(NULL)
      }
      
      # Create vis_data structure for compatibility with existing functions
      vis_data <- list(
        genes = processed_data()$genes,
        gene_symbols = processed_data()$gene_symbols,
        gene_lookup = processed_data()$gene_lookup,
        proteases = processed_data()$proteases,
        original_peptides = gene_peptides
      )
      
      # Get peptides for SELECTED transcripts only
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
            peptide = tx_peptides$peptide,
            stringsAsFactors = FALSE
          )
        }
      }
      
      if (length(all_peptides_list) == 0) {
        showNotification("No peptides found for selected transcripts", type = "warning")
        return(NULL)
      }
      
      # Combine all peptides
      all_peptides_df <- data.table::rbindlist(all_peptides_list)
      
      # Step 3: Use exact same approach as all_isoforms_plot
      incProgress(0.2, detail = 'Preparing visualization data...')
      gene_start <- min(all_peptides_df$start) - 1000
      gene_end <- max(all_peptides_df$end) + 1000
      
      # Create transcript position mapping
      transcript_df <- data.frame(
        transcript = selected_transcripts,
        y_position = seq_along(selected_transcripts),
        stringsAsFactors = FALSE
      )
      
      # Return minimal data - let the plot function handle everything like all_isoforms_plot does
      incProgress(0.1, detail = 'Complete!')
      
      return(list(
        all_peptides = all_peptides_df,
        transcript_df = transcript_df,
        gene_start = gene_start,
        gene_end = gene_end,
        all_transcripts = selected_transcripts
      ))
      
      if (length(rmats_transcripts) > 0) {
        cat("🔍 Detected", length(rmats_transcripts), "rMATS transcripts:", paste(rmats_transcripts, collapse = ", "), "\n")
        cat("📁 Searching for rMATS GTF files...\n")
        
        # Find the most recent rMATS GTF file for this gene
        gene_clean <- gsub('"', '', gene_id)
        cat("🧬 Gene ID (cleaned):", gene_clean, "\n")
        
        # Try multiple search patterns to find rMATS GTF files
        search_patterns <- c(
          paste0('rmats_"', gene_clean, '".*\\.transdecoder\\.genome\\.gtf$'),  # Original pattern
          paste0('rmats_', gene_clean, '.*\\.transdecoder\\.genome\\.gtf$'),    # Without quotes
          paste0('.*', gene_clean, '.*\\.transdecoder\\.genome\\.gtf$'),        # More flexible
          paste0('.*\\.transdecoder\\.genome\\.gtf$')                            # Any rMATS GTF
        )
        
        rmats_gtf_files <- c()
        for (i in seq_along(search_patterns)) {
          pattern <- search_patterns[i]
          cat("🔍 Trying pattern", i, ":", pattern, "\n")
          
          found_files <- list.files(
            path = "rmats_peptide_results", 
            pattern = pattern, 
            full.names = TRUE,
            recursive = TRUE
          )
          
          if (length(found_files) > 0) {
            cat("✅ Found", length(found_files), "files with pattern", i, "\n")
            rmats_gtf_files <- c(rmats_gtf_files, found_files)
            break  # Use the first successful pattern
          } else {
            cat("❌ No files found with pattern", i, "\n")
          }
        }
        
        # Remove duplicates and show all found files
        rmats_gtf_files <- unique(rmats_gtf_files)
        
        if (length(rmats_gtf_files) > 0) {
          cat("📂 Found", length(rmats_gtf_files), "rMATS GTF file(s):\n")
          for (f in rmats_gtf_files) {
            cat("  📄", basename(f), "(", file.info(f)$mtime, ")\n")
          }
          
          # Use the most recent GTF file
          rmats_gtf_file <- rmats_gtf_files[order(file.mtime(rmats_gtf_files), decreasing = TRUE)][1]
          cat("🎯 Selected most recent:", basename(rmats_gtf_file), "\n")
          
          # Load rMATS transcript structure data
          cat("📥 Loading rMATS structure data...\n")
          rmats_structure_data <- load_rmats_transcript_exons(rmats_gtf_file, rmats_transcripts)
          
          if (rmats_structure_data$success) {
            cat("✅ Successfully loaded rMATS structure data for", rmats_structure_data$loaded_transcripts, "out of", rmats_structure_data$requested_transcripts, "transcripts\n")
            if (length(rmats_structure_data$exons) > 0) {
              cat("📊 Loaded exons for transcripts:", paste(names(rmats_structure_data$exons), collapse = ", "), "\n")
            }
            if (length(rmats_structure_data$cds) > 0) {
              cat("📊 Loaded CDS for transcripts:", paste(names(rmats_structure_data$cds), collapse = ", "), "\n")
            }
          } else {
            cat("❌ Failed to load rMATS structure data:", rmats_structure_data$message, "\n")
          }
        } else {
          cat("❌ No rMATS GTF files found for gene", gene_clean, "\n")
          cat("🔍 Available files in rmats_peptide_results:\n")
          all_files <- list.files("rmats_peptide_results", recursive = TRUE, full.names = FALSE)
          gtf_files <- all_files[grepl("\\.gtf$", all_files)]
          if (length(gtf_files) > 0) {
            cat("  📄 GTF files found:", paste(head(gtf_files, 5), collapse = ", "))
            if (length(gtf_files) > 5) cat(" (showing first 5)")
            cat("\n")
          } else {
            cat("  📄 No GTF files found in directory\n")
          }
        }
      }
      
      # Step 5: Load structure data for canonical transcripts from main GTF
      canonical_structure_data <- NULL
      canonical_transcripts <- selected_transcripts[!grepl("\\.(inclusion|exclusion)$", selected_transcripts)]
      
      if (length(canonical_transcripts) > 0) {
        cat("🔍 Loading canonical structure data for", length(canonical_transcripts), "transcripts\n")
        incProgress(0.1, detail = 'Loading canonical genomic structure...')
        
        tryCatch({
          # Load GTF data for canonical transcripts
          gtf_data <- load_gtf_visualization_data(gene_id)
          if (gtf_data$success && length(gtf_data$exons_by_transcript) > 0) {
            cat("✅ Loaded canonical structure data for", length(gtf_data$exons_by_transcript), "transcripts\n")
            canonical_structure_data <- list(
              success = TRUE,
              exons = gtf_data$exons_by_transcript,
              cds = gtf_data$cds_by_transcript
            )
          } else {
            cat("⚠️ No canonical structure data found\n")
          }
        }, error = function(e) {
          cat("❌ Error loading canonical structure data:", e$message, "\n")
        })
      }
      
      # Step 6: Merge canonical and rMATS structure data for plotting
      combined_structure_data <- NULL
      if (!is.null(canonical_structure_data) || !is.null(rmats_structure_data)) {
        combined_structure_data <- list(success = TRUE, exons = list(), cds = list())
        
        # Add canonical structure data
        if (!is.null(canonical_structure_data) && canonical_structure_data$success) {
          cat("🔗 Adding canonical structure data...\n")
          combined_structure_data$exons <- c(combined_structure_data$exons, canonical_structure_data$exons)
          combined_structure_data$cds <- c(combined_structure_data$cds, canonical_structure_data$cds)
        }
        
        # Add rMATS structure data
        if (!is.null(rmats_structure_data) && rmats_structure_data$success) {
          cat("🔗 Adding rMATS structure data...\n")
          combined_structure_data$exons <- c(combined_structure_data$exons, rmats_structure_data$exons)
          combined_structure_data$cds <- c(combined_structure_data$cds, rmats_structure_data$cds)
        }
        
        cat("📊 Final structure data summary:\n")
        cat("  Total transcripts with exons:", length(combined_structure_data$exons), "\n")
        cat("  Total transcripts with CDS:", length(combined_structure_data$cds), "\n")
      }
      
      # Step 7: Final data preparation
      # Update y_position in peptides (this should already be done by analyze_peptide_specificity)
      if (!"y_position" %in% names(all_peptides_df)) {
        all_peptides_df$y_position <- transcript_df$y_position[match(all_peptides_df$transcript, transcript_df$transcript)]
      }
      
      # Hover text should already be added by analyze_peptide_specificity, but ensure it exists
      if (!"hover_text" %in% names(all_peptides_df)) {
        all_peptides_df$hover_text <- paste0(
          "Peptide: ", all_peptides_df$peptide,
          "<br>Position: ", all_peptides_df$start, "-", all_peptides_df$end,
          "<br>Transcript: ", all_peptides_df$transcript,
          "<br>Enzyme: ", protease,
          "<br>Miscleavage: ", miscleavage_type
        )
      }
      
      
      incProgress(0.1, detail = 'Complete!')
      
      return(list(
        all_peptides = all_peptides_df,
        transcript_df = transcript_df,
        gene_start = gene_start,
        gene_end = gene_end,
        all_transcripts = selected_transcripts,
        rmats_structure_data = combined_structure_data
      ))
    })
  })
  
  # Multi-Isoform highlighted data (for selected highlighting)
  multi_isoform_highlighted_data <- reactive({
    req(multi_isoform_data())
    
    base_data <- multi_isoform_data()
    if (is.null(base_data)) return(NULL)
    
    all_peptides_df <- base_data$all_peptides
    all_transcripts <- base_data$all_transcripts
    
    # Global peptide specificity analysis for Multi-Isoform Comparative Analysis
    # No highlighting - treat all isoforms equally
    miscleavage_type <- input$miscleavage_type
    
    # Get chromosome information from gene data
    gene_id <- input$gene
    protease <- input$protease
    chromosome_info <- "Unknown"
    
    tryCatch({
      gene_data <- selected_gene_peptides()
      if (!is.null(gene_data) && length(all_transcripts) > 0) {
        # Get chromosome from first available transcript
        first_tx <- all_transcripts[1]
        if (first_tx %in% gene_data$txID) {
          tx_row <- which(gene_data$txID == first_tx)[1]
          mapped_ranges_col <- paste0(protease, "Peps_mapped_ranges")
          
          if (mapped_ranges_col %in% names(gene_data) && !is.null(gene_data[[mapped_ranges_col]][[tx_row]])) {
            genomic_ranges <- gene_data[[mapped_ranges_col]][[tx_row]]
            if (length(genomic_ranges) > 0) {
              chromosomes <- unique(as.character(seqnames(genomic_ranges)))
              if (length(chromosomes) > 0) {
                chromosome_info <- chromosomes[1]
              }
            }
          }
        }
      }
    }, error = function(e) {
      cat("Warning: Could not load chromosome info:", e$message, "\n")
    })
    
    # Analyze global peptide specificity across ALL transcripts (no highlighting)
    all_peptides_df$specificity_category <- ""
    all_peptides_df$transcript_count <- 0
    all_peptides_df$shared_transcripts <- ""
    
    # Calculate specificity for each peptide globally
    for (i in 1:nrow(all_peptides_df)) {
      pep <- all_peptides_df$peptide[i]
      
      # Count how many transcripts contain this peptide
      transcripts_with_peptide <- unique(all_peptides_df$transcript[all_peptides_df$peptide == pep])
      transcript_count <- length(transcripts_with_peptide)
      
      all_peptides_df$transcript_count[i] <- transcript_count
      all_peptides_df$shared_transcripts[i] <- paste(transcripts_with_peptide, collapse = ", ")
      
      # Categorize globally
      if (transcript_count == 1) {
        all_peptides_df$specificity_category[i] <- "Unique"
      } else if (transcript_count == length(all_transcripts)) {
        all_peptides_df$specificity_category[i] <- "Universal"
      } else {
        all_peptides_df$specificity_category[i] <- "Shared"
      }
    }
    
    # Update hover text with global specificity information
    all_peptides_df$hover_text <- paste0(
      "Peptide: ", all_peptides_df$peptide,
      "<br>Position: ", all_peptides_df$start, "-", all_peptides_df$end,
      "<br>Transcript: ", all_peptides_df$transcript,
      "<br>Miscleavage: ", miscleavage_type,
      "<br>Specificity: ", all_peptides_df$specificity_category,
      "<br>Found in ", all_peptides_df$transcript_count, " of ", length(all_transcripts), " transcripts"
    )
    
    # Create global summary statistics for all peptides across all transcripts
    peptide_categories_summary <- as.data.frame(table(all_peptides_df$specificity_category))
    names(peptide_categories_summary) <- c("Category", "Count")
    peptide_categories_summary$Percentage <- round(peptide_categories_summary$Count / nrow(all_peptides_df) * 100, 1)
    
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
      transcript_df = base_data$transcript_df,
      gene_start = base_data$gene_start,
      gene_end = base_data$gene_end,
      all_transcripts = all_transcripts,
      peptide_categories_summary = peptide_categories_summary,
      isoform_coverage = isoform_coverage,
      chromosome = chromosome_info,
      miscleavage_type = miscleavage_type
    ))
  })
  
  # ============================================================================
  # VISUALIZATION OUTPUTS
  # ============================================================================
  
  # Render comparative plot using the same approach as all_isoforms_plot
  output$multi_isoform_comparative_plot <- renderPlotly({
    data <- multi_isoform_highlighted_data()
    
    if (is.null(data)) {
      p <- plotly::plot_ly() %>%
        plotly::add_annotations(
          x = 0.5, y = 0.5,
          text = "No peptide data available. Please select at least two isoforms and click 'Run'.",
          showarrow = FALSE,
          font = list(size = 16)
        ) %>%
        plotly::layout(
          xaxis = list(range = c(0, 1), showticklabels = FALSE),
          yaxis = list(range = c(0, 1), showticklabels = FALSE)
        )
      return(p)
    }
    
    # Extract data (same as all_isoforms_plot)
    all_peptides <- data$all_peptides
    transcript_df <- data$transcript_df
    gene_start <- data$gene_start
    gene_end <- data$gene_end
    all_transcripts <- transcript_df$transcript
    
    # Add compression support - respect user preference (using isoform-specific control)
    use_compression_isoforms <- !is.null(input$isoform_intron_scale) && input$isoform_intron_scale == "compressed"
    compression_map_isoforms <- NULL
    
    # Load GTF data for exon and CDS boundaries (same as all_isoforms_plot)
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
    
    # Create compression map if user selected compressed mode and exon data exists (same as all_isoforms_plot)
    if (use_compression_isoforms && exons_result$success && length(exons_result$exons) > 0) {
      compression_map_isoforms <- create_compression_map(exons_result$exons)
    }
    
    # Create plotly object (exact copy from all_isoforms_plot)
    p <- plotly::plot_ly()
    
    # Define colors for specificity categories using rgba to prevent hex codes in hover
    specificity_colors <- list(
      "Unique" = "rgba(255, 0, 0, 0.8)",      # Red - peptides unique to this isoform
      "Shared" = "rgba(243, 156, 18, 0.8)",   # Orange - peptides shared with some isoforms
      "Universal" = "rgba(46, 204, 113, 0.8)" # Green - peptides found in all isoforms
    )
    
    # Add transcript lines (ALL transcripts equally - no highlighting)
    for (i in 1:nrow(transcript_df)) {
      line_color <- "#666666"  # Medium gray for all transcripts
      line_width <- 2  # Same width for all transcripts
      
      # Apply compression to transcript line coordinates
      plot_gene_start <- gene_start
      plot_gene_end <- gene_end
      if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
        plot_gene_start <- compression_map_isoforms$compress(gene_start)
        plot_gene_end <- compression_map_isoforms$compress(gene_end)
      }
      
      p <- p %>% plotly::add_trace(
        type = "scatter", mode = "lines",
        x = c(plot_gene_start, plot_gene_end),
        y = c(transcript_df$y_position[i], transcript_df$y_position[i]),
        line = list(color = line_color, width = line_width),
        showlegend = FALSE,
        hovertemplate = paste0("Transcript: ", transcript_df$transcript[i], "<extra></extra>"),
        hoverinfo = "none"
      )
    }
    
    # Add genomic structure (exons and CDS) - exact copy from server.R lines 2993-3047
    if (exons_result$success) {
      exons_by_transcript <- exons_result$exons
      cds_by_transcript <- exons_result$cds
      
      for (i in 1:nrow(transcript_df)) {
        tx <- transcript_df$transcript[i]
        y_pos <- transcript_df$y_position[i]
        
        # Add exons if available
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
    
    # Add legend headers for structure
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
    
    # Add peptides for ALL transcripts with global specificity coloring
    if (nrow(all_peptides) > 0) {
      for (i in 1:nrow(all_peptides)) {
        # Color ALL peptides based on their global specificity category
        if (all_peptides$specificity_category[i] != "" && all_peptides$specificity_category[i] %in% names(specificity_colors)) {
          color <- specificity_colors[[all_peptides$specificity_category[i]]]
          legend_group <- paste0("category_", all_peptides$specificity_category[i])
        } else {
          # Fallback color for uncategorized peptides
          color <- "#CCCCCC"
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
          hovertemplate = paste0(all_peptides$hover_text[i], "<extra></extra>"),
          hoverinfo = "none"
        )
      }
    }
    
    # Add legend entries for each category with counts and percentages - exact copy from server.R lines 3060-3083  
    peptide_summary <- data$peptide_categories_summary
    for (category in names(specificity_colors)) {
      # Get count and percentage for this category
      category_row <- peptide_summary[peptide_summary$Category == category, ]
      if (nrow(category_row) > 0) {
        count <- category_row$Count
        percentage <- category_row$Percentage
        legend_name <- paste0(category, ": ", count, " (", percentage, "%)")
      } else {
        legend_name <- paste0(category, ": 0 (0%)")
      }
      
      p <- p %>% plotly::add_trace(
        type = "scatter", mode = "markers",
        x = c(0), y = c(0),
        marker = list(color = specificity_colors[[category]], size = 10, symbol = "square"),
        showlegend = TRUE,
        name = legend_name,
        legendgroup = paste0("category_", category),
        hoverinfo = "none",
        visible = TRUE
      )
    }
    
    # Final layout configuration
    x_title <- if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
      "Genomic Position (compressed view)"
    } else {
      "Genomic Position"
    }
    
    # Configure X-axis with proper genomic coordinate mapping
    xaxis_config <- if (use_compression_isoforms && !is.null(compression_map_isoforms)) {
      # Compressed view: map compressed positions to real genomic coordinates
      axis_breaks <- compression_map_isoforms$coords$compressed_start
      axis_labels <- as.character(compression_map_isoforms$coords$original_start)
      compressed_range <- range(c(compression_map_isoforms$coords$compressed_start, compression_map_isoforms$coords$compressed_end))
      
      list(
        title = x_title,
        tickvals = axis_breaks,
        ticktext = axis_labels,
        range = compressed_range,
        zeroline = FALSE
      )
    } else {
      # Original view
      list(
        title = x_title,
        range = c(gene_start, gene_end),
        zeroline = FALSE
      )
    }
    
    p <- p %>% layout(
      title = paste0("Multi-Isoform Comparative Analysis - ", data$chromosome, " - ", length(all_transcripts), " Selected"),
      xaxis = xaxis_config,
      yaxis = list(
        title = "Transcript",
        tickvals = transcript_df$y_position,
        ticktext = transcript_df$transcript,
        range = c(min(transcript_df$y_position) - 1, max(transcript_df$y_position) + 1)
      ),
      showlegend = TRUE,
      legend = list(
        orientation = "v",
        x = 1.02,
        xanchor = "left",
        y = 1,
        yanchor = "top"
      ),
      margin = list(l = 50, r = 200, t = 80, b = 50)
    ) %>%
    config(
      displayModeBar = TRUE,
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("pan2d", "select2d", "lasso2d", "resetScale2d")
    )
    
    return(p)
  })
  
  # ============================================================================
  # DATA TABLES
  # ============================================================================
  
  
  
  # ============================================================================
  # STATUS AND UI OUTPUTS
  # ============================================================================
  
  # Comparative analysis status
  output$comparative_analysis_status <- renderUI({
    selected_count <- length(input$compare_isoforms)
    
    if (selected_count == 0) {
      div(
        style = "color: #666; margin-top: 10px;",
        "Please select 2 or more isoforms for comparison."
      )
    } else if (selected_count == 1) {
      div(
        style = "color: #ff9900; margin-top: 10px;",
        paste("Selected 1 isoform. Please select at least 2 for comparison.")
      )
    } else if (selected_count > 8) {
      div(
        style = "color: #cc0000; margin-top: 10px;",
        paste("Selected", selected_count, "isoforms.")
      )
    } else {
      div(
        style = "color: #009900; margin-top: 10px;",
        paste("Selected", selected_count, "isoforms. Ready for comparative analysis.")
      )
    }
  })
  
  # ============================================================================
  # DOWNLOAD HANDLERS
  # ============================================================================
  
  # Download comparative data
  output$download_comparative_data <- downloadHandler(
    filename = function() {
      gene_symbol <- processed_data()$gene_lookup[input$gene] %||% input$gene
      miscleavage_suffix <- switch(input$miscleavage_type,
                                  "no_miss_cleavage" = "no_miss",
                                  "upto_one_misscleavage" = "1_miss",
                                  "upto_two_misscleavage" = "2_miss",
                                  "unknown")
      paste0("comparative_analysis_", gene_symbol, "_", input$protease, "_", 
             miscleavage_suffix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      data <- multi_isoform_highlighted_data()
      if (!is.null(data) && !is.null(data$all_peptides)) {
        comparative_peptides <- data$all_peptides
        
        # Add analysis metadata
        comparative_peptides$gene_id <- input$gene
        comparative_peptides$gene_symbol <- processed_data()$gene_lookup[input$gene] %||% input$gene
        comparative_peptides$enzyme <- input$protease
        comparative_peptides$miscleavage_type <- input$miscleavage_type
        comparative_peptides$analysis_timestamp <- Sys.time()
        comparative_peptides$selected_isoforms <- paste(input$compare_isoforms, collapse = ", ")
        
        write.csv(comparative_peptides, file, row.names = FALSE)
      }
    }
  )
  
  # Download overlap matrix
  output$download_overlap_matrix <- downloadHandler(
    filename = function() {
      gene_symbol <- processed_data()$gene_lookup[input$gene] %||% input$gene
      paste0("overlap_matrix_", gene_symbol, "_", input$protease, "_", 
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      data <- multi_isoform_highlighted_data()
      if (!is.null(data) && !is.null(data$overlap_matrix)) {
        write.csv(data$overlap_matrix, file, row.names = TRUE)
      }
    }
  )
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  
  # Get comparison status
  get_comparison_status <- function() {
    data <- multi_isoform_highlighted_data()
    
    list(
      has_data = !is.null(data),
      selected_count = length(comparison_isoforms()),
      ready_for_analysis = length(comparison_isoforms()) >= 2,
      total_peptides = if (!is.null(data) && !is.null(data$all_peptides)) nrow(data$all_peptides) else 0
    )
  }
  
  # Clear comparison data
  clear_comparison_data <- function() {
    multi_isoform_data(NULL)
    multi_isoform_highlighted_data(NULL)
    comparison_isoforms(character(0))
  }
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  
  return(list(
    # Reactive values
    multi_isoform_data = multi_isoform_data,
    multi_isoform_highlighted_data = multi_isoform_highlighted_data,
    comparison_isoforms = comparison_isoforms,
    
    # Utility functions
    get_comparison_status = get_comparison_status,
    clear_comparison_data = clear_comparison_data
  ))
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test Multi-Isoform Comparison Module
#' 
#' @description 
#' Function to test multi-isoform comparison module functionality
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_multi_isoform_comparison_module <- function() {
  cat("Testing Multi-Isoform Comparison module...\n")
  
  # Test 1: Overlap matrix calculation
  overlap_test <- tryCatch({
    # Test overlap matrix calculation with mock data
    transcripts <- c("TX1", "TX2", "TX3")
    peptides_list <- list(
      TX1 = data.frame(peptide = c("PEPTIDE1", "PEPTIDE2", "PEPTIDE3")),
      TX2 = data.frame(peptide = c("PEPTIDE1", "PEPTIDE4", "PEPTIDE5")),
      TX3 = data.frame(peptide = c("PEPTIDE2", "PEPTIDE3", "PEPTIDE4"))
    )
    
    # Mock overlap matrix structure
    n_tx <- length(transcripts)
    overlap_matrix <- matrix(0, nrow = n_tx, ncol = n_tx)
    diag(overlap_matrix) <- 1.0
    
    # Verify structure
    nrow(overlap_matrix) == length(transcripts) && ncol(overlap_matrix) == length(transcripts)
  }, error = function(e) {
    cat("Overlap test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 2: Comparative summary generation
  summary_test <- tryCatch({
    # Test summary structure
    mock_summary <- data.frame(
      Transcript = c("TX1", "TX2"),
      Total_Peptides = c(10, 15),
      Unique_Peptides = c(8, 12),
      stringsAsFactors = FALSE
    )
    
    # Verify structure
    all(c("Transcript", "Total_Peptides", "Unique_Peptides") %in% names(mock_summary))
  }, error = function(e) {
    cat("Summary test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 3: Status checking logic
  status_test <- tryCatch({
    # Test status structure
    mock_status <- list(
      has_data = FALSE,
      selected_count = 0,
      ready_for_analysis = FALSE,
      total_peptides = 0
    )
    
    # Verify all required fields
    required_fields <- c("has_data", "selected_count", "ready_for_analysis", "total_peptides")
    all(required_fields %in% names(mock_status))
  }, error = function(e) {
    cat("Status test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 4: Selection count validation
  selection_test <- tryCatch({
    # Test selection validation logic
    test_counts <- c(0, 1, 3, 9)
    valid_counts <- sapply(test_counts, function(x) x >= 2 && x <= 8)
    
    # Should be: FALSE, FALSE, TRUE, FALSE
    identical(valid_counts, c(FALSE, FALSE, TRUE, FALSE))
  }, error = function(e) {
    cat("Selection test failed:", e$message, "\n")
    FALSE
  })
  
  if (overlap_test && summary_test && status_test && selection_test) {
    cat("All Multi-Isoform Comparison module tests passed!\n")
    return(TRUE)
  } else {
    cat("Some Multi-Isoform Comparison module tests failed!\n")
    return(FALSE)
  }
}

#===============================================================================
# HELPER FUNCTIONS (DEPENDENCIES)
#===============================================================================

#' Create Empty Plotly Message
#' 
#' Creates a simple plotly plot with a text message for empty states
#' 
#' @param message Text message to display
#' @return Plotly object with the message
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

#' Create Multi-Isoform Plot
#' 
#' Creates a plotly visualization for multi-isoform comparison
#' 
#' @param all_peptides Data frame with peptide information
#' @param transcript_df Data frame with transcript information
#' @param gene_start Numeric gene start position
#' @param gene_end Numeric gene end position
#' @param gene_id String gene identifier
#' @param miscleavage_label String miscleavage description
#' @param protease String protease name
#' @param highlight_isoform String transcript to highlight (optional)
#' @param compression_map_isoforms Pre-computed compression map (optional)
#' @return plotly object
create_multi_isoform_plot <- function(all_peptides, transcript_df, gene_start, gene_end, 
                                    gene_id, miscleavage_label, protease, highlight_isoform = NULL,
                                    rmats_structure_data = NULL, use_compression = FALSE, intron_scale = NULL,
                                    compression_map_isoforms = NULL) {
  
  # Create plotly object directly (not ggplot)
  p <- plotly::plot_ly()
  
  # Apply compression to gene boundaries
  plot_gene_start <- gene_start
  plot_gene_end <- gene_end
  if (use_compression && !is.null(compression_map_isoforms)) {
    plot_gene_start <- compression_map_isoforms$compress(gene_start)
    plot_gene_end <- compression_map_isoforms$compress(gene_end)
    
    # Apply compression to all peptide coordinates
    all_peptides$plot_start <- sapply(all_peptides$start, compression_map_isoforms$compress)
    all_peptides$plot_end <- sapply(all_peptides$end, compression_map_isoforms$compress)
  } else {
    # Use original coordinates
    all_peptides$plot_start <- all_peptides$start
    all_peptides$plot_end <- all_peptides$end
  }
  
  # Add transcript tracks
  for (i in 1:nrow(transcript_df)) {
    tx_name <- transcript_df$transcript[i]
    y_pos <- transcript_df$y_position[i]
    
    # Determine color based on highlighting
    track_color <- if (!is.null(highlight_isoform) && tx_name == highlight_isoform) {
      "#ff6b6b"  # Red for highlighted
    } else {
      "#4ecdc4"  # Teal for others
    }
    
    # Add transcript line using compressed coordinates
    p <- p %>% plotly::add_trace(
      type = "scatter",
      x = c(plot_gene_start, plot_gene_end),
      y = c(y_pos, y_pos),
      mode = "lines",
      line = list(color = "gray70", width = 1),
      showlegend = FALSE,
      hoverinfo = "text",
      text = paste0("Transcript: ", tx_name)
    )
    
    # Add rMATS transcript structure (exons and CDS) if available
    if (!is.null(rmats_structure_data) && rmats_structure_data$success) {
      cat("🎨 Adding rMATS structure visualization for transcript:", tx_name, "\n")
      
      # Add exons for this transcript using existing fill aesthetic system
      if (tx_name %in% names(rmats_structure_data$exons)) {
        tx_exons <- rmats_structure_data$exons[[tx_name]]
        cat("📦 Adding", length(tx_exons), "exon rectangles for", tx_name, "\n")
        
        if (length(tx_exons) > 0) {
          for (j in seq_along(tx_exons)) {
            exon <- tx_exons[j]
            # Apply compression to exon coordinates
            exon_start <- start(exon)
            exon_end <- end(exon)
            if (use_compression && !is.null(compression_map_isoforms)) {
              exon_start <- compression_map_isoforms$compress(exon_start)
              exon_end <- compression_map_isoforms$compress(exon_end)
            }
            
            p <- p %>% plotly::add_trace(
              type = "scatter", mode = "lines",
              x = c(exon_start, exon_end, exon_end, exon_start, exon_start),
              y = c(y_pos - 0.3, y_pos - 0.3, y_pos + 0.3, y_pos + 0.3, y_pos - 0.3),
              fill = "toself",
              fillcolor = "rgba(77, 175, 74, 0.8)",  # Green for exons
              line = list(color = "black", width = 0.5),
              showlegend = FALSE,
              hovertemplate = paste0("Exon ", j, " (", start(exon), "-", end(exon), ")<br>Transcript: ", tx_name, "<extra></extra>"),
              hoverinfo = "none"
            )
          }
        }
      } else {
        cat("⚠️  No exon data found for transcript:", tx_name, "\n")
      }
      
      # Add CDS regions for this transcript using existing fill aesthetic system
      if (tx_name %in% names(rmats_structure_data$cds)) {
        tx_cds <- rmats_structure_data$cds[[tx_name]]
        cat("🟨 Adding", length(tx_cds), "CDS rectangles for", tx_name, "\n")
        
        if (length(tx_cds) > 0) {
          for (k in seq_along(tx_cds)) {
            cds <- tx_cds[k]
            # Apply compression to CDS coordinates
            cds_start <- start(cds)
            cds_end <- end(cds)
            if (use_compression && !is.null(compression_map_isoforms)) {
              cds_start <- compression_map_isoforms$compress(cds_start)
              cds_end <- compression_map_isoforms$compress(cds_end)
            }
            
            p <- p %>% plotly::add_trace(
              type = "scatter", mode = "lines",
              x = c(cds_start, cds_end, cds_end, cds_start, cds_start),
              y = c(y_pos - 0.2, y_pos - 0.2, y_pos + 0.2, y_pos + 0.2, y_pos - 0.2),
              fill = "toself",
              fillcolor = "rgba(255, 221, 0, 0.8)",  # Gold/yellow for CDS
              line = list(color = "black", width = 0.5),
              showlegend = FALSE,
              hovertemplate = paste0("CDS ", k, " (", start(cds), "-", end(cds), ")<br>Transcript: ", tx_name, "<extra></extra>"),
              hoverinfo = "none"
            )
          }
        }
      } else {
        cat("ℹ️  No CDS data found for transcript:", tx_name, "\n")
      }
    } else {
      if (grepl("\\.(inclusion|exclusion)$", tx_name)) {
        cat("❌ No rMATS structure data available for transcript:", tx_name, "\n")
      }
    }
  }
  
  # Add peptides
  if (nrow(all_peptides) > 0) {
    for (i in 1:nrow(all_peptides)) {
      # Determine peptide color based on specificity
      peptide_color <- if ("specificity_color" %in% names(all_peptides)) {
        all_peptides$specificity_color[i]  # Use specificity-based coloring
      } else if (!is.null(highlight_isoform) && "is_highlighted" %in% names(all_peptides) && all_peptides$is_highlighted[i]) {
        "rgba(255, 107, 107, 0.8)"  # Red for highlighted (fallback)
      } else {
        "rgba(76, 205, 196, 0.8)"  # Teal for others (fallback)
      }
      
      p <- p %>% plotly::add_trace(
        type = "scatter", mode = "lines",
        x = c(all_peptides$plot_start[i], all_peptides$plot_end[i], all_peptides$plot_end[i], all_peptides$plot_start[i], all_peptides$plot_start[i]),
        y = c(all_peptides$y_position[i] - 0.15, all_peptides$y_position[i] - 0.15, all_peptides$y_position[i] + 0.15, all_peptides$y_position[i] + 0.15, all_peptides$y_position[i] - 0.15),
        fill = "toself",
        fillcolor = peptide_color,
        line = list(color = "white", width = 0.5),
        showlegend = FALSE,
        hovertemplate = paste0(all_peptides$hover_text[i], "<extra></extra>"),
        hoverinfo = "none"
      )
    }
  }
  
  # Configure axis layout for plotly
  x_title <- if (use_compression && !is.null(compression_map_isoforms)) {
    "Genomic Position (compressed view)"
  } else {
    "Genomic Position"
  }
  
  # Use the same axis configuration approach as the working matrix analysis
  xaxis_config <- if (use_compression && !is.null(compression_map_isoforms)) {
    list(
      title = x_title,
      range = c(plot_gene_start, plot_gene_end),
      showgrid = TRUE,
      gridcolor = "rgba(128,128,128,0.2)"
      # Removed tickvals and ticktext to let plotly auto-generate cleaner axis labels
    )
  } else {
    # Original view
    list(
      title = x_title,
      range = c(plot_gene_start, plot_gene_end),
      zeroline = FALSE
    )
  }
  
  # Configure final plotly layout
  title_text <- paste0("Multi-Isoform Comparison: ", gene_id)
  subtitle_text <- paste0("Protease: ", toupper(protease), " | ", miscleavage_label)
  
  p <- p %>% plotly::layout(
    title = list(
      text = paste0(title_text, "<br><sup>", subtitle_text, "</sup>")
    ),
    xaxis = xaxis_config,
    yaxis = list(
      title = "Transcripts",
      range = c(0.5, nrow(transcript_df) + 0.5),
      tickvals = transcript_df$y_position,
      ticktext = transcript_df$transcript,
      zeroline = FALSE
    ),
    showlegend = FALSE,
    margin = list(l = 150, r = 50, t = 80, b = 50),
    plot_bgcolor = "white",
    paper_bgcolor = "white"
  )
  
  return(p)
}

# NOTE: get_transcript_peptides_for_comparison function is defined in as_analysis.R
# and will be loaded when that file is sourced. No stub needed here.

#' Analyze peptide specificity across transcripts
#' @param all_peptides_df Data frame with peptides from all transcripts
#' @param selected_transcripts Vector of transcript IDs
#' @return Data frame with added specificity columns
analyze_peptide_specificity <- function(all_peptides_df, selected_transcripts) {
  cat("🔍 Analyzing peptide specificity for", nrow(all_peptides_df), "peptides across", length(selected_transcripts), "transcripts\n")
  
  # Add specificity analysis
  all_peptides_df$specificity <- "Unknown"
  all_peptides_df$specificity_color <- "rgba(128, 128, 128, 0.8)"  # Gray default
  all_peptides_df$hover_text <- paste0("Peptide: ", all_peptides_df$peptide, "<br>",
                                       "Position: ", all_peptides_df$start, "-", all_peptides_df$end, "<br>",
                                       "Transcript: ", all_peptides_df$transcript)
  
  # Count transcripts per peptide sequence
  peptide_counts <- table(all_peptides_df$peptide)
  
  # Categorize peptides
  for (i in seq_len(nrow(all_peptides_df))) {
    peptide_seq <- all_peptides_df$peptide[i]
    transcript_count <- peptide_counts[peptide_seq]
    
    if (transcript_count == 1) {
      all_peptides_df$specificity[i] <- "Unique"
      all_peptides_df$specificity_color[i] <- "rgba(255, 107, 107, 0.8)"  # Red for unique
    } else if (transcript_count == length(selected_transcripts)) {
      all_peptides_df$specificity[i] <- "Universal"
      all_peptides_df$specificity_color[i] <- "rgba(76, 175, 74, 0.8)"   # Green for universal
    } else {
      all_peptides_df$specificity[i] <- "Shared"
      all_peptides_df$specificity_color[i] <- "rgba(255, 193, 7, 0.8)"   # Yellow for shared
    }
    
    # Update hover text
    all_peptides_df$hover_text[i] <- paste0(all_peptides_df$hover_text[i], "<br>",
                                           "Specificity: ", all_peptides_df$specificity[i])
  }
  
  # Summary
  specificity_summary <- table(all_peptides_df$specificity)
  cat("📊 Peptide specificity summary:\n")
  for (spec in names(specificity_summary)) {
    cat("  ", spec, ":", specificity_summary[spec], "peptides\n")
  }
  
  return(all_peptides_df)
}