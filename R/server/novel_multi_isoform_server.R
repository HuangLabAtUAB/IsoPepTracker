#============================================================================
# NOVEL MULTI-ISOFORM COMPARATIVE ANALYSIS
#============================================================================

# Extract gene ID from novel data for cache lookup (FIXED VERSION)
extract_gene_id_from_novel_data <- function(novel_data) {
  if (is.null(novel_data) || nrow(novel_data) == 0) {
    return(NULL)
  }
  
  cat("DEBUG: Attempting to extract Ensembl gene ID from novel data...\n")
  cat("DEBUG: Available columns:", paste(names(novel_data), collapse = ", "), "\n")
  
  # Strategy 1: Look for direct gene ID column
  if ("geneID" %in% names(novel_data)) {
    gene_ids <- unique(novel_data$geneID)
    cat("DEBUG: Found geneID column values:", paste(head(gene_ids, 3), collapse = ", "), "\n")
    
    for (gene_id in gene_ids) {
      if (!is.null(gene_id) && !is.na(gene_id) && nchar(as.character(gene_id)) > 0) {
        # Clean up the gene ID 
        clean_gene_id <- gsub('^"+"', '', as.character(gene_id))  # Remove quotes
        clean_gene_id <- gsub('"+$', '', clean_gene_id)
        
        # Check if this looks like an Ensembl gene ID
        if (grepl("^ENSG[0-9]+(\\.[0-9]+)?$", clean_gene_id)) {
          cat("DEBUG: Found valid Ensembl gene ID:", clean_gene_id, "\n")
          return(clean_gene_id)
        }
      }
    }
  }
  
  # Strategy 2: Smart search through GTF cache files using transcript IDs
  if ("txID" %in% names(novel_data) && dir.exists("data/gtf_cache")) {
    transcript_ids <- unique(novel_data$txID)
    cat("DEBUG: Found transcript IDs, searching GTF cache:", paste(head(transcript_ids, 3), collapse = ", "), "\n")
    
    # Get reference transcript IDs (non-novel)
    ref_transcript_ids <- transcript_ids[!grepl("^Novel_sequence_", transcript_ids)]
    if (length(ref_transcript_ids) > 0) {
      cat("DEBUG: Searching for reference transcripts:", paste(head(ref_transcript_ids, 2), collapse = ", "), "\n")
      
      # Get list of cached gene files to search
      cache_files <- list.files("data/gtf_cache", pattern = "\\.rds$", full.names = TRUE)
      
      # Search first 20 cache files for efficiency (common genes first)
      for (cache_file in head(cache_files, 20)) {
        tryCatch({
          gene_details <- readRDS(cache_file)
          if (!is.null(gene_details$transcript_ids)) {
            # Check if any of our reference transcript IDs match this gene
            if (any(ref_transcript_ids %in% gene_details$transcript_ids)) {
              gene_id <- gene_details$gene_id
              matching_txs <- intersect(ref_transcript_ids, gene_details$transcript_ids)
              cat("DEBUG: Found matching transcripts", paste(matching_txs, collapse = ", "), "in gene", gene_id, "\n")
              return(as.character(gene_id))
            }
          }
        }, error = function(e) {
          # Skip files that can't be read
        })
      }
    }
  }
  
  # Strategy 3: Extract from transcript IDs if they contain ENSG (unlikely but try)
  if ("txID" %in% names(novel_data)) {
    transcript_ids <- unique(novel_data$txID)
    for (tx_id in transcript_ids) {
      if (!grepl("^Novel_sequence_", tx_id)) {
        # Try to extract ENSG ID if embedded in transcript ID
        ensembl_match <- regmatches(tx_id, regexpr('ENSG[0-9]+\\.[0-9]+', tx_id))
        if (length(ensembl_match) > 0 && nchar(ensembl_match) > 0) {
          cat("DEBUG: Extracted embedded gene ID:", ensembl_match, "\n")
          return(ensembl_match)
        }
      }
    }
  }
  
  cat("WARNING: Could not extract gene ID from novel data\n")
  cat("DEBUG: Data summary - rows:", nrow(novel_data), "cols:", ncol(novel_data), "\n")
  return(NULL)
}

# Function to load reference transcript exons and CDS from reference GTF (ORIGINAL)
load_reference_transcript_exons <- function(ref_gtf_file, transcript_ids) {
  if (!file.exists(ref_gtf_file)) {
    return(list(success = FALSE, message = "Reference GTF file not found"))
  }
  
  tryCatch({
    cat("Loading reference GTF:", ref_gtf_file, "\n")
    
    # Import exons from reference GTF
    exons <- rtracklayer::import(ref_gtf_file, format = "gtf", feature.type = "exon")
    
    # Import CDS features from reference GTF
    cds <- rtracklayer::import(ref_gtf_file, format = "gtf", feature.type = "CDS")
    
    cat("Found", length(exons), "exons and", length(cds), "CDS features in reference GTF\n")
    
    # Check for transcript IDs
    if (length(exons) > 0 && !"transcript_id" %in% names(mcols(exons))) {
      return(list(success = FALSE, message = "Reference GTF file does not contain transcript_id attribute for exons"))
    }
    
    if (length(cds) > 0 && !"transcript_id" %in% names(mcols(cds))) {
      message("WARNING: Reference GTF file does not contain transcript_id attribute for CDS features")
      cds <- GRanges() # Empty GRanges if no transcript_id
    }
    
    # Extract exons and CDS for each requested transcript
    result <- list(exons = list(), cds = list())
    found_transcripts <- c()
    
    for (tx_id in transcript_ids) {
      # Get exons for this transcript
      tx_exons <- exons[exons$transcript_id == tx_id]
      if (length(tx_exons) > 0) {
        result$exons[[tx_id]] <- tx_exons
        found_transcripts <- c(found_transcripts, tx_id)
        cat("Found", length(tx_exons), "exons for reference transcript", tx_id, "\n")
      }
      
      # Get CDS for this transcript
      if (length(cds) > 0) {
        tx_cds <- cds[cds$transcript_id == tx_id]
        if (length(tx_cds) > 0) {
          result$cds[[tx_id]] <- tx_cds
          cat("Found", length(tx_cds), "CDS regions for reference transcript", tx_id, "\n") 
        }
      }
    }
    
    if (length(found_transcripts) == 0) {
      return(list(success = FALSE, message = paste("No transcripts found in reference GTF for IDs:", paste(transcript_ids, collapse = ", "))))
    }
    
    cat("Successfully loaded reference data for", length(found_transcripts), "transcripts:", paste(found_transcripts, collapse = ", "), "\n")
    
    return(list(
      success = TRUE,
      exons = result$exons,
      cds = result$cds
    ))
  }, error = function(e) {
    return(list(success = FALSE, message = paste("Error processing reference GTF:", e$message)))
  })
}

# Function to load novel transcript exons and CDS from novel GTF (ORIGINAL)
load_novel_transcript_exons <- function(novel_gtf_file, transcript_ids) {
  if (!file.exists(novel_gtf_file)) {
    return(list(success = FALSE, message = "Novel GTF file not found"))
  }
  
  tryCatch({
    cat("=== NOVEL GTF LOADING ===\n")
    cat("Loading novel GTF:", novel_gtf_file, "\n")
    cat("Looking for transcript IDs:", paste(transcript_ids, collapse = ", "), "\n")
    
    # Import exons from novel GTF
    exons <- rtracklayer::import(novel_gtf_file, format = "gtf", feature.type = "exon")
    
    # Import CDS features from novel GTF
    cds <- rtracklayer::import(novel_gtf_file, format = "gtf", feature.type = "CDS")
    
    cat("Found", length(exons), "exons and", length(cds), "CDS features in novel GTF\n")
    
    # Show all available transcript IDs in the GTF
    if (length(exons) > 0 && "transcript_id" %in% names(mcols(exons))) {
      available_tx_ids <- unique(exons$transcript_id)
      cat("Available transcript IDs in novel GTF:", paste(available_tx_ids, collapse = ", "), "\n")
    }
    
    # Check for transcript IDs
    if (length(exons) > 0 && !"transcript_id" %in% names(mcols(exons))) {
      return(list(success = FALSE, message = "Novel GTF file does not contain transcript_id attribute for exons"))
    }
    
    if (length(cds) > 0 && !"transcript_id" %in% names(mcols(cds))) {
      message("WARNING: Novel GTF file does not contain transcript_id attribute for CDS features")
      cds <- GRanges() # Empty GRanges if no transcript_id
    }
    
    # Extract exons and CDS for each requested transcript
    result <- list(exons = list(), cds = list())
    found_transcripts <- c()
    
    for (tx_id in transcript_ids) {
      cat("Processing transcript:", tx_id, "\n")
      
      # Get exons for this transcript (try exact match first)
      tx_exons <- exons[exons$transcript_id == tx_id]
      
      # If no exact match found, try to find close matches (for versioning issues)
      if (length(tx_exons) == 0) {
        cat("No exact match found for", tx_id, ", trying pattern matching...\n")
        
        # Try different patterns
        possible_matches <- c()
        
        # Pattern 1: Try removing version numbers (.1, .2, etc.)
        base_id <- gsub("\\.\\d+$", "", tx_id)
        possible_matches <- c(possible_matches, available_tx_ids[grepl(paste0("^", base_id), available_tx_ids)])
        
        # Pattern 2: Try Novel_sequence_1 variants
        if (grepl("Novel_sequence", tx_id)) {
          possible_matches <- c(possible_matches, available_tx_ids[grepl("Novel_sequence", available_tx_ids)])
        }
        
        possible_matches <- unique(possible_matches)
        
        if (length(possible_matches) > 0) {
          matched_id <- possible_matches[1]  # Use first match
          tx_exons <- exons[exons$transcript_id == matched_id]
          cat("Using matched ID:", matched_id, "for", tx_id, "\n")
        }
      }
      
      if (length(tx_exons) > 0) {
        result$exons[[tx_id]] <- tx_exons
        found_transcripts <- c(found_transcripts, tx_id)
        cat("Found", length(tx_exons), "exons for novel transcript", tx_id, "\n")
      } else {
        cat("No exons found for transcript", tx_id, "\n")
      }
      
      # Get CDS for this transcript (using same logic)
      if (length(cds) > 0) {
        tx_cds <- cds[cds$transcript_id == tx_id]
        
        # Try pattern matching for CDS if exact match failed
        if (length(tx_cds) == 0 && length(possible_matches) > 0) {
          matched_id <- possible_matches[1]
          tx_cds <- cds[cds$transcript_id == matched_id]
        }
        
        if (length(tx_cds) > 0) {
          result$cds[[tx_id]] <- tx_cds
          cat("Found", length(tx_cds), "CDS regions for novel transcript", tx_id, "\n")
        }
      }
    }
    
    if (length(found_transcripts) == 0) {
      return(list(success = FALSE, message = paste("No transcripts found in novel GTF for IDs:", paste(transcript_ids, collapse = ", "))))
    }
    
    cat("Successfully loaded novel data for", length(found_transcripts), "transcripts:", paste(found_transcripts, collapse = ", "), "\n")
    cat("=========================\n")
    
    return(list(
      success = TRUE,
      exons = result$exons,
      cds = result$cds
    ))
  }, error = function(e) {
    return(list(success = FALSE, message = paste("Error processing novel GTF:", e$message)))
  })
}

# Update novel_compare_isoforms choices when novel data is available
observeEvent(novel_merged_data(), {
  merged_data <- novel_merged_data()
  
  if (!is.null(merged_data) && nrow(merged_data) > 0) {
    # Get all available transcripts (novel + reference)
    all_transcripts <- unique(merged_data$txID)
    
    # Separate novel and reference transcripts - Fix pattern to match actual novel transcript naming
    novel_transcripts <- all_transcripts[grepl("^Novel_sequence_", all_transcripts)]
    reference_transcripts <- all_transcripts[!grepl("^Novel_sequence_", all_transcripts)]
    
    # Create choices with labels
    all_choices <- character(0)
    
    # Add novel transcripts
    if (length(novel_transcripts) > 0) {
      novel_choices <- setNames(novel_transcripts, paste0(novel_transcripts, " (Novel)"))
      all_choices <- c(all_choices, novel_choices)
    }
    
    # Add reference transcripts
    if (length(reference_transcripts) > 0) {
      ref_choices <- setNames(reference_transcripts, paste0(reference_transcripts, " (Reference)"))
      all_choices <- c(all_choices, ref_choices)
    }
    
    # Update dropdown
    updateSelectizeInput(session, "novel_compare_isoforms", 
                        choices = all_choices,
                        selected = character(0))
    
    cat("✅ Novel Multi-Isoform: Added", length(novel_transcripts), "novel +", length(reference_transcripts), "reference transcripts\n")
  } else {
    updateSelectizeInput(session, "novel_compare_isoforms", 
                         choices = NULL, selected = character(0))
  }
})

# Enable/disable novel comparative analysis button based on selection count
observeEvent(input$novel_compare_isoforms, {
  selected_count <- length(input$novel_compare_isoforms)
  
  if (selected_count >= 2 && selected_count <= 8) {
    shinyjs::enable("run_novel_comparative_analysis")
    shinyjs::removeClass("run_novel_comparative_analysis", "btn-secondary")
    shinyjs::addClass("run_novel_comparative_analysis", "btn-warning")
  } else {
    shinyjs::disable("run_novel_comparative_analysis")
    shinyjs::removeClass("run_novel_comparative_analysis", "btn-warning") 
    shinyjs::addClass("run_novel_comparative_analysis", "btn-secondary")
  }
})

# Novel Multi-Isoform data loader
novel_multi_isoform_data <- reactive({
  req(input$run_novel_comparative_analysis > 0, input$novel_protease, input$novel_miscleavage_type, 
      novel_merged_data(), input$novel_compare_isoforms, 
      length(input$novel_compare_isoforms) >= 2)
  
  withProgress(message = 'Loading novel multi-isoform analysis...', value = 0, {
    protease <- input$novel_protease
    miscleavage_type <- input$novel_miscleavage_type
    selected_transcripts <- input$novel_compare_isoforms
    merged_data <- novel_merged_data()
    
    cat("Selected transcripts for novel multi-isoform analysis:", paste(selected_transcripts, collapse=", "), "\n")
    cat("Using enzyme:", protease, "and miscleavage type:", miscleavage_type, "\n")
    
    incProgress(0.2, detail = 'Loading peptides for selected transcripts...')
    
    # Create vis_data structure for compatibility with existing functions
    vis_data <- list(
      genes = c(),
      gene_symbols = c(),
      gene_lookup = c(),
      proteases = c("trp", "chymo", "aspn", "lysc", "lysn", "gluc"),
      original_peptides = merged_data
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
    
    incProgress(0.3, detail = 'Combining peptide data...')
    all_peptides_df <- data.table::rbindlist(all_peptides_list)
    
    # Calculate gene boundaries
    gene_start <- min(all_peptides_df$start) - 1000
    gene_end <- max(all_peptides_df$end) + 1000
    
    # Create transcript position mapping
    transcript_df <- data.frame(
      transcript = selected_transcripts,
      y_position = seq_along(selected_transcripts),
      stringsAsFactors = FALSE
    )
    
    # Add hover text for all peptides
    all_peptides_df$hover_text <- paste0(
      "Peptide: ", all_peptides_df$peptide,
      "<br>Position: ", all_peptides_df$start, "-", all_peptides_df$end,
      "<br>Transcript: ", all_peptides_df$transcript,
      "<br>Enzyme: ", protease,
      "<br>Miscleavage: ", miscleavage_type
    )
    
    incProgress(0.3, detail = 'Calculating specificity...')
    
    # Calculate specificity for ALL peptides across selected transcripts
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
    
    incProgress(0.2, detail = 'Completed')
    
    return(list(
      all_peptides = all_peptides_df,
      transcript_df = transcript_df,
      gene_start = gene_start,
      gene_end = gene_end,
      all_transcripts = selected_transcripts
    ))
  })
})

# Novel Multi-Isoform highlighted data
novel_multi_isoform_highlighted_data <- reactive({
  req(novel_multi_isoform_data(), input$novel_highlight_isoform)
  
  base_data <- novel_multi_isoform_data()
  highlight_isoform <- input$novel_highlight_isoform
  
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
    "<br>Enzyme: ", input$novel_protease,
    "<br>Miscleavage: ", input$novel_miscleavage_type
  )
  
  return(list(
    all_peptides = all_peptides_df,
    transcript_df = base_data$transcript_df,
    gene_start = base_data$gene_start,
    gene_end = base_data$gene_end,
    all_transcripts = base_data$all_transcripts
  ))
})

# Render novel comparative plot (RESTORED ORIGINAL WORKING VERSION)
output$novel_comparative_plot <- renderPlotly({
  req(novel_multi_isoform_highlighted_data())
  
  tryCatch({
    data <- novel_multi_isoform_highlighted_data()
    
    if (is.null(data) || is.null(data$all_peptides) || nrow(data$all_peptides) == 0) {
      p <- plotly::plot_ly() %>%
        plotly::add_annotations(
          x = 0.5, y = 0.5,
          text = "No comparative data available. Please select 2 or more isoforms and run analysis.",
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
    
    # Load novel GTF data for exon and CDS boundaries  
    novel_results <- novel_pipeline_results()
    novel_gtf_file <- NULL
    
    # Get the latest novel analysis GTF file
    if (!is.null(novel_results) && novel_results$success) {
      cat("=== NOVEL RESULTS ACCESS ===\n")
      cat("Novel results success:", novel_results$success, "\n")
      cat("Novel results output_dir:", novel_results$output_dir, "\n")
      cat("Novel results gtf_file:", novel_results$gtf_file, "\n")
      
      # Check if GTF file is in results
      if (!is.null(novel_results$gtf_file) && file.exists(novel_results$gtf_file)) {
        novel_gtf_file <- novel_results$gtf_file
        cat("Using GTF file from results\n")
      } else {
        # Look for GTF file in the output directory
        output_dir <- novel_results$output_dir
        potential_gtf <- file.path(output_dir, "novel_transcript_nt.transdecoder.genome.gtf")
        cat("Looking for GTF at:", potential_gtf, "\n")
        if (file.exists(potential_gtf)) {
          novel_gtf_file <- potential_gtf
          cat("Found GTF file in output directory\n")
        } else {
          cat("GTF file not found in output directory\n")
        }
      }
      cat("============================\n")
    } else {
      cat("Novel results not available or failed\n")
    }
    
    cat("Final novel GTF file path:", novel_gtf_file, "\n")
    cat("Novel GTF exists:", !is.null(novel_gtf_file) && file.exists(novel_gtf_file), "\n")
    
    # Separate novel and reference transcripts - Fix pattern to match actual novel transcript naming  
    novel_transcripts <- all_transcripts[grepl("^Novel_sequence_", all_transcripts)]
    reference_transcripts <- all_transcripts[!grepl("^Novel_sequence_", all_transcripts)]
    
    # Initialize combined exon/CDS data
    exons_by_transcript <- list()
    cds_by_transcript <- list()
    
    # Load novel GTF data for novel transcripts (if available)
    if (!is.null(novel_gtf_file) && file.exists(novel_gtf_file) && length(novel_transcripts) > 0) {
      cat("Loading novel GTF for transcripts:", paste(novel_transcripts, collapse = ", "), "\n")
      novel_exons_result <- load_novel_transcript_exons(novel_gtf_file, novel_transcripts)
      
      if (novel_exons_result$success && length(novel_exons_result$exons) > 0) {
        # Add novel exons and CDS to combined data
        exons_by_transcript <- c(exons_by_transcript, novel_exons_result$exons)
        cds_by_transcript <- c(cds_by_transcript, novel_exons_result$cds)
        cat("Loaded", length(novel_exons_result$exons), "novel transcript exon sets\n")
      } else {
        cat("Failed to load novel exons:", novel_exons_result$message, "\n")
      }
    } else {
      if (length(novel_transcripts) > 0) {
        cat("Novel GTF file not available for novel transcripts:", paste(novel_transcripts, collapse = ", "), "\n")
      }
    }
    
    # Load reference GTF data for reference transcripts (if available)
    if (length(reference_transcripts) > 0) {
      cat("Loading reference GTF for transcripts:", paste(reference_transcripts, collapse = ", "), "\n")
      
      # Use cached GTF system if available (FAST - same as rMATS approach)
      if (dir.exists("data/gtf_cache")) {
        # Extract gene ID from novel data to find GTF cache file
        gene_id <- extract_gene_id_from_novel_data(novel_merged_data())
        if (!is.null(gene_id)) {
          
          # Load GTF cache directly using same approach as rMATS system
          cache_file <- file.path("data/gtf_cache", paste0(gene_id, ".rds"))
          cat("Attempting to load reference GTF cache from:", cache_file, "\n")
          
          if (file.exists(cache_file)) {
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
                
                # Store chromosome information for later use
                chromosome <- if (!is.null(gene_details$chromosome)) gene_details$chromosome else "Unknown"
                
                cat("✅ Loaded", length(ref_exons), "reference transcript exon sets from GTF cache\n")
                cat("✅ Loaded", length(ref_cds), "reference transcript CDS sets from GTF cache\n")
                cat("✅ Chromosome information:", chromosome, "\n")
                
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
      
      # Fallback to long GTF parsing if cache not available or failed
      if (length(exons_by_transcript) == 0 && exists("gtf_file") && file.exists(gtf_file)) {
        cat("Falling back to long GTF parsing...\n")
        ref_exons_result <- load_reference_transcript_exons(gtf_file, reference_transcripts)
        
        if (ref_exons_result$success && length(ref_exons_result$exons) > 0) {
          # Add reference exons and CDS to combined data
          exons_by_transcript <- c(exons_by_transcript, ref_exons_result$exons)
          cds_by_transcript <- c(cds_by_transcript, ref_exons_result$cds)
          cat("Loaded", length(ref_exons_result$exons), "reference transcript exon sets from long GTF\n")
        } else {
          cat("Failed to load reference exons:", ref_exons_result$message, "\n")
        }
      }
    }
    
    # Initialize chromosome if not set by cache
    if (!exists("chromosome")) {
      chromosome <- "Unknown"
    }
    
    # Debug: show final combined exon/CDS data
    cat("=== GTF LOADING DEBUG ===\n")
    cat("All transcripts in visualization:", paste(all_transcripts, collapse = ", "), "\n")
    cat("Novel transcripts:", paste(novel_transcripts, collapse = ", "), "\n")
    cat("Reference transcripts:", paste(reference_transcripts, collapse = ", "), "\n")
    cat("Total exon sets loaded:", length(exons_by_transcript), "\n")
    cat("Total CDS sets loaded:", length(cds_by_transcript), "\n")
    cat("Exon transcript IDs:", paste(names(exons_by_transcript), collapse = ", "), "\n")
    cat("CDS transcript IDs:", paste(names(cds_by_transcript), collapse = ", "), "\n")
    
    # Debug: check if novel transcripts have exon/CDS data
    for (novel_tx in novel_transcripts) {
      has_exons <- !is.null(exons_by_transcript[[novel_tx]]) && length(exons_by_transcript[[novel_tx]]) > 0
      has_cds <- !is.null(cds_by_transcript[[novel_tx]]) && length(cds_by_transcript[[novel_tx]]) > 0
      cat("Novel transcript", novel_tx, "has exons:", has_exons, "\n")
      cat("Novel transcript", novel_tx, "has CDS:", has_cds, "\n")
      
      if (!has_exons) {
        cat("ERROR: Novel transcript", novel_tx, "missing exon data - will not show structure!\n")
        cat("ATTEMPTING FALLBACK: Creating mock exon data from peptide coordinates...\n")
        
        # Create mock exon data from peptide coordinates as fallback
        novel_peptides <- all_peptides_df[all_peptides_df$transcript == novel_tx, ]
        if (nrow(novel_peptides) > 0) {
          # Create simple exon covering the peptide range
          exon_start <- min(novel_peptides$start) - 100
          exon_end <- max(novel_peptides$end) + 100
          
          # Create mock exon GRanges object
          mock_exon <- GRanges(
            seqnames = "chr1",  # placeholder
            ranges = IRanges(start = exon_start, end = exon_end),
            strand = "*"
          )
          names(mock_exon) <- novel_tx
          
          # Add to exons_by_transcript
          exons_by_transcript[[novel_tx]] <- mock_exon
          cat("Created mock exon for", novel_tx, ":", exon_start, "-", exon_end, "\n")
        }
      }
    }
    
    # Debug: check if reference transcripts have exon/CDS data  
    for (ref_tx in reference_transcripts) {
      has_exons <- !is.null(exons_by_transcript[[ref_tx]]) && length(exons_by_transcript[[ref_tx]]) > 0
      has_cds <- !is.null(cds_by_transcript[[ref_tx]]) && length(cds_by_transcript[[ref_tx]]) > 0
      cat("Reference transcript", ref_tx, "has exons:", has_exons, "\n")
      cat("Reference transcript", ref_tx, "has CDS:", has_cds, "\n")
      
      if (!has_exons) {
        cat("ERROR: Reference transcript", ref_tx, "missing exon data - will not show structure!\n")
      }
    }
    cat("========================\n")
    
    # ============================================================================
    # COORDINATE COMPRESSION (DEFAULT FOR NOVEL VISUALIZATIONS)
    # ============================================================================
    
    compression_map <- NULL
    use_compression <- TRUE  # Always use compression for novel visualizations
    
    if (use_compression && length(exons_by_transcript) > 0) {
      # Create compression map from exon data
      cat("Creating coordinate compression map for novel visualization...\n")
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
    
    # ============================================================================
    # CREATE PLOT AFTER COMPRESSION IS SET UP
    # ============================================================================
    
    # Create plotly object
    p <- plotly::plot_ly()
    
    # Define colors for specificity categories (exact same as regular analysis)
    specificity_colors <- list(
      "Unique" = "#FF0000",      # Red - peptides unique to this isoform
      "Shared" = "#F39C12",      # Orange - peptides shared with some isoforms
      "Universal" = "#2ECC71"    # Green - peptides found in all isoforms
    )
    
    # Add transcript lines (ALL transcripts, no highlighting) - UPDATED FOR COMPRESSION
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
    
    # Add intron markers (triangle/zigzag lines between exons) - NOVEL FEATURE
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
    
    # Add exon and CDS boundaries if we have any data
    if (length(exons_by_transcript) > 0) {
      # Add exon blocks for each transcript (exactly like regular analysis)
      for (i in 1:nrow(transcript_df)) {
        tx <- transcript_df$transcript[i]
        y_pos <- transcript_df$y_position[i]
        
        # Add exons if available (transparent grey - same as regular analysis) - UPDATED FOR COMPRESSION
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
        
        # Add CDS overlay if available (yellow - same as regular analysis) - UPDATED FOR COMPRESSION
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
    
    # Add legend entries for each category (same as regular analysis)
    for (category in names(specificity_colors)) {
      p <- p %>% plotly::add_trace(
        type = "scatter", mode = "markers",
        x = c(0), y = c(0),
        marker = list(color = specificity_colors[[category]], size = 10),
        showlegend = TRUE,
        name = category,
        legendgroup = paste0("category_", category),
        hoverinfo = "none",
        visible = TRUE
      )
    }
    
    # Add legend entries for gene structure (same as regular analysis)
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
    
    # Add peptides for ALL selected transcripts (same pattern as regular analysis) - UPDATED FOR COMPRESSION
    if (nrow(all_peptides_df) > 0) {
      for (i in 1:nrow(all_peptides_df)) {
        # All peptides get specificity colors
        color <- specificity_colors[[all_peptides_df$specificity_category[i]]]
        legend_group <- paste0("category_", all_peptides_df$specificity_category[i])
        
        # Apply compression to peptide coordinates
        peptide_start <- if (!is.null(compression_map)) compression_map$compress(all_peptides_df$start[i]) else all_peptides_df$start[i]
        peptide_end <- if (!is.null(compression_map)) compression_map$compress(all_peptides_df$end[i]) else all_peptides_df$end[i]
        
        # Add peptide as filled rectangle (same pattern as regular analysis)
        p <- p %>% plotly::add_trace(
          type = "scatter", mode = "lines",
          x = c(peptide_start, peptide_end, peptide_end, peptide_start, peptide_start),
          y = c(all_peptides_df$y_position[i] - 0.15, all_peptides_df$y_position[i] - 0.15, all_peptides_df$y_position[i] + 0.15, all_peptides_df$y_position[i] + 0.15, all_peptides_df$y_position[i] - 0.15),
          fill = "toself",
          fillcolor = color,
          line = list(color = "black", width = 0.5),
          marker = list(opacity = 0),
          legendgroup = legend_group,
          showlegend = FALSE,
          hoverinfo = "text",
          text = all_peptides_df$hover_text[i]
        )
      }
    }
    
    # Set layout (same as regular analysis)
    miscleavage_label <- switch(input$novel_miscleavage_type,
      "no_miss_cleavage" = "No Miscleavage",
      "upto_two_misscleavage" = "Up to 2 Miscleavages"
    )
    
    # Configure X-axis with proper genomic coordinate mapping (like other tabs)
    axis_breaks <- compression_map$coords$compressed_start
    axis_labels <- as.character(compression_map$coords$original_start)
    compressed_range <- range(c(compression_map$coords$compressed_start, compression_map$coords$compressed_end))
    
    p <- p %>% plotly::layout(
      title = paste0('Novel vs Reference Multi-Isoform Comparison (', toupper(input$novel_protease), ', ', miscleavage_label, ') - ', length(all_transcripts), ' Selected - Compressed View - Chromosome ', chromosome),
      xaxis = list(
        title = paste0("Genomic Position (compressed view) - chromosome ", chromosome),
        tickvals = axis_breaks,
        ticktext = axis_labels,
        range = compressed_range,
        zeroline = FALSE
      ),
      yaxis = list(
        title = "Transcripts",
        tickvals = transcript_df$y_position,
        ticktext = transcript_df$transcript,
        range = c(0.5, max(transcript_df$y_position) + 0.5)
      ),
      showlegend = TRUE,
      height = 400,
      hovermode = 'closest'
    )
    
    return(p)
     
  }, error = function(e) {
    cat("Error in novel_comparative_plot:", e$message, "\n")
    return(NULL)
  })
})

# Render novel highlighted isoform table (for multi-isoform comparison)
output$novel_highlighted_isoform_table <- DT::renderDataTable({
  req(novel_multi_isoform_highlighted_data(), input$novel_highlight_isoform)
  
  data <- novel_multi_isoform_highlighted_data()
  highlight_isoform <- input$novel_highlight_isoform
  
  if (is.null(data) || is.null(data$all_peptides) || highlight_isoform == "") {
    return(data.frame(Message = "No data available"))
  }

  all_peptides <- data$all_peptides
  has_genomic_info <- "chromosome" %in% names(all_peptides) && "genomic_positions" %in% names(all_peptides)

  # 1. Group by peptide sequence to see which transcripts they appear in.
  if (has_genomic_info) {
      peptides_grouped <- all_peptides %>%
          dplyr::group_by(peptide, chromosome, genomic_positions) %>%
          dplyr::summarise(transcripts = list(unique(transcript)), .groups = 'drop')
  } else {
      peptides_grouped <- all_peptides %>%
          dplyr::group_by(peptide) %>%
          dplyr::summarise(transcripts = list(unique(transcript)), .groups = 'drop')
  }

  # 2. Filter this grouped list to only include peptides present in our highlighted isoform
  peptides_filtered <- peptides_grouped %>%
    dplyr::filter(sapply(transcripts, function(tx_list) highlight_isoform %in% tx_list))

  if (nrow(peptides_filtered) == 0) {
    return(data.frame(Message = paste("No peptides found for", highlight_isoform)))
  }
    
  # 3. Determine specificity category and other transcript info
  selected_transcripts_for_comparison <- unique(all_peptides$transcript)
  peptides_final <- peptides_filtered %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      specificity_category = dplyr::case_when(
        length(transcripts) == 1 ~ "Unique",
        length(transcripts) == length(selected_transcripts_for_comparison) ~ "Universal",
        TRUE ~ "Shared"
      ),
      other_transcripts = paste(setdiff(transcripts, highlight_isoform), collapse = ", "),
      other_transcript_count = length(setdiff(transcripts, highlight_isoform))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-transcripts) %>%
    dplyr::arrange(factor(specificity_category, levels = c("Unique", "Shared", "Universal")))

  # 4. Prepare for display
  if (has_genomic_info) {
    display_df <- peptides_final[, c("peptide", "chromosome", "genomic_positions", "specificity_category", "other_transcript_count", "other_transcripts")]
    colnames(display_df) <- c("Peptide", "Chromosome", "Genomic Positions", "Specificity", "Other Isoforms Count", "Other Isoforms")
  } else {
    display_df <- peptides_final[, c("peptide", "specificity_category", "other_transcript_count", "other_transcripts")]
    colnames(display_df) <- c("Peptide", "Specificity", "Other Isoforms Count", "Other Isoforms")
  }

  DT::datatable(
    display_df,
    options = list(
      pageLength = 10,
      scrollX = TRUE,
      order = list(list(which(names(display_df) == "Specificity") - 1, 'asc'))
    ),
    rownames = FALSE
  ) %>%
  DT::formatStyle(
    'Specificity',
    backgroundColor = styleEqual(c("Unique", "Shared", "Universal"), c("#FF0000", "#F39C12", "#2ECC71")),
    color = styleEqual(c("Unique", "Shared", "Universal"), c("white", "white", "white"))
  )
})


# Download handlers for novel comparison data
output$download_novel_comparison_data <- downloadHandler(
  filename = function() {
    paste0("novel_comparison_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  },
  content = function(file) {
    data <- novel_multi_isoform_highlighted_data()
    if (!is.null(data) && !is.null(data$all_peptides)) {
      write.csv(data$all_peptides, file, row.names = FALSE)
    }
  }
)

output$download_novel_peptides_analysis <- downloadHandler(
  filename = function() {
    paste0("novel_peptides_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  },
  content = function(file) {
    data <- novel_multi_isoform_highlighted_data()
    if (!is.null(data) && !is.null(data$all_peptides)) {
      # Create detailed analysis
      analysis_data <- data$all_peptides
      analysis_data$enzyme <- input$novel_protease
      analysis_data$miscleavage_type <- input$novel_miscleavage_type
      analysis_data$analysis_timestamp <- Sys.time()
      write.csv(analysis_data, file, row.names = FALSE)
    }
  }
)

              # Multi-Isoform Comparative Analysis - Independent data loading
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
         core_data_module$current_gene_symbol()
   }, error = function(e) {
         input$gene
       })
       if (is.null(gene_symbol)) gene_symbol <- input$gene
     }
   }
   
   # Get novel transcripts if available
   novel_transcripts <- NULL
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
   
   # Combine choices
   all_choices <- character(0)
   
   # Add known gene transcripts
   if (!is.null(known_transcripts) && length(known_transcripts) > 0) {
     known_choices <- setNames(known_transcripts, paste0(known_transcripts, " (", gene_symbol, ")"))
     all_choices <- c(all_choices, known_choices)
   }
   
   # Add novel transcripts
   if (!is.null(novel_transcripts) && length(novel_transcripts) > 0) {
     novel_choices <- setNames(novel_transcripts, paste0(novel_transcripts, " (Novel)"))
     all_choices <- c(all_choices, novel_choices)
   }
   
   # Update dropdown
   if (length(all_choices) == 0) {
   updateSelectizeInput(session, "compare_isoforms", 
                         choices = NULL, selected = character(0))
   } else {
     updateSelectizeInput(session, "compare_isoforms", 
                         choices = all_choices,
                       selected = character(0))
     
     # Show notification about available options
     if (!is.null(known_transcripts) && !is.null(novel_transcripts)) {
       cat("✅ Multiple Isoform Comparison: Added", length(known_transcripts), "known +", length(novel_transcripts), "novel transcripts\n")
     }
   }
 })

# Multi-Isoform data loader (exact copy of all_isoforms_data but triggered by run_comparative_analysis)
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
    # Use the gene-by-gene loaded peptides which already have the correct miscleavage type
    gene_peptides <- selected_gene_peptides()
    if (is.null(gene_peptides) || nrow(gene_peptides) == 0) {
      showNotification(paste("No peptides found for gene", gene_id, "with miscleavage type", miscleavage_type), type = "warning")
      return(NULL)
    }
    
    # Create standardized vis_data structure using core data module  
    vis_data <- core_data_module$create_vis_data_structure()
    
    # Get peptides for SELECTED transcripts only (key difference from all_isoforms_data)
    all_peptides_list <- list()
    for (i in seq_along(selected_transcripts)) {
      tx <- selected_transcripts[i]
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
    
    # Create transcript position mapping for SELECTED transcripts only
    transcript_df <- data.frame(
      transcript = selected_transcripts,
      y_position = seq_along(selected_transcripts),
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
      all_transcripts = selected_transcripts  # Use selected instead of all
    ))
  })
   })

# Multi-Isoform highlighted data (exact copy of highlighted_isoform_data logic but for all selected transcripts)
multi_isoform_highlighted_data <- reactive({
  req(multi_isoform_data())
  
  base_data <- multi_isoform_data()
  if (is.null(base_data)) return(NULL)
  
  all_peptides_df <- base_data$all_peptides
  all_transcripts <- base_data$all_transcripts
  
  # Calculate specificity for ALL peptides across selected transcripts
  all_peptides_df$specificity_category <- ""
  
  for (i in 1:nrow(all_peptides_df)) {
    peptide_seq <- all_peptides_df$peptide[i]
    current_transcript <- all_peptides_df$transcript[i]
    
    # Count transcripts (excluding current) that have this peptide
    other_tx_with_peptide <- unique(all_peptides_df$transcript[
      all_peptides_df$peptide == peptide_seq & all_peptides_df$transcript != current_transcript
    ])
    
    # Classify specificity into three categories: Unique, Shared, Universal
    total_selected_isoforms <- length(all_transcripts)
    other_count <- length(other_tx_with_peptide)
    
    if (other_count == 0) {
      # Peptide found only in this isoform among selected
      all_peptides_df$specificity_category[i] <- "Unique"
    } else if (other_count == (total_selected_isoforms - 1)) {
      # Peptide found in ALL selected isoforms
      all_peptides_df$specificity_category[i] <- "Universal"
    } else {
      # Peptide found in some but not all selected isoforms
      all_peptides_df$specificity_category[i] <- "Shared"
    }
  }
  
  # Update hover text with specificity
  all_peptides_df$hover_text <- paste0(
    "Peptide: ", all_peptides_df$peptide,
    "<br>Position: ", all_peptides_df$start, "-", all_peptides_df$end,
    "<br>Transcript: ", all_peptides_df$transcript,
    "<br>Specificity: ", all_peptides_df$specificity_category,
    "<br>Miscleavage: ", input$miscleavage_type
  )
  
  return(list(
    all_peptides = all_peptides_df,
    transcript_df = base_data$transcript_df,
    gene_start = base_data$gene_start,
    gene_end = base_data$gene_end,
    all_transcripts = all_transcripts
  ))
})

# Enable/disable comparative analysis button based on selection count
observeEvent(input$compare_isoforms, {
  selected_count <- length(input$compare_isoforms)
  
  if (selected_count >= 2 && selected_count <= 8) {
    shinyjs::enable("run_comparative_analysis")
    shinyjs::removeClass("run_comparative_analysis", "btn-secondary")
    shinyjs::addClass("run_comparative_analysis", "btn-warning")
  } else {
    shinyjs::disable("run_comparative_analysis")
    shinyjs::removeClass("run_comparative_analysis", "btn-warning") 
    shinyjs::addClass("run_comparative_analysis", "btn-secondary")
  }
})

# Reactive for comparative analysis data - Using same approach as highlighted_isoform_data
comparative_analysis_data <- eventReactive(input$run_comparative_analysis, {
  req(input$compare_isoforms, length(input$compare_isoforms) >= 2,
      all_isoforms_data())
  
  withProgress(message = 'Running comparative analysis...', value = 0, {
    compare_isoforms <- input$compare_isoforms
    
    incProgress(0.2, detail = 'Loading data from all_isoforms_data...')
    
    # Get base data (same approach as highlighted_isoform_data)
    base_data <- all_isoforms_data()
    if (is.null(base_data)) return(NULL)
    
    # Filter to only selected isoforms (this is the key difference from Isoform-Specific)
    all_peptides_df <- base_data$all_peptides
    selected_peptides <- all_peptides_df[all_peptides_df$transcript %in% compare_isoforms, ]
    
    # Use the same gene boundaries as the base data
    gene_start <- base_data$gene_start
    gene_end <- base_data$gene_end
    
    incProgress(0.3, detail = 'Merging duplicate peptides...')
    
    # Merge duplicate peptides for each transcript to avoid counting the same peptide multiple times
    merged_by_transcript <- list()
    for (tx in compare_isoforms) {
      tx_peptides <- selected_peptides[selected_peptides$transcript == tx, ]
      if (nrow(tx_peptides) > 0) {
        # Merge duplicate peptides (same peptide appears multiple times due to multiple genomic locations)
        merged_tx_peptides <- tx_peptides %>%
          dplyr::group_by(peptide) %>%
          dplyr::summarise(
            transcript = first(transcript),
            start_pos = min(start, na.rm = TRUE),
            end_pos = max(end, na.rm = TRUE),
            y_position = match(tx, compare_isoforms),  # Set y position based on order
            .groups = 'drop'
          )
        merged_by_transcript[[tx]] <- merged_tx_peptides
      }
    }
    
    # Combine merged data
    selected_peptides_merged <- data.table::rbindlist(merged_by_transcript)
    
    incProgress(0.2, detail = 'Calculating peptide overlaps...')
    
    # Calculate peptide overlap matrix using merged data
    overlap_matrix <- matrix(0, nrow = length(compare_isoforms), ncol = length(compare_isoforms))
    rownames(overlap_matrix) <- compare_isoforms
    colnames(overlap_matrix) <- compare_isoforms
    
    for (i in 1:length(compare_isoforms)) {
      for (j in 1:length(compare_isoforms)) {
        peptides_i <- unique(selected_peptides_merged$peptide[selected_peptides_merged$transcript == compare_isoforms[i]])
        peptides_j <- unique(selected_peptides_merged$peptide[selected_peptides_merged$transcript == compare_isoforms[j]])
        
        if (i == j) {
          overlap_matrix[i, j] <- length(peptides_i)
        } else {
          overlap <- length(intersect(peptides_i, peptides_j))
          overlap_matrix[i, j] <- overlap
        }
      }
    }
    
    incProgress(0.3, detail = 'Creating visualization data...')
    
    # Create transcript_df for selected transcripts only (filter from base_data)
    base_transcript_df <- base_data$transcript_df
    transcript_df <- base_transcript_df[base_transcript_df$transcript %in% compare_isoforms, ]
    # Re-assign y positions for the selected transcripts
    transcript_df$y_position <- seq_along(compare_isoforms)
    
    # Update peptide y positions using merged data
    selected_peptides_merged$y_position <- transcript_df$y_position[match(selected_peptides_merged$transcript, transcript_df$transcript)]
    
    # Add comparative specificity categories (3 categories only) using merged data
    selected_peptides_merged$comparative_specificity <- ""
    for (i in 1:nrow(selected_peptides_merged)) {
      peptide_seq <- selected_peptides_merged$peptide[i]
      other_tx_with_peptide <- unique(selected_peptides_merged$transcript[
        selected_peptides_merged$peptide == peptide_seq & selected_peptides_merged$transcript != selected_peptides_merged$transcript[i]
      ])
      
      if (length(other_tx_with_peptide) == 0) {
        selected_peptides_merged$comparative_specificity[i] <- "Unique"
      } else if (length(other_tx_with_peptide) == (length(compare_isoforms) - 1)) {
        selected_peptides_merged$comparative_specificity[i] <- "Universal"
      } else {
        selected_peptides_merged$comparative_specificity[i] <- "Shared"
      }
    }
    
    # Create summary statistics using merged data
    summary_stats <- list()
    for (tx in compare_isoforms) {
      tx_peptides <- selected_peptides_merged[selected_peptides_merged$transcript == tx, ]
      summary_stats[[tx]] <- data.frame(
        Transcript = tx,
        Total_Peptides = nrow(tx_peptides),
        Unique_Peptides = sum(tx_peptides$comparative_specificity == "Unique"),
        Shared_Peptides = sum(tx_peptides$comparative_specificity == "Shared"),
        Universal_Peptides = sum(tx_peptides$comparative_specificity == "Universal"),
        stringsAsFactors = FALSE
      )
    }
    
    summary_df <- data.table::rbindlist(summary_stats)
    
    incProgress(0.2, detail = 'Finalizing...')
    
    return(list(
      peptides = selected_peptides_merged,
      transcript_df = transcript_df,
      gene_start = gene_start,
      gene_end = gene_end,
      overlap_matrix = overlap_matrix,
      summary = summary_df,
      compare_isoforms = compare_isoforms
    ))
  })
})

  # Render comparative plot (exact same pattern as all_isoforms_plot but without highlighting)
output$comparative_plot <- renderPlotly({
  data <- multi_isoform_highlighted_data()
  
  if (is.null(data)) {
    return(empty_plotly_message("No data available for comparative analysis. Please select at least two isoforms and click 'Run'."))
  }
  
  # Define miscleavage_label based on input for the plot title
  miscleavage_label <- switch(input$miscleavage_type,
                              "no_miss_cleavage" = "No Missed Cleavages",
                              "upto_one_misscleavage" = "Up to 1 Missed Cleavage",
                              "upto_two_misscleavage" = "Up to 2 Missed Cleavages",
                              "Unknown")
  
  all_peptides <- data$all_peptides
  transcript_df <- data$transcript_df
  gene_start <- data$gene_start
  gene_end <- data$gene_end
  all_transcripts <- data$all_transcripts
  
  # Load GTF data for exon and CDS boundaries (lightning-fast)
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
  
  # Create plotly object
  p <- plotly::plot_ly()
  
  # Define colors for specificity categories using rgba to prevent hex codes in hover
  specificity_colors <- list(
    "Unique" = "rgba(255, 0, 0, 0.8)",      # Red - peptides unique to this isoform
    "Shared" = "rgba(243, 156, 18, 0.8)",   # Orange - peptides shared with some isoforms
    "Universal" = "rgba(46, 204, 113, 0.8)" # Green - peptides found in all isoforms
  )
    
        # Add transcript lines (ALL transcripts, no highlighting - difference from main plot)
    for (i in 1:nrow(transcript_df)) {
      p <- p %>% plotly::add_trace(
        type = "scatter",
        x = c(gene_start, gene_end),
        y = c(transcript_df$y_position[i], transcript_df$y_position[i]),
        mode = "lines",
        line = list(color = "#666666", width = 2),  # Same color for all (no highlighting)
        showlegend = FALSE,
        hoverinfo = "text",
        text = clean_hover_text(paste0("Transcript: ", transcript_df$transcript[i]))
      )
    }
  
  # Add exon and CDS boundaries if GTF data is available (same as Isoform-Specific analysis)
  if (exons_result$success && length(exons_result$exons) > 0) {
    exons_by_transcript <- exons_result$exons
    cds_by_transcript <- exons_result$cds
    
    # Add exon blocks for each transcript
    for (i in 1:nrow(transcript_df)) {
      tx <- transcript_df$transcript[i]
      y_pos <- transcript_df$y_position[i]
      
      # Add exons if available (transparent grey)
      if (!is.null(exons_by_transcript[[tx]]) && length(exons_by_transcript[[tx]]) > 0) {
        tx_exons <- exons_by_transcript[[tx]]
        for (j in seq_along(tx_exons)) {
          p <- p %>% plotly::add_trace(
            type = "scatter", mode = "lines",
            x = c(start(tx_exons[j]), end(tx_exons[j]), end(tx_exons[j]), start(tx_exons[j]), start(tx_exons[j])),
            y = c(y_pos - 0.35, y_pos - 0.35, y_pos + 0.35, y_pos + 0.35, y_pos - 0.35),
            fill = "toself",
            fillcolor = "rgba(211, 211, 211, 0.3)",  # Light grey for exons (consistent with visualization tab)
            line = list(color = "#D3D3D3", width = 1),  # Light grey border
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
          p <- p %>% plotly::add_trace(
            type = "scatter", mode = "lines",
            x = c(start(tx_cds[j]), end(tx_cds[j]), end(tx_cds[j]), start(tx_cds[j]), start(tx_cds[j])),
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
  
  # Add legend entries for each category
  for (category in names(specificity_colors)) {
    p <- p %>% plotly::add_trace(
      type = "scatter", mode = "markers",
      x = c(0), y = c(0),
      marker = list(color = specificity_colors[[category]], size = 10),
      showlegend = TRUE,
      name = category,
      legendgroup = paste0("category_", category),
      hoverinfo = "none",
      visible = TRUE
    )
  }
  
  # Add legend entries for gene structure (same as Isoform-Specific analysis)
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
  
  # Add peptides for ALL selected transcripts (same pattern as main plot)
  if (nrow(all_peptides) > 0) {
    for (i in 1:nrow(all_peptides)) {
      # All peptides get specificity colors (no highlighting difference)
      color <- specificity_colors[[all_peptides$specificity_category[i]]]
      legend_group <- paste0("category_", all_peptides$specificity_category[i])
      
      # Add peptide as filled rectangle (same pattern as main plot)
      p <- p %>% plotly::add_trace(
        type = "scatter", mode = "lines",
        x = c(all_peptides$start[i], all_peptides$end[i], all_peptides$end[i], all_peptides$start[i], all_peptides$start[i]),
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
  chromosome <- if(exists("gene_details") && !is.null(gene_details) && gene_details$success) {
    gene_details$chromosome
  } else {
    # Try to get from gene loading
    gene_details_local <- load_gene_details(input$gene)
    if(gene_details_local$success) gene_details_local$chromosome else "Unknown"
  }
  
  # Extract strand from peptide data if available
  strand_display <- "Unknown"
  if(nrow(all_peptides) > 0) {
    # Get first transcript's peptides to check strand (assuming all same gene)
    gene_id <- input$gene
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
  
  p <- p %>% plotly::layout(
    title = paste0('Multi-Isoform Comparative Analysis (', miscleavage_label, ') - ', length(all_transcripts), ' Selected'),
    xaxis = list(
      title = paste0("Genomic Position (chromosome ", chromosome, ")"),
      range = c(gene_start, gene_end + 2000),  # Same range as main plot
      zeroline = FALSE
    ),
    yaxis = list(
      title = "",
      range = c(0.3, nrow(transcript_df) + 0.7),
      tickvals = transcript_df$y_position,
      ticktext = transcript_df$transcript,
      zeroline = FALSE
    ),
    legend = list(
      orientation = "v", x = 1.02, y = 0.5,
      title = list(text = "Legend:"),  # Same title as main plot
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

# Render peptide overlap matrix table
output$overlap_matrix_table <- DT::renderDataTable({
  data <- multi_isoform_highlighted_data()
  
  if (is.null(data)) {
    return(data.frame(Message = "No comparative data available"))
  }
  
  all_peptides <- data$all_peptides
  all_transcripts <- data$all_transcripts
  
  # Create overlap matrix
  n_tx <- length(all_transcripts)
  overlap_matrix <- matrix(0, nrow = n_tx, ncol = n_tx)
  rownames(overlap_matrix) <- all_transcripts
  colnames(overlap_matrix) <- all_transcripts
  
  # Fill diagonal with total peptides per transcript
  for (i in 1:n_tx) {
    tx <- all_transcripts[i]
    overlap_matrix[i, i] <- nrow(all_peptides[all_peptides$transcript == tx, ])
  }
  
  # Fill off-diagonal with shared peptides
  for (i in 1:n_tx) {
    for (j in 1:n_tx) {
      if (i != j) {
        tx1 <- all_transcripts[i]
        tx2 <- all_transcripts[j]
        
        # Get peptides from both transcripts
        tx1_peptides <- unique(all_peptides$peptide[all_peptides$transcript == tx1])
        tx2_peptides <- unique(all_peptides$peptide[all_peptides$transcript == tx2])
        
        # Count shared peptides
        shared_count <- length(intersect(tx1_peptides, tx2_peptides))
        overlap_matrix[i, j] <- shared_count
      }
    }
  }
  
  # Convert matrix to data frame for display
  overlap_df <- as.data.frame(overlap_matrix)
  overlap_df <- cbind(Transcript = rownames(overlap_df), overlap_df)
  
  DT::datatable(
    overlap_df,
    options = list(
      pageLength = 10,
      scrollX = TRUE,
      ordering = FALSE
    ),
    rownames = FALSE,
    caption = "Peptide overlap matrix: diagonal shows total peptides per isoform, off-diagonal shows shared peptides"
  ) %>%
  DT::formatStyle(
    columns = 2:ncol(overlap_df),
    backgroundColor = styleInterval(
      c(0, 5, 10),
      c("#FFE6E6", "#FFD1D1", "#FFB3B3", "#FF9999")
    )
  )
})

# Render comparative summary table
output$comparative_summary_table <- DT::renderDataTable({
  data <- multi_isoform_highlighted_data()
  
  if (is.null(data)) {
    return(data.frame(Message = "No comparative data available"))
  }
  
  all_peptides <- data$all_peptides
  all_transcripts <- data$all_transcripts
  
  # Create summary for each transcript
  summary_list <- list()
  for (tx in all_transcripts) {
    tx_peptides <- all_peptides[all_peptides$transcript == tx, ]
    
    # Count by specificity category
    unique_count <- sum(tx_peptides$specificity_category == "Unique")
    shared_count <- sum(tx_peptides$specificity_category == "Shared")
    universal_count <- sum(tx_peptides$specificity_category == "Universal")
    total_count <- nrow(tx_peptides)
    
    summary_list[[tx]] <- data.frame(
      Transcript = tx,
      Total_Peptides = total_count,
      Unique_Peptides = unique_count,
      Shared_Peptides = shared_count,
      Universal_Peptides = universal_count,
      stringsAsFactors = FALSE
    )
  }
  
  summary_df <- data.table::rbindlist(summary_list)
  
  DT::datatable(
    summary_df,
    options = list(
      pageLength = 10,
      scrollX = TRUE,
      order = list(list(1, 'desc'))  # Sort by total peptides
    ),
    rownames = FALSE,
    caption = "Summary of peptide specificity for each selected isoform"
  ) %>%
  DT::formatStyle(
    'Unique_Peptides',
    backgroundColor = styleInterval(
      c(0, 5),
      c("#FFE6E6", "#FFB3B3", "#FF8080")
    )
  )
})

# Status message for comparative analysis
output$comparative_analysis_status <- renderUI({
  selected_count <- length(input$compare_isoforms)
  
  if (selected_count == 0) {
    div(
      style = "color: #999; font-style: italic; margin-top: 10px;",
      "Please select 2 or more isoforms to compare"
  )
  } else if (selected_count == 1) {
    div(
      style = "color: #ff9900; margin-top: 10px;",
      paste("Selected 1 isoform. Please select at least 2 for comparison.")
    )
  } else if (selected_count > 8) {
    div(
      style = "color: #ff0000; margin-top: 10px;",
      paste("Selected", selected_count, "isoforms. Please select no more than 8 for optimal visualization.")
    )
  } else {
    div(
      style = "color: #00aa00; margin-top: 10px;",
      paste("Ready to compare", selected_count, "isoforms. Click 'Compare Selected Isoforms' to proceed.")
    )
  }
})

# Download handler for comparative analysis data
output$download_comparative_data <- downloadHandler(
  filename = function() {
    miscleavage_suffix <- switch(input$miscleavage_type,
      "no_miss_cleavage" = "no_miss",
      "upto_two_misscleavage" = "upto_2_miss"
    )
    paste0("comparative_analysis_", length(input$compare_isoforms), "_isoforms_", 
          miscleavage_suffix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  },
  content = function(file) {
    data <- multi_isoform_highlighted_data()
    if (!is.null(data)) {
      comparative_peptides <- data$all_peptides
      
      # Add metadata
      comparative_peptides$gene <- input$gene
      comparative_peptides$protease <- input$protease
      comparative_peptides$miscleavage_type <- input$miscleavage_type
      comparative_peptides$compared_isoforms <- paste(input$compare_isoforms, collapse = ", ")
      
      write.csv(comparative_peptides, file, row.names = FALSE)
    }
  }
)

      # When all isoforms analysis is loaded, update highlight choices
 observeEvent(all_isoforms_data(), {
   data <- all_isoforms_data()
   req(data)
   
   # Safely get all transcripts from the plot data
   all_transcripts <- data$all_transcripts
  if (is.null(all_transcripts) || length(all_transcripts) == 0) {
     return()
   }
  
  # Create choices with gene symbols
   gene_id <- input$gene
   if (is.null(gene_id)) return()
   
   gene_symbol <- tryCatch({
     core_data_module$current_gene_symbol()
   }, error = function(e) {
     gene_id
   })
   
   if (is.null(gene_symbol)) gene_symbol <- gene_id
  
  choices <- setNames(all_transcripts, paste0(all_transcripts, " (", gene_symbol, ")"))
  
  updateSelectInput(session, "highlight_isoform", 
                    choices = choices,
                     selected = if(length(choices) > 0) choices[1] else NULL)
})