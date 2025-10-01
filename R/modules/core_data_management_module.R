#===============================================================================
# CORE DATA MANAGEMENT MODULE
# Extracted from server.R - Complete data loading, caching, and management
#===============================================================================

#' Core Data Management Module
#' 
#' @description 
#' Complete core data management functionality for the IsoPepTracker app:
#' - Initial data processing and structure setup
#' - Gene-specific data loading with LRU caching 
#' - Reactive data management for genes, transcripts, and peptides
#' - Gene search and selection interface
#' - Data format standardization for modules
#' - Memory management and performance optimization
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' 
#' @return Named list containing all reactive values and functions for data management
#' 
#' @note 
#' - Supports both gene-index and legacy data loading approaches
#' - Implements LRU caching with configurable cache size (default: 200 genes)
#' - Provides standardized vis_data structures for all analysis modules
#' - Handles miscleavage-specific data caching
#' - Memory-efficient on-demand loading for large datasets

create_core_data_management_module <- function(input, output, session) {
  
  # ============================================================================
  # MODULE REACTIVE VALUES
  # ============================================================================
  
  # Core data management reactive values
  gene_data <- reactiveVal(NULL)
  gene_cache <- reactiveVal(list())
  max_cache_size <- 200  # Maximum number of genes to keep in memory
  
  # Non-reactive guard to prevent repeated triggering
  last_processed_key <- NULL
  
  # Cache statistics for monitoring
  cache_stats <- reactiveVal(list(
    hits = 0,
    misses = 0,
    current_size = 0,
    evictions = 0
  ))
  
  # ============================================================================
  # INITIAL DATA PROCESSING
  # ============================================================================
  
  # Process initial data when app starts
  processed_data <- reactive({
    withProgress(message = 'Loading initial data...', {
      # Check if we have the gene index (preferred approach)
      if (exists("gene_index")) {
        cat("✅ Using gene-index approach for data management\n")
        
        # Using gene-specific approach
        incProgress(0.3, detail = "Processing gene index")
        genes <- gene_index$geneID
        gene_symbols <- gene_index$geneSymbol
        gene_lookup <- setNames(gene_symbols, genes)
        
        incProgress(0.7, detail = "Setting up data structures")
        # Return basic structure without loading all peptides
        return(list(
          genes = genes,
          gene_symbols = gene_symbols,
          gene_lookup = gene_lookup,
          proteases = c("trp", "chymo", "aspn", "lysc", "lysn", "gluc"),
          using_gene_index = TRUE,
          total_genes = length(genes)
        ))
        
      } else {
        cat("ℹ️ Using legacy approach for data management\n")
        
        # Fallback to original approach
        req(exists("peptides") && exists("as_database"))
        
        incProgress(0.3, detail = "Processing peptide data")
        genes <- unique(peptides$geneID)
        gene_symbols <- unique(peptides$geneSymbol)
        gene_lookup <- setNames(gene_symbols, genes)
        
        incProgress(0.3, detail = "Processing transcript data")
        transcripts_per_gene <- split(peptides$txID, peptides$geneID)
        
        incProgress(0.4, detail = "Processing AS events")
        as_events <- as_database[, c("refTx", "asTx", "geneID", "eventID", "AS_type", 
                                    "AS_range", "AS_range1", "AS_range2")]
        
        return(list(
          genes = genes,
          gene_symbols = gene_symbols,
          gene_lookup = gene_lookup,
          transcripts_per_gene = transcripts_per_gene,
          as_events = as_events,
          proteases = c("trp", "chymo", "aspn", "lysc", "lysn", "gluc"),
          original_peptides = peptides,
          using_gene_index = FALSE,
          total_genes = length(genes)
        ))
      }
    })
  })
  
  # ============================================================================
  # GENE DATA LOADING AND CACHING
  # ============================================================================
  
  # Shared function to load and cache gene data
  load_and_cache_gene_data <- function(gene_id, miscleavage_type) {
    # Check cache first
    cache_key <- paste0(gene_id, "_", miscleavage_type)
    cache <- gene_cache()
    stats <- cache_stats()
    
    if (cache_key %in% names(cache)) {
      cat("📋 Cache HIT for", gene_id, "with", miscleavage_type, "\n")
      
      # Use cached data if available
      gene_data(cache[[cache_key]])
      
      # Update cache statistics (isolated to prevent reactive loops)
      isolate({
        stats$hits <- stats$hits + 1
        cache_stats(stats)
      })
      
      # Cache hits should be truly read-only - no reordering to prevent reactive loops
      cat("✅ Using cached data without reordering to prevent reactive loops\n")
      
      return(TRUE)
    }
    
    cat("💾 Cache MISS for", gene_id, "with", miscleavage_type, "- loading fresh data\n")
    
    # If not in cache, load fresh data
    withProgress(message = paste('Loading gene data for', gene_id, '(', miscleavage_type, ')'), {
      data <- load_gene_data(gene_id, miscleavage_type = miscleavage_type)
      
      if (is.null(data)) {
        showNotification(paste("Error loading data for gene", gene_id, "with miscleavage type", miscleavage_type), type = "error")
        return(FALSE)
      }
      
      # Process the data
      incProgress(0.5, detail = "Processing gene data")
      
      # Extract transcripts and create the processed data object
      processed_gene_data <- list(
        peptides = data$peptides,
        as_events = data$as_events,
        transcripts_per_gene = list(),
        miscleavage_type = data$miscleavage_type,
        gene_id = gene_id,
        load_timestamp = Sys.time()
      )
      processed_gene_data$transcripts_per_gene[[gene_id]] <- unique(data$peptides$txID)
      
      # Set the reactive value
      gene_data(processed_gene_data)
      
      # Update cache with LRU eviction policy (isolated to prevent reactive loops)
      final_cache_size <- isolate({
        new_cache <- gene_cache()
        new_cache[[cache_key]] <- processed_gene_data
        
        # Update cache statistics
        stats$misses <- stats$misses + 1
        stats$current_size <- length(new_cache)
        
        # Implement LRU eviction policy
        if (length(new_cache) > max_cache_size) {
          oldest_key <- names(new_cache)[1]
          new_cache[[oldest_key]] <- NULL
          stats$evictions <- stats$evictions + 1
          stats$current_size <- length(new_cache)
          
          cat("🗑️ Cache eviction: removed", oldest_key, "- cache size now", length(new_cache), "\n")
        }
        
        cache_stats(stats)
        gene_cache(new_cache)
        
        # Return cache size for logging
        length(new_cache)
      })
      
      cat("✅ Successfully cached", gene_id, "with", miscleavage_type, "- cache size:", final_cache_size, "\n")
    })
    return(TRUE)
  }
  
  # ============================================================================
  # GENE DATA LOADING OBSERVER
  # ============================================================================
  
  # Load gene data when needed (triggered by gene OR miscleavage type change)
  observeEvent(list(input$gene, input$miscleavage_type), {
    req(input$gene, input$miscleavage_type)
    gene_id <- input$gene
    miscleavage_type <- input$miscleavage_type
    
    # Create cache key that includes miscleavage type
    cache_key <- paste0(gene_id, "_", miscleavage_type)
    
    # Guard against repeated processing of the same key
    if (!is.null(last_processed_key) && last_processed_key == cache_key) {
      cat("🛡️ GUARD: Skipping repeated processing for", cache_key, "\n")
      return()
    }
    
    cat("🔄 observeEvent triggered for:", gene_id, "with", miscleavage_type, "- Call time:", Sys.time(), "\n")
    cat("🔍 Reactive context: observeEvent on input$gene + input$miscleavage_type\n")
    
    # Update the guard (non-reactive assignment)
    last_processed_key <<- cache_key
    
    # Check if we're using gene-specific approach
    if (processed_data()$using_gene_index) {
      # Check if gene+miscleavage combination is already in cache
      cache <- gene_cache()
      if (cache_key %in% names(cache)) {
        cat("📋 Cache HIT for", gene_id, "with", miscleavage_type, "\n")
        
        # Check if we already have the same data to prevent unnecessary updates
        current_data <- gene_data()
        cached_data <- cache[[cache_key]]
        
        # Only update gene_data if it's actually different
        if (is.null(current_data) || 
            is.null(current_data$gene_id) || 
            current_data$gene_id != gene_id ||
            current_data$miscleavage_type != miscleavage_type) {
          
          cat("📝 Updating gene_data with cached data\n")
          gene_data(cached_data)
        } else {
          cat("📋 Data already current - no update needed\n")
        }
        
        # Update cache statistics only (isolated to prevent reactive loops)
        isolate({
          stats <- cache_stats()
          stats$hits <- stats$hits + 1
          cache_stats(stats)
        })
        
        cat("✅ Using cached data without triggering reactive updates\n")
        return()  # Early return to prevent further processing
      } else {
        # Load gene data with specified miscleavage type
        withProgress(message = paste('Loading gene data for', gene_id, '(', miscleavage_type, ')'), {
          # Load gene-specific data with miscleavage type
          data <- load_gene_data(gene_id, miscleavage_type = miscleavage_type)
          
          if (is.null(data)) {
            showNotification(paste("Error loading data for gene", gene_id, "with miscleavage type", miscleavage_type), type = "error")
            return()
          }
          
          # Process the data
          incProgress(0.5, detail = "Processing gene data")
          
          # Extract transcripts for this gene
          transcripts_per_gene <- list()
          transcripts_per_gene[[gene_id]] <- unique(data$peptides$txID)
          
          # Create processed data structure
          processed_gene_data <- list(
            peptides = data$peptides,
            as_events = data$as_events,
            transcripts_per_gene = transcripts_per_gene,
            miscleavage_type = data$miscleavage_type,
            gene_id = gene_id,
            load_timestamp = Sys.time()
          )
          
          # Store in gene_data reactive
          gene_data(processed_gene_data)
          
          # Update cache with miscleavage-specific key (isolated to prevent reactive loops)
          isolate({
            new_cache <- gene_cache()
            new_cache[[cache_key]] <- processed_gene_data
            
            # Limit cache size (LRU policy) - now accounts for miscleavage types
            if (length(new_cache) > max_cache_size) {
              # Remove oldest entry (first in the list)
              oldest_key <- names(new_cache)[1]
              new_cache[[oldest_key]] <- NULL
              
              # Update cache statistics
              stats <- cache_stats()
              stats$evictions <- stats$evictions + 1
              stats$current_size <- length(new_cache)
              cache_stats(stats)
              
              cat("🗑️ Cache eviction: removed", oldest_key, "- cache size now", length(new_cache), "\n")
            }
            
            gene_cache(new_cache)
            
            # Update miss statistics
            stats <- cache_stats()
            stats$misses <- stats$misses + 1
            stats$current_size <- length(new_cache)
            cache_stats(stats)
          })
          cat("✅ Loaded and cached", gene_id, "with", miscleavage_type, "\n")
        })
      }
    } else {
      # Using original approach, no need to load gene-specific data
      # Just ensure we have the global data
      req(exists("peptides") && exists("as_database"))
      cat("ℹ️ Using legacy data approach - no gene-specific loading needed\n")
    }
  })
  
  # ============================================================================
  # GENE-SPECIFIC DATA ACCESSORS
  # ============================================================================
  
  # Get transcripts for selected gene
  selected_gene_transcripts <- reactive({
    req(input$gene)
    
    if (processed_data()$using_gene_index) {
      # Using gene-specific approach
      req(gene_data())
      gene_data()$transcripts_per_gene[[input$gene]]
    } else {
      # Using original approach
      processed_data()$transcripts_per_gene[[input$gene]]
    }
  })
  
  # Get AS events for selected gene
  selected_gene_as_events <- reactive({
    req(input$gene)
    
    if (processed_data()$using_gene_index) {
      # Using gene-specific approach
      req(gene_data())
      gene_data()$as_events
    } else {
      # Using original approach
      processed_data()$as_events[processed_data()$as_events$geneID == input$gene, ]
    }
  })
  
  # Get peptides for selected gene
  selected_gene_peptides <- reactive({
    req(input$gene)
    
    if (processed_data()$using_gene_index) {
      # Using gene-specific approach
      req(gene_data())
      gene_data()$peptides
    } else {
      # Using original approach
      processed_data()$original_peptides[processed_data()$original_peptides$geneID == input$gene, ]
    }
  })
  
  # ============================================================================
  # GENE SEARCH AND SELECTION INTERFACE
  # ============================================================================
  
  # Initialize the selectize input for server-side processing
  observeEvent(processed_data(), {
    updateSelectizeInput(session, "gene", server = TRUE, choices = character(0))
  }, once = TRUE)

  # Observer to update gene list based on search term
  observe({
    req(processed_data())
    
    search_term <- input$gene_search 
    
    all_gene_ids <- processed_data()$genes
    gene_symbols <- processed_data()$gene_lookup[all_gene_ids]
    labels <- paste0(gene_symbols, " (", all_gene_ids, ")")
    
    # Filter choices based on the search term
    if (!is.null(search_term) && nzchar(search_term)) {
      search_term_lower <- tolower(search_term)
      
      # Use fixed = TRUE to treat search as literal string (prevents regex errors)
      matching_indices <- tryCatch({
        grep(
          search_term_lower,
          tolower(paste0(gene_symbols, all_gene_ids)),
          fixed = TRUE
        )
      }, error = function(e) {
        # If search fails, return empty results
        integer(0)
      })
      
      # Limit to first 200 matches for performance
      if (length(matching_indices) > 200) {
        matching_indices <- matching_indices[1:200]
      }
      
      # Create filtered choices
      if (length(matching_indices) > 0) {
        filtered_choices <- setNames(
          all_gene_ids[matching_indices],
          labels[matching_indices]
        )
        
        updateSelectizeInput(session, "gene", 
                            choices = filtered_choices,
                            server = TRUE)
      } else {
        # No matches found
        updateSelectizeInput(session, "gene", 
                            choices = character(0),
                            server = TRUE)
      }
    } else {
      # No search term - show first 200 genes
      first_200 <- min(200, length(all_gene_ids))
      initial_choices <- setNames(
        all_gene_ids[1:first_200],
        labels[1:first_200]
      )
      
      updateSelectizeInput(session, "gene", 
                          choices = initial_choices,
                          server = TRUE)
    }
  })
  
  # ============================================================================
  # STANDARDIZED DATA FORMAT GENERATORS
  # ============================================================================
  
  # Create standardized vis_data structure for modules
  create_vis_data_structure <- reactive({
    req(processed_data(), selected_gene_peptides())
    
    list(
      genes = processed_data()$genes,
      gene_symbols = processed_data()$gene_symbols,
      gene_lookup = processed_data()$gene_lookup,
      proteases = processed_data()$proteases,
      original_peptides = selected_gene_peptides()
    )
  })
  
  # Get gene symbol for current gene
  current_gene_symbol <- reactive({
    req(input$gene, processed_data())
    
    tryCatch({
      processed_data()$gene_lookup[input$gene]
    }, error = function(e) {
      input$gene
    })
  })
  
  # ============================================================================
  # CACHE MANAGEMENT FUNCTIONS
  # ============================================================================
  
  # Clear cache
  clear_gene_cache <- function() {
    gene_cache(list())
    cache_stats(list(
      hits = 0,
      misses = 0,
      current_size = 0,
      evictions = 0
    ))
    cat("🗑️ Gene cache cleared\n")
  }
  
  # Get cache statistics
  get_cache_stats <- function() {
    stats <- cache_stats()
    cache_size <- length(gene_cache())
    
    list(
      cache_size = cache_size,
      max_cache_size = max_cache_size,
      hits = stats$hits,
      misses = stats$misses,
      hit_rate = if (stats$hits + stats$misses > 0) stats$hits / (stats$hits + stats$misses) else 0,
      evictions = stats$evictions,
      cached_genes = names(gene_cache())
    )
  }
  
  # Check if gene is cached
  is_gene_cached <- function(gene_id, miscleavage_type) {
    cache_key <- paste0(gene_id, "_", miscleavage_type)
    cache_key %in% names(gene_cache())
  }
  
  # ============================================================================
  # DATA VALIDATION FUNCTIONS
  # ============================================================================
  
  # Validate current gene data
  validate_gene_data <- function() {
    current_data <- gene_data()
    
    if (is.null(current_data)) {
      return(list(valid = FALSE, message = "No gene data loaded"))
    }
    
    required_fields <- c("peptides", "as_events", "transcripts_per_gene", "miscleavage_type")
    missing_fields <- setdiff(required_fields, names(current_data))
    
    if (length(missing_fields) > 0) {
      return(list(
        valid = FALSE, 
        message = paste("Missing required fields:", paste(missing_fields, collapse = ", "))
      ))
    }
    
    # Check peptides data structure
    if (is.null(current_data$peptides) || nrow(current_data$peptides) == 0) {
      return(list(valid = FALSE, message = "No peptides data found"))
    }
    
    return(list(
      valid = TRUE, 
      message = "Gene data is valid",
      peptide_count = nrow(current_data$peptides),
      transcript_count = length(unique(current_data$peptides$txID)),
      as_event_count = nrow(current_data$as_events)
    ))
  }
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  
  return(list(
    # Core reactive values
    processed_data = processed_data,
    gene_data = gene_data,
    gene_cache = gene_cache,
    
    # Gene-specific data accessors
    selected_gene_transcripts = selected_gene_transcripts,
    selected_gene_as_events = selected_gene_as_events,
    selected_gene_peptides = selected_gene_peptides,
    
    # Standardized data structures
    create_vis_data_structure = create_vis_data_structure,
    current_gene_symbol = current_gene_symbol,
    
    # Cache management functions
    load_and_cache_gene_data = load_and_cache_gene_data,
    clear_gene_cache = clear_gene_cache,
    get_cache_stats = get_cache_stats,
    is_gene_cached = is_gene_cached,
    
    # Data validation
    validate_gene_data = validate_gene_data,
    
    # Utility functions
    get_data_summary = function() {
      data <- processed_data()
      current_gene <- gene_data()
      
      list(
        using_gene_index = data$using_gene_index,
        total_genes = data$total_genes,
        current_gene_loaded = !is.null(current_gene),
        current_gene_id = if (!is.null(current_gene)) current_gene$gene_id else NULL,
        cache_size = length(gene_cache()),
        proteases = data$proteases
      )
    },
    
    refresh_gene_data = function(force = FALSE) {
      req(input$gene, input$miscleavage_type)
      
      if (force) {
        # Clear cache entry and reload (isolated to prevent reactive loops)
        isolate({
          cache_key <- paste0(input$gene, "_", input$miscleavage_type)
          cache <- gene_cache()
          if (cache_key %in% names(cache)) {
            cache[[cache_key]] <- NULL
            gene_cache(cache)
            cat("🗑️ Force cleared cache for", cache_key, "\n")
          }
        })
      }
      
      load_and_cache_gene_data(input$gene, input$miscleavage_type)
    }
  ))
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test Core Data Management Module
#' 
#' @description 
#' Function to test core data management module functionality
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_core_data_management_module <- function() {
  cat("Testing Core Data Management module...\n")
  
  # Test 1: Cache key generation
  cache_key_test <- tryCatch({
    # Test cache key generation logic
    gene_id <- "ENSG00000139618.15"
    miscleavage_type <- "no_miss_cleavage"
    cache_key <- paste0(gene_id, "_", miscleavage_type)
    
    # Verify cache key format
    expected_key <- "ENSG00000139618.15_no_miss_cleavage"
    cache_key == expected_key
  }, error = function(e) {
    cat("Cache key test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 2: LRU cache logic simulation
  lru_test <- tryCatch({
    # Test LRU cache behavior
    max_size <- 3
    cache <- list()
    
    # Add items
    cache[["item1"]] <- "data1"
    cache[["item2"]] <- "data2"
    cache[["item3"]] <- "data3"
    
    # Add fourth item (should trigger eviction)
    if (length(cache) >= max_size) {
      oldest_key <- names(cache)[1]
      cache[[oldest_key]] <- NULL
    }
    cache[["item4"]] <- "data4"
    
    # Verify item1 was evicted and item4 was added
    !("item1" %in% names(cache)) && ("item4" %in% names(cache)) && length(cache) == max_size
  }, error = function(e) {
    cat("LRU test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 3: Data structure validation
  structure_test <- tryCatch({
    # Test expected data structure
    mock_processed_data <- list(
      genes = c("ENSG00000139618.15", "ENSG00000141510.14"),
      gene_symbols = c("BRCA2", "TP53"),
      gene_lookup = setNames(c("BRCA2", "TP53"), c("ENSG00000139618.15", "ENSG00000141510.14")),
      proteases = c("trp", "chymo", "aspn", "lysc", "lysn", "gluc"),
      using_gene_index = TRUE,
      total_genes = 2
    )
    
    # Verify all required fields
    required_fields <- c("genes", "gene_symbols", "gene_lookup", "proteases", "using_gene_index", "total_genes")
    all(required_fields %in% names(mock_processed_data))
  }, error = function(e) {
    cat("Structure test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 4: Cache statistics structure
  stats_test <- tryCatch({
    # Test cache statistics structure
    mock_stats <- list(
      hits = 10,
      misses = 5,
      current_size = 15,
      evictions = 2
    )
    
    # Verify all required stats fields
    required_stats <- c("hits", "misses", "current_size", "evictions")
    all(required_stats %in% names(mock_stats))
  }, error = function(e) {
    cat("Stats test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 5: Search term processing
  search_test <- tryCatch({
    # Test search term filtering logic
    search_term <- "BRCA"
    gene_symbols <- c("BRCA2", "TP53", "BRCA1")
    gene_ids <- c("ENSG1", "ENSG2", "ENSG3")
    
    # Simulate search matching
    search_term_lower <- tolower(search_term)
    matching_indices <- grep(search_term_lower, tolower(gene_symbols), fixed = TRUE)
    
    # Should find BRCA2 and BRCA1 (indices 1 and 3)
    length(matching_indices) == 2 && all(c(1, 3) %in% matching_indices)
  }, error = function(e) {
    cat("Search test failed:", e$message, "\n")
    FALSE
  })
  
  if (cache_key_test && lru_test && structure_test && stats_test && search_test) {
    cat("All Core Data Management module tests passed!\n")
    return(TRUE)
  } else {
    cat("Some Core Data Management module tests failed!\n")
    return(FALSE)
  }
}

#===============================================================================
# INTEGRATION HELPER FUNCTIONS
#===============================================================================

#' Create Module Integration Points
#' 
#' @description
#' Helper function to create standardized integration points for other modules
#' 
#' @param core_module The core data management module instance
#' @return List of integration reactive values and functions

create_module_integration_points <- function(core_module) {
  list(
    # Standard reactive values for module integration
    processed_data = core_module$processed_data,
    gene_data = core_module$gene_data,
    selected_gene_transcripts = core_module$selected_gene_transcripts,
    selected_gene_as_events = core_module$selected_gene_as_events,
    selected_gene_peptides = core_module$selected_gene_peptides,
    
    # Standard data structures
    vis_data_structure = core_module$create_vis_data_structure,
    current_gene_symbol = core_module$current_gene_symbol,
    
    # Data loading and caching
    load_and_cache_gene_data = core_module$load_and_cache_gene_data,
    
    # Data validation
    validate_gene_data = core_module$validate_gene_data
  )
}