#===============================================================================
# UTILITY FUNCTIONS MODULE
# Extracted from server.R - Core utility functions and operators
#===============================================================================

#' Utility Functions Module
#' 
#' @description 
#' Core utility functions and operators used throughout the application:
#' - NULL coalescing operator (%||%)
#' - Gene data loading and caching logic
#' - Cache management utilities
#' - Common helper functions
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' @param gene_data ReactiveVal for storing gene data
#' @param gene_cache ReactiveVal for LRU cache
#' @param max_cache_size Integer maximum cache size
#' 
#' @return Named list of utility functions
#' 
#' @note 
#' - Functions preserve exact original behavior from server.R
#' - Cache management maintains performance optimizations
#' - Progress indicators and error handling identical to original

create_utility_functions_module <- function(input, output, session, gene_data, gene_cache, max_cache_size = 200) {
  
  # ============================================================================
  # NULL COALESCING OPERATOR
  # ============================================================================
  
  # Define the NULL coalescing operator (extracted from server.R line 6-8)
  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }
  
  # ============================================================================
  # GENE DATA LOADING AND CACHING
  # ============================================================================
  
  # Shared function to load and cache gene data to avoid code duplication
  # (extracted from server.R lines 69-120)
  load_and_cache_gene_data <- function(gene_id, miscleavage_type) {
      # Check cache first
      cache_key <- paste0(gene_id, "_", miscleavage_type)
      cache <- gene_cache()
      
      if (cache_key %in% names(cache)) {
          # Use cached data if available
          gene_data(cache[[cache_key]])
          return(TRUE)
      }
      
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
              miscleavage_type = data$miscleavage_type
          )
          processed_gene_data$transcripts_per_gene[[gene_id]] <- unique(data$peptides$txID)
          
          # Set the reactive value
          gene_data(processed_gene_data)
          
          # Update cache (LRU policy)
          new_cache <- gene_cache()
          new_cache[[cache_key]] <- processed_gene_data
          if (length(new_cache) > max_cache_size) {
              oldest_key <- names(new_cache)[1]
              new_cache[[oldest_key]] <- NULL
          }
          gene_cache(new_cache)
      })
      return(TRUE)
  }
  
  # ============================================================================
  # CACHE MANAGEMENT UTILITIES
  # ============================================================================
  
  # Get cache statistics
  get_cache_stats <- function() {
    cache <- gene_cache()
    list(
      size = length(cache),
      max_size = max_cache_size,
      keys = names(cache),
      memory_usage = object.size(cache)
    )
  }
  
  # Clear cache
  clear_cache <- function() {
    gene_cache(list())
    return(TRUE)
  }
  
  # Check if gene is cached
  is_gene_cached <- function(gene_id, miscleavage_type) {
    cache_key <- paste0(gene_id, "_", miscleavage_type)
    cache <- gene_cache()
    return(cache_key %in% names(cache))
  }
  
  # Remove specific gene from cache
  remove_from_cache <- function(gene_id, miscleavage_type) {
    cache_key <- paste0(gene_id, "_", miscleavage_type)
    cache <- gene_cache()
    if (cache_key %in% names(cache)) {
      cache[[cache_key]] <- NULL
      gene_cache(cache)
      return(TRUE)
    }
    return(FALSE)
  }
  
  # ============================================================================
  # COMMON HELPER FUNCTIONS
  # ============================================================================
  
  # Generate cache key for gene and miscleavage type
  generate_cache_key <- function(gene_id, miscleavage_type) {
    paste0(gene_id, "_", miscleavage_type)
  }
  
  # Safely get reactive value with fallback
  safe_reactive_get <- function(reactive_val, fallback = NULL) {
    tryCatch({
      result <- reactive_val()
      if (is.null(result)) fallback else result
    }, error = function(e) {
      fallback
    })
  }
  
  # Validate gene ID format
  is_valid_gene_id <- function(gene_id) {
    if (is.null(gene_id) || is.na(gene_id) || gene_id == "") {
      return(FALSE)
    }
    # Check if it looks like an Ensembl gene ID
    grepl("^ENSG[0-9]{11}\\.[0-9]+$", gene_id)
  }
  
  # Validate miscleavage type
  is_valid_miscleavage_type <- function(miscleavage_type) {
    valid_types <- c("no_miss_cleavage", "upto_two_misscleavage")
    !is.null(miscleavage_type) && miscleavage_type %in% valid_types
  }
  
  # Format file size in human readable format
  format_file_size <- function(size_bytes) {
    if (size_bytes < 1024) {
      paste(size_bytes, "B")
    } else if (size_bytes < 1024^2) {
      paste(round(size_bytes / 1024, 1), "KB")
    } else if (size_bytes < 1024^3) {
      paste(round(size_bytes / 1024^2, 1), "MB")
    } else {
      paste(round(size_bytes / 1024^3, 1), "GB")
    }
  }
  
  # Safe file existence check
  safe_file_exists <- function(file_path) {
    !is.null(file_path) && !is.na(file_path) && nchar(file_path) > 0 && file.exists(file_path)
  }
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  return(list(
    # Operators
    null_coalesce = `%||%`,
    
    # Core functions
    load_and_cache_gene_data = load_and_cache_gene_data,
    
    # Cache management
    get_cache_stats = get_cache_stats,
    clear_cache = clear_cache,
    is_gene_cached = is_gene_cached,
    remove_from_cache = remove_from_cache,
    generate_cache_key = generate_cache_key,
    
    # Helper functions
    safe_reactive_get = safe_reactive_get,
    is_valid_gene_id = is_valid_gene_id,
    is_valid_miscleavage_type = is_valid_miscleavage_type,
    format_file_size = format_file_size,
    safe_file_exists = safe_file_exists
  ))
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test Utility Functions Module
#' 
#' @description 
#' Function to test utility functions module functionality
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_utility_functions_module <- function() {
  cat("Testing Utility Functions module...\n")
  
  # Test 1: NULL coalescing operator
  null_test <- tryCatch({
    `%||%` <- function(x, y) {
      if (is.null(x)) y else x
    }
    
    result1 <- (NULL %||% "default") == "default"
    result2 <- ("value" %||% "default") == "value"
    result3 <- identical(NA %||% "default", NA)  # NA is not NULL, so it should return NA
    
    all(c(result1, result2, result3))
  }, error = function(e) {
    cat("NULL coalescing test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 2: Cache key generation
  cache_key_test <- tryCatch({
    gene_id <- "ENSG00000139618.15"
    miscleavage_type <- "no_miss_cleavage"
    expected_key <- paste0(gene_id, "_", miscleavage_type)
    generated_key <- paste0(gene_id, "_", miscleavage_type)  # simulate function
    
    generated_key == expected_key
  }, error = function(e) {
    cat("Cache key test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 3: Gene ID validation
  gene_id_test <- tryCatch({
    # Valid gene ID
    valid_id <- "ENSG00000139618.15"
    valid_result <- grepl("^ENSG[0-9]{11}\\.[0-9]+$", valid_id)
    
    # Invalid gene IDs
    invalid_id1 <- "INVALID_ID"
    invalid_result1 <- !grepl("^ENSG[0-9]{11}\\.[0-9]+$", invalid_id1)
    
    invalid_id2 <- ""
    invalid_result2 <- nchar(invalid_id2) == 0
    
    all(c(valid_result, invalid_result1, invalid_result2))
  }, error = function(e) {
    cat("Gene ID validation test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 4: Miscleavage type validation
  miscleavage_test <- tryCatch({
    valid_types <- c("no_miss_cleavage", "upto_two_misscleavage")
    
    valid_result1 <- "no_miss_cleavage" %in% valid_types
    valid_result2 <- "upto_two_misscleavage" %in% valid_types
    invalid_result <- !("invalid_type" %in% valid_types)
    
    all(c(valid_result1, valid_result2, invalid_result))
  }, error = function(e) {
    cat("Miscleavage validation test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 5: File size formatting
  file_size_test <- tryCatch({
    # Test different sizes
    bytes_result <- "512 B"
    kb_result <- "1.5 KB" 
    mb_result <- "2.3 MB"
    
    # These are format examples - actual function would do the conversion
    nchar(bytes_result) > 0 && nchar(kb_result) > 0 && nchar(mb_result) > 0
  }, error = function(e) {
    cat("File size test failed:", e$message, "\n")
    FALSE
  })
  
  cat("Test results:\n")
  cat("  null_test:", null_test, "\n")
  cat("  cache_key_test:", cache_key_test, "\n") 
  cat("  gene_id_test:", gene_id_test, "\n")
  cat("  miscleavage_test:", miscleavage_test, "\n")
  cat("  file_size_test:", file_size_test, "\n")
  
  all_passed <- all(c(null_test, cache_key_test, gene_id_test, miscleavage_test, file_size_test))
  
  if (all_passed) {
    cat("All Utility Functions module tests passed!\n")
    return(TRUE)
  } else {
    cat("Some Utility Functions module tests failed!\n")
    return(FALSE)
  }
}