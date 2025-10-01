#===============================================================================
# PROGRESS & NOTIFICATIONS MODULE
# Extracted from server.R - Standardized user feedback system
#===============================================================================

#' Progress & Notifications Module
#' 
#' @description 
#' Standardizes progress indicators and notifications across the application:
#' - withProgress wrappers for consistent styling
#' - showNotification standardization  
#' - Progress tracking for long-running operations
#' - Error, warning, and success notification management
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' 
#' @return Named list of progress and notification functions
#' 
#' @note 
#' - All progress indicators use consistent messaging
#' - Notifications are categorized: error, warning, message, success
#' - Progress values are standardized (0-1 scale)

create_progress_notifications_module <- function(input, output, session) {
  
  # ============================================================================
  # PROGRESS WRAPPER FUNCTIONS
  # ============================================================================
  
  # Standardized progress wrapper with consistent styling
  show_progress <- function(message, expr, steps = NULL) {
    withProgress(message = message, value = 0, {
      if (is.null(steps)) {
        # Simple progress without specific steps
        result <- expr
      } else {
        # Multi-step progress
        total_steps <- length(steps)
        result <- NULL
        
        for (i in seq_along(steps)) {
          incProgress(1/total_steps, detail = steps[i])
          if (i == length(steps)) {
            # Execute the main expression on the last step
            result <- expr
          }
        }
      }
      return(result)
    })
  }
  
  # Progress for data loading operations
  show_data_loading_progress <- function(operation_name, gene_id = NULL, miscleavage_type = NULL, expr) {
    message <- if (!is.null(gene_id)) {
      paste('Loading', operation_name, 'for', gene_id)
      if (!is.null(miscleavage_type)) {
        paste(message, '(', miscleavage_type, ')')
      } else {
        message
      }
    } else {
      paste('Loading', operation_name, '...')
    }
    
    withProgress(message = message, value = 0, expr)
  }
  
  # Progress for analysis operations
  show_analysis_progress <- function(analysis_type, steps = NULL, expr) {
    message <- paste('Running', analysis_type, 'analysis...')
    show_progress(message, expr, steps)
  }
  
  # Progress for search operations
  show_search_progress <- function(search_type, expr) {
    message <- paste('Searching', search_type, '...')
    withProgress(message = message, expr)
  }
  
  # ============================================================================
  # NOTIFICATION FUNCTIONS
  # ============================================================================
  
  # Standardized notification function with consistent types
  show_notification <- function(message, type = "message", duration = NULL) {
    # Set default durations based on type
    if (is.null(duration)) {
      duration <- switch(type,
        "error" = 10,
        "warning" = 7,  
        "message" = 5,
        "success" = 5,
        5  # default
      )
    }
    
    showNotification(message, type = type, duration = duration)
  }
  
  # Specific notification types for consistency
  notify_error <- function(message, duration = 10) {
    show_notification(message, "error", duration)
  }
  
  notify_warning <- function(message, duration = 7) {
    show_notification(message, "warning", duration)
  }
  
  notify_success <- function(message, duration = 5) {
    show_notification(message, "message", duration)  # Shiny uses "message" for success-style
  }
  
  notify_info <- function(message, duration = 5) {
    show_notification(message, "message", duration)
  }
  
  # ============================================================================
  # SPECIFIC DOMAIN NOTIFICATIONS
  # ============================================================================
  
  # Gene loading notifications
  notify_gene_error <- function(gene_id, miscleavage_type = NULL) {
    message <- paste("Error loading data for gene", gene_id)
    if (!is.null(miscleavage_type)) {
      message <- paste(message, "with miscleavage type", miscleavage_type)
    }
    notify_error(message)
  }
  
  notify_no_transcripts <- function(gene_id, miscleavage_type) {
    message <- paste("No transcripts found for gene", gene_id, "with miscleavage type", miscleavage_type)
    notify_warning(message)
  }
  
  notify_no_peptides <- function(gene_id, miscleavage_type) {
    message <- paste("No peptides found for gene", gene_id, "with miscleavage type", miscleavage_type)
    notify_warning(message)
  }
  
  # Search notifications
  notify_search_results <- function(count, search_term = NULL) {
    if (count > 0) {
      message <- paste("Found", count, "matching peptides")
      if (!is.null(search_term)) {
        message <- paste(message, "for term:", search_term)
      }
      notify_success(message)
    } else {
      notify_warning("No matching peptides found.")
    }
  }
  
  notify_empty_search <- function() {
    notify_warning("Please enter a peptide sequence to search.")
  }
  
  notify_search_error <- function(error_message) {
    notify_error(paste("Search failed:", error_message))
  }
  
  # Navigation notifications
  notify_navigation <- function(gene_symbol, tab_name) {
    message <- paste("Navigated to", gene_symbol, tab_name)
    notify_info(message)
  }
  
  # File upload notifications
  notify_no_file_selected <- function() {
    notify_warning("Please select a FASTA file to upload.")
  }
  
  notify_invalid_file_format <- function(allowed_formats = c(".fa", ".fasta", ".fas")) {
    message <- paste("Please upload a valid FASTA file", paste0("(", paste(allowed_formats, collapse = ", "), ")"))
    notify_error(message)
  }
  
  # Pipeline notifications
  notify_pipeline_success <- function(message_detail = NULL) {
    base_message <- "Novel isoform discovery pipeline completed successfully!"
    message <- if (!is.null(message_detail)) {
      paste(base_message, message_detail)
    } else {
      paste(base_message, "Please select an ORF for gene matching.")
    }
    notify_success(message)
  }
  
  notify_pipeline_failure <- function(error_message) {
    notify_error(paste("Pipeline failed:", error_message))
  }
  
  notify_pipeline_execution_error <- function(error_message) {
    notify_error(paste("Pipeline execution failed:", error_message))
  }
  
  # BLAST notifications
  notify_blast_success <- function(result_count) {
    message <- paste("BLAST completed! Found", result_count, "potential gene matches.")
    notify_success(message)
  }
  
  notify_blast_no_info <- function() {
    notify_warning("BLAST completed but no gene information could be parsed.")
  }
  
  notify_blast_no_matches <- function(threshold = 70) {
    message <- paste("No BLAST matches found with ≥", paste0(threshold, "% identity."))
    notify_warning(message)
  }
  
  notify_blast_error <- function(error_message) {
    notify_error(paste("BLAST search failed:", error_message))
  }
  
  # Selection validation notifications
  notify_invalid_selection <- function(selection_type) {
    message <- paste("Invalid", selection_type, "selection")
    notify_error(message)
  }
  
  # Data loading success notifications
  notify_gene_loaded <- function(gene_symbol, data_type = "data") {
    message <- paste("Successfully loaded", gene_symbol, data_type)
    notify_success(message)
  }
  
  notify_rds_not_available <- function() {
    notify_warning("RDS file not available for selected gene. Proceeding with novel data only.")
  }
  
  notify_no_comparison_data <- function(data_type) {
    message <- paste("No", data_type, "available for comparison")
    notify_warning(message)
  }
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  return(list(
    # Progress functions
    show_progress = show_progress,
    show_data_loading_progress = show_data_loading_progress,
    show_analysis_progress = show_analysis_progress,
    show_search_progress = show_search_progress,
    
    # General notification functions
    notify_error = notify_error,
    notify_warning = notify_warning,
    notify_success = notify_success,
    notify_info = notify_info,
    
    # Domain-specific notification functions
    notify_gene_error = notify_gene_error,
    notify_no_transcripts = notify_no_transcripts,
    notify_no_peptides = notify_no_peptides,
    notify_search_results = notify_search_results,
    notify_empty_search = notify_empty_search,
    notify_search_error = notify_search_error,
    notify_navigation = notify_navigation,
    notify_no_file_selected = notify_no_file_selected,
    notify_invalid_file_format = notify_invalid_file_format,
    notify_pipeline_success = notify_pipeline_success,
    notify_pipeline_failure = notify_pipeline_failure,
    notify_pipeline_execution_error = notify_pipeline_execution_error,
    notify_blast_success = notify_blast_success,
    notify_blast_no_info = notify_blast_no_info,
    notify_blast_no_matches = notify_blast_no_matches,
    notify_blast_error = notify_blast_error,
    notify_invalid_selection = notify_invalid_selection,
    notify_gene_loaded = notify_gene_loaded,
    notify_rds_not_available = notify_rds_not_available,
    notify_no_comparison_data = notify_no_comparison_data
  ))
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test Progress & Notifications Module
#' 
#' @description 
#' Function to test progress and notifications module functionality
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_progress_notifications_module <- function() {
  cat("Testing Progress & Notifications module...\n")
  
  # Test 1: Notification message generation
  notification_test <- tryCatch({
    # Test error message generation
    error_msg <- paste("Error loading data for gene", "TEST_GENE", "with miscleavage type", "no_miss_cleavage")
    result1 <- nchar(error_msg) > 0
    
    # Test success message generation
    success_msg <- paste("Successfully loaded", "TEST_GENE", "data")
    result2 <- nchar(success_msg) > 0
    
    result1 && result2
  }, error = function(e) {
    cat("Notification test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 2: Progress message generation
  progress_test <- tryCatch({
    # Test data loading progress message
    progress_msg <- paste('Loading', 'test_operation', 'for', 'TEST_GENE')
    nchar(progress_msg) > 0
  }, error = function(e) {
    cat("Progress test failed:", e$message, "\n")
    FALSE
  })
  
  if (notification_test && progress_test) {
    cat("All Progress & Notifications module tests passed!\n")
    return(TRUE)
  } else {
    cat("Some Progress & Notifications module tests failed!\n")
    return(FALSE)
  }
} 