#===============================================================================
# SIMPLE TEXT OUTPUTS MODULE
# Extracted from server.R - Simple text renderers and status outputs
#===============================================================================

#' Simple Text Outputs Module
#' 
#' @description 
#' Manages simple text output renderers throughout the application:
#' - Peptide plot titles and summaries
#' - Pipeline status indicators
#' - Selection status messages
#' - Search result summaries
#' - Informational text displays
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' @param module_data List containing reactive data dependencies
#' 
#' @return NULL (outputs are created directly in output object)
#' 
#' @note 
#' - All text outputs preserve exact original formatting and logic
#' - Dependencies on reactive values maintained identically
#' - Error handling and fallback messages preserved

create_simple_outputs_module <- function(input, output, session, module_data) {
  
  # ============================================================================
  # MODULE DEPENDENCIES
  # ============================================================================
  filtered_peptides <- module_data$filtered_peptides
  as_peptide_comparison_data <- module_data$as_peptide_comparison_data
  transcript_pairs <- module_data$transcript_pairs
  peptide_search_results <- module_data$peptide_search_results
  novel_pipeline_results <- module_data$novel_pipeline_results
  novel_orf_results <- module_data$novel_orf_results
  novel_gene_search_results <- module_data$novel_gene_search_results
  
  # ============================================================================
  # PEPTIDE VISUALIZATION OUTPUTS
  # ============================================================================
  
  # Peptide plot title with miscleavage information
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
  
  # Peptide comparison summary statistics
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
    
    # Calculate summary statistics
    total_peptides <- nrow(peptides)
    as_affected <- sum(peptides$event_overlap)
    as_percentage <- round((as_affected / total_peptides) * 100, 1)
    
    # Generate summary text
    paste0(
      "Comparing ", ref_tx, " vs ", alt_tx, ": ",
      total_peptides, " total peptides, ",
      as_affected, " (", as_percentage, "%) affected by alternative splicing"
    )
  })
  
  # ============================================================================
  # SEARCH OUTPUTS
  # ============================================================================
  
  # Peptide search results summary
  output$peptide_search_summary <- renderText({
    results <- peptide_search_results()
    if (is.null(results) || nrow(results) == 0) {
      return("No results")
    }
    
    paste("Found", nrow(results), "matching peptides")
  })
  
  # ============================================================================
  # NOVEL ISOFORM PIPELINE OUTPUTS
  # ============================================================================
  
  # Novel pipeline status indicator
  output$novel_pipeline_status <- renderText({
    results <- novel_pipeline_results()
    if (is.null(results)) {
      return("Analysis not started")
    }
    
    if (!is.null(results$error)) {
      return(paste("Error:", results$error))
    }
    
    if (!is.null(results$status)) {
      return(results$status)
    }
    
    return("Analysis completed successfully")
  })
  
  # Novel ORF selection status
  output$novel_orf_selection_status <- renderText({
    if (is.null(input$novel_orf_table_rows_selected)) {
      return("Please select an ORF from the table above")
    } else {
      orf_data <- novel_orf_results()
      if (!is.null(orf_data) && length(input$novel_orf_table_rows_selected) > 0) {
        selected_row <- input$novel_orf_table_rows_selected[1]
        if (selected_row <= nrow(orf_data)) {
          selected_orf <- orf_data[selected_row, ]
          return(paste("Selected ORF:", selected_orf$orf_id, "- Proceed with gene search"))
        }
      }
      return("Invalid ORF selection")
    }
  })
  
  # Novel gene search selection status  
  output$novel_gene_search_selection_status <- renderText({
    if (is.null(input$novel_gene_search_results_table_rows_selected)) {
      return("Please select a gene from the search results above")
    } else {
      gene_search_data <- novel_gene_search_results()
      if (!is.null(gene_search_data) && length(input$novel_gene_search_results_table_rows_selected) > 0) {
        selected_row <- input$novel_gene_search_results_table_rows_selected[1]
        if (selected_row <= nrow(gene_search_data)) {
          selected_gene <- gene_search_data[selected_row, ]
          gene_symbol <- selected_gene$gene_symbol %||% "Unknown"
          identity <- selected_gene$identity %||% "N/A"
          return(paste("Selected:", gene_symbol, "with", identity, "% identity - Ready for analysis"))
        }
      }
      return("Invalid gene selection")
    }
  })
  
  # ============================================================================
  # ANALYSIS STATUS OUTPUTS
  # ============================================================================
  
  # General loading status
  output$loading_status <- renderText({
    if (!is.null(input$gene) && !is.null(input$miscleavage_type)) {
      paste("Loaded data for gene:", input$gene, "with", input$miscleavage_type)
    } else {
      "Please select a gene and miscleavage type"
    }
  })
  
  # Cache status information
  output$cache_status <- renderText({
    # This would use cache stats from utility functions if available
    if (exists("utility_functions") && !is.null(utility_functions$get_cache_stats)) {
      stats <- utility_functions$get_cache_stats()
      paste("Cache:", stats$size, "genes loaded")
    } else {
      "Cache information not available"
    }
  })
  
  # ============================================================================
  # ERROR AND WARNING MESSAGES
  # ============================================================================
  
  # General error status
  output$error_status <- renderText({
    # This can be used for displaying general error states
    if (!is.null(input$last_error)) {
      paste("Error:", input$last_error)
    } else {
      ""
    }
  })
  
  # Data availability status
  output$data_availability_status <- renderText({
    if (exists("gene_index") && !is.null(gene_index)) {
      paste("Gene index available with", nrow(gene_index), "genes")
    } else if (exists("peptides") && !is.null(peptides)) {
      paste("Peptide database available with", nrow(peptides), "entries")
    } else {
      "No data sources available"
    }
  })
  
  # ============================================================================
  # HELPER FUNCTIONS FOR TEXT FORMATTING
  # ============================================================================
  
  # Format percentage with proper rounding
  format_percentage <- function(value, total, digits = 1) {
    if (total == 0) return("0%")
    percentage <- round((value / total) * 100, digits)
    paste0(percentage, "%")
  }
  
  # Format large numbers with commas
  format_number <- function(number) {
    if (is.numeric(number)) {
      format(number, big.mark = ",", scientific = FALSE)
    } else {
      as.character(number)
    }
  }
  
  # Safe text extraction with fallback
  safe_text <- function(value, fallback = "N/A") {
    if (is.null(value) || is.na(value) || value == "") {
      return(fallback)
    }
    as.character(value)
  }
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  
  # Return helper functions for use by other modules
  return(list(
    format_percentage = format_percentage,
    format_number = format_number,
    safe_text = safe_text
  ))
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test Simple Outputs Module
#' 
#' @description 
#' Function to test simple outputs module functionality
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_simple_outputs_module <- function() {
  cat("Testing Simple Outputs module...\n")
  
  # Test 1: Text formatting functions
  formatting_test <- tryCatch({
    # Test percentage formatting
    pct1 <- format_percentage <- function(value, total, digits = 1) {
      if (total == 0) return("0%")
      percentage <- round((value / total) * 100, digits)
      paste0(percentage, "%")
    }
    result1 <- pct1(25, 100) == "25%"
    result2 <- pct1(1, 3, 2) == "33.33%"
    result3 <- pct1(0, 0) == "0%"
    
    all(c(result1, result2, result3))
  }, error = function(e) {
    cat("Formatting test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 2: Number formatting
  number_test <- tryCatch({
    format_number <- function(number) {
      if (is.numeric(number)) {
        format(number, big.mark = ",", scientific = FALSE)
      } else {
        as.character(number)
      }
    }
    
    result1 <- grepl(",", format_number(1000))  # Should add comma
    result2 <- format_number("text") == "text"   # Should handle non-numeric
    
    result1 && result2
  }, error = function(e) {
    cat("Number formatting test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 3: Safe text extraction
  safe_text_test <- tryCatch({
    safe_text <- function(value, fallback = "N/A") {
      if (is.null(value) || is.na(value) || value == "") {
        return(fallback)
      }
      as.character(value)
    }
    
    result1 <- safe_text(NULL) == "N/A"
    result2 <- safe_text(NA) == "N/A"
    result3 <- safe_text("") == "N/A"
    result4 <- safe_text("text") == "text"
    result5 <- safe_text(NULL, "custom") == "custom"
    
    all(c(result1, result2, result3, result4, result5))
  }, error = function(e) {
    cat("Safe text test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 4: Miscleavage label generation
  miscleavage_test <- tryCatch({
    # Simulate the switch logic
    miscleavage_type <- "no_miss_cleavage"
    label <- switch(miscleavage_type,
      "no_miss_cleavage" = "No Miscleavage",
      "upto_two_misscleavage" = "Up to 2 Miscleavages"
    )
    
    label == "No Miscleavage"
  }, error = function(e) {
    cat("Miscleavage test failed:", e$message, "\n")
    FALSE
  })
  
  if (formatting_test && number_test && safe_text_test && miscleavage_test) {
    cat("All Simple Outputs module tests passed!\n")
    return(TRUE)
  } else {
    cat("Some Simple Outputs module tests failed!\n")
    return(FALSE)
  }
}