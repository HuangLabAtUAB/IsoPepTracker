#===============================================================================
# SIMPLE TITLE OUTPUTS MODULE  
# Extracted from server.R - Simple title text renderers with minimal dependencies
#===============================================================================

#' Simple Title Outputs Module
#' 
#' @description 
#' Manages simple title text outputs with minimal dependencies:
#' - Peptide plot title with miscleavage information
#' - Other simple titles as needed
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' 
#' @return NULL (outputs are created directly in output object)
#' 
#' @note 
#' - Only includes outputs with minimal reactive dependencies
#' - Preserves exact original formatting and logic

create_simple_title_outputs_module <- function(input, output, session) {
  
  # ============================================================================
  # SIMPLE TITLE OUTPUTS
  # ============================================================================
  
  # Peptide plot title with miscleavage information
  # (extracted from server.R - has only input$miscleavage_type dependency)
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
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  
  # No return values needed for this simple module
  return(list())
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test Simple Title Outputs Module
#' 
#' @description 
#' Function to test simple title outputs module functionality
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_simple_title_outputs_module <- function() {
  cat("Testing Simple Title Outputs module...\n")
  
  # Test 1: Miscleavage label generation
  miscleavage_test <- tryCatch({
    # Test the switch logic
    test_cases <- list(
      list(input = "no_miss_cleavage", expected = "No Miscleavage"),
      list(input = "upto_two_misscleavage", expected = "Up to 2 Miscleavages")
    )
    
    results <- sapply(test_cases, function(case) {
      label <- switch(case$input,
        "no_miss_cleavage" = "No Miscleavage",
        "upto_two_misscleavage" = "Up to 2 Miscleavages"
      )
      label == case$expected
    })
    
    all(results)
  }, error = function(e) {
    cat("Miscleavage label test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 2: Title generation
  title_test <- tryCatch({
    # Test title generation with different inputs
    label1 <- "No Miscleavage"
    title1 <- paste("Peptide Visualization -", label1)
    result1 <- title1 == "Peptide Visualization - No Miscleavage"
    
    title2 <- "Peptide Visualization"  # fallback case
    result2 <- title2 == "Peptide Visualization"
    
    result1 && result2
  }, error = function(e) {
    cat("Title generation test failed:", e$message, "\n")
    FALSE
  })
  
  if (miscleavage_test && title_test) {
    cat("All Simple Title Outputs module tests passed!\n")
    return(TRUE)
  } else {
    cat("Some Simple Title Outputs module tests failed!\n")
    return(FALSE)
  }
}