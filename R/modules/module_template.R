#===============================================================================
# MODULE TEMPLATE - Copy this for all new modules
# Replace MODULE_NAME with your specific module name
#===============================================================================

#' MODULE_NAME Module
#' 
#' @description 
#' Brief description of what this module does
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' @param ... Additional parameters passed from main server
#' 
#' @return Named list of reactive values/functions to be used by main server
#' 
#' @note 
#' - All reactive values should be clearly documented
#' - Dependencies should be listed in the function documentation
#' - Use consistent naming conventions (snake_case for functions, camelCase for inputs)

create_MODULE_NAME_module <- function(input, output, session, ...) {
  
  # ============================================================================
  # MODULE DEPENDENCIES
  # ============================================================================
  # List all dependencies passed in from main server via ...
  # Example:
  # gene_data <- list(...)$gene_data
  # processed_data <- list(...)$processed_data
  
  # ============================================================================
  # MODULE REACTIVE VALUES 
  # ============================================================================
  # Create any reactive values specific to this module
  module_reactive_values <- reactiveValues(
    # example_data = NULL,
    # processing_status = "idle"
  )
  
  # ============================================================================
  # MODULE REACTIVE EXPRESSIONS
  # ============================================================================
  # Create reactive expressions for data processing
  
  # example_reactive <- reactive({
  #   req(input$some_input)
  #   # Processing logic here
  #   return(processed_result)
  # })
  
  # ============================================================================
  # MODULE OBSERVERS
  # ============================================================================
  # Create observers for handling user interactions
  
  # observeEvent(input$some_button, {
  #   # Button click logic here
  # })
  
  # ============================================================================
  # MODULE OUTPUTS
  # ============================================================================
  # Create output renderers
  
  # output$some_output <- renderSomeType({
  #   # Rendering logic here
  # })
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  # Return reactive values/functions that main server needs
  return(list(
    # Reactive values
    # module_data = module_reactive_values,
    
    # Reactive expressions  
    # processed_data = example_reactive,
    
    # Functions
    # utility_function = function() { ... }
  ))
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test MODULE_NAME Module
#' 
#' @description 
#' Function to test module functionality in isolation
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_MODULE_NAME_module <- function() {
  # Add module-specific tests here
  cat("Testing MODULE_NAME module...\n")
  
  # Test 1: Basic functionality
  # Test 2: Edge cases  
  # Test 3: Error handling
  
  cat("All MODULE_NAME module tests passed!\n")
  return(TRUE)
} 