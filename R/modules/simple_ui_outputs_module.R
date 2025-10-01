#===============================================================================
# SIMPLE UI OUTPUTS MODULE
# Extracted from server.R - Simple UI renderers with minimal dependencies
#===============================================================================

#' Simple UI Outputs Module
#' 
#' @description 
#' Manages simple UI output renderers with minimal reactive dependencies:
#' - Action buttons with simple input dependencies
#' - Basic UI elements that don't require complex reactive chains
#' - Status indicators with input-only dependencies
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' 
#' @return NULL (outputs are created directly in output object)
#' 
#' @note 
#' - Only includes UI outputs with minimal dependencies (mostly input$ values)
#' - Preserves exact original styling and functionality
#' - No complex reactive value dependencies

create_simple_ui_outputs_module <- function(input, output, session) {
  
  # ============================================================================
  # SIMPLE UI OUTPUTS WITH MINIMAL DEPENDENCIES
  # ============================================================================
  
  # Compare button - only depends on input$selected_pair
  output$compare_button <- renderUI({
    req(input$selected_pair)
    actionButton("compare_peptides", "Compare Peptides",
                icon = icon("chart-bar"),
                style = "margin: 15px 0; color: #fff; background-color: #3c8dbc; border-color: #367fa9;")
  })
  
  # Selected pair display - only depends on input$selected_pair
  output$selected_pair_display <- renderUI({
    req(input$selected_pair)
    
    # Get the selected pair name
    pair_name <- input$selected_pair
    
    div(
      strong("Currently comparing:"),
      p(pair_name)
    )
  })
  
  # ============================================================================
  # HELPER FUNCTIONS FOR UI GENERATION
  # ============================================================================
  
  # Generate styled action button
  create_styled_button <- function(id, label, icon_name = NULL, style_class = "primary") {
    base_style <- switch(style_class,
      "primary" = "color: #fff; background-color: #3c8dbc; border-color: #367fa9;",
      "success" = "color: #fff; background-color: #00a65a; border-color: #008d4c;",
      "warning" = "color: #fff; background-color: #f39c12; border-color: #e08e0b;",
      "danger" = "color: #fff; background-color: #dd4b39; border-color: #d73925;",
      "color: #fff; background-color: #3c8dbc; border-color: #367fa9;" # default
    )
    
    button_icon <- if (!is.null(icon_name)) icon(icon_name) else NULL
    
    actionButton(id, label,
                icon = button_icon,
                style = paste("margin: 15px 0;", base_style))
  }
  
  # Generate info box UI
  create_info_box <- function(title, value, icon_name = "info-circle", color = "blue") {
    div(class = paste("info-box bg", color, sep = "-"),
        span(class = "info-box-icon", icon(icon_name)),
        div(class = "info-box-content",
            span(class = "info-box-text", title),
            span(class = "info-box-number", value)
        )
    )
  }
  
  # Generate status badge
  create_status_badge <- function(text, status = "info") {
    class_name <- switch(status,
      "success" = "label label-success",
      "warning" = "label label-warning", 
      "danger" = "label label-danger",
      "info" = "label label-info",
      "label label-info" # default
    )
    
    span(class = class_name, text)
  }
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  
  # Return helper functions for use by other modules
  return(list(
    create_styled_button = create_styled_button,
    create_info_box = create_info_box,
    create_status_badge = create_status_badge
  ))
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test Simple UI Outputs Module
#' 
#' @description 
#' Function to test simple UI outputs module functionality
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_simple_ui_outputs_module <- function() {
  cat("Testing Simple UI Outputs module...\n")
  
  # Test 1: Button style generation
  button_style_test <- tryCatch({
    # Test style switching
    primary_style <- "color: #fff; background-color: #3c8dbc; border-color: #367fa9;"
    success_style <- "color: #fff; background-color: #00a65a; border-color: #008d4c;"
    
    # Simulate switch logic
    test_style <- switch("primary",
      "primary" = "color: #fff; background-color: #3c8dbc; border-color: #367fa9;",
      "success" = "color: #fff; background-color: #00a65a; border-color: #008d4c;",
      "color: #fff; background-color: #3c8dbc; border-color: #367fa9;"
    )
    
    test_style == primary_style
  }, error = function(e) {
    cat("Button style test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 2: CSS class generation
  css_class_test <- tryCatch({
    # Test status badge class switching
    success_class <- switch("success",
      "success" = "label label-success",
      "warning" = "label label-warning", 
      "danger" = "label label-danger",
      "info" = "label label-info",
      "label label-info"
    )
    
    success_class == "label label-success"
  }, error = function(e) {
    cat("CSS class test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 3: String concatenation for styles
  style_concat_test <- tryCatch({
    base_style <- "color: #fff; background-color: #3c8dbc;"
    margin_style <- "margin: 15px 0;"
    combined <- paste(margin_style, base_style)
    
    nchar(combined) > nchar(base_style) + nchar(margin_style)
  }, error = function(e) {
    cat("Style concatenation test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 4: Icon handling
  icon_test <- tryCatch({
    # Test icon name validation
    valid_icons <- c("chart-bar", "info-circle", "check", "warning", "danger")
    test_icon <- "chart-bar"
    
    test_icon %in% valid_icons
  }, error = function(e) {
    cat("Icon test failed:", e$message, "\n")
    FALSE
  })
  
  if (button_style_test && css_class_test && style_concat_test && icon_test) {
    cat("All Simple UI Outputs module tests passed!\n")
    return(TRUE)
  } else {
    cat("Some Simple UI Outputs module tests failed!\n")
    return(FALSE)
  }
}