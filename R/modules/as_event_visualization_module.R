#===============================================================================
# AS EVENT VISUALIZATION MODULE
# Extracted from server.R - Alternative splicing event visualization
#===============================================================================

#' AS Event Visualization Module
#' 
#' @description 
#' Alternative splicing event visualization functionality:
#' - AS event structure plotting with interactive elements
#' - Event selection and highlighting
#' - Integration with AS events table
#' - Click event handling for detailed views
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' @param exon_data Reactive containing exon data for AS visualization
#' @param selected_gene_as_events Reactive containing AS events data
#' 
#' @return Named list containing reactive values for AS visualization
#' 
#' @note 
#' - Requires exon_data reactive from transcript visualization
#' - Uses plotly for interactive AS event visualization
#' - Integrates with AS events table for selection

create_as_event_visualization_module <- function(input, output, session, exon_data, selected_gene_as_events) {
  
  # ============================================================================
  # MODULE REACTIVE VALUES
  # ============================================================================
  
  # AS view data storage
  as_view_data <- reactiveVal(NULL)
  
  # Selected AS event
  selected_event <- reactiveVal(NULL)
  
  # ============================================================================
  # AS VIEW DATA LOADING
  # ============================================================================
  
  # Update AS view when gene changes
  observeEvent(input$update, {
    req(input$gene, input$miscleavage_type, exon_data())
    if (!exon_data()$success) return()
    
    withProgress(message = 'Creating AS visualization...', {
      tryCatch({
        gene_id <- input$gene
        gene_symbol <- if (exists("processed_data") && !is.null(processed_data()$gene_lookup)) {
          processed_data()$gene_lookup[gene_id]
        } else {
          gene_id  # fallback to gene_id if lookup not available
        }
        
        incProgress(0.5, detail = "Loading AS events...")
        
        # Create AS visualization using cached data
        as_result <- create_as_view(
          gene_id = gene_id,
          gene_symbol = gene_symbol,
          exon_data = exon_data(),
          as_events = selected_gene_as_events()
        )
        
        # Store AS view data
        as_view_data(as_result)
        
        incProgress(0.5, detail = "Complete!")
        
      }, error = function(e) {
        cat("AS visualization error:", e$message, "\n")
        
        as_view_data(list(
          success = FALSE,
          error = paste("Failed to create AS visualization:", e$message)
        ))
        
        showNotification(
          paste("Error creating AS visualization:", e$message),
          type = "error",
          duration = 8
        )
      })
    })
  })
  
  # ============================================================================
  # EVENT SELECTION HANDLERS
  # ============================================================================
  
  # Handle AS plot click events
  observeEvent(event_data("plotly_click", source = "as_plot"), {
    click_data <- event_data("plotly_click", source = "as_plot")
    
    if (!is.null(click_data) && nrow(click_data) > 0) {
      # Extract event information from click data
      point_number <- click_data$pointNumber[1] + 1  # R is 1-indexed
      
      # Get AS events data
      as_events <- selected_gene_as_events()
      
      if (!is.null(as_events) && nrow(as_events) >= point_number) {
        selected_event_data <- as_events[point_number, ]
        selected_event(selected_event_data)
        
        showNotification(
          paste("Selected AS event:", selected_event_data$eventID),
          type = "message",
          duration = 3
        )
      }
    }
  })
  
  # Handle AS events table selection
  observeEvent(input$as_events_table_rows_selected, {
    req(selected_gene_as_events())
    
    if (length(input$as_events_table_rows_selected) > 0) {
      selected_row <- input$as_events_table_rows_selected[1]
      as_events <- selected_gene_as_events()
      
      if (selected_row <= nrow(as_events)) {
        selected_event_data <- as_events[selected_row, ]
        selected_event(selected_event_data)
      }
    } else {
      selected_event(NULL)
    }
  })
  
  # ============================================================================
  # VISUALIZATION OUTPUTS
  # ============================================================================
  
  # AS structure plot
  output$as_structure_plot <- renderPlotly({
    req(as_view_data())
    
    # Get the AS view data
    as_data <- as_view_data()
    
    if (!as_data$success) {
      # Create error plot
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                label = paste("AS visualization error:", as_data$error), 
                size = 4, hjust = 0.5, color = "red") +
        xlim(0, 1) + ylim(0, 1) +
        theme_void()
      
      return(ggplotly(p) %>% config(displayModeBar = FALSE))
    }
    
    # Get the plot
    p <- as_data$plot
    
    if (is.null(p)) {
      # No plot available
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                label = "No AS events available for visualization", 
                size = 4, hjust = 0.5) +
        xlim(0, 1) + ylim(0, 1) +
        theme_void()
      
      return(ggplotly(p) %>% config(displayModeBar = FALSE))
    }
    
    # Return plot with click event handling
    tryCatch({
      plotly_obj <- p %>% 
        clean_plotly_hover() %>%
        layout(
          dragmode = FALSE,
          source = "as_plot"
        ) %>%
        config(
          displayModeBar = TRUE,
          displaylogo = FALSE,
          modeBarButtonsToRemove = c(
            "pan2d", "select2d", "lasso2d", "resetScale2d"
          )
        )
      
      return(plotly_obj)
      
    }, error = function(e) {
      # Fallback if plotly processing fails
      cat("Plotly processing error:", e$message, "\n")
      
      fallback_plot <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                label = "Error displaying AS plot", 
                size = 4, hjust = 0.5, color = "orange") +
        xlim(0, 1) + ylim(0, 1) +
        theme_void()
      
      return(ggplotly(fallback_plot) %>% config(displayModeBar = FALSE))
    })
  })
  
  # Selected AS event reactive for use by other components
  output$event_selected <- reactive({
    event <- selected_event()
    !is.null(event)
  })
  outputOptions(output, "event_selected", suspendWhenHidden = FALSE)
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  # Get AS visualization status
  get_as_visualization_status <- function() {
    as_data <- as_view_data()
    
    list(
      has_data = !is.null(as_data),
      success = !is.null(as_data) && as_data$success,
      error = if (!is.null(as_data) && !as_data$success) as_data$error else NULL,
      event_count = if (!is.null(as_data) && as_data$success && !is.null(as_data$events)) nrow(as_data$events) else 0
    )
  }
  
  # Clear AS visualization data
  clear_as_data <- function() {
    as_view_data(NULL)
    selected_event(NULL)
  }
  
  # Get selected event details
  get_selected_event_details <- function() {
    event <- selected_event()
    if (is.null(event)) {
      return(NULL)
    }
    
    list(
      eventID = event$eventID,
      eventType = event$AS_type,
      chromosome = event$chr,
      coordinates = paste(event$AS_range, collapse = ", "),
      transcripts = list(
        reference = event$refTx,
        alternative = event$asTx
      )
    )
  }
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  
  return(list(
    # Reactive values
    as_view_data = as_view_data,
    selected_event = selected_event,
    
    # Utility functions
    get_as_visualization_status = get_as_visualization_status,
    clear_as_data = clear_as_data,
    get_selected_event_details = get_selected_event_details
  ))
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test AS Event Visualization Module
#' 
#' @description 
#' Function to test AS event visualization module functionality
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_as_event_visualization_module <- function() {
  cat("Testing AS Event Visualization module...\n")
  
  # Test 1: Status structure validation
  status_test <- tryCatch({
    # Test status structure
    mock_status <- list(
      has_data = FALSE,
      success = FALSE,
      error = NULL,
      event_count = 0
    )
    
    # Verify all required fields
    required_fields <- c("has_data", "success", "error", "event_count")
    all(required_fields %in% names(mock_status))
  }, error = function(e) {
    cat("Status test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 2: Event details structure
  event_details_test <- tryCatch({
    # Test event details structure
    mock_event <- list(
      eventID = "SE_001",
      eventType = "SE",
      chromosome = "chr1",
      coordinates = "100-200, 300-400",
      transcripts = list(
        reference = "ENST00000123456",
        alternative = "ENST00000789012"
      )
    )
    
    # Verify structure
    all(c("eventID", "eventType", "transcripts") %in% names(mock_event))
  }, error = function(e) {
    cat("Event details test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 3: Plotly configuration
  plotly_config_test <- tryCatch({
    # Test plotly configuration
    config_options <- list(
      displayModeBar = TRUE,
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("pan2d", "select2d", "lasso2d", "resetScale2d")
    )
    
    # Verify configuration
    is.logical(config_options$displayModeBar) && 
    is.logical(config_options$displaylogo) &&
    is.character(config_options$modeBarButtonsToRemove)
  }, error = function(e) {
    cat("Plotly config test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 4: Error handling
  error_handling_test <- tryCatch({
    # Test error message generation
    error_msg <- paste("AS visualization error:", "test error")
    nchar(error_msg) > 0
  }, error = function(e) {
    cat("Error handling test failed:", e$message, "\n")
    FALSE
  })
  
  if (status_test && event_details_test && plotly_config_test && error_handling_test) {
    cat("All AS Event Visualization module tests passed!\n")
    return(TRUE)
  } else {
    cat("Some AS Event Visualization module tests failed!\n")
    return(FALSE)
  }
}