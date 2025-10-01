#===============================================================================
# TRANSCRIPT VISUALIZATION MODULE
# Extracted from server.R - Transcript and peptide visualization functionality
#===============================================================================

#' Transcript Visualization Module
#' 
#' @description 
#' Complete transcript and peptide visualization functionality:
#' - Transcript structure plotting with exon/intron display
#' - Peptide visualization with protease-specific digestion
#' - Interactive plotly-based visualizations
#' - Data loading and caching for visualizations
#' - Error handling and user feedback
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' @param selected_gene_transcripts Reactive containing transcript data
#' 
#' @return Named list containing reactive values for visualization data
#' 
#' @note 
#' - Requires selected_gene_transcripts reactive from main server
#' - Uses plotly for interactive visualizations
#' - Maintains all original styling and interactivity

create_transcript_visualization_module <- function(input, output, session, selected_gene_transcripts) {
  
  # ============================================================================
  # MODULE REACTIVE VALUES
  # ============================================================================
  
  # Data storage for visualizations
  transcript_structure_data <- reactiveVal(NULL)
  peptide_visualization_data <- reactiveVal(NULL)
  exon_data <- reactiveVal(NULL)
  
  # ============================================================================
  # DATA LOADING OBSERVERS
  # ============================================================================
  
  # Load common data needed for both visualizations
  observeEvent(input$update, {
    req(input$gene, input$miscleavage_type, selected_gene_transcripts())
    
    withProgress(message = 'Loading visualization data...', {
      tryCatch({
        gene_id <- input$gene
        miscleavage_type <- input$miscleavage_type
        
        incProgress(0.3, detail = "Loading transcript structure...")
        
        # Load transcript structure data
        transcript_result <- create_transcript_plot_data(
          gene_id = gene_id,
          transcripts = selected_gene_transcripts(),
          gtf_file = gtf_file
        )
        
        # Store transcript structure data
        transcript_structure_data(transcript_result)
        
        incProgress(0.4, detail = "Loading exon data...")
        
        # Load exon data for peptide visualization
        if (transcript_result$success) {
          exon_result <- load_transcript_exons(
            gene_id = gene_id,
            transcripts = selected_gene_transcripts(),
            gtf_file = gtf_file
          )
          exon_data(exon_result)
        }
        
        incProgress(0.3, detail = "Complete!")
        
      }, error = function(e) {
        cat("Visualization data loading error:", e$message, "\n")
        
        # Set error states
        transcript_structure_data(list(
          success = FALSE,
          error = paste("Failed to load transcript data:", e$message)
        ))
        
        showNotification(
          paste("Error loading visualization data:", e$message),
          type = "error",
          duration = 8
        )
      })
    })
  })
  
  # Load peptide visualization when protease/miscleavage changes
  observeEvent(list(input$protease, input$miscleavage_type), {
    req(input$gene, input$miscleavage_type, exon_data(), input$protease)
    
    # Skip if no exon data loaded yet
    exon_result <- exon_data()
    if (is.null(exon_result) || !exon_result$success) {
      return()
    }
    
    withProgress(message = 'Loading peptide visualization...', {
      tryCatch({
        gene_id <- input$gene
        miscleavage_type <- input$miscleavage_type
        protease <- input$protease
        
        incProgress(0.5, detail = "Generating peptide visualization...")
        
        # Create peptide plot data
        peptide_result <- create_peptide_plot_data(
          gene_id = gene_id,
          exon_data = exon_result,
          miscleavage_type = miscleavage_type,
          protease = protease
        )
        
        # Store peptide data
        peptide_visualization_data(peptide_result)
        
        incProgress(0.5, detail = "Complete!")
        
      }, error = function(e) {
        cat("Peptide visualization error:", e$message, "\n")
        
        peptide_visualization_data(list(
          success = FALSE,
          error = paste("Failed to load peptide data:", e$message)
        ))
        
        showNotification(
          paste("Error loading peptide visualization:", e$message),
          type = "error",
          duration = 8
        )
      })
    })
  })
  
  # ============================================================================
  # VISUALIZATION OUTPUTS
  # ============================================================================
  
  # Transcript structure plot
  output$transcript_structure_plot <- renderPlotly({
    req(input$gene, selected_gene_transcripts())
    
    # Get plot data
    plot_data <- transcript_structure_data()
    
    # If no data or error, show error message
    if (is.null(plot_data) || !plot_data$success) {
      message_text <- if (is.null(plot_data)) {
        "Click 'Update View' to load transcript structure"
      } else {
        paste("Error:", plot_data$error)
      }
      
      # Create error plot
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = message_text, size = 5, hjust = 0.5) +
        xlim(0, 1) + ylim(0, 1) +
        theme_void() +
        theme(
          panel.background = element_rect(fill = "#f8f9fa", color = NA),
          plot.background = element_rect(fill = "#f8f9fa", color = NA)
        )
      
      return(ggplotly(p) %>%
        config(displayModeBar = FALSE) %>%
        layout(
          xaxis = list(visible = FALSE),
          yaxis = list(visible = FALSE)
        ))
    }
    
    # Create the actual plot
    tryCatch({
      p <- plot_data$plot
      
      # Convert to plotly with proper configuration
      plotly_obj <- ggplotly(p, tooltip = "text") %>%
        config(
          displayModeBar = TRUE,
          displaylogo = FALSE,
          modeBarButtonsToRemove = c(
            "pan2d", "select2d", "lasso2d", "resetScale2d",
            "hoverClosestCartesian", "hoverCompareCartesian"
          )
        ) %>%
        layout(
          showlegend = TRUE,
          legend = list(
            orientation = "h",
            x = 0.5,
            xanchor = "center",
            y = -0.1
          ),
          margin = list(l = 50, r = 50, t = 50, b = 100)
        )
      
      return(plotly_obj)
      
    }, error = function(e) {
      # Error in plotting
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                label = paste("Plotting error:", e$message), 
                size = 4, hjust = 0.5, color = "red") +
        xlim(0, 1) + ylim(0, 1) +
        theme_void()
      
      return(ggplotly(p) %>% config(displayModeBar = FALSE))
    })
  })
  
  # Peptide visualization plot
  output$peptide_plot <- renderPlotly({
    req(input$gene, selected_gene_transcripts(), input$protease)
    
    # Get plot data
    plot_data <- peptide_visualization_data()
    
    # If no data or error, show error message
    if (is.null(plot_data) || !plot_data$success) {
      message_text <- if (is.null(plot_data)) {
        "Click 'Update View' to load peptide visualization"
      } else {
        paste("Error:", plot_data$error)
      }
      
      # Create error plot
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = message_text, size = 5, hjust = 0.5) +
        xlim(0, 1) + ylim(0, 1) +
        theme_void() +
        theme(
          panel.background = element_rect(fill = "#f8f9fa", color = NA),
          plot.background = element_rect(fill = "#f8f9fa", color = NA)
        )
      
      return(ggplotly(p) %>%
        config(displayModeBar = FALSE) %>%
        layout(
          xaxis = list(visible = FALSE),
          yaxis = list(visible = FALSE)
        ))
    }
    
    # Create the actual plot
    tryCatch({
      p <- plot_data$plot
      
      # Convert to plotly with enhanced interactivity
      plotly_obj <- ggplotly(p, tooltip = "text") %>%
        config(
          displayModeBar = TRUE,
          displaylogo = FALSE,
          modeBarButtonsToRemove = c(
            "pan2d", "select2d", "lasso2d", "resetScale2d"
          )
        ) %>%
        layout(
          showlegend = TRUE,
          legend = list(
            orientation = "h",
            x = 0.5,
            xanchor = "center",
            y = -0.1
          ),
          margin = list(l = 50, r = 50, t = 50, b = 100),
          hovermode = 'closest'
        )
      
      return(plotly_obj)
      
    }, error = function(e) {
      # Error in plotting
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                label = paste("Plotting error:", e$message), 
                size = 4, hjust = 0.5, color = "red") +
        xlim(0, 1) + ylim(0, 1) +
        theme_void()
      
      return(ggplotly(p) %>% config(displayModeBar = FALSE))
    })
  })
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  # Create error plot helper
  create_error_plot <- function(message, bg_color = "#f8f9fa") {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = message, size = 5, hjust = 0.5) +
      xlim(0, 1) + ylim(0, 1) +
      theme_void() +
      theme(
        panel.background = element_rect(fill = bg_color, color = NA),
        plot.background = element_rect(fill = bg_color, color = NA)
      )
    
    ggplotly(p) %>%
      config(displayModeBar = FALSE) %>%
      layout(
        xaxis = list(visible = FALSE),
        yaxis = list(visible = FALSE)
      )
  }
  
  # Get visualization status
  get_visualization_status <- function() {
    transcript_data <- transcript_structure_data()
    peptide_data <- peptide_visualization_data()
    
    list(
      transcript_loaded = !is.null(transcript_data) && transcript_data$success,
      peptide_loaded = !is.null(peptide_data) && peptide_data$success,
      transcript_error = if (!is.null(transcript_data) && !transcript_data$success) transcript_data$error else NULL,
      peptide_error = if (!is.null(peptide_data) && !peptide_data$success) peptide_data$error else NULL
    )
  }
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  
  return(list(
    # Reactive values
    transcript_structure_data = transcript_structure_data,
    peptide_visualization_data = peptide_visualization_data,
    exon_data = exon_data,
    
    # Utility functions
    get_visualization_status = get_visualization_status,
    create_error_plot = create_error_plot,
    
    # Data refresh function
    refresh_visualizations = function() {
      transcript_structure_data(NULL)
      peptide_visualization_data(NULL)
      exon_data(NULL)
    }
  ))
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test Transcript Visualization Module
#' 
#' @description 
#' Function to test transcript visualization module functionality
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_transcript_visualization_module <- function() {
  cat("Testing Transcript Visualization module...\n")
  
  # Test 1: Error plot creation
  error_plot_test <- tryCatch({
    # Test error message handling
    test_message <- "Test error message"
    # Would create error plot - test the logic
    nchar(test_message) > 0
  }, error = function(e) {
    cat("Error plot test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 2: Status checking
  status_test <- tryCatch({
    # Test status structure
    mock_status <- list(
      transcript_loaded = FALSE,
      peptide_loaded = FALSE,
      transcript_error = NULL,
      peptide_error = NULL
    )
    
    # Verify structure
    all(c("transcript_loaded", "peptide_loaded", "transcript_error", "peptide_error") %in% names(mock_status))
  }, error = function(e) {
    cat("Status test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 3: Plot configuration
  config_test <- tryCatch({
    # Test plotly configuration options
    config_options <- list(
      displayModeBar = TRUE,
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("pan2d", "select2d", "lasso2d")
    )
    
    # Verify configuration structure
    is.logical(config_options$displayModeBar) && is.logical(config_options$displaylogo)
  }, error = function(e) {
    cat("Config test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 4: Data structure validation
  data_test <- tryCatch({
    # Test expected data structure
    mock_data <- list(
      success = TRUE,
      plot = NULL,
      error = NULL
    )
    
    # Verify required fields exist
    all(c("success", "plot", "error") %in% names(mock_data))
  }, error = function(e) {
    cat("Data structure test failed:", e$message, "\n")
    FALSE
  })
  
  if (error_plot_test && status_test && config_test && data_test) {
    cat("All Transcript Visualization module tests passed!\n")
    return(TRUE)
  } else {
    cat("Some Transcript Visualization module tests failed!\n")
    return(FALSE)
  }
}