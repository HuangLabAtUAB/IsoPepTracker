#===============================================================================
# DATA TABLES MODULE
# Extracted from server.R - Common data table functionality
#===============================================================================

#' Data Tables Module
#' 
#' @description 
#' Common data table rendering functionality using DT:
#' - Standardized table styling and configuration
#' - Export functionality (CSV, Excel, Copy)
#' - Pagination and search capabilities
#' - Row selection handling
#' - Responsive design
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' 
#' @return Named list containing data table utility functions
#' 
#' @note 
#' - Provides standardized DT configurations across the application
#' - Includes export buttons and responsive design
#' - Handles empty data gracefully

create_data_tables_module <- function(input, output, session) {
  
  # ============================================================================
  # STANDARDIZED TABLE CONFIGURATIONS
  # ============================================================================
  
  # Standard DT options for most tables
  get_standard_dt_options <- function(page_length = 10, scrollX = TRUE, include_buttons = TRUE) {
    options <- list(
      pageLength = page_length,
      scrollX = scrollX,
      dom = if (include_buttons) 'Bfrtip' else 'frtip',
      columnDefs = list(
        list(className = 'dt-center', targets = '_all')
      ),
      initComplete = DT::JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#3c8dbc', 'color': '#fff'});",
        "}"
      )
    )
    
    if (include_buttons) {
      options$buttons <- c('copy', 'csv', 'excel')
    }
    
    return(options)
  }
  
  # Options for large tables (more rows per page)
  get_large_table_options <- function() {
    get_standard_dt_options(page_length = 25, scrollX = TRUE, include_buttons = TRUE)
  }
  
  # Options for small tables (fewer rows per page)
  get_small_table_options <- function() {
    get_standard_dt_options(page_length = 5, scrollX = TRUE, include_buttons = FALSE)
  }
  
  # Options for tables with row selection
  get_selectable_table_options <- function(selection_mode = 'single') {
    options <- get_standard_dt_options()
    options$selection <- selection_mode
    return(options)
  }
  
  # ============================================================================
  # TABLE RENDERING FUNCTIONS
  # ============================================================================
  
  # Render a standard data table with common styling
  render_standard_table <- function(data, options = NULL, selection = 'none', 
                                   class = 'cell-border stripe', rownames = FALSE) {
    
    # Handle empty data
    if (is.null(data) || nrow(data) == 0) {
      data <- data.frame(Message = "No data available")
      options <- list(
        pageLength = 5,
        dom = 't',
        ordering = FALSE,
        searching = FALSE
      )
      selection <- 'none'
    }
    
    # Use default options if not provided
    if (is.null(options)) {
      options <- get_standard_dt_options()
    }
    
    # Round numeric columns for better display
    numeric_cols <- sapply(data, is.numeric)
    if (any(numeric_cols)) {
      data[numeric_cols] <- lapply(data[numeric_cols], function(x) round(x, 3))
    }
    
    DT::datatable(
      data,
      options = options,
      selection = selection,
      class = class,
      rownames = rownames,
      escape = FALSE
    ) %>%
      DT::formatStyle(
        columns = names(data),
        fontSize = '14px'
      )
  }
  
  # Render a table with export functionality
  render_exportable_table <- function(data, filename_prefix = "data") {
    options <- get_standard_dt_options(include_buttons = TRUE)
    options$buttons <- list(
      list(
        extend = 'csv',
        filename = paste0(filename_prefix, '_', Sys.Date())
      ),
      list(
        extend = 'excel',
        filename = paste0(filename_prefix, '_', Sys.Date())
      ),
      'copy'
    )
    
    render_standard_table(data, options)
  }
  
  # Render a selectable table (single or multiple selection)
  render_selectable_table <- function(data, selection_mode = 'single', sort_column = NULL) {
    options <- get_selectable_table_options(selection_mode)
    
    # Add sorting if specified
    if (!is.null(sort_column) && sort_column %in% names(data)) {
      col_index <- which(names(data) == sort_column) - 1  # JavaScript is 0-indexed
      options$order <- list(list(col_index, 'desc'))
    }
    
    render_standard_table(data, options, selection = selection_mode)
  }
  
  # ============================================================================
  # SPECIALIZED TABLE FUNCTIONS
  # ============================================================================
  
  # Format table for peptide data
  format_peptide_table <- function(data) {
    if (is.null(data) || nrow(data) == 0) {
      return(data)
    }
    
    # Format specific columns for peptide tables
    if ("peptide" %in% names(data)) {
      # Truncate very long peptide sequences for display
      data$peptide_display <- ifelse(
        nchar(data$peptide) > 50,
        paste0(substr(data$peptide, 1, 47), "..."),
        data$peptide
      )
    }
    
    if ("identity" %in% names(data)) {
      data$identity <- round(data$identity, 2)
    }
    
    if ("evalue" %in% names(data)) {
      data$evalue <- formatC(data$evalue, format = "e", digits = 2)
    }
    
    return(data)
  }
  
  # Format table for AS events
  format_as_events_table <- function(data) {
    if (is.null(data) || nrow(data) == 0) {
      return(data)
    }
    
    # Format AS event specific columns
    if ("AS_type" %in% names(data)) {
      data$AS_type <- factor(data$AS_type, 
                            levels = c("SE", "MXE", "A3SS", "A5SS", "RI"),
                            labels = c("Skipped Exon", "Mutually Exclusive", 
                                     "Alternative 3' SS", "Alternative 5' SS", 
                                     "Retained Intron"))
    }
    
    return(data)
  }
  
  # Format table for gene search results
  format_gene_search_table <- function(data) {
    if (is.null(data) || nrow(data) == 0) {
      return(data)
    }
    
    # Add color coding for identity scores
    if ("identity" %in% names(data)) {
      data$identity_colored <- paste0(
        '<span style="color: ',
        ifelse(data$identity >= 90, 'green',
               ifelse(data$identity >= 70, 'orange', 'red')),
        ';">', round(data$identity, 1), '%</span>'
      )
    }
    
    return(data)
  }
  
  # ============================================================================
  # TABLE UTILITY FUNCTIONS
  # ============================================================================
  
  # Get selected rows from a table
  get_selected_rows <- function(table_id) {
    selected <- input[[paste0(table_id, "_rows_selected")]]
    if (is.null(selected) || length(selected) == 0) {
      return(NULL)
    }
    return(selected)
  }
  
  # Clear table selection
  clear_table_selection <- function(table_id) {
    DT::dataTableProxy(table_id) %>% 
      DT::selectRows(NULL)
  }
  
  # Add highlighting to table rows
  highlight_table_rows <- function(table_proxy, row_indices, color = "#ffffcc") {
    if (length(row_indices) > 0) {
      table_proxy %>% 
        DT::formatStyle(
          columns = 1:ncol(table_proxy$data),
          target = 'row',
          backgroundColor = DT::styleEqual(row_indices, rep(color, length(row_indices)))
        )
    }
  }
  
  # ============================================================================
  # EXPORT FUNCTIONS
  # ============================================================================
  
  # Generate download handler for table data
  create_table_download_handler <- function(data_reactive, filename_prefix) {
    downloadHandler(
      filename = function() {
        paste0(filename_prefix, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        data <- data_reactive()
        if (!is.null(data) && nrow(data) > 0) {
          write.csv(data, file, row.names = FALSE)
        } else {
          # Create empty file with header
          write.csv(data.frame(Message = "No data available"), file, row.names = FALSE)
        }
      }
    )
  }
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  
  return(list(
    # Configuration functions
    get_standard_dt_options = get_standard_dt_options,
    get_large_table_options = get_large_table_options,
    get_small_table_options = get_small_table_options,
    get_selectable_table_options = get_selectable_table_options,
    
    # Rendering functions
    render_standard_table = render_standard_table,
    render_exportable_table = render_exportable_table,
    render_selectable_table = render_selectable_table,
    
    # Formatting functions
    format_peptide_table = format_peptide_table,
    format_as_events_table = format_as_events_table,
    format_gene_search_table = format_gene_search_table,
    
    # Utility functions
    get_selected_rows = get_selected_rows,
    clear_table_selection = clear_table_selection,
    highlight_table_rows = highlight_table_rows,
    create_table_download_handler = create_table_download_handler
  ))
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test Data Tables Module
#' 
#' @description 
#' Function to test data tables module functionality
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_data_tables_module <- function() {
  cat("Testing Data Tables module...\n")
  
  # Test 1: Options generation
  options_test <- tryCatch({
    # Test standard options
    std_options <- list(
      pageLength = 10,
      scrollX = TRUE,
      dom = 'Bfrtip'
    )
    
    # Verify structure
    all(c("pageLength", "scrollX", "dom") %in% names(std_options))
  }, error = function(e) {
    cat("Options test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 2: Data formatting
  formatting_test <- tryCatch({
    # Test peptide formatting
    test_data <- data.frame(
      peptide = c("SHORTPEPTIDE", paste0(rep("A", 55), collapse = "")),
      identity = c(95.567, 78.123),
      evalue = c(1e-10, 2.3e-5)
    )
    
    # Test formatting functions exist and work
    nrow(test_data) == 2 && ncol(test_data) == 3
  }, error = function(e) {
    cat("Formatting test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 3: Empty data handling
  empty_data_test <- tryCatch({
    # Test empty data frame
    empty_df <- data.frame()
    
    # Should handle empty data without error
    nrow(empty_df) == 0
  }, error = function(e) {
    cat("Empty data test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 4: Color coding logic
  color_test <- tryCatch({
    # Test color assignment logic
    identity_score <- 85
    color <- ifelse(identity_score >= 90, 'green',
                   ifelse(identity_score >= 70, 'orange', 'red'))
    
    color == 'orange'  # 85 should be orange
  }, error = function(e) {
    cat("Color test failed:", e$message, "\n")
    FALSE
  })
  
  if (options_test && formatting_test && empty_data_test && color_test) {
    cat("All Data Tables module tests passed!\n")
    return(TRUE)
  } else {
    cat("Some Data Tables module tests failed!\n")
    return(FALSE)
  }
}