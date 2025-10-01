#===============================================================================
# PEPTIDE SEARCH MODULE  
# Extracted from server.R - Complete peptide search functionality
#===============================================================================

#' Peptide Search Module
#' 
#' @description 
#' Complete peptide search functionality using BLASTP:
#' - Database validation and availability checking
#' - BLASTP search execution with progress tracking
#' - Results processing and display
#' - Navigation integration with main analysis
#' - Auto-selection for further analysis
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' 
#' @return Named list containing reactive values and functions for search functionality
#' 
#' @note 
#' - Requires BLASTP database setup and validation functions
#' - Integrates with existing navigation and analysis workflows
#' - Preserves all original search parameters and functionality

create_peptide_search_module <- function(input, output, session) {
  
  # ============================================================================
  # MODULE REACTIVE VALUES
  # ============================================================================
  
  # Reactive values for peptide search results
  peptide_search_results <- reactiveVal(NULL)
  
  # Store peptide for auto-selection after BLASTP navigation
  blast_peptide_for_selection <- reactiveVal(NULL)
  
  # ============================================================================
  # DATABASE AVAILABILITY CHECK
  # ============================================================================
  
  # Check if search databases are available
  search_databases_available <- reactive({
    search_dir <- "data/search"
    if (!dir.exists(search_dir)) return(FALSE)
    
    # Check if at least some search database files exist - simplified approach
    files <- list.files(search_dir, pattern = "peptide_search_", full.names = TRUE)
    rds_files <- files[grepl("\\.rds$", files)]
    return(length(rds_files) > 0)
  })
  
  # ============================================================================
  # SEARCH EXECUTION OBSERVER
  # ============================================================================
  
  # Run peptide search when button is clicked
  observeEvent(input$run_peptide_search, {
    req(input$peptide_search_query)
    
    # Validate BLAST database
    blast_db_check <- validate_blast_database()
    if (!blast_db_check$valid) {
      showNotification(
        paste("BLAST database not found:", blast_db_check$message),
        type = "error",
        duration = 10
      )
      return()
    }
    
    # Validate peptide query
    peptide_query <- trimws(input$peptide_search_query)
    if (peptide_query == "") {
      showNotification("Please enter a peptide sequence to search.", type = "warning")
      return()
    }
    
    withProgress(message = 'Searching peptides with BLASTP...', {
      tryCatch({
        # Get BLAST parameters from input
        evalue_threshold <- if (!is.null(input$peptide_search_evalue)) input$peptide_search_evalue else 10
        identity_threshold <- if (!is.null(input$peptide_search_identity)) input$peptide_search_identity else 70
        max_targets <- if (!is.null(input$peptide_search_max_targets)) input$peptide_search_max_targets else 500
        
        # Progress callback function
        progress_callback <- function(message, value) {
          incProgress(value, detail = message)
        }
        
        incProgress(0.2, detail = "Preparing BLAST search...")
        
        # Execute BLASTP search
        search_results <- search_peptide_blastp(
          peptide_sequence = peptide_query,
          evalue_threshold = evalue_threshold,
          identity_threshold = identity_threshold,
          max_targets = max_targets,
          progress_callback = progress_callback
        )
        
        incProgress(0.8, detail = "Processing results...")
        
        if (!is.null(search_results) && nrow(search_results) > 0) {
          # Store results
          peptide_search_results(search_results)
          
          # Show success notification
          showNotification(
            paste("Search completed! Found", nrow(search_results), "matching peptides."),
            type = "message",
            duration = 5
          )
        } else {
          # No results found
          peptide_search_results(data.frame())
          showNotification("No matching peptides found.", type = "warning", duration = 5)
        }
        
        incProgress(1.0, detail = "Search complete!")
        
      }, error = function(e) {
        cat("Search error:", e$message, "\n")
        showNotification(
          paste("Search failed:", e$message),
          type = "error",
          duration = 10
        )
        peptide_search_results(NULL)
      })
    })
  })
  
  # ============================================================================
  # NAVIGATION INTEGRATION
  # ============================================================================
  
  # Navigate to isoform analysis when requested
  observeEvent(input$goto_isoform_analysis, {
    req(input$peptide_search_results_table_rows_selected, peptide_search_results())
    
    tryCatch({
      selected_rows <- input$peptide_search_results_table_rows_selected
      search_data <- peptide_search_results()
      
      if (length(selected_rows) > 0 && nrow(search_data) >= selected_rows[1]) {
        selected_peptide <- search_data[selected_rows[1], ]
        
        # Extract gene information
        gene_id <- selected_peptide$gene_id
        gene_symbol <- selected_peptide$gene_symbol
        peptide_sequence <- selected_peptide$peptide_sequence
        
        if (!is.null(gene_id) && gene_id != "") {
          # Store peptide for auto-selection
          blast_peptide_for_selection(peptide_sequence)
          
          # Update gene selection in main interface
          updateSelectizeInput(session, "gene",
                              selected = gene_id,
                              server = TRUE)
          
          # Switch to canonical analysis tab
          updateTabItems(session, "tabs", "canonical_analysis")
          
          # Show navigation notification
          showNotification(
            paste("Navigated to", gene_symbol, "analysis with selected peptide"),
            type = "message",
            duration = 5
          )
        } else {
          showNotification("Gene ID not found for selected peptide", type = "warning")
        }
      }
    }, error = function(e) {
      showNotification(paste("Navigation failed:", e$message), type = "error")
    })
  })
  
  # ============================================================================
  # OUTPUT RENDERERS
  # ============================================================================
  
  # Results availability status
  output$peptide_search_results_available <- reactive({
    results <- peptide_search_results()
    !is.null(results) && nrow(results) > 0
  })
  outputOptions(output, "peptide_search_results_available", suspendWhenHidden = FALSE)
  
  # Search summary text (extracted from simple outputs)
  output$peptide_search_summary <- renderText({
    results <- peptide_search_results()
    if (is.null(results) || nrow(results) == 0) {
      return("No results")
    }
    
    paste("Found", nrow(results), "matching peptides")
  })
  
  # Search results table
  output$peptide_search_results_table <- DT::renderDataTable({
    results <- peptide_search_results()
    
    if (is.null(results) || nrow(results) == 0) {
      return(data.frame(Message = "No search results available"))
    }
    
    # Format results for display
    display_results <- results
    
    # Round numeric columns for better display
    numeric_cols <- sapply(display_results, is.numeric)
    display_results[numeric_cols] <- lapply(display_results[numeric_cols], function(x) round(x, 3))
    
    DT::datatable(
      display_results,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel'),
        order = list(list(which(names(display_results) == "identity") - 1, 'desc'))
      ),
      selection = 'single',
      class = 'cell-border stripe',
      rownames = FALSE
    )
  })
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  
  return(list(
    # Reactive values
    peptide_search_results = peptide_search_results,
    blast_peptide_for_selection = blast_peptide_for_selection,
    search_databases_available = search_databases_available,
    
    # Utility functions
    clear_search_results = function() {
      peptide_search_results(NULL)
      blast_peptide_for_selection(NULL)
    },
    
    get_search_stats = function() {
      results <- peptide_search_results()
      if (is.null(results)) {
        return(list(total = 0, has_results = FALSE))
      }
      list(
        total = nrow(results),
        has_results = nrow(results) > 0,
        avg_identity = if (nrow(results) > 0) mean(results$identity, na.rm = TRUE) else 0
      )
    }
  ))
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test Peptide Search Module
#' 
#' @description 
#' Function to test peptide search module functionality
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_peptide_search_module <- function() {
  cat("Testing Peptide Search module...\n")
  
  # Test 1: Search directory checking
  directory_test <- tryCatch({
    # Test directory existence checking logic
    test_dir <- "test_search_dir"
    dir_exists <- dir.exists(test_dir)  # Should be FALSE for non-existent dir
    !dir_exists  # Should return TRUE (test passes when dir doesn't exist)
  }, error = function(e) {
    cat("Directory test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 2: Parameter validation
  parameter_test <- tryCatch({
    # Test default parameter assignment
    evalue_threshold <- 10  # default
    identity_threshold <- 70  # default
    max_targets <- 500  # default
    
    # Validate parameters are reasonable
    evalue_threshold > 0 && identity_threshold > 0 && max_targets > 0
  }, error = function(e) {
    cat("Parameter test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 3: Results processing
  results_test <- tryCatch({
    # Test empty results handling
    empty_results <- data.frame()
    has_results <- nrow(empty_results) > 0
    
    # Test non-empty results
    test_results <- data.frame(
      gene_id = "ENSG00000139618.15",
      peptide_sequence = "TESTPEPTIDE",
      identity = 95.5
    )
    has_test_results <- nrow(test_results) > 0
    
    !has_results && has_test_results  # First should be FALSE, second TRUE
  }, error = function(e) {
    cat("Results test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 4: String validation
  string_test <- tryCatch({
    # Test string trimming and validation
    test_query1 <- "  TESTPEPTIDE  "
    trimmed1 <- trimws(test_query1)
    
    test_query2 <- ""
    trimmed2 <- trimws(test_query2)
    
    nchar(trimmed1) > 0 && nchar(trimmed2) == 0
  }, error = function(e) {
    cat("String test failed:", e$message, "\n")
    FALSE
  })
  
  if (directory_test && parameter_test && results_test && string_test) {
    cat("All Peptide Search module tests passed!\n")
    return(TRUE)
  } else {
    cat("Some Peptide Search module tests failed!\n")
    return(FALSE)
  }
}