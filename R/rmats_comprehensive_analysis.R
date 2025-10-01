#===============================================================================
# rMATS COMPREHENSIVE ANALYSIS
# The missing function that chains coordinate extraction → GTF generation → peptide analysis
#===============================================================================

# Load required dependencies
library(data.table)

# Note: rMATS modules will be loaded by the rMATS server
# The analyze_rmats_event_comprehensive function will be available when called from server context

#===============================================================================
# MAIN COMPREHENSIVE ANALYSIS FUNCTION
#===============================================================================

analyze_rmats_event_comprehensive <- function(selected_event, event_type) {
  
  cat("🔄 Starting comprehensive rMATS analysis...\n")
  cat("Event Type:", event_type, "\n")
  cat("Gene:", selected_event$geneSymbol, "ID:", selected_event$GeneID, "\n")
  
  tryCatch({
    
    # Check that required functions are available
    if (!exists("extract_event_coordinates")) {
      stop("rMATS modules not loaded. This function must be called from rMATS server context.")
    }
    
    # STEP 1: Extract event coordinates using existing rMATS parser
    cat("Step 1: Extracting event coordinates...\n")
    event_coords <- extract_event_coordinates(selected_event, event_type)
    
    # STEP 2: Generate GTF and FASTA files using existing GTF generator
    cat("Step 2: Generating GTF and FASTA files...\n")
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_dir <- file.path("rmats_peptide_results", timestamp)
    
    gtf_results <- generate_rmats_gtf_with_phases(selected_event, event_type, output_dir)
    
    if (!gtf_results$success) {
      return(list(
        success = FALSE,
        error = paste("GTF generation failed:", gtf_results$error)
      ))
    }
    
    # STEP 3: Generate peptide data using existing peptide generator
    cat("Step 3: Generating peptide data...\n")
    gtf_file <- file.path(output_dir, "novel_transcript_nt.transdecoder.genome.gtf")
    fasta_file <- file.path(output_dir, "novel_transcript_nt.fa.transdecoder.pep")
    
    if (!file.exists(gtf_file) || !file.exists(fasta_file)) {
      return(list(
        success = FALSE,
        error = "GTF or FASTA file not found after generation"
      ))
    }
    
    peptide_results <- generate_novel_peptide_data(
      gtf_file = gtf_file,
      fasta_file = fasta_file,
      missedCleavages = 0,
      minLength = 7,
      maxLength = 60
    )
    
    # STEP 4: Save results in expected format
    cat("Step 4: Saving analysis results...\n")
    results_dir <- file.path(output_dir, "results")
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    dataframe_file <- file.path(results_dir, "rmats_isoform_dataframe.rds")
    saveRDS(peptide_results, dataframe_file)
    
    # Create pipeline results structure (matching expected format)
    pipeline_results <- list(
      dataframe_file = dataframe_file,
      gtf_file = gtf_file,
      fasta_file = fasta_file,
      output_dir = output_dir,
      timestamp = timestamp,
      event_info = list(
        event_type = event_type,
        gene_id = selected_event$GeneID,
        gene_symbol = selected_event$geneSymbol,
        coordinates = event_coords
      )
    )
    
    cat("✅ Comprehensive rMATS analysis completed successfully!\n")
    cat("📁 Results directory:", output_dir, "\n")
    cat("📊 Generated", nrow(peptide_results), "protein isoforms\n")
    
    return(list(
      success = TRUE,
      pipeline_results = pipeline_results,
      peptide_data = peptide_results,
      summary = paste("Successfully analyzed", event_type, "event for gene", selected_event$geneSymbol)
    ))
    
  }, error = function(e) {
    cat("❌ Error in comprehensive rMATS analysis:", conditionMessage(e), "\n")
    return(list(
      success = FALSE,
      error = conditionMessage(e)
    ))
  })
}

cat("✅ rMATS comprehensive analysis module loaded\n")