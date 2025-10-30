#===============================================================================
# rMATS COMPREHENSIVE ANALYSIS
# The missing function that chains coordinate extraction → GTF generation → peptide analysis
#===============================================================================

# Load required dependencies
library(data.table)

# Note: Functions will be sourced on-demand to avoid startup errors
# The rMATS GTF generator and novel peptide generator will be loaded when needed

# Note: Other rMATS modules will be loaded by the rMATS server
# The extract_event_coordinates function will be available when called from server context

#===============================================================================
# STREAMLINED ANALYSIS FUNCTION (ELIMINATES REDUNDANT CDS SEARCHES)
#===============================================================================

analyze_rmats_event_streamlined <- function(selected_event, event_type, missedCleavages = 0) {
  
  cat("🚀 Starting STREAMLINED rMATS analysis...\n")
  cat("Event Type:", event_type, "\n")
  cat("Gene:", selected_event$geneSymbol, "ID:", selected_event$GeneID, "\n")
  cat("Target miscleavages:", missedCleavages, "\n")
  
  tryCatch({
    
    # Check that required functions are available
    if (!exists("extract_event_coordinates")) {
      stop("rMATS modules not loaded. This function must be called from rMATS server context.")
    }
    
    # CRITICAL FIX: Load CDS index if not already loaded
    if (!exists("rmats_cds_index")) {
      cat("Loading CDS index for proper phase extraction...\n")
      tryCatch({
        # Source CDS index functions
        source("rmats_peptide/modules/create_cds_index.R")
        # Load the CDS index
        rmats_cds_index <- load_cds_index("rmats_peptide/real_cds_index.RDS")
        # Make it available globally for the flanking_exons functions
        assign("rmats_cds_index", rmats_cds_index, envir = .GlobalEnv)
        cat("✅ CDS index loaded successfully with", nrow(rmats_cds_index), "entries\n")
      }, error = function(e) {
        cat("⚠️  Warning: Could not load CDS index:", e$message, "\n")
        cat("⚠️  Will use simplified phase calculation (may be less accurate)\n")
        assign("rmats_cds_index", NULL, envir = .GlobalEnv)
      })
    } else {
      cat("✅ CDS index already loaded\n")
    }
    
    # STEP 1: Extract event coordinates (same as before)
    cat("Step 1: Extracting event coordinates...\n")
    event_coords <- extract_event_coordinates(selected_event, event_type)
    
    # STEP 2: Build GTF structures (same as before)
    cat("Step 2: Building GTF structures...\n")
    gtf_structures <- build_gtf_structures(event_coords, event_type)
    
    # STEP 3: Identify flanking exons (same as before)
    cat("Step 3: Identifying flanking exons...\n")
    flanking_exons <- identify_flanking_exons(gtf_structures)
    
    # STEP 4: SINGLE CDS SEARCH (the only one we need!)
    cat("Step 4: Performing SINGLE CDS search...\n")
    cds_results <- search_all_exons_in_cds(flanking_exons, 'rmats_peptide/real_cds_index.RDS')
    
    if (!cds_results$ready_for_step5) {
      return(list(
        success = FALSE,
        error = paste("CDS search failed:", cds_results$outcome)
      ))
    }
    
    cat("CDS search successful!\n")
    
    # STEP 5: Extract phase information (same as before)
    cat("Step 5: Extracting phase information...\n")
    phase_results <- extract_phase_information(cds_results)
    
    if (!phase_results$ready_for_step6) {
      return(list(
        success = FALSE,
        error = "No translatable isoforms found"
      ))
    }
    
    cat("Phase extraction successful!\n")
    cat("Inclusion translatable:", phase_results$inclusion_translatable, "\n")
    cat("Exclusion translatable:", phase_results$exclusion_translatable, "\n")
    
    # STEP 6: Generate protein sequences and create GTF/FASTA files
    cat("Step 6: Generating protein sequences and GTF/FASTA files...\n")
    translation_results <- analyze_protein_translation(phase_results)
    
    if (!translation_results$success) {
      return(list(
        success = FALSE,
        error = paste("Translation failed:", translation_results$summary)
      ))
    }
    
    cat("Translation successful!\n")
    cat("Case type:", translation_results$case_type, "\n")
    cat("Functional consequence:", translation_results$functional_consequence, "\n")
    
    # STEP 7: Generate GTF and FASTA files for the existing peptide generator
    cat("Step 7: Creating GTF and FASTA files for peptide generation...\n")
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_dir <- file.path("rmats_peptide_results", timestamp)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Load GTF generator functions if not already available
    if (!exists("generate_complete_gtf") || !exists("generate_protein_fasta")) {
      cat("Loading GTF generator functions...\n")
      source("rmats_exact_gtf_generator.R")
    }
    
    # Generate GTF file with phase information
    gtf_output <- generate_complete_gtf(
      phase_results, 
      selected_event, 
      event_type,
      file.path(output_dir, "novel_transcript_nt.transdecoder.genome.gtf")
    )
    
    # Generate FASTA file with protein sequences
    fasta_output <- generate_protein_fasta(
      translation_results,
      selected_event,
      file.path(output_dir, "novel_transcript_nt.fa.transdecoder.pep")
    )
    
    if (!gtf_output$success || !fasta_output$success) {
      return(list(
        success = FALSE,
        error = "Failed to generate GTF or FASTA files"
      ))
    }
    
    # STEP 8: Use EXISTING novel_peptide_generator (both miscleavage settings)
    cat("Step 8: Using EXISTING peptide generator (both miscleavage settings)...\n")
    gtf_file <- file.path(output_dir, "novel_transcript_nt.transdecoder.genome.gtf")
    fasta_file <- file.path(output_dir, "novel_transcript_nt.fa.transdecoder.pep")
    
    # DEBUG: Check if files actually exist
    cat("🔍 DEBUG: GTF file path:", gtf_file, "\n")
    cat("🔍 DEBUG: GTF file exists:", file.exists(gtf_file), "\n")
    cat("🔍 DEBUG: FASTA file path:", fasta_file, "\n")
    cat("🔍 DEBUG: FASTA file exists:", file.exists(fasta_file), "\n")
    
    # Load novel peptide generator if not already available  
    if (!exists("generate_novel_peptide_data")) {
      cat("Loading novel peptide generator...\n")
      source("novel_peptide_generator.R")
    }
    
    # Generate peptide data using EXISTING generator
    peptide_results_no_miss <- generate_novel_peptide_data(
      gtf_file = gtf_file,
      fasta_file = fasta_file,
      missedCleavages = 0,
      minLength = 7,
      maxLength = 60
    )
    
    peptide_results_2miss <- generate_novel_peptide_data(
      gtf_file = gtf_file,
      fasta_file = fasta_file,
      missedCleavages = 2,
      minLength = 7,
      maxLength = 60
    )
    
    # STEP 9: Save results (maintain compatibility)
    cat("Step 9: Saving results...\n")
    results_dir <- file.path(output_dir, "results")
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save both dataframes
    dataframe_file_no_miss <- file.path(results_dir, "rmats_isoform_dataframe_no_miss.rds")
    dataframe_file_2miss <- file.path(results_dir, "rmats_isoform_dataframe_2miss.rds")
    
    saveRDS(peptide_results_no_miss, dataframe_file_no_miss)
    saveRDS(peptide_results_2miss, dataframe_file_2miss)
    
    # Determine current selection
    if (missedCleavages == 0) {
      current_peptide_results <- peptide_results_no_miss
      current_dataframe_file <- dataframe_file_no_miss
    } else {
      current_peptide_results <- peptide_results_2miss
      current_dataframe_file <- dataframe_file_2miss
    }
    
    # Save current selection for backward compatibility
    original_dataframe_file <- file.path(results_dir, "rmats_isoform_dataframe.rds")
    saveRDS(current_peptide_results, original_dataframe_file)
    
    # Create pipeline results structure
    pipeline_results <- list(
      success = TRUE,  # CRITICAL: Add success field for UI
      gtf_file = gtf_file,  # CRITICAL: Add GTF file path for visualization
      fasta_file = fasta_file,  # CRITICAL: Add FASTA file path for visualization
      dataframe_file = current_dataframe_file,
      dataframe_file_no_miss = dataframe_file_no_miss,
      dataframe_file_2miss = dataframe_file_2miss,
      original_dataframe_file = original_dataframe_file,
      output_dir = output_dir,
      timestamp = timestamp,
      missedCleavages = missedCleavages,
      event_info = list(
        event_type = event_type,
        gene_id = selected_event$GeneID,
        gene_symbol = selected_event$geneSymbol,
        coordinates = event_coords
      )
    )
    
    cat("✅ STREAMLINED rMATS analysis completed successfully!\n")
    cat("📁 Results directory:", output_dir, "\n")
    cat("📊 Generated protein isoforms:\n")
    cat("   - No missed cleavages:", nrow(peptide_results_no_miss), "rows\n")
    cat("   - Up to 2 missed cleavages:", nrow(peptide_results_2miss), "rows\n")
    cat("🚀 SINGLE CDS search eliminated state corruption!\n")
    
    # DEBUG: Confirm final result structure
    cat("🔍 DEBUG: Returning pipeline_results with:\n")
    cat("  - success:", pipeline_results$success, "\n")
    cat("  - gtf_file:", pipeline_results$gtf_file, "\n")
    cat("  - gtf_file exists:", file.exists(pipeline_results$gtf_file), "\n")
    
    # Return flattened structure so server can access gtf_file directly
    return(list(
      success = TRUE,
      gtf_file = pipeline_results$gtf_file,
      fasta_file = pipeline_results$fasta_file,
      dataframe_file = pipeline_results$dataframe_file,
      dataframe_file_no_miss = pipeline_results$dataframe_file_no_miss,
      dataframe_file_2miss = pipeline_results$dataframe_file_2miss,
      output_dir = pipeline_results$output_dir,
      peptide_data = current_peptide_results,
      peptide_data_no_miss = peptide_results_no_miss,
      peptide_data_2miss = peptide_results_2miss,
      event_info = pipeline_results$event_info,
      summary = paste("Successfully analyzed", event_type, "event for gene", selected_event$geneSymbol, "with streamlined pipeline")
    ))
    
  }, error = function(e) {
    cat("❌ Error in streamlined rMATS analysis:", conditionMessage(e), "\n")
    return(list(
      success = FALSE,
      error = conditionMessage(e)
    ))
  })
}

#===============================================================================
# ORIGINAL COMPREHENSIVE ANALYSIS FUNCTION (DEPRECATED - CAUSES CDS CORRUPTION)
#===============================================================================

analyze_rmats_event_with_miscleavage <- function(selected_event, event_type, missedCleavages = 0) {
  
  cat("🔄 Starting comprehensive rMATS analysis...\n")
  cat("Event Type:", event_type, "\n")
  cat("Gene:", selected_event$geneSymbol, "ID:", selected_event$GeneID, "\n")
  
  tryCatch({
    
    # Check that required functions are available and source if needed
    if (!exists("extract_event_coordinates")) {
      stop("rMATS modules not loaded. This function must be called from rMATS server context.")
    }
    
    # Load GTF generator if not already available
    if (!exists("generate_rmats_gtf_with_phases")) {
      cat("Loading rMATS GTF generator...\n")
      source("rmats_exact_gtf_generator.R")
    }
    
    # Load novel peptide generator if not already available  
    if (!exists("generate_novel_peptide_data")) {
      cat("Loading novel peptide generator...\n")
      source("novel_peptide_generator.R")
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
    
    # STEP 3: Generate peptide data for BOTH miscleavage settings (like Novel Isoform)
    cat("Step 3: Generating peptide data for both miscleavage settings...\n")
    gtf_file <- file.path(output_dir, "novel_transcript_nt.transdecoder.genome.gtf")
    fasta_file <- file.path(output_dir, "novel_transcript_nt.fa.transdecoder.pep")
    
    if (!file.exists(gtf_file) || !file.exists(fasta_file)) {
      return(list(
        success = FALSE,
        error = "GTF or FASTA file not found after generation"
      ))
    }
    
    # Generate peptide data for no missed cleavages
    cat("=== Processing 0 Missed Cleavages ===\n")
    peptide_results_no_miss <- generate_novel_peptide_data(
      gtf_file = gtf_file,
      fasta_file = fasta_file,
      missedCleavages = 0,
      minLength = 7,
      maxLength = 60
    )
    
    # Generate peptide data for up to 2 missed cleavages
    cat("=== Processing 2 Missed Cleavages ===\n")
    peptide_results_2miss <- generate_novel_peptide_data(
      gtf_file = gtf_file,
      fasta_file = fasta_file,
      missedCleavages = 2,
      minLength = 7,
      maxLength = 60
    )
    
    # STEP 4: Save BOTH results in expected format
    cat("Step 4: Saving both miscleavage result files...\n")
    results_dir <- file.path(output_dir, "results")
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save both dataframes with specific filenames
    dataframe_file_no_miss <- file.path(results_dir, "rmats_isoform_dataframe_no_miss.rds")
    dataframe_file_2miss <- file.path(results_dir, "rmats_isoform_dataframe_2miss.rds")
    
    saveRDS(peptide_results_no_miss, dataframe_file_no_miss)
    saveRDS(peptide_results_2miss, dataframe_file_2miss)
    
    # Determine which file to use as the primary result based on user's current selection
    if (missedCleavages == 0) {
      current_peptide_results <- peptide_results_no_miss
      current_dataframe_file <- dataframe_file_no_miss
    } else {
      current_peptide_results <- peptide_results_2miss
      current_dataframe_file <- dataframe_file_2miss
    }
    
    # Also save the current selection with the original filename for backward compatibility
    original_dataframe_file <- file.path(results_dir, "rmats_isoform_dataframe.rds")
    saveRDS(current_peptide_results, original_dataframe_file)
    
    # Create pipeline results structure (matching expected format)
    pipeline_results <- list(
      dataframe_file = current_dataframe_file,  # Current selection for backward compatibility
      dataframe_file_no_miss = dataframe_file_no_miss,  # CRITICAL: Add both file paths
      dataframe_file_2miss = dataframe_file_2miss,      # CRITICAL: Add both file paths
      original_dataframe_file = original_dataframe_file,
      gtf_file = gtf_file,
      fasta_file = fasta_file,
      output_dir = output_dir,
      timestamp = timestamp,
      missedCleavages = missedCleavages,
      event_info = list(
        event_type = event_type,
        gene_id = selected_event$GeneID,
        gene_symbol = selected_event$geneSymbol,
        coordinates = event_coords
      )
    )
    
    cat("✅ Comprehensive rMATS analysis completed successfully!\n")
    cat("📁 Results directory:", output_dir, "\n")
    cat("📊 Generated protein isoforms:\n")
    cat("   - No missed cleavages:", nrow(peptide_results_no_miss), "rows\n")
    cat("   - Up to 2 missed cleavages:", nrow(peptide_results_2miss), "rows\n")
    cat("📁 Saved files:\n")
    cat("   - No miss:", dataframe_file_no_miss, "\n")
    cat("   - 2 miss:", dataframe_file_2miss, "\n")
    
    return(list(
      success = TRUE,
      pipeline_results = pipeline_results,
      peptide_data = current_peptide_results,  # Return current selection
      peptide_data_no_miss = peptide_results_no_miss,  # Also return both datasets
      peptide_data_2miss = peptide_results_2miss,      # Also return both datasets
      summary = paste("Successfully analyzed", event_type, "event for gene", selected_event$geneSymbol, "with both miscleavage settings")
    ))
    
  }, error = function(e) {
    cat("❌ Error in comprehensive rMATS analysis:", conditionMessage(e), "\n")
    return(list(
      success = FALSE,
      error = conditionMessage(e)
    ))
  })
}

# NOTE: Removed custom peptide generation function - using existing novel_peptide_generator.R instead

cat("✅ rMATS comprehensive analysis module loaded\n")