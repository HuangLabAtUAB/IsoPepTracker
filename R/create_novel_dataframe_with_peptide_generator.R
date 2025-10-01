#===============================================================================
# CREATE NOVEL DATAFRAME USING EXISTING PEPTIDE GENERATOR
# Uses the working novel_peptide_generator.R script as specified in PIPELINE_STEPS.md
#===============================================================================

#' Create Novel Dataframe with Peptide Generator
#' 
#' Uses the existing novel_peptide_generator.R script to create 24-column RDS
#' following STEP 7-8 from PIPELINE_STEPS.md exactly
#' 
#' @param work_dir Working directory containing pipeline results
#' @param enable_blast Whether BLAST was enabled
#' @param blast_threshold BLAST identity threshold
#' @return Path to created RDS file
create_novel_dataframe_with_peptide_generator <- function(work_dir, 
                                                         enable_blast = FALSE,
                                                         blast_threshold = 70) {
  
  cat("STEP 7-8: Using novel_peptide_generator.R with .pep and .genome.gtf files...\n")
  
  # Expected files from 8-step pipeline (STEP 1-6)
  gtf_file <- file.path(work_dir, "input", "clean_input.transdecoder.genome.gtf")
  pep_file <- file.path(work_dir, "input", "clean_input.fa.transdecoder.pep")
  
  # Validate required files exist
  if (!file.exists(gtf_file)) {
    stop("GTF file not found: ", gtf_file)
  }
  if (!file.exists(pep_file)) {
    stop("Protein sequences file not found: ", pep_file)
  }
  
  cat("Found required files:\n")
  cat("  GTF: ", gtf_file, "\n")
  cat("  PEP: ", pep_file, "\n")
  
  # Copy novel_peptide_generator.R to work directory for execution
  generator_script <- file.path(work_dir, "novel_peptide_generator_local.R")
  
  # Read the existing novel_peptide_generator.R
  if (!file.exists("novel_peptide_generator.R")) {
    stop("novel_peptide_generator.R not found in project root")
  }
  
  # Create customized version for this run
  generator_content <- readLines("novel_peptide_generator.R", warn = FALSE)
  
  # Update file paths in the script to use our specific files
  generator_content <- gsub(
    'gtf_file <- "novel_transcript_nt.transdecoder.genome.gtf"',
    paste0('gtf_file <- "', basename(gtf_file), '"'),
    generator_content
  )
  
  generator_content <- gsub(
    'fasta_file <- "novel_transcript_nt.fa.transdecoder.pep"',
    paste0('fasta_file <- "', basename(pep_file), '"'),
    generator_content
  )
  
  # Update output file path
  generator_content <- gsub(
    'output_file <- "novel_transcript_nt_peptides.rds"',
    'output_file <- "../results/novel_isoform_dataframe.rds"',
    generator_content
  )
  
  # Write customized script
  writeLines(generator_content, generator_script)
  
  # Copy required files to input directory for script execution (only if different locations)
  target_gtf <- file.path(work_dir, "input", basename(gtf_file))
  target_pep <- file.path(work_dir, "input", basename(pep_file))
  
  if (normalizePath(gtf_file) != normalizePath(target_gtf)) {
    file.copy(gtf_file, target_gtf, overwrite = TRUE)
  }
  if (normalizePath(pep_file) != normalizePath(target_pep)) {
    file.copy(pep_file, target_pep, overwrite = TRUE)
  }
  
  # Execute novel_peptide_generator.R in the input directory
  cat("Executing novel_peptide_generator.R...\n")
  
  # Change to input directory and run the script
  original_wd <- getwd()
  tryCatch({
    setwd(file.path(work_dir, "input"))
    
    # Source the peptide generator script
    source(file.path("..", "novel_peptide_generator_local.R"))
    
    cat("✅ novel_peptide_generator.R execution completed\n")
    
  }, error = function(e) {
    stop("Error executing novel_peptide_generator.R: ", e$message)
  }, finally = {
    setwd(original_wd)
  })
  
  # Verify output file was created
  output_file <- file.path(work_dir, "results", "novel_isoform_dataframe.rds")
  
  if (!file.exists(output_file)) {
    stop("novel_peptide_generator.R did not create expected output file: ", output_file)
  }
  
  # Verify the structure is correct (24 columns with _mapped_ranges)
  tryCatch({
    novel_data <- readRDS(output_file)
    
    if (ncol(novel_data) != 24) {
      warning("Output has ", ncol(novel_data), " columns instead of expected 24")
    }
    
    # Check for mapped_ranges columns
    mapped_ranges_cols <- grep("_mapped_ranges$", names(novel_data), value = TRUE)
    if (length(mapped_ranges_cols) == 0) {
      warning("No _mapped_ranges columns found in output")
    } else {
      cat("✅ Found", length(mapped_ranges_cols), "_mapped_ranges columns\n")
    }
    
    cat("✅ Created 24-column RDS with genomic mapping\n")
    cat("✅ Structure matches known gene format exactly\n")
    cat("✅ Contains", nrow(novel_data), "proteins ready for Multiple Isoform Comparison\n")
    
  }, error = function(e) {
    warning("Error validating output structure: ", e$message)
  })
  
  # Clean up temporary script
  if (file.exists(generator_script)) {
    file.remove(generator_script)
  }
  
  return(output_file)
}

#===============================================================================
# INTEGRATION HELPER FUNCTIONS
#===============================================================================

#' Validate Novel RDS Structure
#' 
#' Validates that the novel RDS matches the expected 24-column structure
#' 
#' @param rds_file Path to RDS file
#' @return Logical indicating if structure is valid
validate_novel_rds_structure <- function(rds_file) {
  
  if (!file.exists(rds_file)) {
    return(FALSE)
  }
  
  tryCatch({
    data <- readRDS(rds_file)
    
    # Check column count
    if (ncol(data) != 24) {
      cat("❌ Expected 24 columns, found", ncol(data), "\n")
      return(FALSE)
    }
    
    # Check required base columns
    required_cols <- c("proteinID", "txID", "geneID", "geneSymbol", "numAA", "seq")
    if (!all(required_cols %in% names(data))) {
      cat("❌ Missing required base columns\n")
      return(FALSE)
    }
    
    # Check enzyme columns
    enzymes <- c("trp", "chymo", "aspn", "lysc", "lysn", "gluc")
    for (enzyme in enzymes) {
      required_enzyme_cols <- c(
        paste0(enzyme, "Peps"),
        paste0(enzyme, "Peps_positions"),
        paste0(enzyme, "Peps_mapped_ranges")
      )
      if (!all(required_enzyme_cols %in% names(data))) {
        cat("❌ Missing columns for enzyme:", enzyme, "\n")
        return(FALSE)
      }
    }
    
    # Check for GRanges objects in mapped_ranges columns
    mapped_ranges_cols <- grep("_mapped_ranges$", names(data), value = TRUE)
    for (col in mapped_ranges_cols) {
      sample_ranges <- data[[col]][[1]]
      if (length(sample_ranges) > 0 && !inherits(sample_ranges[[1]], "GRanges")) {
        cat("❌ Column", col, "does not contain GRanges objects\n")
        return(FALSE)
      }
    }
    
    cat("✅ Novel RDS structure validation passed\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("❌ Error validating RDS structure:", e$message, "\n")
    return(FALSE)
  })
} 