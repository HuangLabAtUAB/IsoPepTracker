# rMATS Peptide Analysis System
# Dedicated module for rMATS alternative splicing peptide analysis
# This module creates proper 24-column RDS files for rMATS inclusion/exclusion isoforms

library(stringr)
library(data.table)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(Biostrings)
library(IRanges)
library(cleaver)
library(seqinr)

#===============================================================================
# rMATS PEPTIDE ANALYSIS FUNCTIONS
#===============================================================================

#' Create rMATS peptide dataframe with proper 24-column structure
#' 
#' @param inclusion_protein Character string of inclusion isoform protein sequence
#' @param exclusion_protein Character string of exclusion isoform protein sequence
#' @param gene_id Character string of gene ID
#' @param gene_symbol Character string of gene symbol
#' @param gtf_file Path to rMATS GTF file with CDS coordinates
#' @return data.frame with 24 columns matching novel isoform structure
create_rmats_peptide_dataframe <- function(inclusion_protein, exclusion_protein, gene_id, gene_symbol, gtf_file) {
  
  cat("=== rMATS PEPTIDE DATAFRAME GENERATION (using truly simple generator) ===\n")
  cat("Gene ID:", gene_id, "\n")
  cat("Gene Symbol:", gene_symbol, "\n")
  cat("GTF file:", gtf_file, "\n")
  
  # Use the working generator that we know works!
  # Robust path resolution that works from any subdirectory
  
  # First, try to find the app root directory by looking for characteristic files
  current_dir <- getwd()
  app_root <- NULL
  
  # Look for the generator file going up directory levels
  test_dirs <- c(
    current_dir,
    dirname(current_dir),
    dirname(dirname(current_dir)),
    "/Users/Mahmuda/Desktop/AS_peptides/APP/final_app"  # Fallback absolute path
  )
  
  for (test_dir in test_dirs) {
    potential_path <- file.path(test_dir, "rmats_peptide_generator_truly_simple.R")
    if (file.exists(potential_path)) {
      script_path <- potential_path
      app_root <- test_dir
      break
    }
  }
  
  if (is.null(app_root)) {
    stop("Cannot find rmats_peptide_generator_truly_simple.R in any expected location")
  }
  
  cat("Using generator from:", script_path, "\n")
  source(script_path, local = TRUE)
  
  # Create temporary FASTA file with both proteins
  temp_fasta <- tempfile(fileext = ".fa")
  fasta_content <- paste0(
    '>"', gene_id, '".inclusion GENE."', gene_symbol, '"~~"', gene_id, '".inclusion  ORF type:complete len:', nchar(inclusion_protein), ' (+),score=100.00\n',
    inclusion_protein, '\n',
    '>"', gene_id, '".exclusion GENE."', gene_symbol, '"~~"', gene_id, '".exclusion  ORF type:complete len:', nchar(exclusion_protein), ' (+),score=100.00\n',
    exclusion_protein
  )
  writeLines(fasta_content, temp_fasta)
  
  # Generate peptides using the working generator
  result <- generate_novel_peptide_data(gtf_file, temp_fasta)
  
  # Clean up
  unlink(temp_fasta)
  
  return(result)
}