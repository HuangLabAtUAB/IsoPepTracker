# Wrapper to integrate the working rMATS peptide generator with the server
# This uses the truly simple generator that we know works

source("rmats_peptide_generator_truly_simple.R")

#' Generate rMATS peptide data using the working truly simple generator
#' 
#' @param gtf_file Path to rMATS GTF file  
#' @param pep_file Path to rMATS peptide FASTA file
#' @return data.frame with proper 24-column structure
generate_rmats_peptides_wrapper <- function(gtf_file, pep_file) {
  
  # Use the working generator directly
  result <- generate_novel_peptide_data(gtf_file, pep_file)
  
  return(result)
}