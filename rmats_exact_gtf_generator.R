# EXACT rMATS GTF Generator with Phase Information
# Uses the EXACT implementation from rmats_peptide system
# Supports ALL event types: SE, A3SS, A5SS, MXE, RI

library(data.table)

# Source the EXACT working functions from rmats_peptide
source('rmats_peptide/utils/coordinate_utils.R')
source('rmats_peptide/modules/rmats_parser.R')  
source('rmats_peptide/modules/flanking_exons.R')
source('rmats_peptide/modules/protein_translation.R')

#===============================================================================
# EXACT rMATS GTF GENERATOR FUNCTION
#===============================================================================

generate_rmats_gtf_with_phases <- function(rmats_event, event_type, output_dir = NULL) {
  
  cat("=== EXACT rMATS GTF GENERATOR ===\n")
  cat("Event Type:", event_type, "\n")
  cat("Gene:", rmats_event$geneSymbol, "ID:", rmats_event$GeneID, "\n")
  
  if (is.null(output_dir)) {
    output_dir <- paste0("rmats_gtf_", event_type, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  tryCatch({
    #===========================================================================
    # STEP 1: Extract event coordinates using EXACT rmats_peptide function
    #===========================================================================
    
    cat("Step 1: Extract coordinates for", event_type, "event...\n")
    event_coords <- extract_event_coordinates(rmats_event, event_type)
    
    #===========================================================================
    # STEP 2: Build GTF structures using EXACT rmats_peptide function
    #===========================================================================
    
    cat("Step 2: Build GTF structures...\n") 
    gtf_structures <- build_gtf_structures(event_coords, event_type)
    
    #===========================================================================
    # STEP 3: Identify flanking exons using EXACT rmats_peptide function
    #===========================================================================
    
    cat("Step 3: Identify flanking exons...\n")
    flanking_exons <- identify_flanking_exons(gtf_structures)
    
    #===========================================================================
    # STEP 4: Search CDS index using EXACT rmats_peptide function
    #===========================================================================
    
    cat("Step 4: Search CDS index for exact matches...\n")
    gtf_results <- search_all_exons_in_cds(flanking_exons, 'rmats_peptide/real_cds_index.RDS')
    
    if (!gtf_results$ready_for_step5) {
      stop("CDS search failed: ", gtf_results$outcome)
    }
    
    cat("CDS search successful:", gtf_results$cds_summary, "\n")
    
    #===========================================================================
    # STEP 5: Extract phase information using EXACT rmats_peptide function
    #===========================================================================
    
    cat("Step 5: Extract phase information...\n")
    phase_results <- extract_phase_information(gtf_results)
    
    if (!phase_results$ready_for_step6) {
      cat("WARNING: No translatable isoforms found\n")
      return(list(success = FALSE, message = "No translatable isoforms"))
    }
    
    cat("Phase extraction successful!\n")
    cat("Inclusion translatable:", phase_results$inclusion_translatable, "\n") 
    cat("Exclusion translatable:", phase_results$exclusion_translatable, "\n")
    
    #===========================================================================
    # STEP 6: Generate protein sequences using EXACT rmats_peptide function
    #===========================================================================
    
    cat("Step 6: Generate protein sequences...\n")
    translation_results <- analyze_protein_translation(phase_results)
    
    if (!translation_results$success) {
      cat("WARNING: Translation failed:", translation_results$summary, "\n")
    }
    
    cat("Translation case type:", translation_results$case_type, "\n")
    cat("Functional consequence:", translation_results$functional_consequence, "\n")
    
    #===========================================================================
    # OUTPUT GENERATION: Create GTF and FASTA files with EXACT format
    #===========================================================================
    
    results <- list(
      success = TRUE,
      event_type = event_type,
      gene_id = rmats_event$GeneID,
      gene_symbol = rmats_event$geneSymbol,
      output_dir = output_dir,
      files_created = c()
    )
    
    # Generate GTF file with phase information
    gtf_output <- generate_complete_gtf(
      phase_results, 
      rmats_event, 
      event_type,
      file.path(output_dir, "novel_transcript_nt.transdecoder.genome.gtf")
    )
    
    if (gtf_output$success) {
      results$files_created <- c(results$files_created, gtf_output$filename)
      cat("GTF file created:", gtf_output$filename, "\n")
      cat("Phase values preserved:", gtf_output$phase_count, "CDS entries\n")
    }
    
    # Generate FASTA file with protein sequences
    fasta_output <- generate_protein_fasta(
      translation_results,
      rmats_event,
      file.path(output_dir, "novel_transcript_nt.fa.transdecoder.pep")
    )
    
    if (fasta_output$success) {
      results$files_created <- c(results$files_created, fasta_output$filename)
      cat("FASTA file created:", fasta_output$filename, "\n") 
      cat("Protein sequences:", fasta_output$protein_count, "\n")
    }
    
    # Generate event-specific GTF file
    event_gtf_output <- generate_event_gtf(
      phase_results,
      rmats_event,
      event_type,
      file.path(output_dir, paste0('rmats_', rmats_event$GeneID, '_', event_type, '.transdecoder.genome.gtf'))
    )
    
    if (event_gtf_output$success) {
      results$files_created <- c(results$files_created, event_gtf_output$filename)
      cat("Event-specific GTF created:", event_gtf_output$filename, "\n")
    }
    
    cat("\n=== GTF GENERATION COMPLETE ===\n")
    cat("Output directory:", output_dir, "\n")
    cat("Files created:", length(results$files_created), "\n")
    
    return(results)
    
  }, error = function(e) {
    cat("ERROR in GTF generation:", e$message, "\n")
    return(list(success = FALSE, message = e$message))
  })
}

#===============================================================================
# OUTPUT GENERATION FUNCTIONS
#===============================================================================

generate_complete_gtf <- function(phase_results, rmats_event, event_type, filename) {
  
  cat("Generating complete GTF file with phase information...\n")
  
  tryCatch({
    all_gtf_lines <- c()
    phase_count <- 0
    
    # Generate inclusion isoform GTF
    if (!is.null(phase_results$inclusion_cds_gtf) && nrow(phase_results$inclusion_cds_gtf) > 0) {
      
      inclusion_lines <- generate_isoform_gtf(
        phase_results$inclusion_cds_gtf,
        rmats_event,
        "inclusion",
        event_type
      )
      
      all_gtf_lines <- c(all_gtf_lines, inclusion_lines$gtf_lines)
      phase_count <- phase_count + inclusion_lines$cds_count
    }
    
    # Generate exclusion isoform GTF
    if (!is.null(phase_results$exclusion_cds_gtf) && nrow(phase_results$exclusion_cds_gtf) > 0) {
      
      exclusion_lines <- generate_isoform_gtf(
        phase_results$exclusion_cds_gtf,
        rmats_event,
        "exclusion", 
        event_type
      )
      
      all_gtf_lines <- c(all_gtf_lines, exclusion_lines$gtf_lines)
      phase_count <- phase_count + exclusion_lines$cds_count
    }
    
    if (length(all_gtf_lines) == 0) {
      return(list(success = FALSE, message = "No GTF lines generated"))
    }
    
    # Write GTF file
    writeLines(all_gtf_lines, filename)
    
    return(list(
      success = TRUE,
      filename = filename,
      phase_count = phase_count,
      line_count = length(all_gtf_lines)
    ))
    
  }, error = function(e) {
    return(list(success = FALSE, message = e$message))
  })
}

generate_isoform_gtf <- function(cds_gtf, rmats_event, isoform_type, event_type) {
  
  gtf_lines <- c()
  
  # Create transcript ID and gene ID
  transcript_id <- paste0('"', rmats_event$GeneID, '"', ".", isoform_type)
  gene_id <- paste0('"', rmats_event$geneSymbol, '"', ".1")
  
  # Get chromosome and strand
  chromosome <- rmats_event$chr
  strand <- rmats_event$strand
  
  # Determine transcript boundaries from CDS
  transcript_start <- min(cds_gtf$start)
  transcript_end <- max(cds_gtf$end)
  
  # Add transcript line
  gtf_lines <- c(gtf_lines,
    paste(chromosome, "rmats", "transcript", 
          transcript_start, transcript_end, ".", strand, ".",
          paste0('transcript_id ', transcript_id, '; gene_id ', gene_id),
          sep = "\t")
  )
  
  # Add exon lines (derive from CDS coordinates)
  exon_coords <- unique(data.frame(
    start = cds_gtf$start,
    end = cds_gtf$end
  ))
  
  for (i in seq_len(nrow(exon_coords))) {
    gtf_lines <- c(gtf_lines,
      paste(chromosome, "rmats", "exon",
            exon_coords$start[i], exon_coords$end[i], ".", strand, ".",
            paste0('transcript_id ', transcript_id, '; gene_id ', gene_id, ';'),
            sep = "\t")
    )
  }
  
  # Add CDS lines with EXACT phase information from rmats_peptide
  for (i in seq_len(nrow(cds_gtf))) {
    cds_row <- cds_gtf[i, ]
    
    gtf_lines <- c(gtf_lines,
      paste(cds_row$seqname, "rmats", "CDS",
            cds_row$start, cds_row$end, ".", cds_row$strand, 
            cds_row$frame,  # CRITICAL: This preserves the EXACT phase information
            paste0('transcript_id ', transcript_id, '; gene_id ', gene_id, ';'),
            sep = "\t")
    )
  }
  
  return(list(
    gtf_lines = gtf_lines,
    cds_count = nrow(cds_gtf)
  ))
}

generate_event_gtf <- function(phase_results, rmats_event, event_type, filename) {
  
  cat("Generating event-specific GTF file...\n")
  
  # This creates the rmats_GENEID_EVENTTYPE.gtf file format
  result <- generate_complete_gtf(phase_results, rmats_event, event_type, filename)
  return(result)
}

generate_protein_fasta <- function(translation_results, rmats_event, filename) {
  
  cat("Generating protein FASTA file...\n")
  
  tryCatch({
    fasta_lines <- c()
    protein_count <- 0
    
    # Add inclusion protein
    if (!is.null(translation_results$inclusion_protein) && nchar(translation_results$inclusion_protein) > 0) {
      
      header <- paste0('>', rmats_event$GeneID, '.inclusion GENE.', rmats_event$geneSymbol, 
                      '~~', rmats_event$GeneID, '.inclusion  ORF type:complete len:', 
                      nchar(translation_results$inclusion_protein), ' (+),score=100.00 ', 
                      rmats_event$GeneID, ':1-', nchar(translation_results$inclusion_protein)*3, '(+)')
      
      fasta_lines <- c(fasta_lines, header, translation_results$inclusion_protein)
      protein_count <- protein_count + 1
    }
    
    # Add exclusion protein
    if (!is.null(translation_results$exclusion_protein) && nchar(translation_results$exclusion_protein) > 0) {
      
      header <- paste0('>', rmats_event$GeneID, '.exclusion GENE.', rmats_event$geneSymbol,
                      '~~', rmats_event$GeneID, '.exclusion  ORF type:complete len:',
                      nchar(translation_results$exclusion_protein), ' (+),score=100.00 ',
                      rmats_event$GeneID, ':1-', nchar(translation_results$exclusion_protein)*3, '(+)')
      
      fasta_lines <- c(fasta_lines, header, translation_results$exclusion_protein)
      protein_count <- protein_count + 1
    }
    
    if (length(fasta_lines) == 0) {
      return(list(success = FALSE, message = "No protein sequences to write"))
    }
    
    writeLines(fasta_lines, filename)
    
    return(list(
      success = TRUE,
      filename = filename,
      protein_count = protein_count
    ))
    
  }, error = function(e) {
    return(list(success = FALSE, message = e$message))
  })
}

#===============================================================================
# EXAMPLE USAGE FOR ALL EVENT TYPES
#===============================================================================

if (FALSE) {  # Set to TRUE to run examples
  
  # Example SE event - replace with actual gene parameters
  se_event <- data.frame(
    GeneID = 'YOUR_GENE_ID_HERE',    # e.g. 'ENSG00000099721.15'
    geneSymbol = 'YOUR_GENE_SYMBOL', # e.g. 'AMELY' 
    chr = 'YOUR_CHROMOSOME',         # e.g. 'chrY'
    strand = 'YOUR_STRAND',          # '+' or '-'
    exonStart_0base = 0,             # Actual exon start coordinate
    exonEnd = 0,                     # Actual exon end coordinate
    upstreamES = 0,                  # Upstream exon start
    upstreamEE = 0,                  # Upstream exon end
    downstreamES = 0,                # Downstream exon start
    downstreamEE = 0,                # Downstream exon end
    stringsAsFactors = FALSE
  )
  
  # Generate GTF with phases for SE event
  se_result <- generate_rmats_gtf_with_phases(se_event, "SE")
  
  # Add examples for other event types (A3SS, A5SS, MXE, RI) here...
}

cat("=== rMATS EXACT GTF GENERATOR LOADED ===\n")
cat("Usage: generate_rmats_gtf_with_phases(rmats_event, event_type)\n")
cat("Supports event types: SE, A3SS, A5SS, MXE, RI\n")
cat("Uses EXACT rmats_peptide implementation with proper phase handling\n")