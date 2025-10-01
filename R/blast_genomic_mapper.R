#===============================================================================
# BLAST GENOMIC MAPPER
# Safe integration of novel_peptide_generator mapping logic for BLAST peptides
#===============================================================================

#' Map BLAST Peptide to Genomic Coordinates
#' 
#' Uses proven genomic mapping logic from novel_peptide_generator.R to map
#' BLAST search peptides to precise genomic coordinates on transcripts
#' 
#' @param blast_peptide Peptide sequence from BLAST search
#' @param transcript_id Target transcript ID from BLAST results
#' @param gene_id Gene ID from BLAST results
#' @param transcript_structure Transcript structure from load_transcript_exons()
#' @return List with genomic_ranges (GRanges object) or NULL if mapping fails
map_blast_peptide_to_transcript <- function(blast_peptide, transcript_id, gene_id, transcript_structure) {
  
  cat("DEBUG: map_blast_peptide_to_transcript called\n")
  cat("DEBUG: Peptide:", blast_peptide, "\n")
  cat("DEBUG: Transcript:", transcript_id, "\n")
  cat("DEBUG: Gene:", gene_id, "\n")
  
  tryCatch({
    # Step 1: Get protein sequence for transcript
    protein_sequence <- get_protein_sequence_for_transcript(transcript_id, gene_id)
    
    if (is.null(protein_sequence)) {
      cat("DEBUG: Could not get protein sequence for", transcript_id, "\n")
      return(NULL)
    }
    
    cat("DEBUG: Protein sequence length:", nchar(protein_sequence), "AA\n")
    
    # Step 2: Find peptide position in protein (using novel_peptide_generator logic)
    peptide_positions <- stringr::str_locate_all(protein_sequence, stringr::fixed(blast_peptide))[[1]]
    
    if (nrow(peptide_positions) == 0) {
      cat("DEBUG: Peptide not found in protein sequence:", blast_peptide, "\n")
      return(NULL)
    }
    
    cat("DEBUG: Found", nrow(peptide_positions), "occurrences of peptide in protein\n")
    
    # Step 3: Get CDS data from transcript structure
    if (is.null(transcript_structure) || !transcript_structure$success) {
      cat("DEBUG: Invalid transcript structure\n")
      return(NULL)
    }
    
    cds_data <- transcript_structure$cds[[transcript_id]]
    
    if (is.null(cds_data) || length(cds_data) == 0) {
      cat("DEBUG: No CDS data found for", transcript_id, "\n")
      return(NULL)
    }
    
    cat("DEBUG: Found", length(cds_data), "CDS segments for transcript\n")
    
    # Step 4: Apply novel_peptide_generator mapping pipeline
    cds_info <- retrieve_sorted_cds_safe(transcript_id, cds_data)
    
    if (is.null(cds_info)) {
      cat("DEBUG: Failed to retrieve sorted CDS\n")
      return(NULL)
    }
    
    cds_cumulative <- calculate_cumulative_cds_safe(
      cds_info$cds_tx, 
      mcols(cds_info$cds_tx)$phase, 
      cds_info$strand_tx
    )
    
    if (is.null(cds_cumulative)) {
      cat("DEBUG: Failed to calculate cumulative CDS\n")
      return(NULL)
    }
    
    # Step 5: Map each peptide occurrence to genomic coordinates
    genomic_ranges_list <- list()
    
    for (i in 1:nrow(peptide_positions)) {
      peptide_info <- list(
        peptide = blast_peptide,
        aa_start = peptide_positions[i, "start"],
        aa_end = peptide_positions[i, "end"]
      )
      
      cat("DEBUG: Mapping peptide occurrence", i, "at AA positions", 
          peptide_info$aa_start, "-", peptide_info$aa_end, "\n")
      
      genomic_ranges <- map_peptide_to_genomic_safe(
        peptide_info, cds_cumulative, cds_info$strand_tx, 
        as.character(seqnames(cds_data)[1])
      )
      
      if (!is.null(genomic_ranges)) {
        genomic_ranges_list[[i]] <- genomic_ranges
        cat("DEBUG: Successfully mapped occurrence", i, "to", length(genomic_ranges), "genomic segments\n")
      }
    }
    
    # Step 6: Combine all genomic ranges
    if (length(genomic_ranges_list) > 0) {
      all_ranges <- do.call(c, genomic_ranges_list)
      
      cat("DEBUG: Final mapping result:", length(all_ranges), "total genomic segments\n")
      
      # Create structure compatible with existing visualization system
      return(list(
        genomic_ranges = all_ranges,
        success = TRUE
      ))
    }
    
    cat("DEBUG: No successful mappings\n")
    return(NULL)
    
  }, error = function(e) {
    cat("ERROR in map_blast_peptide_to_transcript:", e$message, "\n")
    return(NULL)
  })
}

#' Get Protein Sequence for Transcript
#' 
#' Retrieves protein sequence for a transcript, trying multiple sources
#' 
#' @param transcript_id Transcript ID
#' @param gene_id Gene ID (for fallback approaches)
#' @return Character string with protein sequence or NULL
get_protein_sequence_for_transcript <- function(transcript_id, gene_id) {
  
  tryCatch({
    cat("DEBUG: Loading protein sequence for transcript", transcript_id, "in gene", gene_id, "\n")
    
    # Method 1: Load gene RDS file directly from data/genes/
    gene_file <- paste0("data/genes/", gene_id, ".rds")
    
    if (file.exists(gene_file)) {
      cat("DEBUG: Loading gene file:", gene_file, "\n")
      gene_data <- readRDS(gene_file)
      
      if (is.data.frame(gene_data) && "txID" %in% names(gene_data) && "seq" %in% names(gene_data)) {
        cat("DEBUG: Loaded gene data with", nrow(gene_data), "transcript entries\n")
        
        # Look for exact transcript match first
        tx_rows <- which(gene_data$txID == transcript_id)
        
        if (length(tx_rows) > 0) {
          protein_seq <- gene_data$seq[tx_rows[1]]
          if (!is.null(protein_seq) && nchar(protein_seq) > 0) {
            cat("DEBUG: Found protein sequence for", transcript_id, "- length:", nchar(protein_seq), "AA\n")
            return(as.character(protein_seq))
          }
        } else {
          # Try transcript ID without version (e.g., ENST00000379409 instead of ENST00000379409.6)
          transcript_base <- gsub("\\..*", "", transcript_id)
          cat("DEBUG: Exact match failed, trying base transcript ID:", transcript_base, "\n")
          
          # Look for transcripts that start with the base ID
          matching_transcripts <- gene_data$txID[grepl(paste0("^", transcript_base), gene_data$txID)]
          
          if (length(matching_transcripts) > 0) {
            cat("DEBUG: Found similar transcripts:", paste(matching_transcripts, collapse=", "), "\n")
            # Use the first matching transcript
            tx_rows <- which(gene_data$txID == matching_transcripts[1])
            protein_seq <- gene_data$seq[tx_rows[1]]
            
            if (!is.null(protein_seq) && nchar(protein_seq) > 0) {
              cat("DEBUG: Found protein sequence via transcript ID matching - length:", nchar(protein_seq), "AA\n")
              return(as.character(protein_seq))
            }
          }
        }
        
        # Debug: show available transcripts in this gene
        available_transcripts <- unique(gene_data$txID)
        cat("DEBUG: Available transcripts in gene", gene_id, ":", paste(head(available_transcripts, 5), collapse=", "), "\n")
      } else {
        cat("DEBUG: Gene file exists but doesn't have expected structure\n")
      }
    } else {
      cat("DEBUG: Gene file not found:", gene_file, "\n")
    }
    
    # Method 2: Try currently loaded gene data as fallback
    if (exists("gene_data") && !is.null(gene_data()) && !is.null(gene_data()$peptides)) {
      gene_peptides <- gene_data()$peptides
      
      tx_rows <- which(gene_peptides$txID == transcript_id)
      if (length(tx_rows) > 0) {
        protein_seq <- gene_peptides$seq[tx_rows[1]]
        if (!is.null(protein_seq) && nchar(protein_seq) > 0) {
          cat("DEBUG: Found protein sequence in currently loaded gene data\n")
          return(as.character(protein_seq))
        }
      }
    }
    
    cat("DEBUG: Could not find protein sequence for transcript", transcript_id, "in gene", gene_id, "\n")
    return(NULL)
    
  }, error = function(e) {
    cat("ERROR in get_protein_sequence_for_transcript:", e$message, "\n")
    return(NULL)
  })
}

#===============================================================================
# SAFE ADAPTATIONS OF NOVEL_PEPTIDE_GENERATOR FUNCTIONS
#===============================================================================

#' Safe version of retrieve_sorted_cds
#' 
#' Adapts novel_peptide_generator's retrieve_sorted_cds with error handling
retrieve_sorted_cds_safe <- function(transcript_id, cds_data) {
  tryCatch({
    # Validate strand consistency
    strands <- unique(as.character(strand(cds_data)))
    if (length(strands) != 1) {
      cat("WARNING: Inconsistent strand information for", transcript_id, "\n")
      return(NULL)
    }
    strand_tx <- strands[1]
    
    # Sort CDS based on strand
    if (strand_tx == "+") {
      cds_tx <- cds_data[order(start(cds_data))]
    } else {
      cds_tx <- cds_data[order(start(cds_data), decreasing = TRUE)]
    }
    
    return(list(cds_tx = cds_tx, strand_tx = strand_tx))
    
  }, error = function(e) {
    cat("ERROR in retrieve_sorted_cds_safe:", e$message, "\n")
    return(NULL)
  })
}

#' Safe version of calculate_cumulative_cds
#' 
#' Adapts novel_peptide_generator's calculate_cumulative_cds with error handling
calculate_cumulative_cds_safe <- function(cds_tx, phases, strand_tx) {
  tryCatch({
    cumulative_cds <- data.frame(
      exon_idx = integer(),
      cds_start_nt = integer(),
      cds_end_nt = integer(),
      genomic_start = integer(),
      genomic_end = integer(),
      stringsAsFactors = FALSE
    )
    
    cumulative_length <- 0
    for (i in seq_along(cds_tx)) {
      exon <- cds_tx[i]
      exon_phase <- phases[i]
      exon_length <- width(exon)
      
      if (i == 1) {
        cds_start_nt <- 1 + exon_phase
      } else {
        cds_start_nt <- cumulative_length + 1 + exon_phase
      }
      
      cds_end_nt <- cds_start_nt + (exon_length - exon_phase) - 1
      
      if (strand_tx == "+") {
        genomic_start <- start(exon) + exon_phase
        genomic_end <- end(exon)
      } else {
        genomic_start <- end(exon) - exon_phase
        genomic_end <- start(exon)
      }
      
      cumulative_cds <- rbind(cumulative_cds, data.frame(
        exon_idx = i,
        cds_start_nt = cds_start_nt,
        cds_end_nt = cds_end_nt,
        genomic_start = genomic_start,
        genomic_end = genomic_end,
        stringsAsFactors = FALSE
      ))
      
      cumulative_length <- cds_end_nt
    }
    
    return(cumulative_cds)
    
  }, error = function(e) {
    cat("ERROR in calculate_cumulative_cds_safe:", e$message, "\n")
    return(NULL)
  })
}

#' Safe version of map_peptide_to_genomic
#' 
#' Adapts novel_peptide_generator's map_peptide_to_genomic with error handling
map_peptide_to_genomic_safe <- function(peptide, cds_cumulative, strand_tx, chrom) {
  tryCatch({
    peptide_name <- peptide$peptide
    aa_start <- peptide$aa_start
    aa_end <- peptide$aa_end
    
    cds_start_pos <- (aa_start - 1) * 3 + 1
    cds_end_pos <- aa_end * 3
    
    # Find exons spanned by peptide
    exons_spanned <- which(
      (cds_cumulative$cds_start_nt <= cds_end_pos) &
        (cds_cumulative$cds_end_nt >= cds_start_pos)
    )
    
    if (length(exons_spanned) == 0) {
      cat("WARNING: Peptide", peptide_name, "has nucleotide positions outside CDS range\n")
      return(NULL)
    }
    
    peptide_gr_list <- list()
    
    for (exon_idx in exons_spanned) {
      exon_info <- cds_cumulative[exon_idx, ]
      
      overlap_start <- max(cds_start_pos, exon_info$cds_start_nt)
      overlap_end <- min(cds_end_pos, exon_info$cds_end_nt)
      
      genomic_start <- map_cds_nt_to_genomic_safe(cds_cumulative, overlap_start, strand_tx)
      genomic_end <- map_cds_nt_to_genomic_safe(cds_cumulative, overlap_end, strand_tx)
      
      if (is.na(genomic_start) || is.na(genomic_end)) {
        cat("WARNING: Peptide", peptide_name, "has mapping issues in exon", exon_idx, "\n")
        next
      }
      
      ir_start <- min(genomic_start, genomic_end)
      ir_end <- max(genomic_start, genomic_end)
      
      # Create GRanges object for this peptide segment
      partial_gr <- GenomicRanges::GRanges(
        seqnames = chrom,
        ranges = IRanges::IRanges(start = ir_start, end = ir_end),
        strand = strand_tx,
        peptide = peptide_name
      )
      
      peptide_gr_list[[length(peptide_gr_list) + 1]] <- partial_gr
    }
    
    if (length(peptide_gr_list) == 0) {
      return(NULL)
    }
    
    peptide_gr <- do.call(c, peptide_gr_list)
    return(peptide_gr)
    
  }, error = function(e) {
    cat("ERROR in map_peptide_to_genomic_safe:", e$message, "\n")
    return(NULL)
  })
}

#' Safe version of map_cds_nt_to_genomic
#' 
#' Maps CDS nucleotide position to genomic position
map_cds_nt_to_genomic_safe <- function(cds_cumulative, nt_pos, strand_tx) {
  tryCatch({
    exon_row <- cds_cumulative[cds_cumulative$cds_start_nt <= nt_pos & cds_cumulative$cds_end_nt >= nt_pos, ]
    
    if (nrow(exon_row) == 0) {
      return(NA)
    }
    
    offset <- nt_pos - exon_row$cds_start_nt
    
    if (strand_tx == "+") {
      genomic_pos <- exon_row$genomic_start + offset
    } else {
      genomic_pos <- exon_row$genomic_start - offset
    }
    
    return(genomic_pos)
    
  }, error = function(e) {
    cat("ERROR in map_cds_nt_to_genomic_safe:", e$message, "\n")
    return(NA)
  })
}

#===============================================================================
# INTEGRATION HELPER FUNCTIONS
#===============================================================================

#' Create Transcript Structure Data for BLAST Visualization
#' 
#' Converts transcript structure to format compatible with existing visualization
create_blast_transcript_data <- function(transcript_id, gene_id, peptide_mapping) {
  
  tryCatch({
    if (is.null(peptide_mapping) || !peptide_mapping$success) {
      return(NULL)
    }
    
    # Create structure compatible with get_transcript_peptides_genomic output
    transcript_peptides_data <- list(
      genomic_ranges = peptide_mapping$genomic_ranges
    )
    
    return(transcript_peptides_data)
    
  }, error = function(e) {
    cat("ERROR in create_blast_transcript_data:", e$message, "\n")
    return(NULL)
  })
}

cat("✅ BLAST genomic mapper loaded successfully\n")