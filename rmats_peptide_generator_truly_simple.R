# rmats_peptide_generator_truly_simple.R
# Script to generate peptide data for rMATS sequences using EXACT original logic

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
# CONFIGURATION
#===============================================================================

# Default file paths (only used when running standalone)
default_gtf_file <- "novel_transcript_nt.transdecoder.genome.gtf"
default_fasta_file <- "novel_transcript_nt.fa.transdecoder.pep"

# All enzymes must be processed for complete output structure
enzymes <- c("trp", "chymo", "aspn", "lysc", "lysn", "gluc")

# Enzyme to cleaver mapping
enzyme_map <- list(
  trp = "trypsin",
  chymo = "chymotrypsin-high",
  aspn = "asp-n endopeptidase",
  lysc = "lysc",
  lysn = "lysn",
  gluc = "glutamyl endopeptidase"
)

#===============================================================================
# FUNCTIONS FOR PROTEIN SEQUENCE EXTRACTION
#===============================================================================

# Extract protein sequences for novel sequences - adapted for novel format
extract_novel_proteins <- function(fasta_file, gtf_file) {
  cat("Loading protein sequences from:", fasta_file, "\n")
  
  # Load all protein sequences
  proteinSeqs = read.fasta(fasta_file, as.string = TRUE, seqtype = 'AA')
  
  # Parse rMATS sequence headers (ONLY change from original)
  # Format: >"ENSG00000099721.15".inclusion GENE."AMELY"~~"ENSG00000099721.15".inclusion  ORF type:complete len:45 (+),score=100.00
  protein_data <- data.table(
    proteinID = sapply(names(proteinSeqs), function(x) {
      # Extract the first part before space and clean quotes
      gsub('"', '', strsplit(x, " ")[[1]][1])
    }),
    txID = sapply(names(proteinSeqs), function(x) {
      # For rMATS sequences, transcript ID is same as protein ID, clean quotes
      gsub('"', '', strsplit(x, " ")[[1]][1])
    }),
    geneID = sapply(names(proteinSeqs), function(x) {
      # Extract gene ID from the GENE. part and clean quotes
      gene_part <- str_extract(x, "GENE\\.[^~]+")
      if (!is.na(gene_part)) {
        clean_gene <- gsub("GENE\\.", "", gene_part)
        clean_gene <- gsub('"', '', clean_gene)
        paste0(clean_gene, ".1")
      } else {
        base_name <- strsplit(strsplit(x, " ")[[1]][1], "\\.")[[1]][1]
        clean_base <- gsub('"', '', base_name)
        paste0(clean_base, ".1")
      }
    }),
    geneSymbol = sapply(names(proteinSeqs), function(x) {
      # Extract base gene symbol from GENE part and clean quotes
      gene_part <- str_extract(x, "GENE\\.[^~]+")
      if (!is.na(gene_part)) {
        clean_gene <- gsub("GENE\\.", "", gene_part)
        gsub('"', '', clean_gene)
      } else {
        base_name <- strsplit(strsplit(x, " ")[[1]][1], "\\.")[[1]][1]
        gsub('"', '', base_name)
      }
    }),
    numAA = sapply(names(proteinSeqs), function(x) {
      # Extract length from "len:XX"
      len_match <- str_extract(x, "len:[0-9]+")
      if (!is.na(len_match)) {
        gsub("len:", "", len_match)
      } else {
        as.character(nchar(proteinSeqs[[x]]))
      }
    }),
    seq = unlist(proteinSeqs)
  )
  
  # Remove stop codons (*) from sequences
  protein_data[, seq := gsub("\\*$", "", seq)]
  
  # Load GTF to get available transcripts
  cat("Loading GTF to check available transcripts...\n")
  gtf_gr <- import(gtf_file, format = "gtf")
  available_transcripts <- unique(gtf_gr$transcript_id)
  
  cat("Available transcripts in GTF:", paste(available_transcripts, collapse = ", "), "\n")
  cat("Proteins in FASTA:", paste(protein_data$proteinID, collapse = ", "), "\n")
  
  # Clean quotes from GTF transcript IDs for comparison with cleaned FASTA IDs
  available_transcripts_clean <- gsub('"', '', available_transcripts)
  
  # Filter proteins to only those with matching transcripts in GTF (compare cleaned IDs)
  protein_data <- protein_data[txID %in% available_transcripts_clean]
  
  if (nrow(protein_data) == 0) {
    stop("No proteins found with matching transcripts in GTF file")
  }
  
  cat("Matched proteins:", nrow(protein_data), "\n")
  cat("Processing proteins:", paste(protein_data$proteinID, collapse = ", "), "\n")
  
  return(protein_data)
}

#===============================================================================
# MAIN GENERATOR FUNCTION
#===============================================================================

generate_novel_peptide_data <- function(gtf_file, fasta_file, missedCleavages = 0, minLength = 7, maxLength = 60) {
  cat("Starting peptide data generation for novel sequences\n")
  
  # Step 1: Extract protein sequences for novel sequences
  cat("Extracting protein sequences...\n")
  protein_seqs <- extract_novel_proteins(fasta_file, gtf_file)
  
  # Ensure no factors in protein_seqs
  for (col in names(protein_seqs)) {
    if (is.factor(protein_seqs[[col]])) {
      protein_seqs[[col]] <- as.character(protein_seqs[[col]])
    }
  }
  
  # Step 2: Process ALL enzymes to match expected output structure  
  # FIX: Process row by row to avoid data.table assignment error with sapply
  for (enzyme in enzymes) {
    cat("Processing enzyme:", enzyme, "\n")
    cleaver_enzyme <- enzyme_map[[enzyme]]
    
    # Create list to store results
    peptides_list <- list()
    
    for (i in 1:nrow(protein_seqs)) {
      seq_i <- protein_seqs$seq[i]
      
      # Digest protein
      peptides <- cleaver::cleave(
        seq_i, 
        enzym = cleaver_enzyme,
        missedCleavages = missedCleavages,
        custom = NULL, 
        unique = TRUE
      )
      
      # Filter by length and convert to simple character vector
      if (length(peptides) > 0) {
        # Extract the actual peptide vector from cleaver's list result
        peptide_vector <- peptides[[1]]  # cleaver returns a named list, get the vector
        aaNo <- sapply(peptide_vector, nchar)
        filtered_peptides <- peptide_vector[aaNo >= minLength & aaNo <= maxLength]
        # Ensure it's a simple character vector
        filtered_peptides <- as.character(filtered_peptides)
      } else {
        filtered_peptides <- character(0)
      }
      
      peptides_list[[i]] <- filtered_peptides
    }
    
    # Assign the simple list structure to the data.table column
    protein_seqs[, paste0(enzyme, "Peps") := peptides_list]
  }
  
  # Step 3: Calculate peptide positions for ALL enzymes
  cat("Calculating peptide positions...\n")
  
  # Function to get peptide positions
  get_peptide_positions <- function(txID_i, peptides) {
    protein_seq <- protein_seqs[txID == txID_i, seq]
    
    if (length(protein_seq) == 0 || is.null(peptides) || length(peptides) == 0 || all(is.na(peptides))) {
      return(NA)
    }
    
    protein_seq <- as.character(protein_seq[[1]])
    peptide_positions_list <- lapply(peptides, function(pep) {
      match_positions <- tryCatch({
        str_locate_all(protein_seq, fixed(pep))[[1]]
      }, error = function(e) {
        warning("Error locating peptide positions: ", e$message)
        return(matrix(nrow=0, ncol=2, dimnames=list(NULL, c("start", "end"))))
      })
      
      if (nrow(match_positions) == 0) {
        return(NULL)
      }
      
      data.frame(
        peptide = pep,
        aa_start = match_positions[, "start"],
        aa_end = match_positions[, "end"],
        stringsAsFactors = FALSE
      )
    })
    
    peptide_positions_df <- do.call(rbind, peptide_positions_list)
    
    if (!is.null(peptide_positions_df) && nrow(peptide_positions_df) > 0) {
      return(peptide_positions_df)
    } else {
      return(NA)
    }
  }
  
  # Apply get_peptide_positions to ALL enzyme peptide types
  for (enzyme in enzymes) {
    enzyme_pep <- paste0(enzyme, "Peps")
    peptide_col <- paste0(enzyme_pep, "_positions")
    
    protein_seqs[, (peptide_col) := lapply(.I, function(i) {
      txID_i <- txID[i]
      peptides_i <- get(enzyme_pep)[[i]]
      
      if (is.null(peptides_i) || length(peptides_i) == 0 || all(is.na(peptides_i))) {
        return(NA)
      }
      
      positions_df <- tryCatch({
        result <- get_peptide_positions(txID_i, peptides_i)
        # Ensure no factors in the result
        if (is.data.frame(result)) {
          for (col in names(result)) {
            if (is.factor(result[[col]])) {
              result[[col]] <- as.character(result[[col]])
            }
          }
        }
        result
      }, error = function(e) {
        warning("Error calculating peptide positions for ", txID_i, ": ", e$message)
        return(NA)
      })
      
      return(positions_df)
    })]
  }
  
  # Step 4: Import GTF and extract CDS information
  cat("Importing GTF and mapping to genomic coordinates...\n")
  gtf_gr <- import(gtf_file, format = "gtf")
  cds_gr <- gtf_gr[gtf_gr$type == "CDS"]
  phases <- mcols(cds_gr)$phase
  cds_by_tx <- split(cds_gr, cds_gr$transcript_id)
  
  # Create mapping from cleaned IDs to original quoted IDs for CDS lookup
  original_transcript_ids <- unique(gtf_gr$transcript_id)
  cleaned_transcript_ids <- gsub('"', '', original_transcript_ids)
  id_mapping <- setNames(original_transcript_ids, cleaned_transcript_ids)
  
  # Step 5: Functions for genomic mapping
  
  # Validate strand consistency
  validate_strand_consistency <- function(cds_tx) {
    strands <- unique(as.character(strand(cds_tx)))
    if (length(strands) != 1) {
      stop("Inconsistent strand information across CDS entries for the transcript.")
    }
    return(strands[1])
  }
  
  # Retrieve and sort CDS information
  retrieve_sorted_cds <- function(txID, cds_by_tx) {
    cds_tx <- cds_by_tx[[txID]]
    if (is.null(cds_tx)) {
      warning(paste("Transcript ID", txID, "not found in CDS entries."))
      return(NULL)
    }
    
    # Validate strand consistency
    strand_tx <- validate_strand_consistency(cds_tx)
    
    # Sort CDS based on strand
    if (strand_tx == "+") {
      cds_tx <- cds_tx[order(start(cds_tx))]
    } else {
      cds_tx <- cds_tx[order(start(cds_tx), decreasing = TRUE)]
    }
    
    return(list(cds_tx = cds_tx, strand_tx = strand_tx))
  }
  
  # Calculate cumulative CDS positions
  calculate_cumulative_cds <- function(cds_tx, phases, strand_tx) {
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
  }
  
  # Map CDS nucleotide to genomic position
  map_cds_nt_to_genomic <- function(cds_cumulative, nt_pos, strand_tx) {
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
  }
  
  # Map peptide to genomic coordinates
  map_peptide_to_genomic <- function(peptide, cds_cumulative, strand_tx, chrom) {
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
      warning(paste("Peptide", peptide_name, "has nucleotide positions outside CDS range."))
      return(NULL)
    }
    
    peptide_gr_list <- list()
    
    for (exon_idx in exons_spanned) {
      exon_info <- cds_cumulative[exon_idx, ]
      
      overlap_start <- max(cds_start_pos, exon_info$cds_start_nt)
      overlap_end <- min(cds_end_pos, exon_info$cds_end_nt)
      
      genomic_start <- map_cds_nt_to_genomic(cds_cumulative, overlap_start, strand_tx)
      genomic_end <- map_cds_nt_to_genomic(cds_cumulative, overlap_end, strand_tx)
      
      if (is.na(genomic_start) || is.na(genomic_end)) {
        warning(paste("Peptide", peptide_name, "has mapping issues in exon", exon_idx))
        next
      }
      
      ir_start <- min(genomic_start, genomic_end)
      ir_end <- max(genomic_start, genomic_end)
      
      # Create GRanges object for this peptide segment
      partial_gr <- GRanges(
        seqnames = chrom,
        ranges = IRanges(start = ir_start, end = ir_end),
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
  }
  
  # Step 6: Initialize mapped_ranges columns for ALL enzymes
  for (enzyme in enzymes) {
    mapped_ranges_column_name <- paste0(enzyme, "Peps_mapped_ranges")
    protein_seqs[[mapped_ranges_column_name]] <- vector("list", nrow(protein_seqs))
  }
  
  # Step 7: Map peptides to genomic coordinates for ALL enzymes
  for (entry_idx in 1:nrow(protein_seqs)) {
    entry <- protein_seqs[entry_idx, ]
    txID_i <- entry$txID
    
    # Retrieve sorted CDS and strand information
    # Map cleaned ID back to original quoted ID for CDS lookup
    original_txID <- id_mapping[[txID_i]]
    if (is.null(original_txID)) {
      warning(paste("No mapping found for transcript", txID_i))
      next
    }
    
    cds_info <- tryCatch(
      retrieve_sorted_cds(original_txID, cds_by_tx),
      error = function(e) {
        warning(paste("Skipping transcript", txID_i, "due to error:", e$message))
        return(NULL)
      }
    )
    
    if (is.null(cds_info)) {
      next
    }
    
    cds_tx <- cds_info$cds_tx
    strand_tx <- cds_info$strand_tx
    chrom <- as.character(seqnames(cds_tx))[1]
    
    # Extract phase information for sorted CDS
    phases_tx <- mcols(cds_tx)$phase
    
    # Calculate cumulative CDS positions considering phase
    cds_cumulative <- calculate_cumulative_cds(cds_tx, phases_tx, strand_tx)
    
    # Loop over ALL enzymes
    for (enzyme in enzymes) {
      # Construct the column name for the peptide positions
      positions_column_name <- paste0(enzyme, "Peps_positions")
      
      # Check if the peptide positions column exists and is not empty
      if (!is.null(entry[[positions_column_name]]) && length(entry[[positions_column_name]]) > 0) {
        # Retrieve peptide positions data frame
        peptides_df <- as.data.frame(entry[[positions_column_name]][[1]])
        
        if (is.data.frame(peptides_df) && nrow(peptides_df) > 0) {
          # Initialize list to store GRanges for this enzyme
          peptide_gr_list <- list()
          
          # Iterate through each peptide
          for (peptide_idx in 1:nrow(peptides_df)) {
            peptide <- as.list(peptides_df[peptide_idx, ])
            
            # Map peptide to genomic ranges
            peptide_gr <- tryCatch({
              map_peptide_to_genomic(peptide, cds_cumulative, strand_tx, chrom)
            }, error = function(e) {
              warning("Error mapping peptide to genomic coordinates: ", e$message)
              return(NULL)
            })
            
            if (!is.null(peptide_gr)) {
              # Tag each peptide GRanges with transcript ID for correct faceting
              mcols(peptide_gr)$txID <- txID_i
              peptide_gr_list[[length(peptide_gr_list) + 1]] <- peptide_gr
            }
          }
          
          # Aggregate mapped peptides for this transcript and enzyme
          if (length(peptide_gr_list) > 0) {
            all_peptides_genomic_gr <- do.call(c, peptide_gr_list)
            # Store in the protein_seqs data frame
            mapped_ranges_column_name <- paste0(enzyme, "Peps_mapped_ranges")
            protein_seqs[[mapped_ranges_column_name]][[entry_idx]] <- all_peptides_genomic_gr
          }
        }
      }
    }
  }
  
  # Step 8: Ensure exact column order matching ENSG*.rds structure (24 columns)
  expected_columns <- c(
    "proteinID", "txID", "geneID", "geneSymbol", "numAA", "seq",
    "trpPeps", "chymoPeps", "aspnPeps", "lyscPeps", "lysnPeps", "glucPeps",
    "trpPeps_positions", "chymoPeps_positions", "aspnPeps_positions", 
    "lyscPeps_positions", "lysnPeps_positions", "glucPeps_positions",
    "trpPeps_mapped_ranges", "aspnPeps_mapped_ranges", "chymoPeps_mapped_ranges",
    "glucPeps_mapped_ranges", "lysnPeps_mapped_ranges", "lyscPeps_mapped_ranges"
  )
  
  # Reorder columns to match expected structure
  protein_seqs <- protein_seqs[, ..expected_columns]
  
  # Convert to data.frame (not data.table) to match expected format
  return_data <- as.data.frame(protein_seqs)
  
  # Explicitly ensure all character columns stay as character
  for (col in names(return_data)) {
    if (is.character(return_data[[col]])) {
      return_data[[col]] <- as.character(return_data[[col]])
    }
  }
  
  cat("Peptide data generation completed\n")
  cat("Generated", nrow(return_data), "rows with", ncol(return_data), "columns\n")
  cat("Column names:", paste(names(return_data), collapse = ", "), "\n")
  
  return(return_data)
}

#===============================================================================
# MAIN EXECUTION
#===============================================================================

# Only run if executed directly (not when sourced by server)
if (!exists(".GlobalEnv") || is.null(.GlobalEnv$.rmats_server_mode)) {
  # Generate peptide data for novel sequences
  cat("=== NOVEL SEQUENCE PEPTIDE GENERATOR ===\n")
  cat("Starting peptide data generation for novel sequences\n")
  
  if (file.exists(default_gtf_file) && file.exists(default_fasta_file)) {
    novel_peptides <- generate_novel_peptide_data(default_gtf_file, default_fasta_file)
    
    # Save to RDS file
    output_file <- "novel_transcript_nt_peptides.rds"
    saveRDS(novel_peptides, output_file)
    cat("Saved results to:", output_file, "\n")
  } else {
    cat("Required files not found:\n")
    cat("  GTF file:", default_gtf_file, "exists:", file.exists(default_gtf_file), "\n")
    cat("  FASTA file:", default_fasta_file, "exists:", file.exists(default_fasta_file), "\n")
    cat("Skipping automatic execution. Files will be generated when rMATS analysis runs.\n")
  }
}

# Display structure summary (only when running standalone)
if (!exists(".GlobalEnv") || is.null(.GlobalEnv$.rmats_server_mode)) {
  cat("\n=== OUTPUT SUMMARY ===\n")
  cat("Structure of first row:\n")
  str(novel_peptides[1,])
  cat("\nDimensions:", nrow(novel_peptides), "rows x", ncol(novel_peptides), "columns\n")
  cat("Expected 24 columns - Actual:", ncol(novel_peptides), "\n")
} 
