# GTF Indexer Module for Step 4
# Creates efficient index for exact coordinate matching and phase extraction

library(data.table)
library(dplyr)

# Load pre-built GTF index (no longer builds from scratch)
load_gtf_index <- function(index_path) {
  
  if (!file.exists(index_path)) {
    stop("GTF index file not found: ", index_path, 
         "\nPlease build the index first using gtf_index_builder.R")
  }
  
  cat("Loading pre-built GTF index from:", index_path, "\n")
  gtf_index <- readRDS(index_path)
  
  cat("GTF index loaded successfully:\n")
  cat("- Total exons:", gtf_index$stats$total_exons, "\n")
  cat("- Genes indexed:", gtf_index$stats$genes_indexed, "\n")
  cat("- Source GTF:", basename(gtf_index$source_file), "\n")
  
  return(gtf_index)
}

# Create mock GTF index for testing (replace with real GTF parsing)
create_mock_gtf_index <- function(gene_id) {
  # Mock GTF entries for testing
  mock_exons <- data.frame(
    seqname = "chr1",
    start = c(1000, 1400, 2000, 2500),
    end = c(1100, 1500, 2200, 2700),
    strand = c("-", "-", "-", "-"),
    feature = "exon",
    gene_id = gene_id,
    transcript_id = c("ENST00000123456", "ENST00000123456", "ENST00000123457", "ENST00000123457"),
    exon_id = c("ENSE00001", "ENSE00002", "ENSE00003", "ENSE00004"),
    phase_start = c(0, 1, 0, 2),
    phase_end = c(1, 0, 2, 0),
    in_cds = c(TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  
  gtf_index <- list(
    by_gene = list(),
    by_coords = mock_exons,
    gtf_path = "mock_gtf",
    indexed_at = Sys.time(),
    total_exons = nrow(mock_exons)
  )
  
  # Index by gene
  gtf_index$by_gene[[gene_id]] <- mock_exons
  
  return(gtf_index)
}

# Search GTF index for exact coordinate matches
search_gtf_exact <- function(gtf_index, chromosome, start, end, strand, gene_id = NULL) {
  
  # Filter by gene if specified
  if (!is.null(gene_id) && gene_id %in% names(gtf_index$by_gene)) {
    search_data <- gtf_index$by_gene[[gene_id]]
  } else {
    search_data <- gtf_index$by_coords
  }
  
  # Exact coordinate matching
  matches <- search_data[
    search_data$seqname == chromosome &
    search_data$start == start &
    search_data$end == end &
    search_data$strand == strand &
    search_data$feature == "exon",
  ]
  
  if (nrow(matches) == 0) {
    return(list(
      found = FALSE,
      matches = data.frame(),
      reason = "No exact coordinate match found"
    ))
  }
  
  return(list(
    found = TRUE,
    matches = matches,
    count = nrow(matches)
  ))
}

# Validate transcript compatibility and extract phase information
validate_transcript_compatibility <- function(upstream_matches, downstream_matches) {
  
  if (!upstream_matches$found || !downstream_matches$found) {
    return(list(
      compatible = FALSE,
      confidence = "none",
      reason = "Missing matches",
      phase_info = NULL
    ))
  }
  
  up_transcripts <- unique(upstream_matches$matches$transcript_id)
  down_transcripts <- unique(downstream_matches$matches$transcript_id)
  
  # Check for same transcript (highest confidence)
  common_transcripts <- intersect(up_transcripts, down_transcripts)
  
  if (length(common_transcripts) > 0) {
    # Same transcript found - highest confidence
    transcript_id <- common_transcripts[1]
    
    up_match <- upstream_matches$matches[upstream_matches$matches$transcript_id == transcript_id, ][1, ]
    down_match <- downstream_matches$matches[downstream_matches$matches$transcript_id == transcript_id, ][1, ]
    
    # Check if both are in CDS
    if (up_match$in_cds && down_match$in_cds) {
      return(list(
        compatible = TRUE,
        confidence = "high",
        reason = "Same transcript, both in CDS",
        transcript_id = transcript_id,
        phase_info = list(
          upstream_end_phase = up_match$phase_end,
          downstream_start_phase = down_match$phase_start,
          phase_compatible = is_phase_compatible(up_match$phase_end, down_match$phase_start)
        ),
        upstream_match = up_match,
        downstream_match = down_match
      ))
    } else {
      return(list(
        compatible = FALSE,
        confidence = "none",
        reason = "Not both in CDS",
        phase_info = NULL
      ))
    }
  }
  
  # Different transcripts - check phase compatibility
  if (length(up_transcripts) > 0 && length(down_transcripts) > 0) {
    up_match <- upstream_matches$matches[1, ]
    down_match <- downstream_matches$matches[1, ]
    
    if (up_match$in_cds && down_match$in_cds) {
      phase_compatible <- is_phase_compatible(up_match$phase_end, down_match$phase_start)
      
      if (phase_compatible) {
        return(list(
          compatible = TRUE,
          confidence = "medium",
          reason = "Different transcripts, compatible phases",
          phase_info = list(
            upstream_end_phase = up_match$phase_end,
            downstream_start_phase = down_match$phase_start,
            phase_compatible = TRUE
          ),
          upstream_match = up_match,
          downstream_match = down_match
        ))
      } else {
        return(list(
          compatible = FALSE,
          confidence = "none",
          reason = "Incompatible phases",
          phase_info = list(
            upstream_end_phase = up_match$phase_end,
            downstream_start_phase = down_match$phase_start,
            phase_compatible = FALSE
          )
        ))
      }
    }
  }
  
  return(list(
    compatible = FALSE,
    confidence = "none",
    reason = "No compatible matches found",
    phase_info = NULL
  ))
}

# Check if phases are compatible across splice junction
is_phase_compatible <- function(upstream_end_phase, downstream_start_phase) {
  # Phase compatibility rules for splice junctions
  # This is a simplified implementation
  
  # Generally, the phases should match for proper splicing
  # Phase 0: complete codon boundary
  # Phase 1: 1 nucleotide into codon
  # Phase 2: 2 nucleotides into codon
  
  # For now, assume phases are compatible if they sum to 0 or 3 (modulo 3)
  phase_sum <- (upstream_end_phase + downstream_start_phase) %% 3
  return(phase_sum == 0)
}

# Extract gene-specific GTF subset for faster searching
extract_gene_gtf <- function(gtf_index, gene_id) {
  if (gene_id %in% names(gtf_index$by_gene)) {
    return(gtf_index$by_gene[[gene_id]])
  } else {
    # Search in main index
    gene_entries <- gtf_index$by_coords[gtf_index$by_coords$gene_id == gene_id, ]
    return(gene_entries)
  }
}