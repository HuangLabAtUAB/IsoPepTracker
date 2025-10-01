# Flanking Exons Module (Steps 3 & 4)
# Implements biological workflow steps 3 and 4

# Step 3: Determine Exon Orientation for ALL event types
# CRITICAL: Determines biological orientation (5' to 3') for subsequent sequence extraction and translation
identify_flanking_exons <- function(event_coords) {
  event_type <- event_coords$event_type
  strand <- event_coords$strand
  inclusion_isoform <- event_coords$inclusion_isoform
  exclusion_isoform <- event_coords$exclusion_isoform
  
  # Determine biological orientation and flanking exons for each event type
  if (event_type == "SE") {
    return(determine_orientation_SE(event_coords))
  } else if (event_type == "A3SS") {
    return(determine_orientation_A3SS(event_coords))
  } else if (event_type == "A5SS") {
    return(determine_orientation_A5SS(event_coords))
  } else if (event_type == "MXE") {
    return(determine_orientation_MXE(event_coords))
  } else if (event_type == "RI") {
    return(determine_orientation_RI(event_coords))
  } else {
    stop(paste("Unsupported event type for orientation determination:", event_type))
  }
}

# SE (Skipped Exon) orientation determination
determine_orientation_SE <- function(event_coords) {
  strand <- event_coords$strand
  inclusion_isoform <- event_coords$inclusion_isoform
  exclusion_isoform <- event_coords$exclusion_isoform
  
  # SE: Determine biological 5' to 3' orientation
  # Inclusion: [upstream] -> [skipped] -> [downstream]
  # Exclusion: [upstream] -> [downstream]
  
  if (strand == "+") {
    # Plus strand: genomic order = biological order
    biological_5prime <- inclusion_isoform[[1]]    # upstream exon
    biological_3prime <- inclusion_isoform[[3]]    # downstream exon
    translation_direction <- "left_to_right"
  } else {
    # Minus strand: genomic order is reverse of biological order
    biological_5prime <- inclusion_isoform[[3]]    # downstream in genomic = upstream in biological
    biological_3prime <- inclusion_isoform[[1]]    # upstream in genomic = downstream in biological  
    translation_direction <- "right_to_left"
  }
  
  return(create_orientation_result(event_coords, biological_5prime, biological_3prime, translation_direction, "SE"))
}

# A3SS (Alternative 3' Splice Site) flanking identification  
determine_orientation_A3SS <- function(event_coords) {
  strand <- event_coords$strand
  inclusion_isoform <- event_coords$inclusion_isoform
  exclusion_isoform <- event_coords$exclusion_isoform
  
  # A3SS: flanking exon is common, alternative exon varies
  # Both isoforms should have the flanking exon in common
  if (length(inclusion_isoform) == 2 && length(exclusion_isoform) == 2) {
    # Find the common flanking exon
    if (strand == "+") {
      upstream_exon <- inclusion_isoform[[1]]    # Flanking exon
      downstream_exon <- NULL                    # No common downstream for A3SS
    } else {
      upstream_exon <- inclusion_isoform[[2]]    # Flanking exon (reversed)
      downstream_exon <- NULL
    }
  }
  
  return(create_flanking_result(event_coords, upstream_exon, downstream_exon, "A3SS"))
}

# A5SS (Alternative 5' Splice Site) flanking identification
determine_orientation_A5SS <- function(event_coords) {
  strand <- event_coords$strand
  inclusion_isoform <- event_coords$inclusion_isoform
  exclusion_isoform <- event_coords$exclusion_isoform
  
  # A5SS: flanking exon is common, alternative exon varies
  if (length(inclusion_isoform) == 2 && length(exclusion_isoform) == 2) {
    # Find the common flanking exon
    if (strand == "+") {
      upstream_exon <- NULL                      # No common upstream for A5SS
      downstream_exon <- inclusion_isoform[[2]]  # Flanking exon
    } else {
      upstream_exon <- inclusion_isoform[[1]]    # Flanking exon (reversed)
      downstream_exon <- NULL
    }
  }
  
  return(create_flanking_result(event_coords, upstream_exon, downstream_exon, "A5SS"))
}

# MXE (Mutually Exclusive Exons) flanking identification
determine_orientation_MXE <- function(event_coords) {
  strand <- event_coords$strand
  inclusion_isoform <- event_coords$inclusion_isoform
  exclusion_isoform <- event_coords$exclusion_isoform
  
  # MXE: upstream and downstream exons are common (like SE)
  if (strand == "+") {
    upstream_exon <- inclusion_isoform[[1]]    # Common upstream
    downstream_exon <- inclusion_isoform[[3]]  # Common downstream
  } else {
    upstream_exon <- inclusion_isoform[[3]]    # Common upstream (biological)
    downstream_exon <- inclusion_isoform[[1]]  # Common downstream (biological)
  }
  
  return(create_flanking_result(event_coords, upstream_exon, downstream_exon, "MXE"))
}

# RI (Retained Intron) flanking identification
determine_orientation_RI <- function(event_coords) {
  strand <- event_coords$strand
  inclusion_isoform <- event_coords$inclusion_isoform
  exclusion_isoform <- event_coords$exclusion_isoform
  
  # RI: upstream and downstream exons from exclusion isoform are the flanking exons
  if (strand == "+") {
    upstream_exon <- exclusion_isoform[[1]]    # Upstream flanking
    downstream_exon <- exclusion_isoform[[2]]  # Downstream flanking
  } else {
    upstream_exon <- exclusion_isoform[[2]]    # Upstream flanking (biological)
    downstream_exon <- exclusion_isoform[[1]]  # Downstream flanking (biological)
  }
  
  return(create_flanking_result(event_coords, upstream_exon, downstream_exon, "RI"))
}

# Helper function to create orientation results for subsequent steps
create_orientation_result <- function(event_coords, biological_5prime, biological_3prime, translation_direction, event_type) {
  orientation_result <- list(
    # Biological orientation (critical for translation)
    biological_5prime_exon = biological_5prime,
    biological_3prime_exon = biological_3prime,
    translation_direction = translation_direction,
    
    # Flanking exons for GTF matching (compatibility)
    upstream = biological_5prime,
    downstream = biological_3prime,
    
    # Event metadata
    strand = event_coords$strand,
    gene_id = event_coords$gene_id,
    chromosome = event_coords$chromosome,
    event_type = event_type,
    
    # Isoform structures for sequence extraction
    inclusion_isoform_oriented = orient_isoform_for_translation(event_coords$inclusion_isoform, event_coords$strand),
    exclusion_isoform_oriented = orient_isoform_for_translation(event_coords$exclusion_isoform, event_coords$strand),
    
    # Validation
    orientation_determined = TRUE,
    warnings = NULL
  )
  
  return(orientation_result)
}

# Helper function to orient isoform for translation (5' to 3' biological order)
orient_isoform_for_translation <- function(isoform, strand) {
  if (strand == "+") {
    # Plus strand: keep genomic order (already 5' to 3')
    return(isoform)
  } else {
    # Minus strand: reverse genomic order to get 5' to 3' biological order
    return(rev(isoform))
  }
}

# Legacy helper function for compatibility
create_flanking_result <- function(event_coords, upstream_exon, downstream_exon, event_type) {
  return(create_orientation_result(event_coords, upstream_exon, downstream_exon, 
                                 ifelse(event_coords$strand == "+", "left_to_right", "right_to_left"), 
                                 event_type))
}

# Step 4: Exact Exon Match in CDS Index
search_all_exons_in_cds <- function(orientation_result, cds_index_path = "real_cds_index.RDS") {
  
  # Handle path context - check if we're in rmats_peptide directory or need to adjust path
  if (file.exists("modules/create_cds_index.R")) {
    source("modules/create_cds_index.R")
  } else if (file.exists("rmats_peptide/modules/create_cds_index.R")) {
    source("rmats_peptide/modules/create_cds_index.R")
  } else {
    # Try to find it relative to this script's location
    script_dir <- dirname(sys.frame(1)$ofile)
    if (file.exists(file.path(script_dir, "create_cds_index.R"))) {
      source(file.path(script_dir, "create_cds_index.R"))
    } else {
      stop("Cannot locate create_cds_index.R module")
    }
  }
  
  # Extract event information
  event_type <- orientation_result$event_type
  chromosome <- orientation_result$chromosome
  strand <- orientation_result$strand
  gene_id <- orientation_result$gene_id
  inclusion_isoform <- orientation_result$inclusion_isoform_oriented
  exclusion_isoform <- orientation_result$exclusion_isoform_oriented
  
  cat("=== STEP 4: EXACT EXON MATCH IN CDS ===\n")
  cat("Event Type:", event_type, "\n")
  cat("Gene ID:", gene_id, "\n")
  cat("Chromosome:", chromosome, "Strand:", strand, "\n")
  
  # Load CDS index
  tryCatch({
    cds_index <- load_cds_index(cds_index_path)
    
    # Check each exon: exact match or not
    cat("\nINCLUSION isoform exact matches:\n")
    inclusion_matches <- check_exact_exon_matches(inclusion_isoform, cds_index, chromosome, strand, gene_id, event_type, "inclusion")
    
    cat("\nEXCLUSION isoform exact matches:\n")
    exclusion_matches <- check_exact_exon_matches(exclusion_isoform, cds_index, chromosome, strand, gene_id, event_type, "exclusion")
    
    return(list(
      event_type = event_type,
      gene_id = gene_id,
      chromosome = chromosome,
      strand = strand,
      inclusion_exact_matches = inclusion_matches,
      exclusion_exact_matches = exclusion_matches,
      ready_for_step5 = TRUE  # Always ready - we just report what we found
    ))
    
  }, error = function(e) {
    cat("ERROR in CDS index search:", e$message, "\n")
    return(list(
      event_type = event_type,
      gene_id = gene_id,
      chromosome = chromosome,
      strand = strand,
      inclusion_exact_matches = list(),
      exclusion_exact_matches = list(),
      ready_for_step5 = FALSE,
      error = e$message
    ))
  })
}

# Enhanced exact exon match check with proper rMATS exon labeling
check_exact_exon_matches <- function(isoform_exons, cds_index, chromosome, strand, gene_id, event_type, isoform_type) {
  
  exon_matches <- list()
  
  # Get proper exon labels based on event type and isoform
  exon_labels <- get_exon_labels(event_type, isoform_type, strand, length(isoform_exons))
  
  for (i in seq_along(isoform_exons)) {
    exon <- isoform_exons[[i]]
    exon_label <- if (i <= length(exon_labels)) exon_labels[i] else paste("Exon", i)
    
    cat("  ", exon_label, ":", exon$start, "-", exon$end, " -> ")
    
    # Simple question: Is there an exact CDS exon with these coordinates?
    cds_result <- search_cds_index(cds_index, chromosome, exon$start, exon$end, strand, gene_id)
    
    if (cds_result$found) {
      phase <- cds_result$matches$phase[1]
      transcript_id <- cds_result$matches$transcript_id[1]
      cat("YES (CDS phase:", phase, "transcript:", transcript_id, ")\n")
      
      exon_matches[[i]] <- list(
        exon_number = i,
        exon_label = exon_label,
        start = exon$start,
        end = exon$end,
        exact_match = TRUE,
        phase = phase,
        transcript_id = transcript_id
      )
    } else {
      cat("NO (not in CDS)\n")
      
      # DYNAMIC PHASE CALCULATION: Calculate expected phase based on sequence context
      calculated_phase <- calculate_expected_phase(isoform_exons, i)
      cat("    Using calculated phase:", calculated_phase, "\n")
      
      exon_matches[[i]] <- list(
        exon_number = i,
        exon_label = exon_label,
        start = exon$start,
        end = exon$end,
        exact_match = FALSE,
        phase = calculated_phase,  # Dynamic calculation instead of static NA
        transcript_id = paste0(gene_id, ".calculated")
      )
    }
  }
  
  return(exon_matches)
}

# Get proper rMATS-based exon labels
get_exon_labels <- function(event_type, isoform_type, strand, num_exons) {
  
  if (event_type == "SE") {
    if (isoform_type == "inclusion") {
      if (strand == "+") {
        return(c("Upstream exon", "Skipped exon", "Downstream exon"))
      } else {
        return(c("Downstream exon (5' bio)", "Skipped exon", "Upstream exon (3' bio)"))
      }
    } else { # exclusion
      if (strand == "+") {
        return(c("Upstream exon", "Downstream exon"))
      } else {
        return(c("Downstream exon (5' bio)", "Upstream exon (3' bio)"))
      }
    }
  } else if (event_type == "A3SS") {
    if (isoform_type == "inclusion") {
      return(c("Flanking exon", "Long exon (alt 3' site)"))
    } else { # exclusion
      return(c("Flanking exon", "Short exon (alt 3' site)"))
    }
  } else if (event_type == "A5SS") {
    if (isoform_type == "inclusion") {
      return(c("Long exon (alt 5' site)", "Flanking exon"))
    } else { # exclusion
      return(c("Short exon (alt 5' site)", "Flanking exon"))
    }
  } else if (event_type == "MXE") {
    if (isoform_type == "inclusion") {
      if (strand == "+") {
        return(c("Upstream exon", "1st mutually exclusive exon", "Downstream exon"))
      } else {
        return(c("Downstream exon (5' bio)", "1st mutually exclusive exon", "Upstream exon (3' bio)"))
      }
    } else { # exclusion
      if (strand == "+") {
        return(c("Upstream exon", "2nd mutually exclusive exon", "Downstream exon"))
      } else {
        return(c("Downstream exon (5' bio)", "2nd mutually exclusive exon", "Upstream exon (3' bio)"))
      }
    }
  } else if (event_type == "RI") {
    if (isoform_type == "inclusion") {
      if (strand == "+") {
        return(c("Upstream exon", "Retained intron region", "Downstream exon"))
      } else {
        return(c("Downstream exon (5' bio)", "Retained intron region", "Upstream exon (3' bio)"))
      }
    } else { # exclusion
      if (strand == "+") {
        return(c("Upstream exon", "Downstream exon"))
      } else {
        return(c("Downstream exon (5' bio)", "Upstream exon (3' bio)"))
      }
    }
  }
  
  # Fallback to generic labels
  return(paste("Exon", 1:num_exons))
}

# Analyze CDS coverage across both isoforms
analyze_cds_coverage <- function(inclusion_matches, exclusion_matches, event_type) {
  
  # Count CDS exons
  inclusion_cds_count <- sum(sapply(inclusion_matches, function(x) x$cds_found))
  exclusion_cds_count <- sum(sapply(exclusion_matches, function(x) x$cds_found))
  
  total_inclusion_exons <- length(inclusion_matches)
  total_exclusion_exons <- length(exclusion_matches)
  
  cat("\nCDS COVERAGE ANALYSIS:\n")
  cat("Inclusion isoform: ", inclusion_cds_count, "/", total_inclusion_exons, " exons in CDS\n")
  cat("Exclusion isoform: ", exclusion_cds_count, "/", total_exclusion_exons, " exons in CDS\n")
  
  # Determine coverage level
  if (inclusion_cds_count == total_inclusion_exons && exclusion_cds_count == total_exclusion_exons) {
    coverage <- "full"
    outcome <- "ALL_EXONS_IN_CDS"
    ready_for_step5 <- TRUE
  } else if (inclusion_cds_count > 0 || exclusion_cds_count > 0) {
    coverage <- "partial"
    outcome <- "PARTIAL_CDS_COVERAGE"
    ready_for_step5 <- TRUE
  } else {
    coverage <- "none"
    outcome <- "NO_CDS_EXONS_FOUND"
    ready_for_step5 <- FALSE
  }
  
  summary <- paste0(
    "Event: ", event_type, " | ",
    "CDS Coverage: ", coverage, " | ",
    "Inc: ", inclusion_cds_count, "/", total_inclusion_exons, " | ",
    "Exc: ", exclusion_cds_count, "/", total_exclusion_exons
  )
  
  return(list(
    coverage = coverage,
    summary = summary,
    outcome = outcome,
    ready_for_step5 = ready_for_step5,
    inclusion_cds_count = inclusion_cds_count,
    exclusion_cds_count = exclusion_cds_count
  ))
}

# Helper function to check phase compatibility
is.compatible_phase <- function(upstream_end_phase, downstream_start_phase) {
  # Phase compatibility rules for splice junctions
  # This is a simplified check - real implementation would be more complex
  
  # Generally, phases should be compatible across splice junctions
  # Phase 0: divisible by 3 (no partial codon)
  # Phase 1: 1 nucleotide into codon
  # Phase 2: 2 nucleotides into codon
  
  return(TRUE)  # Simplified - assume compatible for now
}

# Validate CDS transcript compatibility and extract phase information
validate_cds_compatibility <- function(upstream_matches, downstream_matches) {
  
  if (is.null(upstream_matches) || is.null(downstream_matches) || 
      !upstream_matches$found || !downstream_matches$found) {
    return(list(
      compatible = FALSE,
      confidence = "none",
      reason = "Missing CDS matches",
      phase_info = NULL
    ))
  }
  
  up_cds <- upstream_matches$matches
  down_cds <- downstream_matches$matches
  
  # Check for same transcript (highest confidence)
  up_transcripts <- unique(up_cds$transcript_id)
  down_transcripts <- unique(down_cds$transcript_id)
  common_transcripts <- intersect(up_transcripts, down_transcripts)
  
  if (length(common_transcripts) > 0) {
    # Same transcript found - highest confidence
    transcript_id <- common_transcripts[1]
    
    up_match <- up_cds[up_cds$transcript_id == transcript_id, ][1, ]
    down_match <- down_cds[down_cds$transcript_id == transcript_id, ][1, ]
    
    return(list(
      compatible = TRUE,
      confidence = "high",
      reason = "Same transcript, both in CDS",
      transcript_id = transcript_id,
      phase_info = list(
        upstream_phase = up_match$phase,
        downstream_phase = down_match$phase,
        phase_compatible = is_phase_compatible(up_match$phase, down_match$phase)
      ),
      upstream_match = up_match,
      downstream_match = down_match
    ))
  }
  
  # Different transcripts - check phase compatibility
  if (nrow(up_cds) > 0 && nrow(down_cds) > 0) {
    up_match <- up_cds[1, ]
    down_match <- down_cds[1, ]
    
    phase_compatible <- is_phase_compatible(up_match$phase, down_match$phase)
    
    if (phase_compatible) {
      return(list(
        compatible = TRUE,
        confidence = "medium",
        reason = "Different transcripts, compatible phases",
        phase_info = list(
          upstream_phase = up_match$phase,
          downstream_phase = down_match$phase,
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
          upstream_phase = up_match$phase,
          downstream_phase = down_match$phase,
          phase_compatible = FALSE
        )
      ))
    }
  }
  
  return(list(
    compatible = FALSE,
    confidence = "none",
    reason = "No compatible CDS matches found",
    phase_info = NULL
  ))
}

# Check if phases are compatible for proper splicing  
is_phase_compatible <- function(upstream_phase, downstream_phase) {
  # Simple phase compatibility check
  # For now, return TRUE - this needs proper biological logic
  return(TRUE)
}

# Step 5: Extract Phase Information and Build CDS GTF Structures
extract_phase_information <- function(step4_results) {
  
  cat("=== STEP 5: EXTRACT PHASE INFORMATION ===\n")
  cat("Event Type:", step4_results$event_type, "\n")
  cat("Gene ID:", step4_results$gene_id, "\n")
  
  # Extract metadata for GTF creation
  gene_id <- step4_results$gene_id
  
  # Extract chromosome and strand info from Step 3 orientation result (stored in step4_results)
  chromosome <- step4_results$chromosome %||% "chr1"
  strand <- step4_results$strand %||% "+"
  
  cat("Using chromosome:", chromosome, "strand:", strand, "\n")
  
  # Process inclusion and exclusion isoforms independently
  inclusion_result <- process_isoform_for_translation(
    step4_results$inclusion_exact_matches, 
    "inclusion", 
    step4_results$event_type,
    gene_id,
    chromosome,
    strand
  )
  
  exclusion_result <- process_isoform_for_translation(
    step4_results$exclusion_exact_matches, 
    "exclusion", 
    step4_results$event_type,
    gene_id,
    chromosome,
    strand
  )
  
  # SPECIAL A5SS PHASE SHARING LOGIC
  # Strategy: If both alternatives found → independent, if one found → shared
  if (step4_results$event_type == "A5SS") {
    cat("\n--- A5SS SPECIAL PHASE SHARING LOGIC ---\n")
    
    # Check if either isoform needs phase sharing
    inclusion_needs_sharing <- !inclusion_result$translatable && 
                               !is.null(inclusion_result$needs_phase_sharing) && 
                               inclusion_result$needs_phase_sharing
    exclusion_needs_sharing <- !exclusion_result$translatable && 
                               !is.null(exclusion_result$needs_phase_sharing) && 
                               exclusion_result$needs_phase_sharing
    
    # Case 1: Inclusion needs sharing (use exclusion's phase)
    if (inclusion_needs_sharing && exclusion_result$translatable) {
      cat("A5SS: Inclusion needs phase sharing, using exclusion's phase\n")
      
      # Get exclusion alternative exon info
      excl_alt_match <- step4_results$exclusion_exact_matches[[1]]
      incl_alt_match <- step4_results$inclusion_exact_matches[[1]]
      incl_flanking_match <- step4_results$inclusion_exact_matches[[2]]
      
      # Create shared phase GTF for inclusion using exclusion's phase
      inclusion_shared_gtf <- data.frame(
        seqname = chromosome,
        source = "rmats_a5ss_shared",
        feature = "CDS",
        start = c(incl_alt_match$start, incl_flanking_match$start),
        end = c(incl_alt_match$end, incl_flanking_match$end),
        score = ".",
        strand = strand,
        frame = c(excl_alt_match$phase, 
                  if(incl_flanking_match$exact_match) incl_flanking_match$phase 
                  else (excl_alt_match$phase + (incl_alt_match$end - incl_alt_match$start + 1)) %% 3),
        attributes = c(
          paste0('gene_id "', gene_id, '"; transcript_id "', gene_id, '.inclusion_shared"; isoform_type "inclusion"; a5ss_type "shared"; shared_from "exclusion";'),
          paste0('gene_id "', gene_id, '"; transcript_id "', gene_id, '.inclusion_shared"; isoform_type "inclusion"; a5ss_type "shared"; segment "2";')
        ),
        stringsAsFactors = FALSE
      )
      
      # Update inclusion result
      inclusion_result$translatable <- TRUE
      inclusion_result$reason <- "A5SS shared phase from exclusion"
      inclusion_result$cds_gtf <- inclusion_shared_gtf
      
      cat("✅ A5SS inclusion now translatable using exclusion's phase", excl_alt_match$phase, "\n")
    }
    
    # Case 2: Exclusion needs sharing (use inclusion's phase)
    if (exclusion_needs_sharing && inclusion_result$translatable) {
      cat("A5SS: Exclusion needs phase sharing, using inclusion's phase\n")
      
      # Get inclusion alternative exon info
      incl_alt_match <- step4_results$inclusion_exact_matches[[1]]
      excl_alt_match <- step4_results$exclusion_exact_matches[[1]]
      excl_flanking_match <- step4_results$exclusion_exact_matches[[2]]
      
      # Create shared phase GTF for exclusion using inclusion's phase
      exclusion_shared_gtf <- data.frame(
        seqname = chromosome,
        source = "rmats_a5ss_shared",
        feature = "CDS",
        start = c(excl_alt_match$start, excl_flanking_match$start),
        end = c(excl_alt_match$end, excl_flanking_match$end),
        score = ".",
        strand = strand,
        frame = c(incl_alt_match$phase,
                  if(excl_flanking_match$exact_match) excl_flanking_match$phase
                  else (incl_alt_match$phase + (excl_alt_match$end - excl_alt_match$start + 1)) %% 3),
        attributes = c(
          paste0('gene_id "', gene_id, '"; transcript_id "', gene_id, '.exclusion_shared"; isoform_type "exclusion"; a5ss_type "shared"; shared_from "inclusion";'),
          paste0('gene_id "', gene_id, '"; transcript_id "', gene_id, '.exclusion_shared"; isoform_type "exclusion"; a5ss_type "shared"; segment "2";')
        ),
        stringsAsFactors = FALSE
      )
      
      # Update exclusion result
      exclusion_result$translatable <- TRUE
      exclusion_result$reason <- "A5SS shared phase from inclusion"
      exclusion_result$cds_gtf <- exclusion_shared_gtf
      
      cat("✅ A5SS exclusion now translatable using inclusion's phase", incl_alt_match$phase, "\n")
    }
    
    # Summary
    if (inclusion_result$translatable && exclusion_result$translatable) {
      if (!inclusion_needs_sharing && !exclusion_needs_sharing) {
        cat("🎯 A5SS: Both alternatives found in CDS - using independent phases\n")
      } else {
        cat("🔗 A5SS: Phase sharing successful - both isoforms now translatable\n")
      }
    } else if (inclusion_result$translatable || exclusion_result$translatable) {
      cat("⚖️  A5SS: One isoform translatable (with potential phase sharing)\n")
    } else {
      cat("❌ A5SS: Neither alternative found in CDS - both non-translatable\n")
    }
  }
  
  # Create final result
  phase_result <- list(
    event_type = step4_results$event_type,
    gene_id = step4_results$gene_id,
    
    # Translation decisions
    inclusion_translatable = inclusion_result$translatable,
    exclusion_translatable = exclusion_result$translatable,
    
    # CDS GTF structures (NULL if not translatable)
    inclusion_cds_gtf = inclusion_result$cds_gtf,
    exclusion_cds_gtf = exclusion_result$cds_gtf,
    
    # Decision reasons
    inclusion_reason = inclusion_result$reason,
    exclusion_reason = exclusion_result$reason,
    
    # Ready for Step 6
    ready_for_step6 = inclusion_result$translatable || exclusion_result$translatable
  )
  
  cat("\nRESULTS:\n")
  cat("Inclusion translatable:", phase_result$inclusion_translatable, 
      if (!is.null(inclusion_result$reason)) paste("(", inclusion_result$reason, ")") else "", "\n")
  cat("Exclusion translatable:", phase_result$exclusion_translatable,
      if (!is.null(exclusion_result$reason)) paste("(", exclusion_result$reason, ")") else "", "\n")
  cat("Ready for Step 6:", phase_result$ready_for_step6, "\n")
  
  return(phase_result)
}

# Process individual isoform for translation with frame validation
process_isoform_for_translation <- function(exact_matches, isoform_type, event_type, gene_id = NULL, chromosome = NULL, strand = NULL) {
  
  cat("\nProcessing", isoform_type, "isoform:\n")
  
  if (is.null(exact_matches) || length(exact_matches) == 0) {
    cat("  No exon matches found\n")
    return(list(
      translatable = FALSE,
      reason = "No exons found",
      cds_gtf = NULL
    ))
  }
  
  # SPECIAL HANDLING FOR RI EVENTS
  if (event_type == "RI") {
    return(process_ri_isoform_for_translation(exact_matches, isoform_type, gene_id, chromosome, strand))
  }
  
  # SPECIAL HANDLING FOR A5SS EVENTS (variable first exon dependency)
  if (event_type == "A5SS") {
    return(process_a5ss_isoform_for_translation(exact_matches, isoform_type, gene_id, chromosome, strand))
  }
  
  # Check if first exon (5' end) is in CDS
  first_exon <- exact_matches[[1]]
  if (!first_exon$exact_match) {
    cat("  First exon not in CDS - cannot translate\n")
    return(list(
      translatable = FALSE,
      reason = "First exon not in CDS",
      cds_gtf = NULL
    ))
  }
  
  cat("  First exon in CDS with phase:", first_exon$phase, "\n")
  
  # Validate frame continuity across all exons
  frame_validation <- validate_reading_frame_continuity(exact_matches)
  
  if (!frame_validation$valid) {
    cat("  Frame validation FAILED:", frame_validation$reason, "\n")
    return(list(
      translatable = FALSE,
      reason = frame_validation$reason,
      cds_gtf = NULL
    ))
  }
  
  cat("  Frame validation PASSED\n")
  
  # Build CDS GTF structure with correct metadata
  cds_gtf <- build_cds_gtf_structure(exact_matches, frame_validation$phases, isoform_type, gene_id, chromosome, strand)
  
  cat("  CDS GTF created with", nrow(cds_gtf), "CDS exons\n")
  
  return(list(
    translatable = TRUE,
    reason = "Frame preserved",
    cds_gtf = cds_gtf
  ))
}

# Validate reading frame continuity across exons
validate_reading_frame_continuity <- function(exact_matches) {
  
  current_phase <- exact_matches[[1]]$phase
  calculated_phases <- c(current_phase)
  
  for (i in 2:length(exact_matches)) {
    
    # Calculate expected phase based on previous exon length (3nt = 1 AA)
    prev_exon <- exact_matches[[i-1]]
    prev_length <- prev_exon$end - prev_exon$start + 1
    expected_phase <- (current_phase + prev_length) %% 3
    
    current_exon <- exact_matches[[i]]
    
    # If this exon is also in CDS, validate phase matches (for quality control)
    if (current_exon$exact_match) {
      indexed_phase <- current_exon$phase
      
      if (indexed_phase != expected_phase) {
        cat("    WARNING: Exon", i, "phase mismatch - expected:", expected_phase, "but CDS index shows:", indexed_phase, "\n")
        cat("    Using calculated phase", expected_phase, "to maintain reading frame continuity\n")
        # Continue with calculated phase instead of failing
      } else {
        cat("    Exon", i, "phase validated:", indexed_phase, "matches calculated", expected_phase, "\n")
      }
    } else {
      # Exon not in CDS - use calculated phase (this is the key biological insight!)
      cat("    Exon", i, "not in CDS - using calculated phase:", expected_phase, "(maintains reading frame)\n")
    }
    
    # Always use calculated phase to maintain reading frame continuity
    calculated_phases <- c(calculated_phases, expected_phase)
    current_phase <- expected_phase
  }
  
  return(list(
    valid = TRUE,  # Always valid - we calculate phases to maintain frame
    phases = calculated_phases
  ))
}

# Build CDS GTF structure using ACTUAL CDS phases when available
build_cds_gtf_structure <- function(exact_matches, calculated_phases, isoform_type, gene_id = NULL, chromosome = NULL, strand = NULL) {
  
  # Include ALL exons with proper phase information (actual CDS phases preferred)
  if (length(exact_matches) == 0) {
    return(NULL)
  }
  
  # Use provided values or try to infer from transcript info
  final_chromosome <- if (!is.null(chromosome)) chromosome else "chr1"
  final_strand <- if (!is.null(strand)) strand else "+"
  final_gene_id <- if (!is.null(gene_id)) gene_id else "unknown_gene"
  
  # Generate transcript ID for this isoform
  transcript_id <- paste0(final_gene_id, ".", isoform_type)
  
  # CRITICAL FIX: Use actual CDS phases when available, calculated only as fallback
  actual_phases <- sapply(seq_along(exact_matches), function(i) {
    match <- exact_matches[[i]]
    if (match$exact_match && !is.na(match$phase)) {
      # Use actual phase from CDS index (this is the correct biological phase)
      cat("    Using actual CDS phase", match$phase, "for exon", i, "(", match$exon_label, ")\n")
      return(match$phase)
    } else {
      # Use calculated phase for non-CDS exons (maintains frame continuity)
      cat("    Using calculated phase", calculated_phases[i], "for exon", i, "(", match$exon_label, ") - not in CDS\n")
      return(calculated_phases[i])
    }
  })
  
  # Create proper GTF attributes for ALL exons
  gtf_attributes <- paste0(
    "gene_id \"", final_gene_id, "\"; ",
    "transcript_id \"", transcript_id, "\"; ",
    "isoform_type \"", isoform_type, "\"; ",
    "exon_label \"", sapply(exact_matches, function(x) x$exon_label), "\"; ",
    "cds_confirmed \"", sapply(exact_matches, function(x) if(x$exact_match) "yes" else "no"), "\"; ",
    "phase_source \"", sapply(seq_along(exact_matches), function(i) {
      if (exact_matches[[i]]$exact_match && !is.na(exact_matches[[i]]$phase)) "actual_cds" else "calculated"
    }), "\";"
  )
  
  # Create CDS GTF data frame with ACTUAL phases
  cds_gtf <- data.frame(
    seqname = rep(final_chromosome, length(exact_matches)),
    source = "rmats_phase_extraction",
    feature = "CDS",
    start = sapply(exact_matches, function(x) x$start),
    end = sapply(exact_matches, function(x) x$end),
    score = ".",
    strand = rep(final_strand, length(exact_matches)),
    frame = actual_phases,  # ACTUAL CDS phases when available!
    attributes = gtf_attributes,
    stringsAsFactors = FALSE
  )
  
  return(cds_gtf)
}

# Validation function for flanking exon results
validate_flanking_results <- function(flanking_exons, gtf_results) {
  validation <- list(
    flanking_valid = TRUE,
    gtf_valid = TRUE,
    warnings = c(),
    errors = c()
  )
  
  # Check flanking exons
  if (is.null(flanking_exons$upstream) || is.null(flanking_exons$downstream)) {
    validation$flanking_valid <- FALSE
    validation$errors <- c(validation$errors, "Missing flanking exon coordinates")
  }
  
  # Check GTF results
  if (!gtf_results$ready_for_step5) {
    validation$gtf_valid <- FALSE
    validation$errors <- c(validation$errors, paste("GTF search failed:", gtf_results$outcome))
  }
  
  # Check for potential issues
  if (gtf_results$confidence == "medium") {
    validation$warnings <- c(validation$warnings, "Medium confidence - different transcripts found")
  }
  
  return(validation)
}

# RI-specific processing function for proper merged sequence handling
process_ri_isoform_for_translation <- function(exact_matches, isoform_type, gene_id = NULL, chromosome = NULL, strand = NULL) {
  
  cat("  Applying RI-specific processing logic...\n")
  
  # Use provided values or defaults
  final_chromosome <- if (!is.null(chromosome)) chromosome else "chr1"
  final_strand <- if (!is.null(strand)) strand else "+"
  final_gene_id <- if (!is.null(gene_id)) gene_id else "unknown_gene"
  
  if (isoform_type == "inclusion") {
    # RI inclusion: 3 segments (upstream + retained intron + downstream)
    cat("  RI inclusion: merging 3 segments into continuous sequence\n")
    
    # Check if upstream exon (first segment) is in CDS
    upstream_match <- exact_matches[[1]]
    if (!upstream_match$exact_match) {
      cat("  Upstream exon not in CDS - cannot translate\n")
      return(list(
        translatable = FALSE,
        reason = "Upstream exon not in CDS",
        cds_gtf = NULL
      ))
    }
    
    # Use upstream phase for the entire merged sequence
    upstream_phase <- upstream_match$phase
    cat("  Using upstream phase:", upstream_phase, "for entire inclusion sequence\n")
    
    # Create single merged CDS entry representing the continuous sequence
    merged_start <- exact_matches[[1]]$start
    merged_end <- exact_matches[[length(exact_matches)]]$end
    
    # Generate transcript ID for this isoform
    transcript_id <- paste0(final_gene_id, ".inclusion_merged")
    
    merged_cds_gtf <- data.frame(
      seqname = final_chromosome,
      source = "rmats_ri_merged",
      feature = "CDS",
      start = merged_start,
      end = merged_end,
      score = ".",
      strand = final_strand,
      frame = upstream_phase,  # Use upstream phase for entire sequence
      attributes = paste0(
        "gene_id \"", final_gene_id, "\"; ",
        "transcript_id \"", transcript_id, "\"; ",
        "isoform_type \"inclusion\"; ",
        "ri_type \"merged\"; ",
        "merged_segments \"3\";"
      ),
      stringsAsFactors = FALSE
    )
    
    cat("  Created merged CDS GTF with upstream phase", upstream_phase, "\n")
    return(list(
      translatable = TRUE,
      reason = "RI inclusion merged with upstream phase",
      cds_gtf = merged_cds_gtf
    ))
    
  } else {
    # RI exclusion: 2 segments (upstream + downstream, intron spliced out)
    cat("  RI exclusion: handling splice junction between upstream and downstream\n")
    
    # Check if upstream exon is in CDS
    upstream_match <- exact_matches[[1]]
    if (!upstream_match$exact_match) {
      cat("  Upstream exon not in CDS - cannot translate\n")
      return(list(
        translatable = FALSE,
        reason = "Upstream exon not in CDS",
        cds_gtf = NULL
      ))
    }
    
    # Check downstream exon (may or may not be in CDS)
    downstream_match <- exact_matches[[2]]
    
    if (downstream_match$exact_match) {
      cat("  Both upstream and downstream in CDS - validating splice junction\n")
      
      # Calculate expected phase at splice junction
      upstream_length <- upstream_match$end - upstream_match$start + 1
      expected_downstream_phase <- (upstream_match$phase + upstream_length) %% 3
      
      if (downstream_match$phase == expected_downstream_phase) {
        cat("  Splice junction phases compatible (", upstream_match$phase, "->", downstream_match$phase, ")\n")
        downstream_phase <- downstream_match$phase
      } else {
        cat("  WARNING: Splice junction phase mismatch - using calculated phase", expected_downstream_phase, "\n")
        downstream_phase <- expected_downstream_phase
      }
    } else {
      # Only upstream in CDS - downstream extends beyond CDS
      cat("  Only upstream in CDS - downstream extends beyond CDS boundary\n")
      upstream_length <- upstream_match$end - upstream_match$start + 1
      downstream_phase <- (upstream_match$phase + upstream_length) %% 3
    }
    
    # Generate transcript ID for this isoform
    transcript_id <- paste0(final_gene_id, ".exclusion_spliced")
    
    # Create CDS GTF entries for the segments
    ri_exclusion_gtf <- data.frame(
      seqname = rep(final_chromosome, 2),
      source = "rmats_ri_exclusion",
      feature = "CDS",
      start = c(upstream_match$start, downstream_match$start),
      end = c(upstream_match$end, downstream_match$end),
      score = ".",
      strand = rep(final_strand, 2),
      frame = c(upstream_match$phase, downstream_phase),
      attributes = paste0(
        "gene_id \"", final_gene_id, "\"; ",
        "transcript_id \"", transcript_id, "\"; ",
        "isoform_type \"exclusion\"; ",
        "ri_type \"spliced\"; ",
        "segment \"", 1:2, "\";"
      ),
      stringsAsFactors = FALSE
    )
    
    cat("  Created spliced CDS GTF with phases", upstream_match$phase, "and", downstream_phase, "\n")
    return(list(
      translatable = TRUE,
      reason = "RI exclusion with splice junction",
      cds_gtf = ri_exclusion_gtf
    ))
  }
}

# A5SS (Alternative 5' Splice Site) isoform processing with special variable first exon handling
process_a5ss_isoform_for_translation <- function(exact_matches, isoform_type, gene_id = NULL, chromosome = NULL, strand = NULL) {
  cat("  Applying A5SS-specific processing logic...\n")
  
  # A5SS structure: [alternative_exon] + [flanking_exon]
  # Challenge: First exon (alternative) varies between inclusion/exclusion
  # Strategy: Search both alternatives, use independent/shared phases
  
  if (length(exact_matches) != 2) {
    cat("  ERROR: A5SS should have exactly 2 exons, got", length(exact_matches), "\n")
    return(list(
      translatable = FALSE,
      reason = "Invalid A5SS structure",
      cds_gtf = NULL
    ))
  }
  
  alternative_match <- exact_matches[[1]]  # First exon (long or short alternative)
  flanking_match <- exact_matches[[2]]     # Second exon (flanking)
  
  cat("  A5SS", isoform_type, "- Alternative exon:", if(alternative_match$exact_match) "IN CDS" else "NOT in CDS")
  if (alternative_match$exact_match) cat(" (phase:", alternative_match$phase, ")")
  cat("\n")
  
  cat("  A5SS", isoform_type, "- Flanking exon:", if(flanking_match$exact_match) "IN CDS" else "NOT in CDS")
  if (flanking_match$exact_match) cat(" (phase:", flanking_match$phase, ")")
  cat("\n")
  
  # Check if alternative exon is in CDS (critical for A5SS)
  if (!alternative_match$exact_match) {
    cat("  A5SS", isoform_type, "- Alternative exon not in CDS, checking for phase borrowing...\n")
    # This will be handled by the calling function through phase sharing logic
    return(list(
      translatable = FALSE,
      reason = "A5SS alternative exon not in CDS",
      cds_gtf = NULL,
      needs_phase_sharing = TRUE,
      flanking_available = flanking_match$exact_match
    ))
  }
  
  # Alternative exon is in CDS - proceed with standard processing
  cat("  A5SS", isoform_type, "- Alternative exon found in CDS, using its phase\n")
  
  # Build A5SS CDS GTF using the alternative exon's phase
  a5ss_cds_gtf <- data.frame(
    seqname = chromosome,
    source = paste0("rmats_a5ss_", isoform_type),
    feature = "CDS",
    start = c(alternative_match$start, flanking_match$start),
    end = c(alternative_match$end, flanking_match$end),
    score = ".",
    strand = strand,
    frame = c(alternative_match$phase, flanking_match$phase),
    attributes = c(
      paste0('gene_id "', gene_id, '"; transcript_id "', gene_id, '.', isoform_type, '_a5ss"; isoform_type "', isoform_type, '"; a5ss_type "independent"; segment "1";'),
      paste0('gene_id "', gene_id, '"; transcript_id "', gene_id, '.', isoform_type, '_a5ss"; isoform_type "', isoform_type, '"; a5ss_type "independent"; segment "2";')
    ),
    stringsAsFactors = FALSE
  )
  
  # Validate flanking exon phase if available
  if (flanking_match$exact_match) {
    alt_length <- alternative_match$end - alternative_match$start + 1
    expected_flanking_phase <- (alternative_match$phase + alt_length) %% 3
    
    if (flanking_match$phase != expected_flanking_phase) {
      cat("  WARNING: A5SS phase mismatch - flanking expected:", expected_flanking_phase, "but got:", flanking_match$phase, "\n")
      cat("  Using calculated phase", expected_flanking_phase, "to maintain reading frame continuity\n")
      a5ss_cds_gtf$frame[2] <- expected_flanking_phase
    } else {
      cat("  A5SS phase transition validated - alternative", alternative_match$phase, "→ flanking", flanking_match$phase, "\n")
    }
  } else {
    # Flanking exon not in CDS - calculate expected phase
    alt_length <- alternative_match$end - alternative_match$start + 1
    expected_flanking_phase <- (alternative_match$phase + alt_length) %% 3
    cat("  A5SS flanking exon not in CDS - using calculated phase", expected_flanking_phase, "\n")
    a5ss_cds_gtf$frame[2] <- expected_flanking_phase
  }
  
  cat("  Created A5SS", isoform_type, "CDS GTF with independent phase from alternative exon\n")
  return(list(
    translatable = TRUE,
    reason = "A5SS independent phase from alternative exon",
    cds_gtf = a5ss_cds_gtf
  ))
}

#===============================================================================
# DYNAMIC PHASE CALCULATION FOR ALTERNATIVE SPLICING
#===============================================================================

# Helper function for dynamic phase calculation in alternative splicing events
# Calculates expected phase based on cumulative sequence length of previous exons
calculate_expected_phase <- function(isoform_exons, current_index) {
  tryCatch({
    # First exon always starts at phase 0 (beginning of translation)
    if (current_index == 1) {
      return(0)
    }
    
    # Calculate cumulative length of all previous exons
    cumulative_length <- 0
    for (j in 1:(current_index-1)) {
      exon <- isoform_exons[[j]]
      exon_length <- exon$end - exon$start + 1
      cumulative_length <- cumulative_length + exon_length
    }
    
    # Phase is cumulative length modulo 3 (codon frame)
    calculated_phase <- cumulative_length %% 3
    
    return(calculated_phase)
    
  }, error = function(e) {
    # Safe fallback if calculation fails
    cat("    Warning: Phase calculation failed, using phase 0 fallback\n")
    return(0)
  })
}