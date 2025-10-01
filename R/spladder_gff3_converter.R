#===============================================================================
# SplAdder GFF3 to rMATS Format Converter
# Converts hierarchical SplAdder GFF3 (gene → mRNA → exon) to rMATS tabular format
#===============================================================================

library(data.table)

#===============================================================================
# CORE PARSING FUNCTIONS
#===============================================================================

parse_spladder_gff3 <- function(gff3_file) {
  
  cat("🔄 Parsing SplAdder GFF3 file:", basename(gff3_file), "\n")
  
  # Validate file exists
  if (!file.exists(gff3_file)) {
    stop("GFF3 file not found: ", gff3_file)
  }
  
  # Read GFF3 file
  gff3_lines <- readLines(gff3_file)
  cat("DEBUG: Read", length(gff3_lines), "lines from file\n")
  
  # Remove comments and empty lines
  gff3_lines <- gff3_lines[!startsWith(gff3_lines, "#") & nzchar(gff3_lines)]
  cat("DEBUG: After filtering,", length(gff3_lines), "data lines remain\n")
  
  if (length(gff3_lines) == 0) {
    stop("No data lines found in GFF3 file after filtering")
  }
  
  # Parse into data frame
  gff3_data <- data.frame(
    seqid = character(),
    source = character(),
    type = character(),
    start = integer(),
    end = integer(),
    score = character(),
    strand = character(),
    phase = character(),
    attributes = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(gff3_lines)) {
    line <- gff3_lines[i]
    
    if (nzchar(line)) {  # Skip empty lines
      parts <- strsplit(line, "\t")[[1]]
      
      if (length(parts) == 9) {
        tryCatch({
          gff3_data <- rbind(gff3_data, data.frame(
            seqid = parts[1],
            source = parts[2],
            type = parts[3],
            start = as.integer(parts[4]),
            end = as.integer(parts[5]),
            score = parts[6],
            strand = parts[7],
            phase = parts[8],
            attributes = parts[9],
            stringsAsFactors = FALSE
          ))
        }, error = function(e) {
          cat("ERROR parsing line", i, ":", line, "\n")
          cat("Error:", e$message, "\n")
          stop("Failed to parse GFF3 line ", i, ": ", e$message)
        })
      } else {
        cat("WARNING: Line", i, "has", length(parts), "parts instead of 9:", substr(line, 1, 100), "\n")
      }
    }
  }
  
  cat("Parsed", nrow(gff3_data), "GFF3 entries\n")
  
  # Convert to rMATS format by event type
  event_types <- unique(gff3_data$source)
  cat("Found event types:", paste(event_types, collapse = ", "), "\n")
  
  rmats_events <- list()
  
  for (event_type in event_types) {
    type_events <- convert_event_type_to_rmats(gff3_data, event_type)
    if (nrow(type_events) > 0) {
      # Map to rMATS event type names
      if (event_type == "exon_skip") {
        rmats_events[["SE"]] <- type_events
      } else if (event_type == "alt_3prime") {
        rmats_events[["A3SS"]] <- type_events
      } else if (event_type == "alt_5prime") {
        rmats_events[["A5SS"]] <- type_events
      } else if (event_type == "mutex_exons") {
        rmats_events[["MXE"]] <- type_events
      } else if (event_type == "intron_retention") {
        rmats_events[["RI"]] <- type_events
      }
    }
  }
  
  total_events <- sum(sapply(rmats_events, nrow))
  cat("✅ Converted", total_events, "events to rMATS format\n")
  
  return(rmats_events)
}

#===============================================================================
# EVENT TYPE CONVERSION FUNCTIONS
#===============================================================================

convert_event_type_to_rmats <- function(gff3_data, event_type) {
  
  cat("Converting", event_type, "events...\n")
  
  # Filter to this event type
  event_data <- gff3_data[gff3_data$source == event_type, ]
  
  if (event_type == "exon_skip") {
    return(convert_exon_skip_to_rmats_se(event_data))
  } else if (event_type == "alt_3prime") {
    return(convert_alt_3prime_to_rmats_a3ss(event_data))
  } else if (event_type == "alt_5prime") {
    return(convert_alt_5prime_to_rmats_a5ss(event_data))
  } else if (event_type == "mutex_exons") {
    return(convert_mutex_exons_to_rmats_mxe(event_data))
  } else if (event_type == "intron_retention") {
    return(convert_intron_retention_to_rmats_ri(event_data))
  } else {
    cat("Warning: Unknown event type", event_type, "\n")
    return(data.frame())
  }
}

#===============================================================================
# EXON SKIP (SE) CONVERSION
#===============================================================================

convert_exon_skip_to_rmats_se <- function(event_data) {
  
  # Get all gene-level events
  genes <- event_data[event_data$type == "gene", ]
  
  rmats_events <- data.frame()
  
  for (i in 1:nrow(genes)) {
    gene_entry <- genes[i, ]
    
    # Extract gene info from attributes
    gene_id <- extract_gene_id(gene_entry$attributes)
    event_id <- extract_event_id(gene_entry$attributes)
    
    # Get all mRNA isoforms for this gene
    gene_mrnas <- event_data[event_data$type == "mRNA" & 
                           grepl(paste0("Parent=", event_id), event_data$attributes), ]
    
    if (nrow(gene_mrnas) < 2) {
      next  # Need at least 2 isoforms to compare
    }
    
    # Get exons for each isoform
    isoform1_exons <- get_exons_for_isoform(event_data, gene_mrnas[1, ])
    isoform2_exons <- get_exons_for_isoform(event_data, gene_mrnas[2, ])
    
    # Find the skipped exon (present in one isoform but not the other)
    skipped_exon <- find_skipped_exon(isoform1_exons, isoform2_exons)
    
    if (is.null(skipped_exon)) {
      next
    }
    
    # Find upstream and downstream exons
    all_exons <- rbind(isoform1_exons, isoform2_exons)
    all_exons <- all_exons[!duplicated(paste(all_exons$start, all_exons$end)), ]
    all_exons <- all_exons[order(all_exons$start), ]
    
    # Find exons flanking the skipped exon
    flanking <- find_flanking_exons(all_exons, skipped_exon)
    
    if (is.null(flanking$upstream) || is.null(flanking$downstream)) {
      next
    }
    
    # Create rMATS SE format entry
    rmats_event <- data.frame(
      ID = as.character(extract_event_number(event_id)),
      GeneID = gene_id,
      geneSymbol = gene_id,  # Use gene_id as symbol for now
      chr = paste0("chr", gene_entry$seqid),
      strand = gene_entry$strand,
      exonStart_0base = skipped_exon$start - 1,  # Convert to 0-based
      exonEnd = skipped_exon$end,               # 1-based
      upstreamES = flanking$upstream$start - 1,  # Convert to 0-based
      upstreamEE = flanking$upstream$end,       # 1-based
      downstreamES = flanking$downstream$start - 1,  # Convert to 0-based
      downstreamEE = flanking$downstream$end,   # 1-based
      # Add rMATS count columns (using dummy values for compatibility)
      IJC_SAMPLE_1 = "10,10,10",  # Inclusion junction counts
      SJC_SAMPLE_1 = "10,10,10",  # Skipping junction counts  
      IJC_SAMPLE_2 = "10,10,10",  # Inclusion junction counts
      SJC_SAMPLE_2 = "10,10,10",  # Skipping junction counts
      IncFormLen = as.integer(skipped_exon$end - skipped_exon$start + 1),  # Inclusion form length
      SkipFormLen = as.integer(0),  # Skip form length (no extra exon)
      stringsAsFactors = FALSE
    )
    
    rmats_events <- rbind(rmats_events, rmats_event)
  }
  
  return(rmats_events)
}

#===============================================================================
# HELPER FUNCTIONS
#===============================================================================

extract_gene_id <- function(attributes) {
  # Extract GeneName from attributes string like: GeneName="ENSG00000173614"
  gene_match <- regexpr('GeneName="([^"]+)"', attributes, perl = TRUE)
  if (gene_match > 0) {
    start_pos <- attr(gene_match, "capture.start")[1]
    length_val <- attr(gene_match, "capture.length")[1]
    return(substr(attributes, start_pos, start_pos + length_val - 1))
  }
  return("Unknown")
}

extract_event_id <- function(attributes) {
  # Extract ID from attributes string like: ID=exon_skip.92
  id_match <- regexpr('ID=([^;]+)', attributes, perl = TRUE)
  if (id_match > 0) {
    start_pos <- attr(id_match, "capture.start")[1]
    length_val <- attr(id_match, "capture.length")[1]
    return(substr(attributes, start_pos, start_pos + length_val - 1))
  }
  return("Unknown")
}

extract_event_number <- function(event_id) {
  # Extract number from event_id like "exon_skip.92" -> "92"
  num_match <- regexpr('\\.([0-9]+)', event_id, perl = TRUE)
  if (num_match > 0) {
    start_pos <- attr(num_match, "capture.start")[1]
    length_val <- attr(num_match, "capture.length")[1]
    return(substr(event_id, start_pos, start_pos + length_val - 1))
  }
  return("0")
}

get_exons_for_isoform <- function(event_data, mrna_entry) {
  # Get all exons that belong to this mRNA isoform
  mrna_id <- extract_event_id(mrna_entry$attributes)
  exons <- event_data[event_data$type == "exon" & 
                     grepl(paste0("Parent=", mrna_id), event_data$attributes), ]
  return(exons)
}

find_skipped_exon <- function(isoform1_exons, isoform2_exons) {
  # Find exon present in one isoform but not the other
  
  # Check if isoform2 has an extra exon (skipped in isoform1)
  for (i in 1:nrow(isoform2_exons)) {
    exon2 <- isoform2_exons[i, ]
    found_in_iso1 <- any(isoform1_exons$start == exon2$start & isoform1_exons$end == exon2$end)
    if (!found_in_iso1) {
      return(exon2)  # This is the skipped exon
    }
  }
  
  # Check if isoform1 has an extra exon (skipped in isoform2)
  for (i in 1:nrow(isoform1_exons)) {
    exon1 <- isoform1_exons[i, ]
    found_in_iso2 <- any(isoform2_exons$start == exon1$start & isoform2_exons$end == exon1$end)
    if (!found_in_iso2) {
      return(exon1)  # This is the skipped exon
    }
  }
  
  return(NULL)
}

find_flanking_exons <- function(all_exons, skipped_exon) {
  
  # Sort exons by position
  all_exons <- all_exons[order(all_exons$start), ]
  
  # Find position of skipped exon
  skipped_idx <- which(all_exons$start == skipped_exon$start & all_exons$end == skipped_exon$end)
  
  if (length(skipped_idx) == 0) {
    return(list(upstream = NULL, downstream = NULL))
  }
  
  skipped_idx <- skipped_idx[1]
  
  # Get upstream and downstream exons
  upstream <- if (skipped_idx > 1) all_exons[skipped_idx - 1, ] else NULL
  downstream <- if (skipped_idx < nrow(all_exons)) all_exons[skipped_idx + 1, ] else NULL
  
  return(list(upstream = upstream, downstream = downstream))
}

#===============================================================================
# ADDITIONAL HELPER FUNCTIONS FOR ALL EVENT TYPES
#===============================================================================

# Find alternative 3' splice site differences
find_alt_3prime_difference <- function(isoform1_exons, isoform2_exons) {
  
  # Look for exons with same start but different end (alternative 3' splice site)
  for (i in 1:nrow(isoform1_exons)) {
    exon1 <- isoform1_exons[i, ]
    for (j in 1:nrow(isoform2_exons)) {
      exon2 <- isoform2_exons[j, ]
      
      # Same start, different end = alternative 3' splice site
      if (exon1$start == exon2$start && exon1$end != exon2$end) {
        return(list(
          shared_start = exon1$start,
          end1 = exon1$end,
          end2 = exon2$end
        ))
      }
    }
  }
  
  return(NULL)
}

# Find alternative 5' splice site differences
find_alt_5prime_difference <- function(isoform1_exons, isoform2_exons) {
  
  # Look for exons with same end but different start (alternative 5' splice site)
  for (i in 1:nrow(isoform1_exons)) {
    exon1 <- isoform1_exons[i, ]
    for (j in 1:nrow(isoform2_exons)) {
      exon2 <- isoform2_exons[j, ]
      
      # Same end, different start = alternative 5' splice site
      if (exon1$end == exon2$end && exon1$start != exon2$start) {
        return(list(
          shared_end = exon1$end,
          start1 = exon1$start,
          start2 = exon2$start
        ))
      }
    }
  }
  
  return(NULL)
}

# Find flanking exon for alternative splice sites
find_flanking_exon_for_alt_splice <- function(all_exons, reference_coord) {
  
  # Find an exon that doesn't contain the reference coordinate
  for (i in 1:nrow(all_exons)) {
    exon <- all_exons[i, ]
    if (exon$start > reference_coord || exon$end < reference_coord) {
      return(exon)
    }
  }
  
  return(NULL)
}

# Find mutually exclusive exons
find_mutex_exons <- function(isoform1_exons, isoform2_exons) {
  
  # Find exons that exist in one isoform but not the other (and vice versa)
  mutex_exon1 <- NULL
  mutex_exon2 <- NULL
  
  # Find exon unique to isoform1
  for (i in 1:nrow(isoform1_exons)) {
    exon1 <- isoform1_exons[i, ]
    found_in_iso2 <- any(isoform2_exons$start == exon1$start & isoform2_exons$end == exon1$end)
    if (!found_in_iso2) {
      mutex_exon1 <- exon1
      break
    }
  }
  
  # Find exon unique to isoform2
  for (i in 1:nrow(isoform2_exons)) {
    exon2 <- isoform2_exons[i, ]
    found_in_iso1 <- any(isoform1_exons$start == exon2$start & isoform1_exons$end == exon2$end)
    if (!found_in_iso1) {
      mutex_exon2 <- exon2
      break
    }
  }
  
  if (is.null(mutex_exon1) || is.null(mutex_exon2)) {
    return(NULL)
  }
  
  return(list(exon1 = mutex_exon1, exon2 = mutex_exon2))
}

# Find flanking exons for mutually exclusive exons
find_flanking_exons_for_mutex <- function(all_exons, mutex_exon1, mutex_exon2) {
  
  # Sort exons by position
  all_exons <- all_exons[order(all_exons$start), ]
  
  # Find positions of mutex exons - FIX DATA TYPES
  mutex1_idx <- which(all_exons$start == mutex_exon1$start[1] & all_exons$end == mutex_exon1$end[1])
  mutex2_idx <- which(all_exons$start == mutex_exon2$start[1] & all_exons$end == mutex_exon2$end[1])
  
  if (length(mutex1_idx) == 0 || length(mutex2_idx) == 0) {
    return(list(upstream = NULL, downstream = NULL))
  }
  
  # Find shared flanking exons (exist before both mutex exons and after both)
  min_idx <- min(mutex1_idx[1], mutex2_idx[1])
  max_idx <- max(mutex1_idx[1], mutex2_idx[1])
  
  upstream <- if (min_idx > 1) all_exons[min_idx - 1, ] else NULL
  downstream <- if (max_idx < nrow(all_exons)) all_exons[max_idx + 1, ] else NULL
  
  return(list(upstream = upstream, downstream = downstream))
}

# Find retained intron
find_retained_intron <- function(isoform1_exons, isoform2_exons) {
  
  # Look for case where one isoform has a single exon spanning what's split in another
  
  # Check if isoform1 has a long exon that spans multiple exons in isoform2
  for (i in 1:nrow(isoform1_exons)) {
    long_exon <- isoform1_exons[i, ]
    
    # Find consecutive exons in isoform2 that fall within this long exon
    spanning_exons <- isoform2_exons[isoform2_exons$start >= long_exon$start & 
                                    isoform2_exons$end <= long_exon$end, ]
    
    if (nrow(spanning_exons) >= 2) {
      # This is a retained intron case
      spanning_exons <- spanning_exons[order(spanning_exons$start), ]
      
      return(list(
        retained_start = long_exon$start,
        retained_end = long_exon$end,
        upstream_start = spanning_exons[1, ]$start,
        upstream_end = spanning_exons[1, ]$end,
        downstream_start = spanning_exons[nrow(spanning_exons), ]$start,
        downstream_end = spanning_exons[nrow(spanning_exons), ]$end
      ))
    }
  }
  
  # Check the reverse case (isoform2 has the retained intron)
  for (i in 1:nrow(isoform2_exons)) {
    long_exon <- isoform2_exons[i, ]
    
    spanning_exons <- isoform1_exons[isoform1_exons$start >= long_exon$start & 
                                    isoform1_exons$end <= long_exon$end, ]
    
    if (nrow(spanning_exons) >= 2) {
      spanning_exons <- spanning_exons[order(spanning_exons$start), ]
      
      return(list(
        retained_start = long_exon$start,
        retained_end = long_exon$end,
        upstream_start = spanning_exons[1, ]$start,
        upstream_end = spanning_exons[1, ]$end,
        downstream_start = spanning_exons[nrow(spanning_exons), ]$start,
        downstream_end = spanning_exons[nrow(spanning_exons), ]$end
      ))
    }
  }
  
  return(NULL)
}

#===============================================================================
# ALTERNATIVE 3' SPLICE SITE (A3SS) CONVERSION
#===============================================================================

convert_alt_3prime_to_rmats_a3ss <- function(event_data) {
  
  # Get all gene-level events
  genes <- event_data[event_data$type == "gene", ]
  
  rmats_events <- data.frame()
  
  for (i in 1:nrow(genes)) {
    gene_entry <- genes[i, ]
    
    # Extract gene info from attributes
    gene_id <- extract_gene_id(gene_entry$attributes)
    event_id <- extract_event_id(gene_entry$attributes)
    
    # Get all mRNA isoforms for this gene
    gene_mrnas <- event_data[event_data$type == "mRNA" & 
                           grepl(paste0("Parent=", event_id), event_data$attributes), ]
    
    if (nrow(gene_mrnas) < 2) {
      next  # Need at least 2 isoforms to compare
    }
    
    # Get exons for each isoform
    isoform1_exons <- get_exons_for_isoform(event_data, gene_mrnas[1, ])
    isoform2_exons <- get_exons_for_isoform(event_data, gene_mrnas[2, ])
    
    # Find the alternative 3' splice site (different 3' end of same exon)
    alt_3ss_info <- find_alt_3prime_difference(isoform1_exons, isoform2_exons)
    
    if (is.null(alt_3ss_info)) {
      next
    }
    
    # Find flanking exon
    all_exons <- rbind(isoform1_exons, isoform2_exons)
    all_exons <- all_exons[!duplicated(paste(all_exons$start, all_exons$end)), ]
    all_exons <- all_exons[order(all_exons$start), ]
    
    flanking <- find_flanking_exon_for_alt_splice(all_exons, alt_3ss_info$shared_start)
    
    if (is.null(flanking)) {
      next
    }
    
    # Create rMATS A3SS format entry
    rmats_event <- data.frame(
      ID = as.character(extract_event_number(event_id)),
      GeneID = gene_id,
      geneSymbol = gene_id,
      chr = paste0("chr", gene_entry$seqid),
      strand = gene_entry$strand,
      longExonStart_0base = alt_3ss_info$shared_start - 1,  # 0-based
      longExonEnd = max(alt_3ss_info$end1, alt_3ss_info$end2),  # Long isoform end
      shortES = alt_3ss_info$shared_start - 1,  # 0-based 
      shortEE = min(alt_3ss_info$end1, alt_3ss_info$end2),  # Short isoform end
      flankingES = flanking$start - 1,  # 0-based
      flankingEE = flanking$end,  # 1-based
      # Add rMATS count columns (using dummy values for compatibility)
      IJC_SAMPLE_1 = "10,10,10",  # Inclusion junction counts
      SJC_SAMPLE_1 = "10,10,10",  # Skipping junction counts  
      IJC_SAMPLE_2 = "10,10,10",  # Inclusion junction counts
      SJC_SAMPLE_2 = "10,10,10",  # Skipping junction counts
      IncFormLen = as.integer(max(alt_3ss_info$end1, alt_3ss_info$end2) - alt_3ss_info$shared_start + 1),  # Long form length
      SkipFormLen = as.integer(min(alt_3ss_info$end1, alt_3ss_info$end2) - alt_3ss_info$shared_start + 1),  # Short form length
      stringsAsFactors = FALSE
    )
    
    rmats_events <- rbind(rmats_events, rmats_event)
  }
  
  return(rmats_events)
}

#===============================================================================
# ALTERNATIVE 5' SPLICE SITE (A5SS) CONVERSION
#===============================================================================

convert_alt_5prime_to_rmats_a5ss <- function(event_data) {
  
  # Get all gene-level events
  genes <- event_data[event_data$type == "gene", ]
  
  rmats_events <- data.frame()
  
  for (i in 1:nrow(genes)) {
    gene_entry <- genes[i, ]
    
    # Extract gene info from attributes
    gene_id <- extract_gene_id(gene_entry$attributes)
    event_id <- extract_event_id(gene_entry$attributes)
    
    # Get all mRNA isoforms for this gene
    gene_mrnas <- event_data[event_data$type == "mRNA" & 
                           grepl(paste0("Parent=", event_id), event_data$attributes), ]
    
    if (nrow(gene_mrnas) < 2) {
      next
    }
    
    # Get exons for each isoform
    isoform1_exons <- get_exons_for_isoform(event_data, gene_mrnas[1, ])
    isoform2_exons <- get_exons_for_isoform(event_data, gene_mrnas[2, ])
    
    # Find the alternative 5' splice site (different 5' start of same exon)
    alt_5ss_info <- find_alt_5prime_difference(isoform1_exons, isoform2_exons)
    
    if (is.null(alt_5ss_info)) {
      next
    }
    
    # Find flanking exon
    all_exons <- rbind(isoform1_exons, isoform2_exons)
    all_exons <- all_exons[!duplicated(paste(all_exons$start, all_exons$end)), ]
    all_exons <- all_exons[order(all_exons$start), ]
    
    flanking <- find_flanking_exon_for_alt_splice(all_exons, alt_5ss_info$shared_end)
    
    if (is.null(flanking)) {
      next
    }
    
    # Create rMATS A5SS format entry
    rmats_event <- data.frame(
      ID = as.character(extract_event_number(event_id)),
      GeneID = gene_id,
      geneSymbol = gene_id,
      chr = paste0("chr", gene_entry$seqid),
      strand = gene_entry$strand,
      longExonStart_0base = min(alt_5ss_info$start1, alt_5ss_info$start2) - 1,  # Long isoform start (0-based)
      longExonEnd = alt_5ss_info$shared_end,  # 1-based
      shortES = max(alt_5ss_info$start1, alt_5ss_info$start2) - 1,  # Short isoform start (0-based)
      shortEE = alt_5ss_info$shared_end,  # 1-based
      flankingES = flanking$start - 1,  # 0-based
      flankingEE = flanking$end,  # 1-based
      # Add rMATS count columns (using dummy values for compatibility)
      IJC_SAMPLE_1 = "10,10,10",  # Inclusion junction counts
      SJC_SAMPLE_1 = "10,10,10",  # Skipping junction counts  
      IJC_SAMPLE_2 = "10,10,10",  # Inclusion junction counts
      SJC_SAMPLE_2 = "10,10,10",  # Skipping junction counts
      IncFormLen = as.integer(alt_5ss_info$shared_end - min(alt_5ss_info$start1, alt_5ss_info$start2) + 1),  # Long form length
      SkipFormLen = as.integer(alt_5ss_info$shared_end - max(alt_5ss_info$start1, alt_5ss_info$start2) + 1),  # Short form length
      stringsAsFactors = FALSE
    )
    
    rmats_events <- rbind(rmats_events, rmats_event)
  }
  
  return(rmats_events)
}

#===============================================================================
# MUTUALLY EXCLUSIVE EXONS (MXE) CONVERSION
#===============================================================================

convert_mutex_exons_to_rmats_mxe <- function(event_data) {
  
  cat("=== CONVERTING MUTEX_EXONS TO rMATS MXE ===\n")
  
  # Get all gene-level events
  genes <- event_data[event_data$type == "gene", ]
  cat("Found", nrow(genes), "mutex_exons gene entries\n")
  
  rmats_events <- data.frame()
  
  for (i in 1:nrow(genes)) {
    gene_entry <- genes[i, ]
    
    # Extract gene info from attributes
    gene_id <- extract_gene_id(gene_entry$attributes)
    event_id <- extract_event_id(gene_entry$attributes)
    
    # Get all mRNA isoforms for this gene
    gene_mrnas <- event_data[event_data$type == "mRNA" & 
                           grepl(paste0("Parent=", event_id), event_data$attributes), ]
    
    if (nrow(gene_mrnas) < 2) {
      next
    }
    
    # Get exons for each isoform
    isoform1_exons <- get_exons_for_isoform(event_data, gene_mrnas[1, ])
    isoform2_exons <- get_exons_for_isoform(event_data, gene_mrnas[2, ])
    
    # Find mutually exclusive exons
    cat("Processing gene", i, ":", gene_id, "\n")
    cat("  Isoform1 exons:", nrow(isoform1_exons), "\n")
    cat("  Isoform2 exons:", nrow(isoform2_exons), "\n")
    
    mutex_info <- find_mutex_exons(isoform1_exons, isoform2_exons)
    
    if (is.null(mutex_info)) {
      cat("  No mutex exons found, skipping\n")
      next
    }
    
    cat("  Found mutex exons - exon1:", mutex_info$exon1$start[1], "-", mutex_info$exon1$end[1], 
        "exon2:", mutex_info$exon2$start[1], "-", mutex_info$exon2$end[1], "\n")
    
    # Find upstream and downstream exons
    all_exons <- rbind(isoform1_exons, isoform2_exons)
    all_exons <- all_exons[!duplicated(paste(all_exons$start, all_exons$end)), ]
    all_exons <- all_exons[order(all_exons$start), ]
    
    flanking <- find_flanking_exons_for_mutex(all_exons, mutex_info$exon1, mutex_info$exon2)
    
    if (is.null(flanking$upstream) || is.null(flanking$downstream)) {
      next
    }
    
    # Create rMATS MXE format entry - FIX DATA TYPES
    rmats_event <- data.frame(
      ID = as.character(extract_event_number(event_id)),
      GeneID = as.character(gene_id),
      geneSymbol = as.character(gene_id),
      chr = as.character(paste0("chr", gene_entry$seqid)),
      strand = as.character(gene_entry$strand),
      "X1stExonStart_0base" = as.integer(mutex_info$exon1$start[1] - 1),  # 0-based, force single value
      "X1stExonEnd" = as.integer(mutex_info$exon1$end[1]),  # 1-based, force single value
      "X2ndExonStart_0base" = as.integer(mutex_info$exon2$start[1] - 1),  # 0-based, force single value
      "X2ndExonEnd" = as.integer(mutex_info$exon2$end[1]),  # 1-based, force single value
      upstreamES = as.integer(flanking$upstream$start[1] - 1),  # 0-based, force single value
      upstreamEE = as.integer(flanking$upstream$end[1]),  # 1-based, force single value
      downstreamES = as.integer(flanking$downstream$start[1] - 1),  # 0-based, force single value
      downstreamEE = as.integer(flanking$downstream$end[1]),  # 1-based, force single value
      # Add rMATS count columns (using dummy values for compatibility)
      IJC_SAMPLE_1 = "10,10,10",  # Inclusion junction counts
      SJC_SAMPLE_1 = "10,10,10",  # Skipping junction counts  
      IJC_SAMPLE_2 = "10,10,10",  # Inclusion junction counts
      SJC_SAMPLE_2 = "10,10,10",  # Skipping junction counts
      IncFormLen = as.integer(mutex_info$exon1$end[1] - mutex_info$exon1$start[1] + 1),  # First exon length
      SkipFormLen = as.integer(mutex_info$exon2$end[1] - mutex_info$exon2$start[1] + 1),  # Second exon length
      stringsAsFactors = FALSE
    )
    
    # Ensure all columns are proper data types before binding
    rmats_event[] <- lapply(rmats_event, function(x) {
      if (is.list(x) && length(x) == 1) {
        return(x[[1]])  # Extract single values from lists
      }
      return(x)
    })
    cat("  Successfully created MXE event data.frame\n")
    
    rmats_events <- rbind(rmats_events, rmats_event)
  }
  
  return(rmats_events)
}

#===============================================================================
# RETAINED INTRON (RI) CONVERSION
#===============================================================================

convert_intron_retention_to_rmats_ri <- function(event_data) {
  
  # Get all gene-level events
  genes <- event_data[event_data$type == "gene", ]
  
  rmats_events <- data.frame()
  
  for (i in 1:nrow(genes)) {
    gene_entry <- genes[i, ]
    
    # Extract gene info from attributes
    gene_id <- extract_gene_id(gene_entry$attributes)
    event_id <- extract_event_id(gene_entry$attributes)
    
    # Get all mRNA isoforms for this gene
    gene_mrnas <- event_data[event_data$type == "mRNA" & 
                           grepl(paste0("Parent=", event_id), event_data$attributes), ]
    
    if (nrow(gene_mrnas) < 2) {
      next
    }
    
    # Get exons for each isoform
    isoform1_exons <- get_exons_for_isoform(event_data, gene_mrnas[1, ])
    isoform2_exons <- get_exons_for_isoform(event_data, gene_mrnas[2, ])
    
    # Find retained intron (one isoform spans across what's split in the other)
    ri_info <- find_retained_intron(isoform1_exons, isoform2_exons)
    
    if (is.null(ri_info)) {
      next
    }
    
    # Create rMATS RI format entry
    rmats_event <- data.frame(
      ID = as.character(extract_event_number(event_id)),
      GeneID = gene_id,
      geneSymbol = gene_id,
      chr = paste0("chr", gene_entry$seqid),
      strand = gene_entry$strand,
      riExonStart_0base = ri_info$retained_start - 1,  # 0-based
      riExonEnd = ri_info$retained_end,  # 1-based
      upstreamES = ri_info$upstream_start - 1,  # 0-based
      upstreamEE = ri_info$upstream_end,  # 1-based
      downstreamES = ri_info$downstream_start - 1,  # 0-based
      downstreamEE = ri_info$downstream_end,  # 1-based
      # Add rMATS count columns (using dummy values for compatibility)
      IJC_SAMPLE_1 = "10,10,10",  # Inclusion junction counts
      SJC_SAMPLE_1 = "10,10,10",  # Skipping junction counts  
      IJC_SAMPLE_2 = "10,10,10",  # Inclusion junction counts
      SJC_SAMPLE_2 = "10,10,10",  # Skipping junction counts
      IncFormLen = as.integer(ri_info$retained_end - ri_info$retained_start + 1),  # Retained intron length
      SkipFormLen = as.integer(0),  # Skip form has no intron (length = 0)
      stringsAsFactors = FALSE
    )
    
    rmats_events <- rbind(rmats_events, rmats_event)
  }
  
  return(rmats_events)
}

#===============================================================================
# HELPER FUNCTIONS FOR EVENT DETECTION
#===============================================================================

#' Find alternative 3' splice site differences between isoforms
find_alt_3prime_difference <- function(isoform1_exons, isoform2_exons) {
  # Look for exons with same start but different ends
  for (i in 1:nrow(isoform1_exons)) {
    exon1 <- isoform1_exons[i, ]
    
    # Find exon in isoform2 with same start
    matching_exon <- isoform2_exons[isoform2_exons$start == exon1$start, ]
    
    if (nrow(matching_exon) == 1 && matching_exon$end != exon1$end) {
      return(list(
        shared_start = exon1$start,
        end1 = exon1$end,
        end2 = matching_exon$end
      ))
    }
  }
  
  return(NULL)
}

#' Find alternative 5' splice site differences between isoforms  
find_alt_5prime_difference <- function(isoform1_exons, isoform2_exons) {
  # Look for exons with same end but different starts
  for (i in 1:nrow(isoform1_exons)) {
    exon1 <- isoform1_exons[i, ]
    
    # Find exon in isoform2 with same end
    matching_exon <- isoform2_exons[isoform2_exons$end == exon1$end, ]
    
    if (nrow(matching_exon) == 1 && matching_exon$start != exon1$start) {
      return(list(
        shared_end = exon1$end,
        start1 = exon1$start,
        start2 = matching_exon$start
      ))
    }
  }
  
  return(NULL)
}

#' Find mutually exclusive exons between isoforms
find_mutex_exons <- function(isoform1_exons, isoform2_exons) {
  # Find exons that exist in one isoform but not the other
  iso1_coords <- paste(isoform1_exons$start, isoform1_exons$end)
  iso2_coords <- paste(isoform2_exons$start, isoform2_exons$end)
  
  # Exons unique to isoform 1
  iso1_unique_idx <- which(!iso1_coords %in% iso2_coords)
  # Exons unique to isoform 2  
  iso2_unique_idx <- which(!iso2_coords %in% iso1_coords)
  
  if (length(iso1_unique_idx) >= 1 && length(iso2_unique_idx) >= 1) {
    # Take first unique exon from each
    exon1 <- isoform1_exons[iso1_unique_idx[1], ]
    exon2 <- isoform2_exons[iso2_unique_idx[1], ]
    
    return(list(
      exon1 = exon1,
      exon2 = exon2
    ))
  }
  
  return(NULL)
}

#' Find retained intron between isoforms
find_retained_intron <- function(isoform1_exons, isoform2_exons) {
  # Look for case where one isoform has a single large exon 
  # and the other has multiple smaller exons that would span it
  
  # Check if isoform1 has fewer exons (potential retention)
  if (nrow(isoform1_exons) < nrow(isoform2_exons)) {
    long_exon <- isoform1_exons[1, ]  # Assume single retained exon
    spanning_exons <- isoform2_exons[
      isoform2_exons$start >= long_exon$start & 
      isoform2_exons$end <= long_exon$end, ]
    
    if (nrow(spanning_exons) >= 2) {
      spanning_exons <- spanning_exons[order(spanning_exons$start), ]
      
      return(list(
        retained_start = long_exon$start,
        retained_end = long_exon$end,
        upstream_start = spanning_exons[1, ]$start,
        upstream_end = spanning_exons[1, ]$end,
        downstream_start = spanning_exons[nrow(spanning_exons), ]$start,
        downstream_end = spanning_exons[nrow(spanning_exons), ]$end
      ))
    }
  }
  
  # Check if isoform2 has fewer exons (potential retention)
  if (nrow(isoform2_exons) < nrow(isoform1_exons)) {
    long_exon <- isoform2_exons[1, ]  # Assume single retained exon
    spanning_exons <- isoform1_exons[
      isoform1_exons$start >= long_exon$start & 
      isoform1_exons$end <= long_exon$end, ]
    
    if (nrow(spanning_exons) >= 2) {
      spanning_exons <- spanning_exons[order(spanning_exons$start), ]
      
      return(list(
        retained_start = long_exon$start,
        retained_end = long_exon$end,
        upstream_start = spanning_exons[1, ]$start,
        upstream_end = spanning_exons[1, ]$end,
        downstream_start = spanning_exons[nrow(spanning_exons), ]$start,
        downstream_end = spanning_exons[nrow(spanning_exons), ]$end
      ))
    }
  }
  
  return(NULL)
}

#' Find flanking exon for alternative splice site events
find_flanking_exon_for_alt_splice <- function(all_exons, reference_coord) {
  # Find the closest neighboring exon
  all_exons <- all_exons[order(all_exons$start), ]
  
  # Find exons that don't contain the reference coordinate
  candidates <- all_exons[!(all_exons$start <= reference_coord & all_exons$end >= reference_coord), ]
  
  if (nrow(candidates) > 0) {
    # Return closest one
    distances <- abs((candidates$start + candidates$end) / 2 - reference_coord)
    closest_idx <- which.min(distances)
    return(candidates[closest_idx, ])
  }
  
  return(NULL)
}

#' Find upstream and downstream flanking exons for mutex events
find_mutex_flanking_exons <- function(all_exons, mutex_info) {
  all_exons <- all_exons[order(all_exons$start), ]
  
  # Find exons flanking the mutex region
  mutex_start <- min(mutex_info$exon1$start, mutex_info$exon2$start)
  mutex_end <- max(mutex_info$exon1$end, mutex_info$exon2$end)
  
  # Upstream exon (before mutex region)
  upstream_candidates <- all_exons[all_exons$end < mutex_start, ]
  upstream <- if (nrow(upstream_candidates) > 0) {
    upstream_candidates[nrow(upstream_candidates), ]  # Closest upstream
  } else {
    NULL
  }
  
  # Downstream exon (after mutex region) 
  downstream_candidates <- all_exons[all_exons$start > mutex_end, ]
  downstream <- if (nrow(downstream_candidates) > 0) {
    downstream_candidates[1, ]  # Closest downstream
  } else {
    NULL
  }
  
  if (is.null(upstream) || is.null(downstream)) {
    return(NULL)
  }
  
  return(list(
    upstream = upstream,
    downstream = downstream
  ))
}

cat("✅ SplAdder GFF3 converter loaded\n")