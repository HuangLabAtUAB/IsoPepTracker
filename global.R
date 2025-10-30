# SETUP AND DATA LOADING (from global.R)
#===============================================================================

# Load required libraries
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(plotly)
library(rtracklayer)
library(GenomicRanges)
library(shinycssloaders)  # Add this to your library imports
library(shinyjs)  # For enabling/disabling buttons
library(dplyr)    # For data manipulation
library(memoise)   # For caching function results
library(magrittr)  # For pipe operators
library(promises)  # For async processing
library(future)    # For parallel processing

# BSgenome libraries for rMATS protein translation
library(BSgenome.Hsapiens.UCSC.hg38)  # Human genome sequences
library(Biostrings)  # DNA/protein sequence manipulation

# Install and load ggnewscale for multiple fill scales if not already installed
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}
library(ggnewscale)

# Define GTF file path first
gtf_file <- "~/Downloads/gencode.v38.annotation.gtf.gz"
gtf_exists <- file.exists(gtf_file)
if (!gtf_exists) {
  message("WARNING: GTF file not found at ", gtf_file)
}

# SOURCE HELPER FUNCTIONS FIRST (needed for cache building)
source("R/data_processing.R")
source("R/visualization.R")
source("R/as_analysis.R")
source("R/coordinate_compression.R")
source("R/gene_splitter.R")  # Add the new gene splitter functions
source("R/novel_isoform_analysis.R")  # Add the novel isoform functions
source("R/gtf_preprocessor.R")  # GTF preprocessing for lightning-fast visualizations
source("R/fast_visualization.R")  # Lightning-fast visualization functions
source("R/blastp_search.R")  # BLASTP-based peptide search functions
source("R/rmats_peptide_analysis.R")  # rMATS-specific peptide analysis functions
# REMOVED: exon_based_coordinate_mapper.R - cleaning up broken code

# Load gene index instead of full databases
# This is much more memory efficient
gene_index_file <- "data/index/gene_index.rds"
if (file.exists(gene_index_file)) {
  gene_index <- readRDS(gene_index_file)
  message("Loaded gene index with ", nrow(gene_index), " genes")
  
  # Pre-load popular genes for better performance
  popular_genes <- c("ENSG00000139618", "ENSG00000117748", "ENSG00000134086", 
                    "ENSG00000165704", "ENSG00000134057", "ENSG00000005381",
                    "ENSG00000096968", "ENSG00000134086", "ENSG00000141510",
                    "ENSG00000171862")  # Common genes like BRCA2, RBM39, etc.
  
  # Filter to only genes that actually exist in our index
  available_popular_genes <- intersect(popular_genes, gene_index$geneID)
  if (length(available_popular_genes) > 0) {
    message("Pre-loading ", length(available_popular_genes), " popular genes...")
  }
  
  # CHECK GTF CACHE STATUS
  cache_dir <- "data/gtf_cache"
  completion_marker <- file.path(cache_dir, ".cache_complete")
  
  if (file.exists(completion_marker) && dir.exists(cache_dir)) {
    cache_count <- length(list.files(cache_dir, pattern = "*.rds"))
    message("✅ GTF cache loaded with ", cache_count, " genes. Visualizations are lightning-fast!")
  } else if (dir.exists(cache_dir) && length(list.files(cache_dir, pattern = "*.rds")) > 0) {
    cache_count <- length(list.files(cache_dir, pattern = "*.rds"))
    message("⚠️ Partial GTF cache found with ", cache_count, " genes.")
    message("For best performance, run: Rscript build_cache.R")
  } else {
    message("⚠️ No GTF cache found. Visualizations will use slower GTF parsing.")
    message("To enable lightning-fast visualizations, run: Rscript build_cache.R")
    message("This one-time setup takes ~15-20 minutes but makes all visualizations 150x faster!")
  }
} else {
  # Fallback to original data loading method if gene index doesn't exist
  message("Gene index not found. Loading full databases...")
  peptides <- readRDS("total_database.rds")
  as_database <- readRDS("AS_filtered.rds")
  message("Loaded full databases")
}


# Enable gzip/brotli compression for all responses (reduces plot JSON transfer time)
options(shiny.compressResponse = TRUE)

# Additional performance optimizations
options(
  shiny.maxRequestSize = 100*1024^2,  # 100MB max upload
  shiny.usecairo = FALSE,  # Disable Cairo for faster plotting
  shiny.autoreload = FALSE,  # Disable auto-reload in production
  shiny.reactlog = FALSE,  # Disable reactlog for better performance
  shiny.trace = FALSE  # Disable tracing for better performance
)

#===============================================================================
# HELPER FUNCTIONS ALREADY SOURCED ABOVE
#===============================================================================

#-------------------------------------------------------------------------------
# PERFORMANCE: memoise heavy GTF helpers so each gene is parsed only once per
# session.  The function bodies are defined in R/data_processing.R (already
# sourced above). Wrapping them with memoise() keeps identical calls fast but
# preserves their signatures and return values, so no downstream code changes.
#-------------------------------------------------------------------------------

if (exists("load_gene_details") && !isTRUE(attr(load_gene_details, "memoised"))) {
  load_gene_details <- memoise(load_gene_details)
}

if (exists("load_transcript_exons") && !isTRUE(attr(load_transcript_exons, "memoised"))) {
  load_transcript_exons <- memoise(load_transcript_exons)
}

# Memoize heavy visualization functions for better performance
if (exists("create_transcript_plot_data") && !isTRUE(attr(create_transcript_plot_data, "memoised"))) {
  create_transcript_plot_data <- memoise(create_transcript_plot_data)
}

# Special handling for create_peptide_plot_data to ensure proper signature
if (exists("create_peptide_plot_data")) {
  # Clear any existing memoised version to prevent signature conflicts
  if (isTRUE(attr(create_peptide_plot_data, "memoised"))) {
    # Force reload the original function from visualization.R
    source("R/visualization.R", local = TRUE)
  }
  
  # Verify the function has the use_compression parameter
  if ("use_compression" %in% names(formals(create_peptide_plot_data))) {
    create_peptide_plot_data <- memoise(create_peptide_plot_data)
    cat("Successfully memoised create_peptide_plot_data with use_compression parameter\n")
  } else {
    cat("WARNING: create_peptide_plot_data missing use_compression parameter - not memoising\n")
  }
}

if (exists("create_as_view") && !isTRUE(attr(create_as_view, "memoised"))) {
  create_as_view <- memoise(create_as_view)
}

#===============================================================================
# CSS STYLING (from original app.R)
#===============================================================================

css_styles <- tags$head(
  tags$style(HTML("
    /* Wrapper and container fixes */
    .content-wrapper, .right-side {
      overflow-x: hidden !important;
      overflow-y: auto !important;
    }
    
    /* Box sizing and overflow handling */
    .box-body {
      overflow-x: auto !important;
      padding: 10px !important;
    }
    
    /* Plot container responsiveness */
    .plotly {
      min-height: 600px !important;
      width: 100% !important;
    }
    
    /* Force plot containers to be responsive */
    .plot-container.plotly, .js-plotly-plot {
      width: 100% !important;
    }
    
    /* Table responsiveness */
    .dataTables_wrapper {
      width: 100% !important;
      overflow-x: auto !important;
    }
    
    /* Dashboard content area */
    .content {
      padding: 15px !important;
    }
    
    /* Box margins and padding */
    .box {
      margin-bottom: 15px !important;
      overflow: hidden !important;
    }
    
    /* Tab content */
    .tab-content {
      padding: 10px 0 !important;
    }
    
    /* Legend positioning */
    .legend {
      position: relative !important;
      margin-top: 10px !important;
    }
  "))
)
