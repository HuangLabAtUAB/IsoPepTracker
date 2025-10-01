#===============================================================================
# PURE CARD-BASED UI - NO SIDEBAR
#===============================================================================
library("shiny")
library("shinydashboard") 
library("shinyjs")
library("shinycssloaders")
library("DT")
library("plotly")
# Enhanced styling with card-based navigation
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
    
    /* Analysis view styling */
    .analysis-container {
      padding: 20px;
      margin: 0;
      min-height: 100vh;
      background: #f4f4f4;
    }
    
    .view-header {
      background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
      color: white;
      padding: 20px 30px;
      margin: -20px -20px 20px -20px;
      border-radius: 0;
      display: flex;
      justify-content: space-between;
      align-items: center;
    }
    
    .view-title {
      font-size: 1.8rem;
      font-weight: 300;
      margin: 0;
      text-shadow: 0 1px 3px rgba(0,0,0,0.3);
    }
    
    .back-btn {
      background: rgba(255, 255, 255, 0.2) !important;
      border: 1px solid rgba(255, 255, 255, 0.3) !important;
      color: white !important;
      padding: 8px 16px !important;
      border-radius: 20px !important;
      font-size: 0.9rem !important;
      transition: all 0.3s ease !important;
    }
    
    .back-btn:hover {
      background: rgba(255, 255, 255, 0.3) !important;
      color: white !important;
      text-decoration: none !important;
    }
    
    /* Tab navigation styling */
    .tab-btn {
      background: #ffffff !important;
      color: #666 !important;
      border: 1px solid #ddd !important;
    }
    
    .tab-btn:hover {
      background: #f8f9fa !important;
      color: #333 !important;
    }
    
    .tab-btn.active {
      background: #667eea !important;
      color: white !important;
      border-color: #667eea !important;
    }
    
    /* Mobile responsiveness for card grid */
    @media (max-width: 768px) {
      .card-grid {
        grid-template-columns: 1fr !important;
        grid-template-rows: auto !important;
        max-width: 400px !important;
        gap: 20px !important;
      }
      
      .card-grid > div {
        width: 100% !important;
        max-width: 350px !important;
      }
    }
  ")),
  
  # JavaScript for tab navigation
  tags$script(HTML("
    $(document).ready(function() {
      // Initialize active tab
      $('.tab-overview').addClass('active');
      
      // Handle tab clicks
      $('.tab-btn').on('click', function() {
        // Remove active class from all tabs
        $('.tab-btn').removeClass('active');
        
        // Add active class to clicked tab
        $(this).addClass('active');
      });
      
      // Listen for programmatic tab changes
      $(document).on('shiny:inputchanged', function(event) {
        if (event.name === 'canonical_current_tab') {
          // Remove active class from all tabs
          $('.tab-btn').removeClass('active');
          
          // Add active class to the corresponding tab
          if (event.value === 'overview') {
            $('.tab-overview').addClass('active');
          } else if (event.value === 'events') {
            $('.tab-events').addClass('active');
          } else if (event.value === 'isoforms') {
            $('.tab-isoforms').addClass('active');
          }
        } else if (event.name === 'alt_splicing_current_tab') {
          // Remove active class from all tabs
          $('.tab-btn').removeClass('active');
          
          // Add active class to the corresponding tab
          if (event.value === 'rmats') {
            $('.tab-rmats').addClass('active');
          } else if (event.value === 'spladder') {
            $('.tab-spladder').addClass('active');
          }
        }
      });
    });
  "))
)
#===============================================================================
# CARD-BASED USER INTERFACE - PURE FLUIDPAGE
#===============================================================================
ui <- fluidPage(
  css_styles,
  shinyjs::useShinyjs(),
  
  # Current view state - hidden from user
  div(style = "display: none;",
    selectInput("current_view", NULL, 
                choices = c("home", "canonical", "novel", "peptide", "alt_splicing"), 
                selected = "home")
  ),
  
  # Home View - Professional Landing Page
  conditionalPanel(
    condition = "input.current_view == 'home'",
    div(
      style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
               min-height: 100vh; 
               margin: 0; 
               padding: 0; 
               position: relative;
               overflow: hidden;",
      
      # Header Section
      div(
        style = "background: rgba(255, 255, 255, 0.1); 
                 backdrop-filter: blur(10px); 
                 padding: 20px 40px; 
                 border-bottom: 1px solid rgba(255, 255, 255, 0.2);",
        div(
          style = "display: flex; justify-content: space-between; align-items: center;",
          div(
            style = "display: flex; align-items: center; gap: 15px; color: white;",
            div(
              style = "background: rgba(255, 255, 255, 0.2); 
                       padding: 12px; 
                       border-radius: 12px; 
                       font-size: 24px;",
              icon("dna")
            ),
            div(
              h2("IsoPepTracker", style = "margin: 0; font-weight: 300; font-size: 28px;"),
              p("Proteomics Analysis Platform", style = "margin: 0; opacity: 0.8; font-size: 14px;")
            )
          )
        )
      ),
      
      # Main Content Section
      div(
        style = "padding: 60px 40px; text-align: center;",
        
        # Title and Subtitle
        h1("Proteomics Analysis Suite", 
           style = "color: white; 
                    font-size: 3.5rem; 
                    font-weight: 300; 
                    margin-bottom: 20px; 
                    text-shadow: 0 2px 4px rgba(0,0,0,0.3);"),
        
        p("Advanced tools for analyzing protein isoforms, discovering novel transcripts, and performing comprehensive peptide sequence analysis with state-of-the-art algorithms.",
          style = "color: rgba(255, 255, 255, 0.9); 
                   font-size: 1.3rem; 
                   font-weight: 300; 
                   max-width: 800px; 
                   margin: 0 auto 60px auto; 
                   line-height: 1.6;"),
        
        # Feature Cards Grid - 2x2 layout
        div(
          class = "card-grid",
          style = "display: grid; 
                   grid-template-columns: repeat(2, 1fr); 
                   grid-template-rows: repeat(2, 1fr);
                   gap: 25px; 
                   max-width: 900px; 
                   margin: 0 auto;
                   justify-items: center;
                   align-items: stretch;",
          
          # Card 1: Canonical Isoform Analysis
          div(
            onclick = "Shiny.setInputValue('current_view', 'canonical');",
            style = "background: rgba(255, 255, 255, 0.15); 
                     backdrop-filter: blur(10px); 
                     border: 1px solid rgba(255, 255, 255, 0.2); 
                     border-radius: 20px; 
                     padding: 40px; 
                     cursor: pointer; 
                     transition: all 0.3s ease; 
                     color: white; 
                     text-align: center;
                     width: 380px;
                     height: 420px;
                     display: flex;
                     flex-direction: column;
                     justify-content: center;",
            onmouseover = "this.style.transform='translateY(-8px)'; 
                           this.style.background='rgba(255, 255, 255, 0.25)'; 
                           this.style.borderColor='rgba(255, 255, 255, 0.4)'; 
                           this.style.boxShadow='0 15px 35px rgba(0,0,0,0.2)';",
            onmouseout = "this.style.transform='translateY(0px)'; 
                          this.style.background='rgba(255, 255, 255, 0.15)'; 
                          this.style.borderColor='rgba(255, 255, 255, 0.2)'; 
                          this.style.boxShadow='none';",
            
            div(
              style = "background: rgba(102, 126, 234, 0.3); 
                       border-radius: 50%; 
                       width: 80px; 
                       height: 120px; 
                       display: flex; 
                       align-items: center; 
                       justify-content: center; 
                       margin: 0 auto 20px auto;",
              icon("dna", style = "font-size: 28px;")
            ),
            h3("Canonical Isoform Analysis", 
               style = "font-size: 1.4rem; 
                        margin-bottom: 15px; 
                        font-weight: 500;"),
            p("Analyze known gene structures and alternative splicing patterns with comprehensive protease digestion predictions.",
              style = "font-size: 1rem; 
                       line-height: 1.6; 
                       opacity: 0.9; 
                       margin: 0;")
          ),
          
          # Card 2: Novel Isoform Discovery
          div(
            onclick = "Shiny.setInputValue('current_view', 'novel');",
            style = "background: rgba(255, 255, 255, 0.15); 
                     backdrop-filter: blur(10px); 
                     border: 1px solid rgba(255, 255, 255, 0.2); 
                     border-radius: 20px; 
                     padding: 40px; 
                     cursor: pointer; 
                     transition: all 0.3s ease; 
                     color: white; 
                     text-align: center;
                     width: 380px;
                     height: 420px;
                     display: flex;
                     flex-direction: column;
                     justify-content: center;",
            onmouseover = "this.style.transform='translateY(-8px)'; 
                           this.style.background='rgba(255, 255, 255, 0.25)'; 
                           this.style.borderColor='rgba(255, 255, 255, 0.4)'; 
                           this.style.boxShadow='0 15px 35px rgba(0,0,0,0.2)';",
            onmouseout = "this.style.transform='translateY(0px)'; 
                          this.style.background='rgba(255, 255, 255, 0.15)'; 
                          this.style.borderColor='rgba(255, 255, 255, 0.2)'; 
                          this.style.boxShadow='none';",
            
            div(
              style = "background: rgba(118, 75, 162, 0.3); 
                       border-radius: 50%; 
                       width: 80px; 
                       height: 120px; 
                       display: flex; 
                       align-items: center; 
                       justify-content: center; 
                       margin: 0 auto 20px auto;",
              icon("rocket", style = "font-size: 28px;")
            ),
            h3("Novel Isoform Discovery", 
               style = "font-size: 1.4rem; 
                        margin-bottom: 15px; 
                        font-weight: 500;"),
            p("Upload transcript sequences to identify novel isoforms and compare against known proteomes with advanced algorithms.",
              style = "font-size: 1rem; 
                       line-height: 1.6; 
                       opacity: 0.9; 
                       margin: 0;")
          ),
          
          # Card 3: Peptide Sequence Search
          div(
            onclick = "Shiny.setInputValue('current_view', 'peptide');",
            style = "background: rgba(255, 255, 255, 0.15); 
                     backdrop-filter: blur(10px); 
                     border: 1px solid rgba(255, 255, 255, 0.2); 
                     border-radius: 20px; 
                     padding: 40px; 
                     cursor: pointer; 
                     transition: all 0.3s ease; 
                     color: white; 
                     text-align: center;
                     width: 380px;
                     height: 420px;
                     display: flex;
                     flex-direction: column;
                     justify-content: center;",
            onmouseover = "this.style.transform='translateY(-8px)'; 
                           this.style.background='rgba(255, 255, 255, 0.25)'; 
                           this.style.borderColor='rgba(255, 255, 255, 0.4)'; 
                           this.style.boxShadow='0 15px 35px rgba(0,0,0,0.2)';",
            onmouseout = "this.style.transform='translateY(0px)'; 
                          this.style.background='rgba(255, 255, 255, 0.15)'; 
                          this.style.borderColor='rgba(255, 255, 255, 0.2)'; 
                          this.style.boxShadow='none';",
            
            div(
              style = "background: rgba(102, 126, 234, 0.3); 
                       border-radius: 50%; 
                       width: 80px; 
                       height: 120px; 
                       display: flex; 
                       align-items: center; 
                       justify-content: center; 
                       margin: 0 auto 20px auto;",
              icon("search", style = "font-size: 28px;")
            ),
            h3("Peptide Sequence Search", 
               style = "font-size: 1.4rem; 
                        margin-bottom: 15px; 
                        font-weight: 500;"),
            p("Perform BLASTP searches against comprehensive protein databases with statistical analysis and advanced filtering.",
              style = "font-size: 1rem; 
                       line-height: 1.6; 
                       opacity: 0.9; 
                       margin: 0;")
          ),
          
          # Card 4: Alternative Splicing Analysis  
          div(
            onclick = "Shiny.setInputValue('current_view', 'alt_splicing');",
            style = "background: rgba(255, 255, 255, 0.15); 
                     backdrop-filter: blur(10px); 
                     border: 1px solid rgba(255, 255, 255, 0.2); 
                     border-radius: 20px; 
                     padding: 40px; 
                     cursor: pointer; 
                     transition: all 0.3s ease; 
                     color: white; 
                     text-align: center;
                     width: 380px;
                     height: 420px;
                     display: flex;
                     flex-direction: column;
                     justify-content: center;",
            onmouseover = "this.style.transform='translateY(-8px)'; 
                           this.style.background='rgba(255, 255, 255, 0.25)'; 
                           this.style.borderColor='rgba(255, 255, 255, 0.4)'; 
                           this.style.boxShadow='0 15px 35px rgba(0,0,0,0.2)';",
            onmouseout = "this.style.transform='translateY(0px)'; 
                          this.style.background='rgba(255, 255, 255, 0.15)'; 
                          this.style.borderColor='rgba(255, 255, 255, 0.2)'; 
                          this.style.boxShadow='none';",
            
            div(
              style = "background: rgba(118, 75, 162, 0.3); 
                       border-radius: 50%; 
                       width: 80px; 
                       height: 120px; 
                       display: flex; 
                       align-items: center; 
                       justify-content: center; 
                       margin: 0 auto 20px auto;",
              icon("project-diagram", style = "font-size: 28px;")
            ),
            h3("Alternative Splicing Analysis", 
               style = "font-size: 1.4rem; 
                        margin-bottom: 15px; 
                        font-weight: 500;"),
            p("Comprehensive analysis of alternative splicing events using rMATS with functional consequence prediction.",
              style = "font-size: 1rem; 
                       line-height: 1.6; 
                       opacity: 0.9; 
                       margin: 0;")
          )
        )
      )
    )
  ), # End Home conditionalPanel
  
  # Canonical Analysis View
  conditionalPanel(
    condition = "input.current_view == 'canonical'",
    div(class = "analysis-container",
      div(class = "view-header",
        h2("Canonical Isoform Analysis", class = "view-title"),
        actionButton("back_to_home_canonical", "← Back to Home", class = "back-btn", 
                    onclick = "Shiny.setInputValue('current_view', 'home');")
      ),
      
      # Hidden tab navigation for canonical analysis
      div(style = "display: none;",
        selectInput("canonical_current_tab", NULL, 
                    choices = c("overview", "events", "isoforms"), 
                    selected = "overview")
      ),
      
      # Analysis Controls Section
      fluidRow(
        box(
          width = 12,
          title = "Analysis Parameters",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          fluidRow(
            column(3,
              textInput("gene_search", "Search Gene:", value = "EGFR", placeholder = "Type gene symbol or ID...")
            ),
            column(3,
              selectizeInput("gene", "Select Gene:", choices = NULL, options = list(
                placeholder = "Type to search genes",
                maxOptions = 20
              ))
            ),
            column(2,
              selectizeInput("protease", "Select Enzyme:", 
                    choices = c("Trypsin" = "trp", 
                               "Chymotrypsin" = "chymo",
                               "AspN" = "aspn",
                               "LysC" = "lysc",
                               "LysN" = "lysn", 
                               "GluC" = "gluc"),
                    selected = "trp")
            ),
            column(2,
              selectInput("miscleavage_type", "Miscleavage Allowance:", 
                 choices = c("No miscleavage" = "no_miss_cleavage",
                            "Up to 2 miscleavages" = "upto_two_misscleavage"),
                 selected = "no_miss_cleavage")
            ),
            column(2,
              radioButtons("intron_scale", "Intron Display:",
                          choices = list(
                            "True scale" = "true",
                            "Compressed" = "compressed"
                          ),
                          selected = "compressed",
                          inline = FALSE)
            ),
            column(2,
              div(style = "margin-top: 25px;",
                actionButton("update", "Update View", icon = icon("refresh"), width = "100%")
              )
            )
          )
        )
      ),
      
      # Tab Navigation Bar
      div(style = "margin: 20px 0;",
        div(style = "display: flex; background: #f4f4f4; border-radius: 6px; padding: 4px;",
          div(
            onclick = "Shiny.setInputValue('canonical_current_tab', 'overview');",
            style = "flex: 1; text-align: center; padding: 12px; cursor: pointer; border-radius: 4px; font-weight: 500; transition: all 0.3s ease;",
            class = "tab-btn tab-overview",
            "Overview"
          ),
          div(
            onclick = "Shiny.setInputValue('canonical_current_tab', 'events');",
            style = "flex: 1; text-align: center; padding: 12px; cursor: pointer; border-radius: 4px; font-weight: 500; transition: all 0.3s ease;",
            class = "tab-btn tab-events",
            "Events-centric view"
          ),
          div(
            onclick = "Shiny.setInputValue('canonical_current_tab', 'isoforms');",
            style = "flex: 1; text-align: center; padding: 12px; cursor: pointer; border-radius: 4px; font-weight: 500; transition: all 0.3s ease;",
            class = "tab-btn tab-isoforms",
            "Isoform-centric view"
          )
        )
      ),
      
      # Tab 1: Overview
      conditionalPanel(
        condition = "input.canonical_current_tab == 'overview'",
        fluidRow(
          box(width = 12,
            title = "Gene Visualization Overview",
            status = "primary",
            solidHeader = TRUE,
            div(
              style = "background-color: #f0f8ff; padding: 15px; border-left: 4px solid #3c8dbc; margin-bottom: 15px;",
              h4("What You Can Do Here", style = "margin-top: 0; color: #3c8dbc;"),
              p("Explore how your gene looks and see where peptides come from when you digest proteins with different enzymes. You can see your gene's structure with all different transcript variants, find exactly which peptides you'll get from each transcript, choose your enzyme to compare how different proteases cut your protein, and interactively explore the data by hovering, zooming, and clicking to dig deeper.", style = "margin-bottom: 10px;")
            )
          )
        ),
        conditionalPanel(
          condition = "input.gene && input.update > 0",
          fluidRow(
            box(width = 12,
              title = textOutput("peptide_plot_title"),
              status = "primary",
              solidHeader = TRUE,
              withSpinner(
                plotlyOutput("peptide_plot", height = "600px"),
                type = 8,
                color = "#3c8dbc"
              )
            )
          )
        )
      ),
      
      # Tab 2: Events-centric view
      conditionalPanel(
        condition = "input.canonical_current_tab == 'events'",
        fluidRow(
          box(width = 12,
            title = "Alternative Splicing Analysis Overview",
            status = "primary", 
            solidHeader = TRUE,
            div(
              style = "background-color: #f0f8ff; padding: 15px; border-left: 4px solid #3c8dbc; margin-bottom: 15px;",
              h4("What You Can Do Here", style = "margin-top: 0; color: #3c8dbc;"),
              p("Discover how alternative splicing creates different peptides. Find out which peptides are unique to specific splice variants! Browse all splicing events to see every way your gene gets spliced, investigate specific events by clicking to dive deep into individual splicing changes, compare splice variants with side-by-side peptide comparison between versions, and see the big picture by visualizing all splicing differences at once.", style = "margin-bottom: 10px;")
            )
          )
        ),
        conditionalPanel(
          condition = "input.gene && input.update > 0",
          fluidRow(
            box(width = 12,
              title = "Alternative Splicing Structure",
              status = "primary",
              solidHeader = TRUE,
              withSpinner(
                plotlyOutput("as_structure_plot", height = "300px"),
                type = 8,
                color = "#3c8dbc"
              )
            )
          ),
          fluidRow(
            box(width = 12,
              title = "Alternative Splicing Events",
              status = "primary",
              solidHeader = TRUE,
              withSpinner(
                DT::dataTableOutput("as_events_table"),
                type = 8,
                color = "#3c8dbc"
              )
            )
          )
        ),
        conditionalPanel(
          condition = "output.event_selected",
          fluidRow(
            box(width = 12,
              title = "Selected Event Information",
              status = "primary",
              solidHeader = TRUE,
              uiOutput("selected_as_event_info")
            )
          ),
          fluidRow(
            box(
              width = 12,
              title = "Peptide Comparison",
              status = "primary",
              solidHeader = TRUE,
              fluidRow(
                column(8, uiOutput("comparison_pair_selector")),
                column(2, 
                  div(style = "margin-top: 25px;",
                    checkboxInput("show_all_peptides", "Show all peptides", value = FALSE)
                  )
                ),
                column(2, 
                  div(style = "margin-top: 25px;",
                    downloadButton("download_peptide_data", "Download Peptide Data")
                  )
                )
              ),
              actionButton("compare_peptides", "Compare Peptides", icon = icon("exchange-alt")),
              hr(),
              withSpinner(
                tabsetPanel(
                  id = "comparison_tabs",
                  type = "tabs",
                  tabPanel("Visualization", 
                    withSpinner(
                      plotlyOutput("peptide_comparison_plot"),
                      type = 8, color = "#3c8dbc"
                    )
                  ),
                  tabPanel("Peptide Table", 
                    h4("Detailed Peptide Data"),
                    p("Browse the full list of peptides for the selected transcript pair."),
                    DT::dataTableOutput("peptide_comparison_table")
                  ),
                  tabPanel("Summary", 
                    h4("Comparison Summary"),
                    p("A text summary highlighting the key differences in peptide content."),
                    verbatimTextOutput("peptide_comparison_summary")
                  )
                ),
                type = 8, color = "#3c8dbc"
              )
            )
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "All Events Peptide Comparison",
            status = "primary",
            solidHeader = TRUE,
            actionButton("load_all_events", "Load All Event Peptides", icon = icon("play"),
                         style = "color: #fff; background-color: #3c8dbc; border-color: #367fa9; margin-bottom: 10px;"),
            checkboxInput("show_event_labels", "Show Event ID Labels", value = TRUE, width = "100%"),
            conditionalPanel(
              condition = "input.load_all_events > 0",
              withSpinner(
                tabsetPanel(
                  id = "all_events_tabs",
                  type = "tabs",
                  tabPanel("Visualization", 
                    withSpinner(
                      plotlyOutput("as_all_peptide_comparison", height = "600px"),
                      type = 8, color = "#3c8dbc"
                    )
                  ),
                  tabPanel("Peptide Table", 
                    h4("All Events Peptide Summary"),
                    p("Detailed breakdown of all peptides and their AS event associations. Junction-spanning peptides show multiple coordinate ranges."),
                    DT::dataTableOutput("as_all_events_table"),
                    br(),
                    downloadButton("download_all_events_table", "Download Peptide Data", 
                                  class = "btn-primary")
                  )
                ),
                type = 8, color = "#3c8dbc"
              )
            )
          )
        )
      ),
      
      # Tab 3: Isoform-centric view
      conditionalPanel(
        condition = "input.canonical_current_tab == 'isoforms'",
        fluidRow(
          box(width = 12,
            title = "Isoform Analysis Overview",
            status = "primary",
            solidHeader = TRUE,
            div(
              style = "background-color: #f0f8ff; padding: 15px; border-left: 4px solid #3c8dbc; margin-bottom: 15px;",
              h4("What You Can Do Here", style = "margin-top: 0; color: #3c8dbc;"),
              p("Analyze which transcript isoforms you can detect in your mass spec experiments by comparing peptides across all transcript variants! See all isoforms at once to compare peptides across every transcript variant, categorize peptides as universal, shared, and unique, visualize gene structure with exon and CDS boundaries with peptide mapping, and find isoform-specific signatures by discovering which peptides are unique to specific isoforms.", style = "margin-bottom: 10px;")
            )
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "Multi-Isoform Comparison",
            status = "primary",
            solidHeader = TRUE,
            actionButton("load_all_isoforms", "Load All Isoforms", icon = icon("dna"),
                        style = "color: #fff; background-color: #3c8dbc; border-color: #367fa9; margin-bottom: 10px;")
          )
        ),
        conditionalPanel(
          condition = "input.load_all_isoforms > 0",
          fluidRow(
            box(
              title = "Isoform-Specific Peptide Analysis", 
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              div(
                style = "background-color: #f0f8ff; padding: 15px; border-left: 4px solid #3c8dbc; margin-bottom: 15px;",
                h4("Focus on One Isoform", style = "margin-top: 0; color: #3c8dbc;"),
                p("Select a specific isoform to highlight and analyze its unique, shared, and universal peptides. Focus on one isoform to highlight and analyze its peptides with color-coded specificity for unique, shared, and universal peptides, get detailed statistics with comprehensive specificity analysis, and download your analysis results.", style = "margin-bottom: 10px;")
              ),
              fluidRow(
                column(4,
                  selectInput("highlight_isoform", 
                            "Choose Isoform to Highlight:",
                            choices = NULL, 
                            width = "100%")
                ),
                column(3,
                  radioButtons("isoform_intron_scale", "Intron Display:",
                              choices = list(
                                "True scale" = "true",
                                "Compressed" = "compressed"
                              ),
                              selected = "compressed",
                              inline = FALSE)
                ),
                column(5,
                  div(style = "margin-top: 25px;",
                    downloadButton("download_isoform_analysis", 
                                 "Download Analysis Data",
                                 style = "color: #fff; background-color: #3c8dbc; border-color: #367fa9; width: 100%;")
                  )
                )
              ),
              hr(),
              h4("All Isoforms Visualization"),
              withSpinner(
                plotlyOutput("all_isoforms_plot", height = "600px"),
                type = 8, color = "#3c8dbc"
              )
            )
          ),
          fluidRow(
            column(12,
              box(
                title = "Highlighted Isoform Peptides", 
                status = "primary",
                solidHeader = TRUE,
                width = NULL,
                DT::dataTableOutput("highlighted_isoform_table")
              )
            )
          )
        ),
        conditionalPanel(
          condition = "input.load_all_isoforms > 0",
          fluidRow(
            box(
              title = "Multi-Isoform Comparative Analysis", 
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              collapsible = TRUE,
              div(
                style = "background-color: #f0f8ff; padding: 15px; border-left: 4px solid #3c8dbc; margin-bottom: 15px;",
                h4("Compare Multiple Isoforms", style = "margin-top: 0; color: #3c8dbc;"),
                p("Select 2 to 8 isoforms to visualize their structures side by side and analyze shared, unique, and universal peptides. Compare multiple isoforms by analyzing 2-8 transcript variants with side-by-side visualization to see structural differences, analyze shared peptides to identify universal, shared, and unique peptides, and get detailed overlap statistics with comparative summary.", style = "margin-bottom: 10px;")
              ),
              fluidRow(
                column(8,
                  selectizeInput("compare_isoforms", 
                               "Select Isoforms to Compare (2-8):",
                               choices = NULL,
                               multiple = TRUE,
                               options = list(
                                 maxItems = 8,
                                 placeholder = "Select at least 2 isoforms to compare..."
                               ),
                               width = "100%")
                ),
                column(4,
                  br(),
                  actionButton("run_comparative_analysis", 
                              "Compare Selected Isoforms",
                              icon = icon("chart-line"),
                              style = "color: #fff; background-color: #3c8dbc; border-color: #367fa9; width: 100%;")
                )
              ),
              conditionalPanel(
                condition = "input.run_comparative_analysis > 0",
                hr(),
                h4("Comparative Visualization"),
                withSpinner(
                  plotlyOutput("multi_isoform_comparative_plot", height = "600px"),
                  type = 8, color = "#3c8dbc"
                )
              )
            )
          ),
          
          # Single Isoform Multi-Enzyme Coverage Analysis
          fluidRow(
            box(
              title = "Single Isoform Multi-Enzyme Coverage Analysis",
              status = "info",
              solidHeader = TRUE,
              width = 12,
              collapsible = TRUE,
              div(
                style = "background-color: #f0fff0; padding: 15px; border-left: 4px solid #28a745; margin-bottom: 15px;",
                h4("Multi-Enzyme Peptide Coverage Comparison", style = "margin-top: 0; color: #28a745;"),
                p("Compare peptide coverage across multiple enzymes for a single isoform. Identify the best enzyme combination for detecting your target isoform! Focus on one isoform by selecting a specific transcript variant for detailed analysis, see coverage differences across up to 6 proteases with multi-enzyme comparison, compare peptide maps and coverage statistics with side-by-side visualization, and optimize your experiments by finding the best enzymes for your proteomics workflow.", style = "margin-bottom: 10px;")
              ),
              
              # Control Panel
              fluidRow(
                column(4,
                  selectInput("enzyme_coverage_isoform", 
                            "Select Isoform for Coverage Analysis:",
                            choices = NULL,
                            width = "100%")
                ),
                column(3,
                  selectInput("enzyme_miscleavage_global", 
                            "Miscleavage Setting:",
                            choices = c("No miscleavage" = "no_miss_cleavage",
                                      "Up to 2 miscleavages" = "upto_two_misscleavage"),
                            selected = "no_miss_cleavage",
                            width = "100%")
                ),
                column(3,
                  br(),
                  actionButton("run_enzyme_coverage", 
                             "Generate Coverage Analysis",
                             icon = icon("chart-line"),
                             style = "color: #fff; background-color: #28a745; border-color: #1e7e34; width: 100%; margin-top: 5px;")
                ),
                column(2,
                  br(),
                  downloadButton("download_enzyme_coverage", 
                               "Download Results",
                               style = "width: 100%; margin-top: 5px;")
                )
              ),
              
              hr(),
              
              # Enzyme Selection Grid
              h5("Select Enzymes to Compare:", style = "margin-bottom: 10px;"),
              fluidRow(
                column(2,
                  div(style = "text-align: center; padding: 8px; border: 1px solid #ddd; border-radius: 4px; background-color: #f8f9fa;",
                    checkboxInput("use_trp", strong("Trypsin"), value = TRUE),
                    helpText("Cuts after K, R", style = "font-size: 10px; margin: 0;")
                  )
                ),
                column(2,
                  div(style = "text-align: center; padding: 8px; border: 1px solid #ddd; border-radius: 4px; background-color: #f8f9fa;",
                    checkboxInput("use_chymo", strong("Chymotrypsin"), value = TRUE),
                    helpText("Cuts after F, W, Y, L", style = "font-size: 10px; margin: 0;")
                  )
                ),
                column(2,
                  div(style = "text-align: center; padding: 8px; border: 1px solid #ddd; border-radius: 4px; background-color: #f8f9fa;",
                    checkboxInput("use_aspn", strong("Asp-N"), value = FALSE),
                    helpText("Cuts before D", style = "font-size: 10px; margin: 0;")
                  )
                ),
                column(2,
                  div(style = "text-align: center; padding: 8px; border: 1px solid #ddd; border-radius: 4px; background-color: #f8f9fa;",
                    checkboxInput("use_lysc", strong("Lys-C"), value = FALSE),
                    helpText("Cuts after K", style = "font-size: 10px; margin: 0;")
                  )
                ),
                column(2,
                  div(style = "text-align: center; padding: 8px; border: 1px solid #ddd; border-radius: 4px; background-color: #f8f9fa;",
                    checkboxInput("use_lysn", strong("Lys-N"), value = FALSE),
                    helpText("Cuts before K", style = "font-size: 10px; margin: 0;")
                  )
                ),
                column(2,
                  div(style = "text-align: center; padding: 8px; border: 1px solid #ddd; border-radius: 4px; background-color: #f8f9fa;",
                    checkboxInput("use_gluc", strong("Glu-C"), value = FALSE),
                    helpText("Cuts after E", style = "font-size: 10px; margin: 0;")
                  )
                )
              ),
              
              # Results Section
              conditionalPanel(
                condition = "input.run_enzyme_coverage > 0",
                hr(),
                h4("Multi-Enzyme Coverage Results"),
                withSpinner(
                  plotlyOutput("enzyme_coverage_plot", height = "600px"),
                  type = 8, color = "#28a745"
                ),
                br(),
                fluidRow(
                  column(12,
                    box(
                      title = "Coverage Statistics",
                      status = "success",
                      solidHeader = TRUE,
                      width = NULL,
                      DT::dataTableOutput("enzyme_coverage_stats_table")
                    )
                  )
                )
              )
            )
          ),
          
          # Advanced Multi-Isoform Multi-Enzyme Matrix Analysis
          fluidRow(
            box(
              title = "Advanced Multi-Isoform Multi-Enzyme Matrix Analysis",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              collapsible = TRUE,
              div(
                style = "background-color: #fff3cd; padding: 15px; border-left: 4px solid #ffc107; margin-bottom: 15px;",
                h4("Matrix Analysis: Multiple Isoforms × Multiple Enzymes", style = "margin-top: 0; color: #856404;"),
                p("Combine the power of multi-isoform and multi-enzyme analysis! Generate comprehensive peptide coverage matrices.", style = "margin-bottom: 10px;"),
                p("Select 2-8 isoforms AND 2-6 enzymes for comprehensive combinations, see all isoform×enzyme combinations with matrix visualization, use advanced filtering for universal, shared, and unique peptides, and get comprehensive analysis by dimension with coverage statistics.", style = "margin-bottom: 10px;")
              ),
              
              # Control Panel
              fluidRow(
                column(6,
                  selectizeInput("matrix_isoforms", 
                                "Select Isoforms for Matrix Analysis (2-8):",
                                choices = NULL,
                                multiple = TRUE,
                                options = list(
                                  maxItems = 8,
                                  placeholder = "Select 2-8 isoforms for matrix analysis..."
                                ),
                                width = "100%")
                ),
                column(6,
                  div(
                    h5("Select Enzymes for Matrix Analysis (2-6):", style = "margin-bottom: 10px;"),
                    fluidRow(
                      column(4,
                        div(style = "text-align: center; padding: 8px; border: 1px solid #ddd; border-radius: 4px; background-color: #f8f9fa;",
                          checkboxInput("matrix_trp", strong("Trypsin"), value = TRUE),
                          helpText("Cuts after K, R", style = "font-size: 10px; margin: 0;")
                        )
                      ),
                      column(4,
                        div(style = "text-align: center; padding: 8px; border: 1px solid #ddd; border-radius: 4px; background-color: #f8f9fa;",
                          checkboxInput("matrix_chymo", strong("Chymotrypsin"), value = TRUE),
                          helpText("Cuts after F, W, Y, L", style = "font-size: 10px; margin: 0;")
                        )
                      ),
                      column(4,
                        div(style = "text-align: center; padding: 8px; border: 1px solid #ddd; border-radius: 4px; background-color: #f8f9fa;",
                          checkboxInput("matrix_aspn", strong("Asp-N"), value = FALSE),
                          helpText("Cuts before D", style = "font-size: 10px; margin: 0;")
                        )
                      )
                    ),
                    fluidRow(
                      column(4,
                        div(style = "text-align: center; padding: 8px; border: 1px solid #ddd; border-radius: 4px; background-color: #f8f9fa;",
                          checkboxInput("matrix_lysc", strong("Lys-C"), value = FALSE),
                          helpText("Cuts after K", style = "font-size: 10px; margin: 0;")
                        )
                      ),
                      column(4,
                        div(style = "text-align: center; padding: 8px; border: 1px solid #ddd; border-radius: 4px; background-color: #f8f9fa;",
                          checkboxInput("matrix_lysn", strong("Lys-N"), value = FALSE),
                          helpText("Cuts before K", style = "font-size: 10px; margin: 0;")
                        )
                      ),
                      column(4,
                        div(style = "text-align: center; padding: 8px; border: 1px solid #ddd; border-radius: 4px; background-color: #f8f9fa;",
                          checkboxInput("matrix_gluc", strong("Glu-C"), value = FALSE),
                          helpText("Cuts after E", style = "font-size: 10px; margin: 0;")
                        )
                      )
                    )
                  )
                )
              ),
              
              # Analysis Options
              fluidRow(
                column(3,
                  selectInput("matrix_miscleavage", "Miscleavage Setting:",
                             choices = c("No miscleavage" = "no_miss_cleavage",
                                       "Up to 2 miscleavages" = "upto_two_misscleavage"),
                             selected = "no_miss_cleavage",
                             width = "100%")
                ),
                column(3,
                  div(
                    h5("Genomic Visualization", style = "margin-bottom: 5px; color: #856404;"),
                    p("Shows exon/CDS structure with peptide overlays", style = "font-size: 12px; color: #666; margin: 0;")
                  )
                ),
                column(3,
                  br(),
                  checkboxInput("matrix_show_statistics", "Show Statistics Panel", value = TRUE)
                ),
                column(3,
                  br(),
                  actionButton("run_matrix_analysis", "Run Analysis",
                              icon = icon("play"),
                              style = "width: 100%; color: #fff; background-color: #ffc107; border-color: #ffc107; margin-top: 5px;")
                )
              ),
              
              # Results Section
              conditionalPanel(
                condition = "input.run_matrix_analysis > 0",
                hr(),
                h4("Multi-Dimensional Analysis Results"),
                
                # Main Visualization
                withSpinner(
                  plotlyOutput("matrix_analysis_plot", height = "750px"),
                  type = 8, color = "#ffc107"
                ),
                
                br(),
                
                # Statistics and Tables
                conditionalPanel(
                  condition = "input.matrix_show_statistics",
                  fluidRow(
                    column(6,
                      box(
                        title = "Coverage Matrix",
                        status = "warning",
                        solidHeader = TRUE,
                        width = NULL,
                        DT::dataTableOutput("coverage_matrix_table")
                      )
                    ),
                    column(6,
                      box(
                        title = "Peptide Specificity Analysis",
                        status = "warning",
                        solidHeader = TRUE,
                        width = NULL,
                        DT::dataTableOutput("specificity_matrix_table")
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
      )  
  ), # End Canonical Analysis View conditionalPanel
  
  # Novel Isoform Discovery - Complete Implementation
  conditionalPanel(
    condition = "input.current_view == 'novel'",
    div(class = "analysis-container",
      div(class = "view-header",
        h2("Novel Isoform Discovery", class = "view-title"),
        actionButton("back_to_home_novel", "← Back to Home", class = "back-btn", 
                    onclick = "Shiny.setInputValue('current_view', 'home');")
      ),
      
      # Overview Section
      fluidRow(
        box(width = 12,
          title = "Novel Isoform Discovery Overview",
          status = "primary",
          solidHeader = TRUE,
          div(
            style = "background-color: #f0f8ff; padding: 10px; border-left: 4px solid #3c8dbc; margin-bottom: 10px;",
            h4("🚀 How to Use This Tool", style = "margin-top: 0; color: #3c8dbc;"),
            p("Analyze your novel transcript sequences in 3 simple steps:", style = "margin-bottom: 8px;"),
            div(style = "display: flex; justify-content: space-between; margin-bottom: 10px;",
              div(style = "flex: 1; text-align: center; padding: 8px; background: #fff; margin: 2px; border-radius: 4px;",
                div(style = "font-size: 20px; color: #3c8dbc;", "📤"),
                div(style = "font-weight: bold; font-size: 12px;", "1. UPLOAD"),
                div(style = "font-size: 11px; color: #666;", "FASTA nucleotide sequence")
              ),
              div(style = "flex: 1; text-align: center; padding: 8px; background: #fff; margin: 2px; border-radius: 4px;",
                div(style = "font-size: 20px; color: #3c8dbc;", "🔍"),
                div(style = "font-weight: bold; font-size: 12px;", "2. ANALYZE"),
                div(style = "font-size: 11px; color: #666;", "Find ORFs & proteins")
              ),
              div(style = "flex: 1; text-align: center; padding: 8px; background: #fff; margin: 2px; border-radius: 4px;",
                div(style = "font-size: 20px; color: #3c8dbc;", "📊"),
                div(style = "font-weight: bold; font-size: 12px;", "3. COMPARE"),
                div(style = "font-size: 11px; color: #666;", "Against known genes")
              )
            ),
            p("Upload RNA-seq transcripts, long-read sequences, or predicted isoforms to discover novel proteins and compare them with known gene variants.", style = "font-size: 13px; margin: 5px 0; color: #666;")
          )
        )
      ),
      
      # File Upload & Analysis Pipeline
      fluidRow(
        box(
          width = 12,
          title = "Novel Transcript Analysis Pipeline",
          status = "primary",
          solidHeader = TRUE,
          h5("Step 1: Upload FASTA File"),
          fileInput("novel_fasta_file", "Choose FASTA File",
                   accept = c(".fa", ".fasta", ".fas", ".txt"),
                   width = "100%"),
          helpText("Upload nucleotide sequences (DNA/RNA). Min. 300 nt recommended."),
          h5("Step 2: Analysis Parameters"),
          fluidRow(
            column(3,
              numericInput("min_protein_length", "Min Protein Length (AA):",
                         value = 30, min = 10, max = 500, step = 5, width = "100%")
            ),
            column(3,
              selectInput("genetic_code", "Genetic Code:",
                         choices = list("Standard (1)" = 1, "Vertebrate Mito (2)" = 2, 
                                      "Yeast Mito (3)" = 3, "Other" = 4),
                         selected = 1, width = "100%")
            ),
            column(3,
              selectInput("strand_specific", "Strand:",
                         choices = list("Both (+/-)" = "both", "Plus (+)" = "plus", "Minus (-)" = "minus"),
                         selected = "both", width = "100%")
            ),
            column(3,
              selectInput("min_overlap_percent", "Gene Overlap Threshold:",
                         choices = list("High (≥50%)" = 50, "Moderate (≥20%)" = 20, "Low (≥10%)" = 10),
                         selected = 10, width = "100%")
            )
          ),
          fluidRow(
            column(4,
              checkboxInput("require_start_codon", "Require Start Codon", value = FALSE)
            ),
            column(4,
              checkboxInput("single_best_orf", "Best ORF Only", value = FALSE)
            ),
            column(4,
              checkboxInput("enable_gene_search", "Enable Gene Search", value = TRUE)
            )
          ),
          h5("Step 3: Execute Analysis"),
          actionButton("run_novel_pipeline", "Start Analysis", 
                      icon = icon("dna"),
                      style = "color: #fff; background-color: #3c8dbc; border-color: #367fa9; padding: 8px 16px; margin-bottom: 10px;"),
          conditionalPanel(
            condition = "input.run_novel_pipeline > 0",
            withSpinner(
              verbatimTextOutput("novel_pipeline_status"),
              type = 8, color = "#3c8dbc"
            )
          )
        )
      ),
      
      # ORF Selection Interface
      conditionalPanel(
        condition = "output.novel_pipeline_completed",
        fluidRow(
          box(
            width = 12,
            title = "Step 2: Select ORF for Analysis",
            status = "primary",
            solidHeader = TRUE,
            helpText("Select the best ORF (Open Reading Frame) for analysis:"),
            DT::dataTableOutput("novel_orf_table"),
            br(),
            actionButton("run_gene_search", "Search Candidate Genes", 
                        icon = icon("search"),
                        style = "color: #fff; background-color: #3c8dbc; border-color: #367fa9; margin-right: 10px;"),
            verbatimTextOutput("novel_orf_selection_status")
          )
        )
      ),
      
      # Gene Selection Interface
      conditionalPanel(
        condition = "output.novel_gene_search_completed",
        fluidRow(
          box(
            width = 12,
            title = "Step 3: Select Gene for Comparison",
            status = "primary",
            solidHeader = TRUE,
            helpText("Select a gene to compare your novel isoform against:"),
            DT::dataTableOutput("novel_gene_search_results_table"),
            br(),
            actionButton("load_selected_gene", "Load Gene for Comparison", 
                        icon = icon("download"),
                        style = "color: #fff; background-color: #3c8dbc; border-color: #367fa9; margin-right: 10px;"),
            verbatimTextOutput("novel_gene_search_selection_status")
          )
        )
      ),
      
      # Novel Isoform Analysis Interface
      conditionalPanel(
        condition = "output.novel_analysis_ready",
        fluidRow(
          box(
            width = 12,
            title = "Novel Isoform Analysis",
            status = "primary",
            solidHeader = TRUE,
            fluidRow(
              column(4, 
                selectInput("novel_highlight_isoform", "Highlight:", 
                           choices = NULL, width = "100%")
              ),
              column(3,
                selectInput("novel_protease", "Enzyme:", 
                           choices = c("Trypsin" = "trp", "Chymo" = "chymo", "AspN" = "aspn", 
                                     "LysC" = "lysc", "LysN" = "lysn", "GluC" = "gluc"),
                           selected = "trp", width = "100%")
              ),
              column(3,
                selectInput("novel_miscleavage_type", "Miscleavage:",
                           choices = c("None" = "no_miss_cleavage", "Up to 2" = "upto_two_misscleavage"),
                           selected = "no_miss_cleavage", width = "100%")
              ),
              column(2,
                actionButton("load_novel_isoforms", "Load Analysis", 
                            icon = icon("play"),
                            style = "color: #fff; background-color: #3c8dbc; border-color: #367fa9; width: 100%; margin-top: 25px;")
              )
            )
          )
        ),
        
        # Comparative Visualization
        conditionalPanel(
          condition = "input.load_novel_isoforms > 0",
          fluidRow(
            box(
              title = "Multi-Isoform Comparative Analysis", 
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              helpText("Compare your novel sequence with known isoforms:"),
              fluidRow(
                column(8,
                  selectizeInput("novel_compare_isoforms", 
                               "Select Isoforms to Compare:",
                               choices = NULL,
                               multiple = TRUE,
                               options = list(maxItems = 8, placeholder = "Select 2-8 isoforms..."),
                               width = "100%")
                ),
                column(4,
                  actionButton("run_novel_comparative_analysis", 
                              "Compare Isoforms",
                              icon = icon("chart-line"),
                              style = "color: #fff; background-color: #3c8dbc; border-color: #367fa9; width: 100%; margin-top: 25px;")
                )
              ),
              conditionalPanel(
                condition = "input.run_novel_comparative_analysis > 0",
                hr(),
                withSpinner(
                  plotlyOutput("novel_comparative_plot", height = "600px"),
                  type = 8, color = "#3c8dbc"
                )
              )
            )
          )
        ),
        
        # Peptide Results Table
        conditionalPanel(
          condition = "input.run_novel_comparative_analysis > 0",
          fluidRow(
            column(12,
              box(
                title = "Highlighted Isoform Peptides", 
                status = "primary",
                solidHeader = TRUE,
                width = NULL,
                DT::dataTableOutput("novel_highlighted_isoform_table")
              )
            )
          )
        ),
        
        # Export Section
        fluidRow(
          box(
            width = 12,
            title = "Export Novel Isoform Results",
            status = "primary",
            solidHeader = TRUE,
            conditionalPanel(
              condition = "input.load_novel_isoforms > 0",
              fluidRow(
                column(3,
                  downloadButton("download_novel_dataframe", "Download Novel Dataframe (RDS)",
                               style = "width: 100%; margin-bottom: 10px;")
                ),
                column(3,
                  downloadButton("download_novel_gtf", "Download GTF with Phases",
                               style = "width: 100%; margin-bottom: 10px;")
                ),
                column(3,
                  downloadButton("download_novel_peptides", "Download Peptide Analysis",
                               style = "width: 100%; margin-bottom: 10px;")
                ),
                column(3,
                  downloadButton("download_novel_pipeline_log", "Download Pipeline Log",
                               style = "width: 100%; margin-bottom: 10px;")
                )
              )
            )
          )
        )
      )
    )
  ), # End Novel Isoform Discovery conditionalPanel
  
  # Peptide Search - Complete BLASTP Implementation with Module Jumping
  conditionalPanel(
    condition = "input.current_view == 'peptide'",
    div(class = "analysis-container",
      div(class = "view-header",
        h2("Peptide Sequence Search", class = "view-title"),
        actionButton("back_to_home_peptide", "← Back to Home", class = "back-btn", 
                    onclick = "Shiny.setInputValue('current_view', 'home');")
      ),
      
      # Peptide Search Overview
      fluidRow(
        box(width = 12,
          title = "Peptide Search Overview",
          status = "primary",
          solidHeader = TRUE,
          div(
            style = "background-color: #f0f8ff; padding: 15px; border-left: 4px solid #3c8dbc; margin-bottom: 15px;",
            h4("What You Can Do Here", style = "margin-top: 0; color: #3c8dbc;"),
            p("Got a specific peptide sequence? Use BLASTP to search through 20,000+ genes and find homologous matches with detailed statistics! This tool is BLASTP powered with industry-standard protein homology search, provides rich statistics including E-values, bit scores, and identity percentages, offers flexible matching to find exact matches and similar sequences, has customizable parameters to adjust sensitivity and specificity, and lets you jump right into full gene analysis by clicking results.", style = "margin-bottom: 10px;")
          )
        )
      ),
      
      # BLASTP Search Interface
      fluidRow(
        box(
          title = "Search Peptides Across All Genes & Isoforms",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          fluidRow(
            column(12,
              textInput("peptide_search_query", "Peptide Sequence:", 
                       placeholder = "Enter peptide sequence (e.g., MTEYKLVVVGAGGVGK)",
                       width = "100%")
            )
          ),
          h5("BLASTP Search Parameters"),
          fluidRow(
            column(3,
              numericInput("peptide_search_evalue", "E-value Threshold:",
                         value = 0.01, min = 0.0001, max = 10, step = 0.001,
                         width = "100%"),
              helpText("Lower values = more stringent (0.01 recommended)")
            ),
            column(3,
              numericInput("peptide_search_identity", "Min Identity %:",
                         value = 50, min = 1, max = 100, step = 1,
                         width = "100%"),
              helpText("Minimum sequence identity (50% for discovery)")
            ),
            column(3,
              numericInput("peptide_search_max_targets", "Max Targets:",
                         value = 500, min = 10, max = 5000, step = 50,
                         width = "100%"),
              helpText("Maximum number of hits")
            ),
            column(3,
              br(),
              actionButton("run_peptide_search", "Search with BLASTP",
                         icon = icon("search"),
                         style = "color: #fff; background-color: #3c8dbc; border-color: #367fa9; width: 100%; margin-top: 5px;")
            )
          )
        )
      ),
      
      # Search Results Display
      fluidRow(
        conditionalPanel(
          condition = "output.peptide_search_results_available",
          box(
            title = "Search Results",
            status = "primary", 
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(6,
                h4(textOutput("peptide_search_summary"))
              ),
              column(6,
                downloadButton("download_peptide_search", "Download Results",
                             style = "float: right; margin-top: 10px;")
              )
            ),
            br(),
            DT::dataTableOutput("peptide_search_results_table")
          )
        )
      ),
      
      # CRITICAL: Module Jumping Navigation
      fluidRow(
        conditionalPanel(
          condition = "output.peptide_search_results_available && input.peptide_search_results_table_rows_selected.length > 0",
          box(
            title = "Navigate to Gene Analysis", 
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            p("Click a result above, then navigate to detailed isoform analysis:"),
            fluidRow(
              column(12,
                actionButton("goto_isoform_analysis", "Go to Isoform Analysis",
                           icon = icon("project-diagram"),
                           onclick = "Shiny.setInputValue('current_view', 'canonical');",
                           style = "width: 100%; color: #fff; background-color: #3c8dbc; border-color: #367fa9;")
              )
            )
          )
        )
      ),
      
      # Perfect Match Info Panel
      fluidRow(
        conditionalPanel(
          condition = "output.blast_perfect_match_selected",
          box(
            title = "100% Perfect Match Found", 
            status = "success",
            solidHeader = TRUE,
            width = 12,
            div(
              style = "background-color: #d4edda; padding: 15px; margin-bottom: 15px; border-left: 4px solid #28a745; border-radius: 5px;",
              h5("🎯 100% Perfect Match Selected", style = "margin-top: 0; color: #155724;"),
              textOutput("blast_selected_info"),
              br(),
              actionButton("show_blast_visualization", 
                         "Show Transcript Visualization",
                         icon = icon("dna"),
                         style = "color: #fff; background-color: #28a745; border-color: #1e7e34; font-weight: bold;")
            )
          )
        )
      ),
      
      # Advanced Transcript Visualization
      fluidRow(
        conditionalPanel(
          condition = "output.blast_visualization_visible",
          box(
            title = "Perfect Match Transcript Visualization", 
            status = "success",
            solidHeader = TRUE,
            width = 12,
            div(
              style = "background-color: #f8f9fa; padding: 10px; margin-bottom: 15px; border-radius: 5px;",
              actionButton("hide_blast_visualization", "Hide Visualization", 
                         icon = icon("eye-slash"),
                         style = "color: #6c757d; background-color: #e9ecef; border-color: #ced4da; float: right;")
            ),
            withSpinner(
              plotlyOutput("blast_transcript_plot", height = "750px"),
              type = 8, color = "#28a745"
            ),
            br(),
            downloadButton("download_blast_plot", "Download Transcript Plot",
                         icon = icon("download"),
                         style = "color: #fff; background-color: #28a745; border-color: #1e7e34;")
          )
        )
      )
    )
  ), # End Peptide Search conditionalPanel
  
  # Alternative Splicing Analysis - rMATS Implementation
  conditionalPanel(
    condition = "input.current_view == 'alt_splicing'",
    div(class = "analysis-container",
      div(class = "view-header",
        h2("Alternative Splicing Analysis", class = "view-title"),
        actionButton("back_to_home_alt_splicing", "← Back to Home", class = "back-btn", 
                    onclick = "Shiny.setInputValue('current_view', 'home');")
      ),
      
      # Tab Navigation Bar
      div(style = "margin: 20px 0;",
        div(style = "display: flex; background: #f4f4f4; border-radius: 6px; padding: 4px;",
          div(
            onclick = "Shiny.setInputValue('alt_splicing_current_tab', 'rmats');",
            style = "flex: 1; text-align: center; padding: 12px; cursor: pointer; border-radius: 4px; font-weight: 500; transition: all 0.3s ease;",
            class = "tab-btn tab-rmats",
            "rMATS Analysis"
          ),
          div(
            onclick = "Shiny.setInputValue('alt_splicing_current_tab', 'spladder');",
            style = "flex: 1; text-align: center; padding: 12px; cursor: pointer; border-radius: 4px; font-weight: 500; transition: all 0.3s ease;",
            class = "tab-btn tab-spladder",
            "SplAdder Analysis"
          )
        )
      ),
      
      # Tab 1: rMATS Analysis
      conditionalPanel(
        condition = "input.alt_splicing_current_tab == 'rmats'",
      fluidRow(
        box(width = 12,
          title = "rMATS Peptide Analysis Overview",
          status = "primary",
          solidHeader = TRUE,
          div(
            style = "background-color: #f0f8ff; padding: 15px; border-left: 4px solid #3c8dbc; margin-bottom: 15px;",
            h4("rMATS Functional Analysis", style = "margin-top: 0; color: #3c8dbc;"),
            p("Analyze alternative splicing events from rMATS output files. Supports SE, A3SS, A5SS, MXE, and RI event types with protein translation and functional prediction.", style = "margin-bottom: 10px;")
          )
        )
      ),
      
      # File Upload Section
      fluidRow(
        box(
          title = "Upload rMATS File", 
          status = "info", 
          solidHeader = TRUE,
          width = 12,
          
          fluidRow(
            column(4,
              h5("Select Event Type:"),
              selectInput("rmats_event_type", NULL,
                choices = list(
                  "Select Event Type" = "",
                  "SE (Skipped Exon)" = "SE",
                  "A3SS (Alternative 3' Splice Site)" = "A3SS", 
                  "A5SS (Alternative 5' Splice Site)" = "A5SS",
                  "MXE (Mutually Exclusive Exons)" = "MXE",
                  "RI (Retained Intron)" = "RI"
                ),
                selected = "")
              
            ),
            
            column(4,
              conditionalPanel(
                condition = "input.rmats_event_type != ''",
                h5("Upload File:"),
                fileInput("rmats_file", NULL,
                         accept = c(".txt", ".JC.txt"))
              )
            ),
            
            column(4,
              conditionalPanel(
                condition = "output.rmats_file_uploaded && input.rmats_event_type != ''",
                div(class = "alert alert-success",
                    icon("check"), " File uploaded successfully!"),
                verbatimTextOutput("rmats_file_info")
              )
            )
          ),
          
          # Expected columns info
          conditionalPanel(
            condition = "input.rmats_event_type != ''",
            hr(),
            h5("Expected Columns for ", textOutput("rmats_event_type_display2", inline = TRUE), ":"),
            verbatimTextOutput("rmats_expected_columns")
          )
        )
      ),
      
      # 8-Step Analysis Pipeline
      fluidRow(
        box(
          title = "Analysis Pipeline", 
          status = "primary", 
          solidHeader = TRUE,
          width = 12,
          
          conditionalPanel(
            condition = "output.rmats_file_uploaded",
            
            # Step 1: Parse rMATS Events
            div(
              style = "background-color: #f8f9fa; padding: 10px; margin: 10px 0; border-radius: 5px;",
              actionButton("rmats_parse_events", "Parse rMATS Events", 
                class = "btn-primary", icon = icon("cogs")),
              br(), br(),
              conditionalPanel(
                condition = "output.rmats_step1_completed",
                div(class = "alert alert-success", icon("check"), " Step 1 completed successfully!"),
                verbatimTextOutput("rmats_step1_results")
              )
            ),
            
            # Event Selection Table
            conditionalPanel(
              condition = "output.rmats_step1_completed",
              div(
                style = "background-color: #fff8dc; padding: 15px; margin: 15px 0; border-radius: 5px; border-left: 4px solid #f39c12;",
                
                fluidRow(
                  column(12,
                    DT::dataTableOutput("rmats_events_table"),
                    br(),
                    conditionalPanel(
                      condition = "input.rmats_events_table_rows_selected.length > 0",
                      verbatimTextOutput("rmats_selected_event_info")
                    )
                  )
                )
              )
            ),
            
            # Processing Status (minimal)
            conditionalPanel(
              condition = "output.rmats_step1_completed && input.rmats_events_table_rows_selected.length > 0",
              conditionalPanel(
                condition = "output.rmats_auto_processing",
                withSpinner(
                  span(""),
                  type = 8, color = "#667eea", size = 0.5
                )
              )
            )
          )
        )
      ),
      
      # Comprehensive Analysis Section
      fluidRow(
        box(
          width = 12,
          title = "Generate Isoforms",
          status = "success",
          solidHeader = TRUE,
          collapsible = TRUE,
          
          # Run Analysis
          fluidRow(
            column(12,
              actionButton("run_rmats_comprehensive_analysis", 
                         "Generate Inclusion/Exclusion Isoforms",
                         icon = icon("dna"),
                         style = "color: #fff; background-color: #667eea; border-color: #667eea; margin-bottom: 10px;"),
              verbatimTextOutput("rmats_comprehensive_status")
            )
          )
        )
      ),
      
      # Multi-Isoform Visualization (conditional on analysis completion)
      conditionalPanel(
        condition = "output.rmats_comprehensive_completed",
        fluidRow(
          box(
            width = 12,
            title = "Multi-Isoform Comparative Visualization",
            status = "success", 
            solidHeader = TRUE,
            
            # Analysis controls
            fluidRow(
              column(3,
                selectInput("rmats_highlight_isoform", "Highlight:", 
                           choices = NULL, width = "100%")
              ),
              column(3,
                selectInput("rmats_protease", "Enzyme:", 
                           choices = c("Trypsin" = "trp", "Chymo" = "chymo", "AspN" = "aspn", 
                                     "LysC" = "lysc", "LysN" = "lysn", "GluC" = "gluc"),
                           selected = "trp", width = "100%")
              ),
              column(3,
                selectInput("rmats_miscleavage_type", "Miscleavage:",
                           choices = c("None" = "no_miss_cleavage", "Up to 2" = "upto_two_misscleavage"),
                           selected = "no_miss_cleavage", width = "100%")
              ),
              column(3,
                actionButton("load_rmats_gene", "Load Gene Data", 
                            icon = icon("play"),
                            style = "color: #fff; background-color: #28a745; border-color: #28a745; width: 100%; margin-top: 25px;")
              )
            )
          )
        ),
        
        # Multi-isoform comparison section
        conditionalPanel(
          condition = "input.load_rmats_gene > 0",
          fluidRow(
            box(
              title = "Multi-Isoform Comparative Analysis", 
              status = "success",
              solidHeader = TRUE,
              width = 12,
              helpText("Compare rMATS inclusion/exclusion isoforms with known gene transcripts:"),
              fluidRow(
                column(8,
                  selectizeInput("rmats_compare_isoforms", 
                               "Select Isoforms to Compare:",
                               choices = NULL,
                               multiple = TRUE,
                               options = list(maxItems = 8, placeholder = "Select 2-8 isoforms (include rMATS variants)..."),
                               width = "100%")
                ),
                column(4,
                  actionButton("run_rmats_comparative_analysis", 
                              "Compare Isoforms",
                              icon = icon("chart-line"),
                              style = "color: #fff; background-color: #28a745; border-color: #28a745; width: 100%; margin-top: 25px;")
                )
              ),
              conditionalPanel(
                condition = "input.run_rmats_comparative_analysis > 0",
                hr(),
                withSpinner(
                  plotlyOutput("rmats_comparative_plot", height = "600px"),
                  type = 8, color = "#28a745"
                )
              )
            )
          ),
          
          # Highlighted isoform table
          conditionalPanel(
            condition = "input.run_rmats_comparative_analysis > 0",
            fluidRow(
              column(12,
                box(
                  title = "Highlighted Isoform Peptides", 
                  status = "success",
                  solidHeader = TRUE,
                  width = NULL,
                  DT::dataTableOutput("rmats_highlighted_isoform_table")
                )
              )
            )
          )
        )
      ) # End rMATS Analysis
      ), # End rMATS conditionalPanel
      
      # Tab 2: SplAdder Analysis  
      conditionalPanel(
        condition = "input.alt_splicing_current_tab == 'spladder'",
      
      # SplAdder Analysis (copy of rMATS with modifications)
      fluidRow(
        box(width = 12,
          title = "SplAdder Peptide Analysis Overview",
          status = "primary",
          solidHeader = TRUE,
          div(
            style = "background-color: #f0f8ff; padding: 15px; border-left: 4px solid #3c8dbc; margin-bottom: 15px;",
            h4("SplAdder Functional Analysis", style = "margin-top: 0; color: #3c8dbc;"),
            p("Analyze alternative splicing events from SplAdder GFF3 output files. Supports SE, A3SS, A5SS, MXE, and RI event types with protein translation and functional prediction.", style = "margin-bottom: 10px;")
          )
        )
      ),
      
      # File Upload Section (modified for GFF3)
      fluidRow(
        box(
          title = "Upload SplAdder GFF3 File", 
          status = "info", 
          solidHeader = TRUE,
          width = 12,
          
          fluidRow(
            column(4,
              h5("Select Event Type:"),
              selectInput("spladder_event_type", NULL,
                choices = list(
                  "Select Event Type" = "",
                  "SE (Skipped Exon)" = "exon_skip",
                  "A3SS (Alternative 3' Splice Site)" = "alt_3prime", 
                  "A5SS (Alternative 5' Splice Site)" = "alt_5prime",
                  "MXE (Mutually Exclusive Exons)" = "mutex_exons",
                  "RI (Retained Intron)" = "intron_retention"
                ),
                selected = "")
              
            ),
            
            column(4,
              conditionalPanel(
                condition = "input.spladder_event_type != ''",
                h5("Upload GFF3 File:"),
                fileInput("spladder_file", NULL,
                  accept = c(".gff3", ".gff"),
                  width = "100%")
              )
            ),
            
            column(4,
              conditionalPanel(
                condition = "input.spladder_file",
                br(),
                div(class = "alert alert-success", 
                  icon("check-circle"), " GFF3 file uploaded!")
              )
            )
          ),
          
          conditionalPanel(
            condition = "input.spladder_event_type != ''",
            hr(),
            h5("Event Type: ", textOutput("spladder_event_type_display2", inline = TRUE)),
            verbatimTextOutput("spladder_expected_format")
          )
        )
      ),
      
      # 8-Step Analysis Pipeline (identical to rMATS)
      fluidRow(
        box(
          title = "Analysis Pipeline", 
          status = "primary", 
          solidHeader = TRUE,
          width = 12,
          
          conditionalPanel(
            condition = "output.spladder_file_uploaded",
            
            # Step 1: Parse SplAdder Events
            div(
              style = "background-color: #f8f9fa; padding: 10px; margin: 10px 0; border-radius: 5px;",
              actionButton("spladder_parse_events", "Parse SplAdder Events", 
                class = "btn-primary", icon = icon("cogs")),
              br(), br(),
              conditionalPanel(
                condition = "output.spladder_step1_completed",
                div(class = "alert alert-success", icon("check"), " Step 1 completed successfully!"),
                verbatimTextOutput("spladder_step1_results")
              )
            ),
            
            # Event Selection Table
            conditionalPanel(
              condition = "output.spladder_step1_completed",
              div(
                style = "background-color: #fff8dc; padding: 15px; margin: 15px 0; border-radius: 5px; border-left: 4px solid #f39c12;",
                
                fluidRow(
                  column(12,
                    DT::dataTableOutput("spladder_events_table"),
                    br(),
                    conditionalPanel(
                      condition = "input.spladder_events_table_rows_selected.length > 0",
                      verbatimTextOutput("spladder_selected_event_info")
                    )
                  )
                )
              )
            )
          )
        )
      ),
      
      # Comprehensive Analysis Section (identical to rMATS)
      fluidRow(
        box(
          width = 12,
          title = "Generate Isoforms",
          status = "success",
          solidHeader = TRUE,
          collapsible = TRUE,
          
          # Run Analysis
          fluidRow(
            column(12,
              actionButton("run_spladder_comprehensive_analysis", 
                         "Generate Inclusion/Exclusion Isoforms",
                         icon = icon("dna"),
                         style = "color: #fff; background-color: #667eea; border-color: #667eea; margin-bottom: 10px;"),
              verbatimTextOutput("spladder_comprehensive_status")
            )
          )
        )
      ),
      
      # Multi-Isoform Visualization (conditional on analysis completion)
      conditionalPanel(
        condition = "output.spladder_comprehensive_completed",
        fluidRow(
          box(
            width = 12,
            title = "Multi-Isoform Comparative Visualization",
            status = "success", 
            solidHeader = TRUE,
            
            # Analysis controls
            fluidRow(
              column(3,
                selectInput("spladder_highlight_isoform", "Highlight:", 
                           choices = NULL, width = "100%")
              ),
              column(3,
                selectInput("spladder_protease", "Enzyme:", 
                           choices = c("Trypsin" = "trp", "Chymo" = "chymo", "AspN" = "aspn", 
                                     "LysC" = "lysc", "LysN" = "lysn", "GluC" = "gluc"),
                           selected = "trp", width = "100%")
              ),
              column(3,
                selectInput("spladder_miscleavage_type", "Miscleavage:",
                           choices = c("None" = "no_miss_cleavage", "Up to 2" = "upto_two_misscleavage"),
                           selected = "no_miss_cleavage", width = "100%")
              ),
              column(3,
                actionButton("load_spladder_gene", "Load Gene Data", 
                            icon = icon("play"),
                            style = "color: #fff; background-color: #28a745; border-color: #28a745; width: 100%; margin-top: 25px;")
              )
            )
          )
        ),
        
        # Multi-isoform comparison section
        conditionalPanel(
          condition = "input.load_spladder_gene > 0",
          fluidRow(
            box(
              title = "Multi-Isoform Comparative Analysis", 
              status = "success",
              solidHeader = TRUE,
              width = 12,
              
              # Status indicator
              conditionalPanel(
                condition = "!output.spladder_step1_completed",
                div(
                  class = "alert alert-info",
                  HTML("<strong>Step 1:</strong> Please complete comprehensive analysis first by selecting an event and clicking 'Generate Inclusion/Exclusion Isoforms'")
                )
              ),
              
              # Multi-isoform selection (only show when step 1 is complete)
              conditionalPanel(
                condition = "output.spladder_step1_completed",
                fluidRow(
                  column(6,
                    selectizeInput("spladder_compare_isoforms", "Select Isoforms to Compare:",
                                  choices = NULL,
                                  multiple = TRUE,
                                  options = list(maxItems = 10),
                                  width = "100%")
                  ),
                  column(6,
                    actionButton("run_spladder_comparative_analysis", 
                                "Generate Comparative Analysis",
                                icon = icon("chart-line"),
                                style = "color: #fff; background-color: #28a745; border-color: #28a745; width: 100%; margin-top: 25px;")
                  )
                )
              ),
              
              # Results - USE SPLADDER VISUALIZATION
              conditionalPanel(
                condition = "input.run_spladder_comparative_analysis > 0",
                hr(),
                withSpinner(
                  plotlyOutput("spladder_comparative_plot", height = "600px"),
                  type = 8, color = "#28a745"
                )
              )
            )
          ),
          
          # Highlighted isoform table
          conditionalPanel(
            condition = "input.run_spladder_comparative_analysis > 0",
            fluidRow(
              column(12,
                box(
                  title = "Highlighted Isoform Peptides", 
                  status = "success",
                  solidHeader = TRUE,
                  width = NULL,
                  DT::dataTableOutput("spladder_highlighted_isoform_table")
                )
              )
            )
          )
        )
      )
      ) # End SplAdder conditionalPanel
      
    ) # End div analysis-container  
  ) # End Alternative Splicing Analysis conditionalPanel
) # Close fluidPage and ui assignment
