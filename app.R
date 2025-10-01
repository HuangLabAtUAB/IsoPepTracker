#===============================================================================
# IsoPepTracker - Shiny Application  
# Refactored into modular structure: global.R, ui.R, and server.R
#===============================================================================
# Load required libraries first
library(shiny)
library(shinydashboard)

cat("IsoPepTracker - Lightning-Fast Visualizations\n")
cat("===============================================\n")
cat("Starting app with automatic performance optimizations...\n\n")

source("global.R")  # Load global configurations and data
source("ui.R")      # Load UI definition  
source("server.R")

cat("\n🎉 App ready! Visualizations are lightning-fast!\n")
cat("Open your browser and start exploring genes!\n\n")

# Run the application
shinyApp(ui = ui, server = server)
