## ---- load packages ----
suppressPackageStartupMessages({
  # Core UI
  library(shiny)
  library(bs4Dash)
  
  # Data
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  
  # Plots / widgets
  library(ggplot2)
  library(ggiraph)
  library(ggupset)
  library(shinyWidgets)
  library(shinycssloaders)
  
  # UI helpers
  library(bslib)
  library(htmltools)
  library(base64enc)
  
  # BioC + GitHub
  library(Biostrings)
  library(ALJEbinf)
})


# ---------------- UI ----------------
ui <- dashboardPage(
  title = "RIFRMap",
  help  = NULL,
  
  header = dashboardHeader(
    title = dashboardBrand(
      title = "RIFRMap",
      color = "primary"
    ),
    skin           = "light",
    status         = "white",
    border         = TRUE,
    sidebarIcon    = icon("bars"),
    controlbarIcon = icon("cog")
  ),
  
  sidebar = dashboardSidebar(
    skin      = "light",
    status    = "primary",
    elevation = 3,
    fixed = TRUE,
    bs4SidebarMenu(
      id = "sidebar",
      bs4SidebarMenuItem("Home",    tabName = "home",    icon = icon("home")),
      bs4SidebarMenuItem("Explore", tabName = "explore", icon = icon("compass"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
      tags$script(src = "script.js")
    ),
    
    bs4TabItems(
      ## ---------------- Home ----------------
      bs4TabItem(
        tabName = "home",
        h2("Welcome to RIFRMap")
      ),
      
      ## ---------------- Explore ----------------
      bs4TabItem(
        tabName = "explore",
        
        tabsetPanel(
          id = "explore_steps",
          
          ## ---- Step 1: Select data ----
          tabPanel(
            title = "Select data",
            value = "filters",
            tags$p(
              class = "text-muted",
              HTML("Select data sources and filtering criteria to define the mutation set used in downstream analyses.")
            ),
      

            fluidRow(
              # CARD
              column(
                width = 12, lg = 3, xl = 2,
                bs4Card(
                  title = "CARD",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  checkboxInput("card_enabled", "Include data", FALSE),
                  selectInput(
                    "card_species", "Species",
                    choices = NULL, selected = "All species", multiple = TRUE
                  )
                )
              ),
              
              # CRyPTIC
              column(
                width = 12, lg = 3, xl = 2,
                bs4Card(
                  title = "CRyPTIC",
                  width = 12,
                  status = "success",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  checkboxInput("cryptic_enabled", "Include data", FALSE),
                  selectInput(
                    "cryptic_species", "Species",
                    choices = NULL, selected = "Mycobacterium tuberculosis", multiple = TRUE
                  ),
                  sliderInput(
                    "cryptic_coef", "Effect size",
                    min = -4.0, max = 9.6, value = c(-4, 9.6), step = 0.1
                  ),
                  sliderInput(
                    "cryptic_counts", "Counts",
                    min = 0, max = 2326, value = c(0, 2326), step = 1
                  ),
                  tags$div(
                    class = "cryptic-bhp-wrap",
                    radioButtons(
                      "cryptic_bhp",
                      "BH-corrected p-value (≤)",
                      choices  = c("0.05" = 0.05, "0.01" = 0.01, "0.005" = 0.005, "0.001" = 0.001),
                      selected = 0.05,
                      inline   = TRUE
                    )
                  ),
                  selectInput(
                    "cryptic_homoplasy", "Homoplasy",
                    choices = c("", "TRUE", "FALSE"),
                    selected = "",
                    width = "140px"
                  )
                )
              ),
              
              # DRAGdb
              column(
                width = 12, lg = 3, xl = 2,
                bs4Card(
                  title = "DRAGdb",
                  width = 12,
                  status = "warning",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  checkboxInput("dragdb_enabled", "Include data", FALSE),
                  selectInput(
                    "dragdb_species", "Species",
                    choices = NULL, selected = "All species", multiple = TRUE
                  ),
                  sliderInput(
                    "dragdb_provean", "PROVEAN",
                    min = -10.7, max = -2.6, value = c(-10.7, -2.6), step = 0.5
                  )
                )
              ),
              
              # Yang
              column(
                width = 12, lg = 3, xl = 2,
                bs4Card(
                  title = "Yang",
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  checkboxInput("yang_enabled", "Include data", FALSE),
                  selectInput(
                    "yang_species", "Species",
                    choices = NULL, selected = "Escherichia coli", multiple = TRUE
                  ),
                  sliderInput(
                    "yang_log2fc", "log2FC",
                    min = 0, max = 2.7, value = c(0, 2.7), step = 0.1
                  )
                )
              ),
              
              # Curated data
              column(
                width = 12, lg = 3, xl = 2,
                bs4Card(
                  title = "Curated data",
                  width = 12,
                  status = "secondary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  checkboxInput("curated_enabled", "Include data", FALSE),
                  selectInput(
                    "curated_species", "Species",
                    choices = NULL, selected = "All species", multiple = TRUE
                  ),
                  sliderInput(
                    "curated_nstudies", "Number of studies (per mutation)",
                    min = 1, max = 59, value = c(1, 59), step = 1
                  ),
                  sliderInput(
                    "curated_nspecies", "Number of species (per mutation)",
                    min = 1, max = 29, value = c(1, 29), step = 1
                  )
                )
              )
            ),
            
            br(),
            
            fluidRow(
              column(
                width = 12,
                actionButton("apply_filters", "Apply filters", icon = icon("play"), class = "btn-primary"),
                actionButton("reset_filters", "Reset", icon = icon("rotate-left"))
              )
            )
          ),
          
          ## ---- Step 2: Species patterns ----
          tabPanel(
            title = "Species patterns",
            value = "cluster",
            
            div(
              class = "step-intro",
              tags$p(
                class = "text-muted",
                HTML("Species-centric mutation cluster view. Sources are ignored; patterns reflect <b>species</b> and <b>cluster</b> structure.")
              )
            ),
            
            fluidRow(
              # Controls
              column(
                width = 2,
                card(
                  card_body(
                    div(
                      class = "clusters-buttons",
                      shinyWidgets::checkboxGroupButtons(
                        inputId   = "clustersToPlot",
                        label     = "Clusters",
                        choices   = c("N", "I", "II", "III"),
                        selected  = c("N", "I", "II", "III"),
                        direction = "horizontal",
                        justified = TRUE
                      )
                    ),
                    uiOutput("legend_panel")
                  )
                )
              ),
              
              # Plot
              column(
                width = 10,
                card(
                  width = 12,
                  collapsible = FALSE,
                  solidHeader = FALSE,
                  style = "height: 90vh;",
                  uiOutput("clusterPlotUI")
                )
              )
            )
          ),
          
          ## ---- Step 3: Evidence support ----
          tabPanel(
            title = "Sources Evidence",
            value = "sources",
            
            div(
              class = "step-intro",
              tags$p(
                class = "text-muted",
                HTML("Source-centric view of mutation evidence. Plots highlight overlap among <b>databases</b> (CARD, CRyPTIC, DRAGdb, Yang, Curated data), independent of species identity.")
              )
            ),
            
            column(
              width = 12,
              bs4Card(
                title = "Upset plot",
                width = 12,
                status = "info",
                solidHeader = TRUE,
                collapsible = TRUE,
                maximizable = TRUE,
                style = "height: 50vh; overflow-y: auto;",
                div(class = "plot-center", uiOutput("upset_container"))
              )
            ),
            
            column(
              width = 12,
              bs4Card(
                title = "Heatmap",
                width = 12,
                status = "info",
                solidHeader = TRUE,
                collapsible = TRUE,
                maximizable = TRUE,
                class = "heatmap-card",  
                style = "height: 90vh; overflow-y: auto;",
                
                div(
                  class = "heatmap-topbar",
                  actionButton(
                    "reset_heatmap_legend",
                    "Reset",
                    icon  = icon("rotate-left"),
                    class = "btn-outline-secondary btn-sm"
                  )
                ),
                
                div(class = "plot-center", uiOutput("heatmap_container"))
              )
            )
          ),
          
          ## ---- Step 4: Summary ----
          tabPanel(
            title = "Summary",
            value = "summary",
            uiOutput("summary_table_ui")
          )
        )
      )
    )
  )
)
