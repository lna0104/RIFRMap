## ---- load packages ----
suppressPackageStartupMessages({
  # Core UI
  library(shiny)
  library(bs4Dash)
  
  # Data
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  library(gtsummary)
  
  
  # Plots / widgets
  library(ggplot2)
  library(ggiraph)
  library(plotly)
  library(ggupset)
  library(shinyWidgets)
  library(shinycssloaders)
  library(gt)
  library(wesanderson)
  library(RColorBrewer)
  
  # UI helpers
  library(bslib)
  library(htmltools)

  
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
      bs4SidebarMenuItem("Documentation", tabName = "doc", icon = icon("book-open")),
      bs4SidebarMenuItem("Explorer", tabName = "explorer", icon = icon("compass"))
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
      
      ## ---------------- Documentation ----------------
      bs4TabItem(
        tabName = "doc",
        h2("Information")
      ),
      
      ## ---------------- Explorer ----------------
      bs4TabItem(
        tabName = "explorer",
        
        tabsetPanel(
          id = "explore_steps",
          
          ## ---- Tab 1: Data Filters ----
          tabPanel(
            title = "Data Filters",
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
                    choices  = "Mycobacterium tuberculosis",
                    selected = "Mycobacterium tuberculosis",
                    multiple = FALSE
                  ),
                  tags$script(HTML("$('#cryptic_species').prop('disabled', true);")),
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
                  status = "secondary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  checkboxInput("yang_enabled", "Include data", FALSE),
                  selectInput(
                    "yang_species", "Species",
                    choices  = "Escherichia coli",
                    selected = "Escherichia coli",
                    multiple = FALSE
                  ),
                  tags$script(HTML("$('#yang_species').prop('disabled', true);")),
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
                  status = "info",
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
                actionButton("apply_filters", "Apply", icon = icon("play"), class = "btn-primary"),
                actionButton("reset_filters", "Reset", icon = icon("rotate-left"))
              )
            )
          ),
          
          ## ---- Tab 2: Mutation Map ----
          tabPanel(
            title = "Mutation Map",
            value = "cluster",
            
            div(
              class = "step-intro",
              tags$p(
                class = "text-muted",
                HTML("Species-centric mutation cluster view. Data sources are not considered; the plot visualises the distribution of resistance mutations within each species across <i>rpoB</i> clusters.")
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
                bs4Card(
                  title = "Species-level distribution of resistance mutations",
                  width = 12,
                  status = "info",
                  collapsible = FALSE,
                  solidHeader = TRUE,
                  style = "height: 90vh;",
                  uiOutput("clusterPlotUI")
                )
              )
            )
          ),
          
          ## ---- Tab 3: Mutation Frequency ----
          tabPanel(
            title = "Mutation Frequency",
            value = "frequency",
            
            div(
              class = "step-intro",
              tags$p(
                class = "text-muted",
                HTML("Heatmap summarises how widely each mutation is distributed across species. Each tile represents an amino-acid position (<i>E. coli</i> coordinates) and substituted residue, coloured by the number of species carrying that mutation. 
                     Data sources are ignored.")
              )
            ),
            
            fluidRow(
              # --- Heatmap card ---
              column(
                width = 8,
                bs4Card(
                  width = 12,
                  title = "Species count per amino acid substitution",
                  status = "info",
                  solidHeader = TRUE,
                  div(
                    style = "height: 850px; overflow-y: auto; border: 1px solid #eee;",
                    ggiraph::girafeOutput("heatmap_pos_mut", height = "1200px")
                  )
                )
              ),
              
              # --- Species list card ---
              column(
                width = 4,
                bs4Card(
                  width = 12,
                  title = "Species in selected tile",
                  status = "info",
                  solidHeader = TRUE,
                  
                  uiOutput("selected_tile_info"),
                  
                  DT::DTOutput("selected_tile_species_table")
                )
              )
            )
          ),
          
          
          ## ---- Tab 4: Source Overlap ----
          tabPanel(
            title = "Source Overlap",
            value = "sources",
            
            div(
              class = "step-intro",
              tags$p(
                class = "text-muted",
                HTML("Source-centric view of mutation evidence. Plots highlight overlap among databases (CARD, CRyPTIC, DRAGdb, Yang, Curated data), independent of species identity.")
              )
            ),
            
            column(
              width = 12,
              bs4Card(
                title = "Cross-database mutation overlap",
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
                title = "Cross-database map of resistance mutations",
                width = 12,
                status = "info",
                solidHeader = TRUE,
                collapsible = TRUE,
                maximizable = TRUE,
                class = "heatmap-card",  
                style = "height: 90vh; overflow-y: auto;",
                div(
                  class = "step-intro",
                  tags$p(
                    class = "text-muted",
                    HTML("Double-click legend items to hide/show source groups. Use Reset to restore the full heatmap view.")
                  )
                ),
                
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
          
          
          ## ---- Tab 5: Summary ----
          tabPanel(
            title = "Summary",
            tags$p(
              class = "text-muted",
              HTML("Summary table and histograms update dynamically based on selected filters.")
            ),
            
            
            fluidRow(
              # left: summary table 
              column(
                width = 4,
                bs4Card(
                  width = 12,
                  title = "Summary",
                  status = "info",
                  solidHeader = TRUE,
                  gt_output("summary_table")
                )
              ),
              
              # right: two histograms stacked 
              column(
                width = 8,
                
                bs4Card(
                  width = 12,
                  title = "Top AA mutations across species",
                  status = "info",
                  solidHeader = TRUE,
                  plotOutput("bar_aa_mut_name", height = "320px")
                ),
          
                
                bs4Card(
                  width = 12,
                  title = "AA_Pos histogram",
                  status = "info",
                  solidHeader = TRUE,
                  plotOutput("hist_aa_pos", height = "280px")
                ),
                
                bs4Card(
                  width = 12,
                  title = "AA_Pos_Ecoli histogram",
                  status = "info",
                  solidHeader = TRUE,
                  plotOutput("hist_aa_pos_ecoli", height = "280px")
                )
              )
            )
          ),
          
          ## ---- Tab 6: Download ----
          tabPanel(
            title = "Download",
            
            fluidRow(
              column(
                width = 12,
                bs4Card(
                  title = "Data table",
                  status = "info",
                  width = 12,
                  div(
                    style = "display:flex; gap:10px; justify-content:flex-end; margin-bottom:8px;",
                    downloadButton("download_table", "Download CSV", class = "btn btn-outline-primary btn-sm")
                  ),
                  DT::DTOutput("table_data")
                )
              )
              
            )
          )
        )
      )
    )
  )
)
