## ---- load packages ----

suppressPackageStartupMessages({
  
  ## CRAN packages
  library(bs4Dash)
  library(shiny)
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  library(ggiraph)
  library(bslib)
  library(ggplot2)
  library(ggupset)
  library(base64enc)
  library(htmltools)
  library(shinyWidgets)

  ## Bioconductor
  library(Biostrings)
  
  ## GitHub
  library(ALJEbinf)
})


# ---- ui ----

ui <- dashboardPage(
  title = "RIFRMap",
  help = NULL,  
  
  header = dashboardHeader(
    title = dashboardBrand(
      title = "RIFRMap",
      color = "primary",
    ),
    skin = "light",
    status = "white",
    border = TRUE,
    sidebarIcon = icon("bars"),
    controlbarIcon = icon("cog")
  ),
  
  sidebar = dashboardSidebar(
    skin = "light",
    status = "primary",
    elevation = 3,
    bs4SidebarMenu(
      id = "sidebar",
      bs4SidebarMenuItem("Home", tabName = "home", icon = icon("home")),
      bs4SidebarMenuItem("Filters", tabName = "filters", icon = icon("sliders")),
      bs4SidebarMenuItem("View", tabName = "view", icon = icon("chart-simple"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$link(
        rel  = 'stylesheet',
        type = 'text/css',
        href = 'styles.css'
      ),
      tags$script(src = "script.js"),
    ),
    
    
    bs4TabItems(
      bs4TabItem(
        tabName = "home",
        h2("Welcome to RIFRMap")
          ),
      # ---------------- Filters ----------------
      bs4TabItem(
        tabName = "filters",
        h2("RIF Resistant Data"),
        
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
              selectInput("card_species", "Species", choices = NULL, selected = "All species", multiple = TRUE)
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
              selectInput("cryptic_species", "Species", choices = NULL, selected = "Mycobacterium tuberculosis", multiple = TRUE),
              sliderInput("cryptic_coef", "Effect size", min = -4.0, max = 9.6, value = c(-4, 9.6), step = 0.1),
              sliderInput("cryptic_counts", "Counts", min = 0, max = 2326, value = c(0, 2326), step = 1),
              tags$div(
                style = "width: 50vw;",
                radioButtons(
                  "cryptic_bhp",
                  "BH-corrected p-value (≤)",
                  choices = c(
                    "0.05"  = 0.05,
                    "0.01"  = 0.01,
                    "0.005" = 0.005,
                    "0.001" = 0.001
                  ),
                  selected = 0.05,
                  inline = TRUE
                )
              ),
              selectInput("cryptic_homoplasy", "Homoplasy", choices = c("", "TRUE", "FALSE"), selected = "", width = "140px")
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
              selectInput("dragdb_species", "Species", choices = NULL, selected = "All species", multiple = TRUE),
              sliderInput("dragdb_provean", "PROVEAN", min = -10.7, max = -2.6, value = c(-10.7, -2.6), step = 0.5)
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
              selectInput("yang_species", "Species", choices = NULL, selected = "Escherichia coli", multiple = TRUE),
              sliderInput("yang_log2fc", "log2FC", min = 0, max = 2.7, value = c(0, 2.7), step = 0.1)
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
              selectInput("curated_species", "Species", choices = NULL, selected = "All species", multiple = TRUE),
              sliderInput("curated_nstudies", "Number of studies (per mutation)", min = 1, max = 59, value = c(1, 59), step = 1),
              sliderInput("curated_nspecies", "Number of species (per mutation)", min = 1, max = 29, value = c(1, 29), step = 1)
            )
          ),
          
          br(),
          
          # Action buttons in a row
          fluidRow(
            column(
              width = 12,
              actionButton("apply_filters", "Apply filters", icon = icon("play"), class = "btn-primary"),
              actionButton("reset_filters", "Reset", icon = icon("rotate-left"))
            )
          )
          
        )
      ),
      
      # ---------------- VIEW ----------------
      bs4TabItem(
        tabName = "view",
        tabsetPanel(
          id = "view_tabs",
          
          # ---- TAB 1: Cluster plot ----
          tabPanel(
            title = "Species–Cluster",
            value = "cluster",
            
            div(
              style = "padding: 12px 16px;",
              tags$p(
                class = "text-muted",
                style = "margin-bottom: 12px; ",
                HTML("Species-centric mutation cluster view. Sources are ignored; patterns reflect <b>species</b> and <b>cluster</b> structure.")
              )
            ),
            
            
            fluidRow(
              ## ---- Controls ----
              column(
                width = 3,
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
                      ),
                    ),

                    uiOutput("legend_panel")
                  )
                )
              ),
              
              ## ---- Plot ----
              
              column(
                width = 9,
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
          
          # ---- TAB 2: Source Plot ----
          tabPanel(
            title = "Sources",
            value = "sources", 
            
            div(
              style = "padding: 12px 16px;",
              tags$p(
                class = "text-muted",
                style = "margin-bottom: 12px; ",
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
                div(
                  style = "display:flex; justify-content:center;",
                  uiOutput("upset_container")
                )
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
                style = "height: 90vh; overflow-y: auto;",
                div(
                  style = "display: flex; justify-content: center;",
                  uiOutput("heatmap_container")
                )
              )
            )),
          
          # ---- TAB 3: Summary Data ----
          
          tabPanel(
            "Summary Data",
            uiOutput("summary_table_ui")
          )
        )
      )
  )
))



