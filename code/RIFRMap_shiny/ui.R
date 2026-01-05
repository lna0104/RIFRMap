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
      bs4SidebarMenuItem("Home", tabName = "home", icon = icon("sliders")),
      bs4SidebarMenuItem("View", tabName = "view", icon = icon("chart-simple"))
    )
  ),
  
  body = dashboardBody(
    tags$head(
      tags$style(HTML("
        .filter-form .form-header { font-weight: 700; }
        .filter-form { min-width: 320px; }
      "))
    ),
    
    bs4TabItems(
      # ---------------- HOME ----------------
      bs4TabItem(
        tabName = "home",
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
              selectInput("cryptic_species", "Species", choices = NULL, selected = "All species", multiple = TRUE),
              sliderInput("cryptic_coef", "Effect size", min = -4.0, max = 9.6, value = c(-4, 9.6), step = 0.1),
              sliderInput("cryptic_counts", "Counts", min = 0, max = 2326, value = c(0, 2326), step = 1),
              sliderInput("cryptic_bhp", "BH p-value", min = 0, max = 1, value = c(0, 1), step = 0.01),
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
              selectInput("yang_species", "Species", choices = NULL, selected = "All species", multiple = TRUE),
              sliderInput("yang_log2fc", "log2FC", min = 0, max = 2.7, value = c(0, 2.7), step = 0.1)
            )
          ),
          
          # Negin
          column(
            width = 12, lg = 3, xl = 2,
            bs4Card(
              title = "Negin",
              width = 12,
              status = "secondary",
              solidHeader = TRUE,
              collapsible = TRUE,
              checkboxInput("negin_enabled", "Include data", FALSE),
              selectInput("negin_species", "Species", choices = NULL, selected = "All species", multiple = TRUE),
              sliderInput("negin_nstudies", "Number of studies (per mutation)", min = 1, max = 59, value = c(1, 59), step = 1),
              sliderInput("negin_nspecies", "Number of species (per mutation)", min = 1, max = 29, value = c(1, 29), step = 1)
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
            
            div(
              style = "padding: 12px 16px;",
              tags$p(
                class = "text-muted",
                style = "margin-bottom: 12px; font-size: 0.9rem;",
                HTML("Species-centric mutation cluster view. Sources are ignored; patterns reflect <b>species</b> and <b>cluster</b> structure.")
              )
            ),
            
            
            fluidRow(
              ## ---- Controls ----
              column(
                width = 2,
                card(
                  selectInput(
                    inputId = "clustersToPlot",
                    label   = "Clusters",
                    choices = c("N", "I", "II", "III"),
                    selected = "",
                    multiple = TRUE
                  ),
                  
                  div(
                    uiOutput("legend_panel")
                  )
                ), 
              ),
              
              ## ---- Plot ----
              
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
          
          # ---- TAB 2: Source Plot ----
          tabPanel(
            title = "Sources",
            
            div(
              style = "padding: 12px 16px;",
              tags$p(
                class = "text-muted",
                style = "margin-bottom: 12px; font-size: 0.9rem;",
                HTML("Source-centric view of mutation evidence. Plots highlight overlap among <b>databases</b> (CARD, CRyPTIC, DRAGdb, Yang, Negin), independent of species identity.")
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
                uiOutput("upset_container")
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
                uiOutput("heatmap_container")
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



