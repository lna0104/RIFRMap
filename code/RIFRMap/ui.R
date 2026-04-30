library(shiny)
library(shinyjs)
library(shinyBS)
library(plotly)
library(leaflet)
library(shinyWidgets)
library(ggiraph)
library(DT)
library(gt)

shinyUI(
  navbarPage(
    title = "RIFRMap",
    id = "main_nav",
    theme = "www/styles.css",
    footer = includeHTML("www/footer.html"),
    fluid = TRUE,
    collapsible = TRUE,
    
    
    header = tagList(
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
        tags$script(src = "script.js")
      ),
      shinyjs::useShinyjs()
    ),
    
    ## ---------------- Home ----------------
    tabPanel(
      title = tagList(icon("home"), "Home"),
      includeHTML("www/home.html")
      
    ),

    ## ---------------- Explore ----------------
    tabPanel(
      title = tagList(icon("chart-line"), "Explore"),
      value = "explore",
      fluidPage(
        class = "filters-row",
        
        tabsetPanel(
          id = "explore_steps",
          
          ## ---- Data Filters ----
          tabPanel(
            "Data Filters",
            value = "filters",
            div(
              class = "text-muted-explore",
              p("Select data sources and filtering criteria to define the mutation set used in downstream analyses."),
              p(tags$i("Last updated: 24/03/2026"))
            ),
            fluidRow(
              column(
                width = 5,
                
                div(
                  class = "filter-box",
                  div(
                    class = "filter-box-header",
                    checkboxInput("card_enabled", NULL, FALSE),
                    h4("CARD")
                  ),
                  selectInput(
                    "card_species", "Species",
                    choices = NULL,
                    selected = "All species",
                    multiple = TRUE,
                    width = "100%"
                  )
                ),
                
                div(
                  class = "filter-box",
                  div(
                    class = "filter-box-header",
                    checkboxInput("yang_enabled", NULL, FALSE),
                    h4("Yang-RIF")
                  ),
                  div(
                    class = "form-group shiny-input-container",
                    style = "width: 100%;",
                    tags$label("Species", `for` = "yang_species"),
                    tags$input(
                      id = "yang_species",
                      type = "text",
                      value = "Escherichia coli",
                      class = "form-control",
                      readonly = "readonly"
                    )
                  ), 
                  sliderInput(
                    "yang_log2fc", "log2FC",
                    min = 0, max = 2.7,
                    value = c(0, 2.7),
                    step = 0.1,
                    width = "100%"
                  )
                ),
                
                div(
                  class = "filter-box",
                  div(
                    class = "filter-box-header",
                    checkboxInput("dragdb_enabled", NULL, FALSE),
                    h4("DRAGdb")
                  ),
                  
                  selectInput(
                    "dragdb_species", "Species",
                    choices = NULL,
                    selected = "All species",
                    multiple = TRUE,
                    width = "100%"
                  ),
                  div(
                    style = "display:flex; align-items:center; gap:6px;",
                    
                    tags$label("PROVEAN score"),
                    
                    tags$i(
                      class = "fa fa-question-circle info-icon",
                      title = "PROVEAN score predicts the functional impact of amino acid substitutions. Variants with scores ≤ −2.5 are considered deleterious, while scores > −2.5 are considered neutral.",
                      `data-toggle` = "tooltip"
                    )
                  ),
                  sliderInput(
                    "dragdb_provean", NULL,
                    min = -10.7, max = -2.6,
                    value = c(-10.7, -2.6),
                    step = 0.5,
                    width = "100%"
                  )
                )
              ),
              
              column(
                width = 7,
                
                div(
                  class = "filter-box",
                  div(
                    class = "filter-box-header",
                    checkboxInput("cryptic_enabled", NULL, FALSE),
                    h4("CRyPTIC")
                  ),
                  
                  div(
                    class = "form-group shiny-input-container",
                    style = "width: 100%;",
                    tags$label("Species", `for` = "cryptic_species"),
                    tags$input(
                      id = "cryptic_species",
                      type = "text",
                      value = "Mycobacterium tuberculosis",
                      class = "form-control",
                      readonly = "readonly"
                    )
                  ),
                  
                  div(
                    style = "display:flex; align-items:center; gap:6px;",
                    tags$label("Effect size"),
                    tags$i(
                      class = "fa fa-question-circle info-icon",
                      title = "Effect size represents the association between a mutation and changes in MIC. Positive values indicate increased resistance, while negative values suggest increased susceptibility.",
                      `data-toggle` = "tooltip"
                    )
                  ),
                  
                  sliderInput(
                    "cryptic_coef", NULL,
                    min = -4.0, max = 9.6,
                    value = c(-4, 9.6),
                    step = 0.1,
                    width = "100%"
                  ),
                  
                  div(
                    style = "display:flex; align-items:center; gap:6px;",
                    
                    tags$label("Counts"),
                    
                    tags$i(
                      class = "fa fa-question-circle info-icon",
                      title = "Counts represent the number of isolates in which a given mutation was observed in the CRyPTIC dataset.",
                      `data-toggle` = "tooltip"
                      )
                  ),
                  
                  sliderInput(
                    "cryptic_counts", NULL,
                    min = 0, max = 2326,
                    value = c(0, 2326),
                    step = 1,
                    width = "100%"
                  ),
                  
                  div(
                    style = "display:flex; align-items:center; gap:6px;",
                    
                    tags$label("BH-corrected p-value (≤)"),
                    
                    tags$i(
                      class = "fa fa-question-circle info-icon",
                      title = "Benjamini–Hochberg corrected p-value controlling for multiple testing. Lower values indicate stronger statistical support for association with resistance.",
                      `data-toggle` = "tooltip"
                    )
                  ), 
                  
                  radioButtons(
                    "cryptic_bhp",
                    NULL,
                    choices = c("0.05" = 0.05, "0.01" = 0.01, "0.005" = 0.005, "0.001" = 0.001),
                    selected = 0.05,
                    inline = TRUE
                  ),
                ),
                
                div(
                  class = "filter-box",
                  div(
                    class = "filter-box-header",
                    checkboxInput("curated_enabled", NULL, FALSE),
                    h4("Curated Data")
                  ),
                  selectInput(
                    "curated_species", "Species",
                    choices = NULL,
                    selected = "All species",
                    multiple = TRUE,
                    width = "100%"
                  ),
                  sliderInput(
                    "curated_nstudies", "Number of studies",
                    min = 1, max = 59,
                    value = c(1, 59),
                    step = 1,
                    width = "100%"
                  ),
                  sliderInput(
                    "curated_nspecies", "Number of species",
                    min = 1, max = 29,
                    value = c(1, 29),
                    step = 1,
                    width = "100%"
                  )
                )
              )
            ),
            br(),
            actionButton("apply_filters", "APPLY", icon = icon("play"), class = "btn-primary"),
            actionButton("reset_filters", "RESET", icon = icon("rotate-left"))
          ),
          
          ## ---- Mutation Map ----
          tabPanel(
            "Mutation Map",
            value = "cluster",
            div(
              class = "text-muted-explore",
              p(
                "Species-centric mutation cluster view. Data sources are not considered; the plot visualises the distribution of resistance mutations within each species across ",
                tags$i("rpoB"),
                " clusters."
              )
            ),
            div(
              class = "plot-container",
              
              div(
                class = "plot-header",
                h3("Species-level distribution of resistance mutations")
              ),
              
              sidebarLayout(
                sidebarPanel(
                  div(
                    class = "plot-controls",
                    
                    selectInput(
                      inputId = "colour_by",
                      label = "Colour by",
                      choices = c(
                        "None" = "none",
                        "Distribution across species" = "distribution",
                        "Amino acid type" = "aa_type",
                        "Amino acid distance" = "aa_distance"
                      ),
                      selected = "none",
                      width = "100%"
                    ),
                    
                    uiOutput("legend_panel"),
                    
                    shinyWidgets::checkboxGroupButtons(
                      inputId = "clustersToPlot",
                      label = "Clusters",
                      choices = c("N", "I", "II", "III"),
                      selected = c("N", "I", "II", "III"),
                      direction = "horizontal",
                      justified = TRUE
                    ),
                    
                    selectInput(
                      "mutation_map_format",
                      "Download format",
                      choices = c("png", "pdf", "jpg", "jpeg", "tiff"),
                      selected = "png"
                    ),
                    
                    downloadButton("download_mutation_map", "Download plot")
                  ),
                  width = 3
                ),
                
                mainPanel(
                  div(
                    class = "plot-main",
                    uiOutput("clusterPlotNote"),
                    uiOutput("clusterPlotUI")
                  ),
                  width = 9
                )
              )
            )
          ),
          
          ## ---- Mutation Frequency ----
          tabPanel(
            "Mutation Frequency",
            div(
              class = "text-muted-explore",
              p(
                "The heatmap summarises the distribution of mutations across species. Each tile represents an amino acid position (",
                tags$i("E. coli"),
                " coordinates) and its substituted residue, with colour indicating the number of species harbouring the mutation."
              )
            ),
            div(
              class = "plot-container",
              div(
                class = "plot-header",
                h3("Species count per amino acid substitution")
              ),
              
              fluidRow(
                column(
                  width = 8,
                  div(
                    class = "plot-main", 
                    style = "height: 850px; overflow-y: auto; border: 1px solid #eee;",
                    ggiraph::girafeOutput("heatmap_pos_mut", height = "1200px")
                  )
                ),
                column(
                  width = 4,
                  div(
                    style = "margin-bottom: 20px;",
                    selectInput(
                      "mutation_freq_format",
                      "Download format",
                      choices = c("png", "pdf", "jpg", "jpeg", "tiff"),
                      selected = "png"
                    ),
                    
                    downloadButton("download_mutation_freq", "Download plot")
                  ),
                  
                  div(
                    uiOutput("selected_tile_info"),
                    DT::DTOutput("selected_tile_species_table")
                  )
                )
              )
            )
          ),
          
          ## ---- Source Overlap ----
          tabPanel(
            "Source Overlap",
            div(
              class = "text-muted-explore",
              p("Source-centric view of mutation evidence. Plots highlight overlap among databases independent of species identity.")
            ),
            fluidRow(
              column(
                width = 12,
                div(
                  class = "plot-main",
                  
                  div(
                    class = "plot-header",
                    h3("Cross-database mutation overlap")
                  ),
                  
                  div(
                    class = "plot-controls-bottom-left",
                    
                    selectInput(
                      "sources_upset_format",
                      "Download format",
                      choices = c("png", "pdf", "jpg", "jpeg", "tiff"),
                      selected = "png",
                      width = "200px"
                    ),
                    
                    downloadButton("download_sources_upset", "Download plot", style = "width:130px;")
                  ),
                  
                  uiOutput("upset_container")
                )
              ),
              
              column(
                width = 12,
                div(
                  class = "plot-main",
                  
                  div(
                    class = "plot-header plot-header-between",
                    
                    div(
                      h3("Cross-database map of resistance mutations"),
                      p(
                        class = "text-muted-explore",
                        "Click legend items to show or hide sources. Use Reset to restore the full heatmap view."
                      )
                    ),
                    
                    actionButton(
                      "reset_heatmap_legend",
                      "Reset",
                      icon = icon("rotate-left"),
                      class = "btn-outline-secondary btn-sm"
                    )
                  ),
                  
                  div(
                    class = "plot-controls-bottom-left",
                    
                    selectInput(
                      "sources_hm_format",
                      "Download format",
                      choices = c("png", "pdf", "jpg", "jpeg", "tiff"),
                      selected = "png",
                      width = "200px"
                    ),
                    
                    downloadButton("download_sources_hm", "Download plot", style = "width:130px;")
                  ),
                  
                  uiOutput("heatmap_container")
                )
              )
            )
          ),
          
          ## ---- Summary ----
          tabPanel(
            "Summary",
            div(
              class = "text-muted-explore",
              p("Summary table and histograms are generated from all rows in the filtered dataset; mutations reported by multiple sources may appear more than once.")
            ),
            fluidRow(
              column(
                width = 4,
                div(
                  class = "plot-main",
                  
                  div(
                    class = "plot-header",
                    h3("Summary")
                  ),
                  
                  gt_output("summary_table")
                )
              ),
              
              column(
                width = 8,
                
                div(
                  class = "plot-main",
                  
                  div(
                    class = "plot-header",
                    h3("Top amino acid mutations")
                  ),
                  
                  plotOutput("bar_aa_mut_name", height = "320px")
                ),
                
                div(
                  class = "plot-main",
                  
                  div(
                    class = "plot-header",
                    h3("AA position (E. coli) distribution")
                  ),
                  
                  plotOutput("hist_aa_pos_ecoli", height = "280px")
                )
              )
            )
          )
        )
      )
    ),
    
    ## ---------------- Data ----------------
    tabPanel(
      title = tagList(icon("database"), "Data"),
      
      div(
        class = "doc-hero",
        div(
          class = "doc-hero-content",
          h2("Data")
        )
      ),
      
      fluidRow(
        column(
          width = 12,
          div(
            class = "table-header",
            h3("Mutation Data")
          ), 
          div(
            style = "display:flex; gap:10px; justify-content:flex-end; margin-bottom:8px;",
            downloadButton(
              "download_mutations",
              "Download mutations (CSV)",
              class = "btn btn-outline-primary btn-sm"
            )
          ),
          DT::DTOutput("mutation_table"),
          br(),
          div(
            class = "table-header",
            h3("Reference Data")
          ), 
          div(
            style = "display:flex; gap:10px; justify-content:flex-end; margin-bottom:8px;",
            downloadButton(
              "download_refs",
              "Download references (CSV)",
              class = "btn btn-outline-secondary btn-sm"
            )
          ),
          DT::DTOutput("ref_table")
        )
      )
    ),
    
    ## ---------------- Documentation ----------------
    tabPanel(
      title = tagList(icon("book"), "Docs"),
      value = "doc",
      includeHTML("www/doc.html")
    )
  )
)