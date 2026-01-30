#shiny::runApp('code/RIFRMap_shiny')
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
  library(fresh)

  
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
      bs4SidebarMenuItem("Home", tabName = "home", icon = icon("home")),
      bs4SidebarMenuItem("Documentation", tabName = "doc", icon = icon("book-open")),
      bs4SidebarMenuItem("Explorer", tabName = "explorer", icon = icon("compass"))
    )
  ),
  
  dashboardBody(
    
    use_theme(create_theme(
      bs4dash_status(
        primary = "#0A9F9D",
        secondary = "#FBA72A",
        success = "#5785C1",
        info = "#CEB175",
        warning = "#AF4E24",
        danger = "#dc3545",
        light = "#f8f9fa",
        dark = "#343a40"
      )
    )),
   
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
      tags$script(src = "script.js")
    ),
    
    
    bs4TabItems(
      ## ---------------- Home ----------------
      bs4TabItem(
        tabName = "home",
        div(
          class = "text-center my-5",
          
          h1(
            class = "display-3 font-weight-bold",
            "Welcome to RIFRMap"
          ),
          
          p(
            class = "lead text-muted mt-3",
            "An interactive platform for exploring rifampicin resistance mutations across bacterial species and curated datasets."
          )
        ),
        
        br(), div(
          class = "container",
          
          h3("About Rifampicin Resistance"),
          
          p(
            "Rifampicin (RIF) is a frontline antibiotic widely used in the treatment of tuberculosis and other bacterial infections. 
       Resistance most commonly arises through mutations in the ",
            tags$i("rpoB"),
            " gene, which encodes the beta subunit of RNA polymerase. 
       These mutations alter the drug-binding site, reducing antibiotic effectiveness and contributing to the global burden 
       of antimicrobial resistance."
          ),
          
          p(
            "Notably, almost all RIF resistance-associated mutations are concentrated within four relatively small regions 
       of", tags$i("rpoB")," commonly referred to as clusters I, II, III, and N. These clusters represent mutational hotspots where 
       resistance repeatedly evolves across diverse bacterial taxa."
          ),
          
          p(
            tags$b("RIFRMap"), " enables systematic exploration of these resistance-associated mutations across diverse species and multiple curated data sources, including CARD, CRyPTIC, DRAGdb, Yang, and additional curated datasets. 
          By integrating these resources into a unified interface, RIFRMap enables comparative exploration, dynamic filtering, and deeper insights into the rifampicin resistant patterns."
          )
        ),
        
        br(),
        
        div(
          class = "text-center my-4",
          
          actionButton(
            "go_docs",
            "View Documentation",
            class = "btn btn-outline-secondary btn-lg mr-3"
          ),
          
          actionButton(
            "go_explorer",
            "Launch Explorer",
            class = "btn btn-primary btn-lg"
          )
        )
        
        
      ),
      
      ## ---------------- Documentation ----------------
      bs4TabItem(
        tabName = "doc",
        h2("Data Sources and Curation"),
        
        br(),
        
        h4("Comprehensive Antibiotic Resistance Database (CARD)"),
        
        p(
          "The Comprehensive Antibiotic Resistance Database (CARD) is a hand-curated resource developed to cover the full spectrum of antimicrobial resistance genes (ARGs). ",
          "To be included, an AMR determinant must meet three criteria: it must be reported in a peer-reviewed publication, have a corresponding sequence deposited in GenBank, ",
          "and demonstrate clear experimental evidence of increased minimum inhibitory concentration (MIC) compared with control strains ", tags$a("[1]", href = "#ref1", class = "ref-link"), "."
        ),
        
        p(
          "For ", tags$i("Mycobacterium tuberculosis"), 
          ", additional curated mutation data were integrated from several major resistance catalogs, including CRyPTIC, ReSeqTB, and the WHO 2023 mutation catalog."
        ),
        
        br(),
        
        h4("Comprehensive Resistance Prediction for Tuberculosis: an International Consortium (CRyPTIC)"),
        
        p(
          "The CRyPTIC Consortium generated a large-scale dataset by measuring quantitative MICs for 13 anti-tuberculosis drugs—including rifampicin—in over 15,000 globally collected ",
          tags$i("M. tuberculosis"), " isolates. Using multivariable linear mixed models alongside whole-genome sequencing, the study identified 492 variants associated with increased MICs ",
          "(resistance) and 91 variants linked to decreased MICs (hypersusceptibility)",
          tags$a("[2]", href = "#ref2", class = "ref-link"), ". ",
          "RIFRMap enables interactive exploration of this dataset with filters for effect size, variant count,
          Benjamini–Hochberg adjusted p-value and homoplasy, allowing users to prioritize high-confidence and convergently evolved mutations."
        ),
        
        br(),
        
        h4("Drug Resistance Associated Genes Database (DRAGdb)"),
        
        p("DRAGdb is a manually curated repository of mutational data for drug resistance–associated genes (DRAGs) across ESKAPE pathogens (i.e. ", tags$i("Enterococcus faecium"), ", ", tags$i("Staphylococcus aureus"), ", ", tags$i("Klebsiella pneumoniae"), ", ", tags$i("Acinetobacter baumannii"), ", ", tags$i("Pseudomonas aeruginosa"), ", and ", tags$i("Enterobacter"), " spp.), as well as other bacteria, with a particular focus on ", tags$i("Mycobacterium tuberculosis"), " (MTB). 
          The database integrates PROVEAN scores to predict the functional impact of gene mutations ", tags$a("[3]", href = "#ref3", class = "ref-link"), ". "),
        
        br(),
        
        h4("Yang-RIF (2023)"),
        
        p(
          "The \"Yang-RIF\" dataset is based on a study by Yang and colleagues, who performed saturation mutagenesis of the rifampicin-binding site in ",
          tags$i("E. coli"), " RNA polymerase. They constructed 760 single-residue ", tags$i("rpoB"), " mutants spanning positions 510–537 and 563–572 ",
          "to systematically assess resistance phenotypes ", tags$a("[4]", href = "#ref4", class = "ref-link"), ". ", "Within RIFRMap, users can explore this dataset by filtering mutations based on their enrichment after rifampicin exposure. ",
          "A customizable threshold on log₂(fold-change) allows identification of variants associated with increased resistance."
        ),
        
        br(),
        
        h4("Curated Data"),
        
        p(
          "This curated dataset was developed by our lab's members through a comprehensive literature review of rifampicin resistance mutations in ",
          tags$i("rpoB"),
          ". Relevant studies were identified using targeted Google Scholar searches with key terms such as ",
          tags$i("\"rifampicin resistance\""), ", ",
          tags$i("\"RNA polymerase subunit beta\""), ", ",
          tags$i("\"point mutation\""), ", and ",
          tags$i("\"amino acid substitution\""), ". ",
          "From these studies, they manually extracted single amino acid substitutions reported to confer rifampicin resistance ", tags$a("[5]", href = "#ref5", class = "ref-link"), ". ",
          "Within RIFRMap, users can filter this dataset based on the number of supporting studies and the diversity of bacterial species reporting each mutation, allowing tailored exploration of high-confidence variants."
        ),
        
        h4("References"),
        
        tags$ol(
          tags$li(
            id = "ref1",
            tags$a(
              href = "https://doi.org/10.1093/nar/gkac920",
              "Alcock, B. P., Huynh, W., Chalil, R., et al. (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 51(D1), D690–D699.",
              target = "_blank"
            )
          ),
          tags$li(
            id = "ref2",
            tags$a(
              href = "https://doi.org/10.1038/s41467-023-44325-5",
              "The CRyPTIC Consortium. (2024, January). Quantitative measurement of antibiotic resistance in Mycobacterium tuberculosis reveals genetic determinants of resistance and susceptibility in a target gene approach. Nature Communications, 15, 488.",
              target = "_blank"
            )
          ),
          tags$li(
            id = "ref3",
            tags$a(
              href = "https://doi.org/10.1038/s41598-020-65766-8",
              "Ghosh, A., N, S., & Saha, S. (2020). Survey of drug resistance associated gene mutations in Mycobacterium tuberculosis, ESKAPE and other bacterial species. Scientific Reports, 10(1), 8957.",
              target = "_blank"
            )
          ),
          tags$li(
            id = "ref4",
            tags$a(
              href = "https://doi.org/10.1038/s41586-023-06495-6",
              "Yang, K.B., Cameranesi, M., Gowder, M., et al. (2023). High-resolution landscape of an antibiotic binding site. Nature, 622, 180–187.",
              target = "_blank"
            )
          ),
          tags$li(
            id = "ref5",
            tags$a(
              href = "https://doi.org/10.1073/pnas.2424307122",
              "Bolourchi, N., Brown, C.R.P., Letten, A.D., & Engelstädter, J. (2025). Evolution and evolvability of rifampicin resistance across the bacterial tree of life. Proceedings of the National Academy of Sciences of the United States of America, 122(31), e2424307122.",
              target = "_blank"
            )
          )
          
        )
        
        
        
        
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
                  status = "info",
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
                  title = "Yang-RIF",
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
                  title = "Curated Data",
                  width = 12,
                  status = "primary",
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
                  status = "primary",
                  collapsible = FALSE,
                  solidHeader = TRUE,
                  style = "height: 90vh;",
                  tags$div(
                    class = "text-muted",
                    style = "font-size: 0.85em; margin-bottom: 6px;",
                    "Note: Mutations outside the defined clusters are not shown."
                  ),
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
                  status = "primary",
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
                  status = "primary",
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
                status = "primary",
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
                status = "primary",
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
                  status = "primary",
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
                  status = "primary",
                  solidHeader = TRUE,
                  plotOutput("bar_aa_mut_name", height = "320px")
                ),
          
                
                bs4Card(
                  width = 12,
                  title = "AA_Pos histogram",
                  status = "primary",
                  solidHeader = TRUE,
                  plotOutput("hist_aa_pos", height = "280px")
                ),
                
                bs4Card(
                  width = 12,
                  title = "AA_Pos_Ecoli histogram",
                  status = "primary",
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
                  status = "primary",
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
