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
      "Home",
      includeHTML("www/home.html")
      
    ),
    
    ## ---------------- Documentation ----------------
    tabPanel(
      "Documentation",
      fluidPage(
        h2("Data Sources and Curation"),
        br(),
        
        h4("Comprehensive Antibiotic Resistance Database (CARD)"),
        p(
          "The Comprehensive Antibiotic Resistance Database (CARD) is a hand-curated resource developed to cover the full spectrum of antimicrobial resistance genes (ARGs). ",
          "To be included, an AMR determinant must be reported in a peer-reviewed publication, have a corresponding sequence deposited in GenBank, ",
          "and show clear experimental evidence of increased minimum inhibitory concentration (MIC) compared with control strains ",
          tags$a("[1]", href = "#ref1", class = "ref-link"), "."
        ),
        p(
          "For ", tags$i("Mycobacterium tuberculosis"),
          ", additional curated mutation data were integrated from several major resistance catalogs, including CRyPTIC, ReSeqTB, and the WHO 2023 mutation catalog."
        ),
        
        br(),
        
        h4("Comprehensive Resistance Prediction for Tuberculosis: an International Consortium (CRyPTIC)"),
        p(
          "The CRyPTIC Consortium generated a large-scale dataset by measuring quantitative MICs for 13 anti-tuberculosis drugs—including rifampicin—in over 15,000 globally collected ",
          tags$i("M. tuberculosis"),
          " isolates. Using multivariable linear mixed models alongside whole-genome sequencing, the study identified variants associated with increased MICs ",
          "(resistance) and decreased MICs (hypersusceptibility) ",
          tags$a("[2]", href = "#ref2", class = "ref-link"), ". ",
          "RIFRMap enables interactive exploration of this dataset with filters for effect size, variant count, adjusted p-value, and homoplasy."
        ),
        
        br(),
        
        h4("Drug Resistance Associated Genes Database (DRAGdb)"),
        p(
          "DRAGdb is a manually curated repository of mutational data for drug resistance-associated genes across ESKAPE pathogens and other bacteria, ",
          "with a particular focus on ", tags$i("Mycobacterium tuberculosis"), ". ",
          "The database integrates PROVEAN scores to predict the functional impact of gene mutations ",
          tags$a("[3]", href = "#ref3", class = "ref-link"), "."
        ),
        
        br(),
        
        h4("Yang-RIF (2023)"),
        p(
          "The Yang-RIF dataset is based on a study by Yang and colleagues, who performed saturation mutagenesis of the rifampicin-binding site in ",
          tags$i("E. coli"), " RNA polymerase. They constructed single-residue ", tags$i("rpoB"),
          " mutants spanning positions 510–537 and 563–572 to systematically assess resistance phenotypes ",
          tags$a("[4]", href = "#ref4", class = "ref-link"), ". ",
          "Within RIFRMap, users can explore this dataset by filtering mutations based on enrichment after rifampicin exposure."
        ),
        
        br(),
        
        h4("Curated Data"),
        p(
          "This curated dataset was developed through a literature review of rifampicin resistance mutations in ",
          tags$i("rpoB"), ". Relevant studies were identified using targeted Google Scholar searches, and single amino acid substitutions reported to confer rifampicin resistance were manually extracted ",
          tags$a("[5]", href = "#ref5", class = "ref-link"), "."
        ),
        
        br(),
        
        h4("References"),
        tags$ol(
          tags$li(
            id = "ref1",
            tags$a(
              href = "https://doi.org/10.1093/nar/gkac920",
              target = "_blank",
              "Alcock, B. P., Huynh, W., Chalil, R., et al. (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 51(D1), D690–D699."
            )
          ),
          tags$li(
            id = "ref2",
            tags$a(
              href = "https://doi.org/10.1038/s41467-023-44325-5",
              target = "_blank",
              "The CRyPTIC Consortium. (2024). Quantitative measurement of antibiotic resistance in Mycobacterium tuberculosis reveals genetic determinants of resistance and susceptibility in a target gene approach. Nature Communications, 15, 488."
            )
          ),
          tags$li(
            id = "ref3",
            tags$a(
              href = "https://doi.org/10.1038/s41598-020-65766-8",
              target = "_blank",
              "Ghosh, A., N, S., & Saha, S. (2020). Survey of drug resistance associated gene mutations in Mycobacterium tuberculosis, ESKAPE and other bacterial species. Scientific Reports, 10(1), 8957."
            )
          ),
          tags$li(
            id = "ref4",
            tags$a(
              href = "https://doi.org/10.1038/s41586-023-06495-6",
              target = "_blank",
              "Yang, K.B., Cameranesi, M., Gowder, M., et al. (2023). High-resolution landscape of an antibiotic binding site. Nature, 622, 180–187."
            )
          ),
          tags$li(
            id = "ref5",
            tags$a(
              href = "https://doi.org/10.1073/pnas.2424307122",
              target = "_blank",
              "Bolourchi, N., Brown, C.R.P., Letten, A.D., & Engelstädter, J. (2025). Evolution and evolvability of rifampicin resistance across the bacterial tree of life. Proceedings of the National Academy of Sciences of the United States of America, 122(31), e2424307122."
            )
          )
        )
      )
    ),
    
    ## ---------------- Explorer ----------------
    tabPanel(
      "Explorer",
      fluidPage(
        tabsetPanel(
          id = "explore_steps",
          
          ## ---- Data Filters ----
          tabPanel(
            "Data Filters",
            div(
              class = "text-muted",
              p("Select data sources and filtering criteria to define the mutation set used in downstream analyses."),
              p(tags$i("Last updated: 24/03/2026"))
            ),
            fluidRow(
              column(
                width = 2,
                wellPanel(
                  h4("CARD"),
                  checkboxInput("card_enabled", "Include data", FALSE),
                  selectInput("card_species", "Species", choices = NULL, selected = "All species", multiple = TRUE)
                )
              ),
              column(
                width = 2,
                wellPanel(
                  h4("CRyPTIC"),
                  checkboxInput("cryptic_enabled", "Include data", FALSE),
                  selectInput(
                    "cryptic_species", "Species",
                    choices = "Mycobacterium tuberculosis",
                    selected = "Mycobacterium tuberculosis",
                    multiple = FALSE
                  ),
                  sliderInput("cryptic_coef", "Effect size", min = -4.0, max = 9.6, value = c(-4, 9.6), step = 0.1),
                  sliderInput("cryptic_counts", "Counts", min = 0, max = 2326, value = c(0, 2326), step = 1),
                  radioButtons(
                    "cryptic_bhp",
                    "BH-corrected p-value (≤)",
                    choices = c("0.05" = 0.05, "0.01" = 0.01, "0.005" = 0.005, "0.001" = 0.001),
                    selected = 0.05,
                    inline = TRUE
                  ),
                  selectInput("cryptic_homoplasy", "Homoplasy", choices = c("", "TRUE", "FALSE"), selected = "")
                )
              ),
              column(
                width = 2,
                wellPanel(
                  h4("DRAGdb"),
                  checkboxInput("dragdb_enabled", "Include data", FALSE),
                  selectInput("dragdb_species", "Species", choices = NULL, selected = "All species", multiple = TRUE),
                  sliderInput("dragdb_provean", "PROVEAN", min = -10.7, max = -2.6, value = c(-10.7, -2.6), step = 0.5)
                )
              ),
              column(
                width = 2,
                wellPanel(
                  h4("Yang-RIF"),
                  checkboxInput("yang_enabled", "Include data", FALSE),
                  selectInput(
                    "yang_species", "Species",
                    choices = "Escherichia coli",
                    selected = "Escherichia coli",
                    multiple = FALSE
                  ),
                  sliderInput("yang_log2fc", "log2FC", min = 0, max = 2.7, value = c(0, 2.7), step = 0.1)
                )
              ),
              column(
                width = 2,
                wellPanel(
                  h4("Curated Data"),
                  checkboxInput("curated_enabled", "Include data", FALSE),
                  selectInput("curated_species", "Species", choices = NULL, selected = "All species", multiple = TRUE),
                  sliderInput("curated_nstudies", "Number of studies", min = 1, max = 59, value = c(1, 59), step = 1),
                  sliderInput("curated_nspecies", "Number of species", min = 1, max = 29, value = c(1, 29), step = 1)
                )
              )
            ),
            br(),
            actionButton("apply_filters", "Apply", icon = icon("play"), class = "btn-primary"),
            actionButton("reset_filters", "Reset", icon = icon("rotate-left"))
          ),
          
          ## ---- Mutation Map ----
          tabPanel(
            "Mutation Map",
            div(
              class = "text-muted",
              p(
                "Species-centric mutation cluster view. Data sources are not considered; the plot visualises the distribution of resistance mutations within each species across ",
                tags$i("rpoB"),
                " clusters."
              )
            ),
            sidebarLayout(
              sidebarPanel(
                selectInput(
                  inputId = "colour_by",
                  label = "Colour by",
                  choices = c(
                    "None" = "none",
                    "Distribution across species" = "distribution",
                    "Amino acid type" = "aa_type",
                    "Amino acid distance" = "aa_distance"
                  ),
                  selected = "none"
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
                width = 3
              ),
              mainPanel(
                h4("Species-level distribution of resistance mutations"),
                uiOutput("clusterPlotNote"),
                uiOutput("clusterPlotUI"),
                width = 9
              )
            )
          ),
          
          ## ---- Mutation Frequency ----
          tabPanel(
            "Mutation Frequency",
            div(
              class = "text-muted",
              p(
                "The heatmap summarises the distribution of mutations across species. Each tile represents an amino acid position (",
                tags$i("E. coli"),
                " coordinates) and its substituted residue, with colour indicating the number of species harbouring the mutation."
              )
            ),
            fluidRow(
              column(
                width = 8,
                h4("Species count per amino acid substitution"),
                div(
                  style = "height: 850px; overflow-y: auto; border: 1px solid #eee;",
                  ggiraph::girafeOutput("heatmap_pos_mut", height = "1200px")
                )
              ),
              column(
                width = 4,
                h4("Species in selected tile"),
                uiOutput("selected_tile_info"),
                DT::DTOutput("selected_tile_species_table")
              )
            )
          ),
          
          ## ---- Source Overlap ----
          tabPanel(
            "Source Overlap",
            div(
              class = "text-muted",
              p("Source-centric view of mutation evidence. Plots highlight overlap among databases independent of species identity.")
            ),
            fluidRow(
              column(
                width = 12,
                h4("Cross-database mutation overlap"),
                uiOutput("upset_container")
              ),
              column(
                width = 12,
                br(),
                h4("Cross-database map of resistance mutations"),
                p(class = "text-muted", "Click legend items to show/hide sources. Use Reset to restore the full heatmap view."),
                actionButton(
                  "reset_heatmap_legend",
                  "Reset",
                  icon = icon("rotate-left"),
                  class = "btn-outline-secondary btn-sm"
                ),
                uiOutput("heatmap_container")
              )
            )
          ),
          
          ## ---- Summary ----
          tabPanel(
            "Summary",
            div(
              class = "text-muted",
              p("Summary table and histograms are generated from all rows in the filtered dataset; mutations reported by multiple sources may appear more than once.")
            ),
            fluidRow(
              column(
                width = 4,
                h4("Summary"),
                gt_output("summary_table")
              ),
              column(
                width = 8,
                h4("Top amino acid mutations"),
                plotOutput("bar_aa_mut_name", height = "320px"),
                br(),
                h4("AA_Pos_Ecoli histogram"),
                plotOutput("hist_aa_pos_ecoli", height = "280px")
              )
            )
          ),
          
          ## ---- Download ----
          tabPanel(
            "Download",
            fluidRow(
              column(
                width = 12,
                h4("Mutation Data"),
                div(
                  style = "display:flex; gap:10px; justify-content:flex-end; margin-bottom:8px;",
                  downloadButton("download_mutations", "Download mutations (CSV)", class = "btn btn-outline-primary btn-sm")
                ),
                DT::DTOutput("mutation_table"),
                br(),
                h4("Reference Data"),
                div(
                  style = "display:flex; gap:10px; justify-content:flex-end; margin-bottom:8px;",
                  downloadButton("download_refs", "Download references (CSV)", class = "btn btn-outline-secondary btn-sm")
                ),
                DT::DTOutput("ref_table")
              )
            )
          )
        )
      )
    )
  )
)