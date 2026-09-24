shinyUI(
  navbarPage(
    title = "RIFRMap",
    id = "main_nav",
    theme = "www/styles.css",
    # footer = includeHTML("www/footer.html"),
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
              p(tags$i(paste0("Last updated: ", DATA_LAST_UPDATED)))
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
                    tags$label("Effect size (log₂ MIC)"),
                    tags$i(
                      class = "fa fa-question-circle info-icon",
                      title = "Effect size is the GWAS coefficient for the association between a mutation and the rifampicin MIC, in log₂(MIC) units (1 = a two-fold change in MIC). Positive values indicate resistance, negative values susceptibility.",
                      `data-toggle` = "tooltip"
                    )
                  ),
                  sliderInput(
                    "cryptic_coef", NULL,
                    min = SLIDER_RANGES$cryptic_coef[1], max = SLIDER_RANGES$cryptic_coef[2],
                    value = SLIDER_RANGES$cryptic_coef,
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
                    min = SLIDER_RANGES$cryptic_counts[1], max = SLIDER_RANGES$cryptic_counts[2],
                    value = SLIDER_RANGES$cryptic_counts,
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
                    checkboxInput("rifreport_enabled", NULL, FALSE),
                    h4("RifReport")
                  ),
                  selectInput(
                    "rifreport_species", "Species",
                    choices = NULL,
                    selected = "All species",
                    multiple = TRUE,
                    width = "100%"
                  ),
                  div(
                    style = "display:flex; align-items:center; gap:6px;",
                    tags$label("Number of studies"),
                    tags$i(
                      class = "fa fa-question-circle info-icon",
                      title = "Number of independent published studies in which the mutation was reported to confer resistance. Higher values indicate stronger literature support.",
                      `data-toggle` = "tooltip"
                    )
                  ),
                  sliderInput(
                    "rifreport_nstudies", NULL,
                    min = SLIDER_RANGES$rifreport_nstudies[1], max = SLIDER_RANGES$rifreport_nstudies[2],
                    value = SLIDER_RANGES$rifreport_nstudies,
                    step = 1,
                    width = "100%"
                  ),
                  div(
                    style = "display:flex; align-items:center; gap:6px;",
                    tags$label("Number of species"),
                    tags$i(
                      class = "fa fa-question-circle info-icon",
                      title = "Number of bacterial species in which the mutation has been reported across the reviewed literature. Higher values indicate a mutation shared by more taxa.",
                      `data-toggle` = "tooltip"
                    )
                  ),
                  sliderInput(
                    "rifreport_nspecies", NULL,
                    min = SLIDER_RANGES$rifreport_nspecies[1], max = SLIDER_RANGES$rifreport_nspecies[2],
                    value = SLIDER_RANGES$rifreport_nspecies,
                    step = 1,
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
                    min = SLIDER_RANGES$dragdb_provean[1], max = SLIDER_RANGES$dragdb_provean[2],
                    value = SLIDER_RANGES$dragdb_provean,
                    step = 0.5,
                    width = "100%"
                  )
                ),
                div(
                  class = "filter-box",
                  div(
                    class = "filter-box-header",
                    checkboxInput("pointfinder_enabled", NULL, FALSE),
                    h4("PointFinder")
                  ),
                  selectInput(
                    "pointfinder_species", "Species",
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
                    checkboxInput("refgenes_enabled", NULL, FALSE),
                    h4("Reference Gene Catalog")
                  ),
                  selectInput(
                    "refgenes_species", "Species",
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
                  div(
                    style = "display:flex; align-items:center; gap:6px;",
                    tags$label("Effect size (log₂FC)"),
                    tags$i(
                      class = "fa fa-question-circle info-icon",
                      title = "Effect size is the log₂ fold-change in mutant frequency under rifampicin selection.",
                      `data-toggle` = "tooltip"
                    )
                  ),
                  sliderInput(
                    "yang_log2fc", NULL,
                    min = SLIDER_RANGES$yang_log2fc[1], max = SLIDER_RANGES$yang_log2fc[2],
                    value = SLIDER_RANGES$yang_log2fc,
                    step = 0.1,
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
                " clusters. The optional \"Show all mutations\" toggle replaces the cluster view with every mapped mutation across the whole gene."
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
                    checkboxInput(
                      inputId = "show_all_mutations",
                      label = "Show all mutations",
                      value = FALSE
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
                  width = 10,
                  div(
                    class = "plot-main",
                    style = "width: 100%; overflow-x: auto; border: 1px solid #eee;",
                    div(
                      style = "min-width: 1500px; margin-top: 40px",
                      ggiraph::girafeOutput(
                        "heatmap_pos_mut",
                        width = "100%",
                        height = "300px"
                      )
                    )
                  )
                ),
                column(
                  width = 2,
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

          ## ---- Dataset Comparison ----
          tabPanel(
            "Dataset Comparison",
            div(
              class = "text-muted-explore",
              p("Source-centric view of mutation evidence. Plots highlight overlap among datasets independent of species identity.")
            ),
            fluidRow(
              column(
                width = 12,
                div(
                  class = "plot-main",
                  div(
                    class = "plot-header",
                    h3("Cross-dataset mutation overlap")
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
                  div(
                    style = "margin-top: 120px;",
                    uiOutput("upset_container")
                  )
                )
              ),
              column(
                width = 12,
                div(
                  class = "plot-main",
                  div(
                    class = "plot-header plot-header-between",
                    div(
                      h3("Cross-dataset map of resistance mutations"),
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
              p("Summary table and histograms are generated from all rows in the filtered dataset; mutations reported by multiple sources may appear more than once. Use \"Download filtered data (CSV)\" to export these rows, with the same columns as the full mutation table on the Data tab.")
            ),
            div(
              style = "display:flex; gap:10px; justify-content:flex-end; margin-bottom:8px;",
              downloadButton(
                "download_filtered_mutations",
                "Download filtered data (CSV)",
                class = "btn btn-outline-primary btn-sm"
              )
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
          uiOutput("changelog_summary"),
          div(
            class = "table-header",
            h3("Mutation Data")
          ),
          p(
            class = "text-muted",
            tags$em(paste0("Version: ", APP_VERSION, " | Last updated: ", DATA_LAST_UPDATED))
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
          p(
            class = "text-muted",
            tags$em(paste0("Version: ", APP_VERSION, " | Last updated: ", DATA_LAST_UPDATED))
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
      uiOutput("docs_content")
    ),

    ## ---------------- Team ----------------
    tabPanel(
      title = tagList(icon("users"), "Team"),
      includeHTML("www/team.html")
    )
  )
)
