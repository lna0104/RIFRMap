library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(shiny)

# Load data
combined_data <- read.csv("data/processed_data/combined_data.csv", stringsAsFactors = FALSE)

# Keep fix a and y axis 
all_AAs_sorted <- str_split_1("HAKRGILPVFWYDESTCMNQ", "")
all_pos <- unique(combined_data$AA_pos_Ecoli) |> sort()

source_filters <- list(
  "CARD"    = c("species"),
  "CRyPTIC" = c("species","coef","se","z","raw_p","counts","ci_hi","ci_lo","bh_corrected_p","homoplasy"),
  "DRAGdb"  = c("species","provean"),
  "Negin"   = c("species","n_studies"),
  "Yang"    = c("species","log2_fc")
)

# shiny app -----

# ui ----
ui <- fluidPage(
  fluidRow(
    column(12, h2("RIFMap"))
  ),
  fluidRow(
    column(
      3,
      checkboxGroupInput(
        "sources",
        "Data sources",
        choices = sort(unique(combined_data$source)),
        selected = "CARD"
      ),
      uiOutput("filters_ui"),
      actionButton("submit", "Generate plot", class = "btn-primary")
    ),
    column(
      9,
      div(
        style = "height: 85vh; overflow-y: auto;",
        plotOutput("mutationPlot", height = "2000px", width = "100%")
      )
    )
  ),
  fluidRow(
    column(12, 
           verbatimTextOutput("debug_filtered")
    )
  )
)




# server -----
server <- function(input, output, session) {
  ## dynamic list of filters to show -----
  selected_filters <- reactive({
    req(input$sources)
    sort(unique(unlist(source_filters[input$sources])))
  })
  
  all_species_rv <- reactiveVal(character(0))
  
  ## dynamic filter UI -----
  output$filters_ui <- renderUI({
    f <- selected_filters()
    
    tagList(
      # shared
      if ("species" %in% f) {
        selectInput("species", "Species", choices = NULL, multiple = TRUE)},
      
      # Negin
      if ("n_studies" %in% f) {
        sliderInput("n_studies", "Number of studies", min = 0, max = 60, value = c(0, 60), step = 1)},
      
      # DRAGdb
      if ("provean" %in% f) {
        sliderInput("provean", "PROVEAN score (DRAGdb)",
                    min = -10.630, max = -2.698, value = c(-10.630, -2.698), step = 0.5)},
      
      # Yang
      if ("log2_fc" %in% f) {
        sliderInput("log2_fc", "log2 FC (Yang)",
                    min = 0.029, max = 2.700, value = c(0.029, 2.700), step = 0.1)},
      
      # CRyPTIC
      if ("coef" %in% f) {
        sliderInput("coef", "COEF (CRyPTIC)", min = -4.000, max = 9.600, value = c(-4.000, 9.600), step = 0.1)},
      
      if ("se" %in% f) {
        sliderInput("se", "SE (CRyPTIC)", min = 0.090, max = 3.000, value = c(0.090, 3.000), step = 0.1)},
      
      if ("z" %in% f) {
        sliderInput("z", "Z (CRyPTIC)", min = -3.650, max = 92.500, value = c(-3.650, 92.500), step = 1)},
      
      if ("raw_p" %in% f) {
        sliderInput("raw_p", "RAW P (CRyPTIC)", min = 0, max = 0.005, value = c(0, 0.005), step = 0.0001)},
      
      if ("counts" %in% f) {
        sliderInput("counts", "COUNTS (CRyPTIC)", min = 0, max = 2326, value = c(0, 2326), step = 1)},
      
      if ("ci_hi" %in% f) {
        sliderInput("ci_hi", "CI_HI (CRyPTIC)", min = -1.900, max = 15.500, value = c(-1.900, 15.500), step = 0.5)},
      
      if ("ci_lo" %in% f) {
        sliderInput("ci_lo", "CI_LO (CRyPTIC)", min = -6.080, max = 8.250, value = c(-6.080, 8.250), step = 0.5)},
      
      if ("bh_corrected_p" %in% f) {
        sliderInput("bh_corrected_p", "BH Corrected P (CRyPTIC)",
                    min = 0, max = 0.040, value = c(0, 0.040), step = 0.0001)},
      
      if ("homoplasy" %in% f) {
        checkboxInput("homoplasy", "HOMOPLASY only (CRyPTIC)", value = FALSE)}
    )
  })
  
  ## update species choices based on selected sources -----
  observeEvent(input$sources, {
    req(input$sources)
    
    all_species <- combined_data |>
      dplyr::filter(source %in% input$sources) |>
      dplyr::pull(species) |>
      unique() |>
      sort()
    
    all_species_rv(all_species)
    
    # only update if species UI exists for the currently selected sources
    req("species" %in% selected_filters())
    
    updateSelectInput(
      session,
      "species",
      choices  = c("All", all_species),
      selected = "All"
    )
  }, ignoreInit = FALSE)
  
  
  species_selected <- reactive({
    req(all_species_rv())
    if (is.null(input$species) || length(input$species) == 0 || "All" %in% input$species) {
      all_species_rv()
    } else {
      input$species
    }
  })
  
  
  ## filter each selected source and combine -----
  
  filtered <- eventReactive(input$submit, {
    req(input$sources)
    
    out <- lapply(input$sources, function(s) {
      df <- combined_data |> dplyr::filter(source == s)
      
      # shared species filter
      df <- df |> dplyr::filter(species %in% species_selected())
      
      # source-specific filters
      if (s == "DRAGdb" && !is.null(input$provean))
        df <- df |> dplyr::filter(dplyr::between(provean, input$provean[1], input$provean[2]))
      
      if (s == "Yang" && !is.null(input$log2_fc))
        df <- df |> dplyr::filter(dplyr::between(log2_fc, input$log2_fc[1], input$log2_fc[2]))
      
      if (s == "Negin" && !is.null(input$n_studies))
        df <- df |> dplyr::filter(dplyr::between(n_studies, input$n_studies[1], input$n_studies[2]))
      
      if (s == "CRyPTIC") {
        if (!is.null(input$coef))   df <- df |> dplyr::filter(dplyr::between(coef, input$coef[1], input$coef[2]))
        if (!is.null(input$se))     df <- df |> dplyr::filter(dplyr::between(se, input$se[1], input$se[2]))
        if (!is.null(input$z))      df <- df |> dplyr::filter(dplyr::between(z, input$z[1], input$z[2]))
        if (!is.null(input$raw_p))  df <- df |> dplyr::filter(dplyr::between(raw_p, input$raw_p[1], input$raw_p[2]))
        if (!is.null(input$counts)) df <- df |> dplyr::filter(dplyr::between(counts, input$counts[1], input$counts[2]))
        if (!is.null(input$ci_hi))  df <- df |> dplyr::filter(dplyr::between(ci_hi, input$ci_hi[1], input$ci_hi[2]))
        if (!is.null(input$ci_lo))  df <- df |> dplyr::filter(dplyr::between(ci_lo, input$ci_lo[1], input$ci_lo[2]))
        if (!is.null(input$bh_corrected_p))
          df <- df |> dplyr::filter(dplyr::between(bh_corrected_p, input$bh_corrected_p[1], input$bh_corrected_p[2]))
        if (isTRUE(input$homoplasy))
          df <- df |> dplyr::filter(homoplasy)
      }
      
      df
    })
    
    dplyr::bind_rows(out)
  }, ignoreInit = TRUE)
  
  
  
  
  ## build fixed-grid presence data for the heatmap -----
  plotting_df <- reactive({
    df <- filtered()
    req(nrow(df) > 0)
    
    base_grid <- tidyr::expand_grid(
      AA_pos_Ecoli = all_pos,
      AA_mutation = all_AAs_sorted
    )
    
    obs <- df |>
      dplyr::distinct(AA_pos_Ecoli, AA_mutation, source) |>
      dplyr::group_by(AA_pos_Ecoli, AA_mutation) |>
      dplyr::summarise(
        source_combo = paste(sort(unique(source)), collapse = "+"),
        .groups = "drop"
      )
    
    base_grid |>
      dplyr::left_join(obs, by = c("AA_pos_Ecoli", "AA_mutation")) |>
      dplyr::mutate(
        source_combo = tidyr::replace_na(source_combo, "None"),
        AA_mutation  = factor(AA_mutation, levels = all_AAs_sorted),
        AA_pos_Ecoli = factor(AA_pos_Ecoli, levels = all_pos)
      )
  })
  
  ## plot -----
  source_palette <- c(
    "CARD"    = "#F4A7A3",  
    "CRyPTIC" = "#8ECEDB",  
    "DRAGdb"  = "#B3B7E6",  
    "Negin"   = "#E6B3E6",  
    "Yang"    = "#F6B3C6",  
    
    "CARD+CRyPTIC"  = "#E8B97D",
    "CARD+DRAGdb"   = "#C4D88A",
    "CARD+Negin"    = "#9ED6B8",
    "CARD+Yang"     = "#D8A7C4",
    
    "CRyPTIC+Negin" = "#9FC6E8",
    "CRyPTIC+Yang"  = "#B6C6F0",
    "DRAGdb+Negin"  = "#C6B6EB",
    "DRAGdb+Yang"   = "#D6B6EB",
    "Negin+Yang"    = "#F2B6D8",
    
    "CARD+CRyPTIC+Negin" = "#D6D88A",
    "CARD+CRyPTIC+Yang"  = "#E0C6A6",
    "CARD+DRAGdb+Negin"  = "#A6D6B8",
    "CARD+DRAGdb+Yang"   = "#C6D6D6",
    "CARD+Negin+Yang"    = "#E8B6D6",
    "CRyPTIC+Negin+Yang" = "#B6C6E8",
    
    "CARD+CRyPTIC+DRAGdb+Negin"      = "#D6C68A",
    "CARD+CRyPTIC+Negin+Yang"        = "#E6D0A6",
    "CARD+CRyPTIC+DRAGdb+Negin+Yang" = "#E6C6A6",
    
    "None" = "#D9D9D9"   
  )
  
  output$mutationPlot <- renderPlot({
    req(input$submit > 0)
    ggplot(plotting_df()) +
      geom_tile(
        aes(x = AA_mutation, y = AA_pos_Ecoli, fill = source_combo),
        color = "grey70",
        linewidth = 0.2
      ) +
      scale_fill_manual(
        values = source_palette,
        drop = FALSE
      ) +
      labs(x = "Amino acid", y = "Position", fill = "Source(s)") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(size = 6))
  }, res = 96)
  
  output$debug_filtered <- renderPrint({
    df <- filtered()
    req(nrow(df) > 0)
    summary(df)
  })
  
}

shinyApp(ui, server)



