function(input, output, session) {
  # ---- Go to Documentation tab ----
  observeEvent(input$go_docs, {
    updateNavbarPage(session, "main_nav", selected = "Documentation")
  })

  # ---- Go to Explorer tab ----
  observeEvent(input$go_explorer, {
    updateNavbarPage(session, "main_nav", selected = "Explorer")
  })

  ## ---- Populate species dropdowns ----
  observeEvent(TRUE,
    {
      updateSelectInput(
        session, "card_species",
        choices = c("All species", sort(unique(card$species))),
        selected = "All species"
      )

      updateSelectInput(
        session, "cryptic_species",
        choices = c(sort(unique(cryptic$species))),
        selected = "Mycobacterium tuberculosis"
      )

      updateSelectInput(
        session, "dragdb_species",
        choices = c("All species", sort(unique(dragdb$species))),
        selected = "All species"
      )

      updateSelectInput(
        session, "yang_species",
        choices = c(sort(unique(yang$species))),
        selected = "Escherichia coli"
      )

      updateSelectInput(
        session, "curated_species",
        choices = c("All species", sort(unique(curated$species))),
        selected = "All species"
      )
    },
    once = TRUE
  )

  ## ---- Combine selected datasets ----
  combined_rows <- reactive({
    dfs <- list()

    ## CARD
    if (isTRUE(input$card_enabled)) {
      d <- card |>
        filter_by_species(input$card_species, species_col = "species")

      dfs[["CARD"]] <- d
    }

    ## CRyPTIC
    if (isTRUE(input$cryptic_enabled)) {
      d <- cryptic |>
        filter_by_species(input$cryptic_species, species_col = "species")

      # numeric filters
      d <- filter_range_if(d, "coef", input$cryptic_coef)
      d <- filter_range_if(d, "counts", input$cryptic_counts)
      d <- d |>
        filter(bh_corrected_p <= as.numeric(input$cryptic_bhp))

      # categorical filter
      d <- filter_in_if(d, "homoplasy", input$cryptic_homoplasy)

      dfs[["CRyPTIC"]] <- d
    }

    ## DRAGdb
    if (isTRUE(input$dragdb_enabled)) {
      d <- dragdb |>
        filter_by_species(input$dragdb_species, species_col = "species")

      d <- filter_range_if(d, "provean", input$dragdb_provean)

      dfs[["DRAGdb"]] <- d
    }

    ## Yang
    if (isTRUE(input$yang_enabled)) {
      d <- yang |>
        filter_by_species(input$yang_species, species_col = "species")

      d <- filter_range_if(d, "log2_fc", input$yang_log2fc)

      dfs[["Yang"]] <- d
    }

    ## Curated data
    if (isTRUE(input$curated_enabled)) {
      d <- curated |>
        filter_by_species(input$curated_species, species_col = "species")

      d <- filter_range_if(d, "n_studies", input$curated_nstudies)
      d <- filter_range_if(d, "n_species", input$curated_nspecies)

      dfs[["Curated"]] <- d
    }

    if (length(dfs) == 0) {
      return(NULL)
    }

    bind_rows(dfs, .id = "source")
  })


  # ---- Redirect to View tab on Apply ----
  observeEvent(input$apply_filters, {
    updateTabItems(session, "sidebar", "explorer")
    updateTabsetPanel(session, "explore_steps", selected = "cluster")
  })

  # ---- Reset button: clear inputs + go back Home
  observeEvent(input$reset_filters, {
    updateCheckboxInput(session, "card_enabled", value = FALSE)
    updateCheckboxInput(session, "cryptic_enabled", value = FALSE)
    updateCheckboxInput(session, "dragdb_enabled", value = FALSE)
    updateCheckboxInput(session, "yang_enabled", value = FALSE)
    updateCheckboxInput(session, "curated_enabled", value = FALSE)

    updateSelectInput(session, "card_species", selected = "All species")
    updateSelectInput(session, "cryptic_species", selected = "Mycobacterium tuberculosis")
    updateSelectInput(session, "dragdb_species", selected = "All species")
    updateSelectInput(session, "yang_species", selected = "Escherichia coli")

    updateSliderInput(session, "cryptic_coef", value = c(-4, 9.6))
    updateSliderInput(session, "cryptic_counts", value = c(0, 2326))
    updateRadioButtons(session, "cryptic_bhp", selected = 0.05)
    updateSelectInput(session, "cryptic_homoplasy", selected = "")

    updateSliderInput(session, "dragdb_provean", value = c(-10.7, -2.6))
    updateSliderInput(session, "yang_log2fc", value = c(0, 2.7))
    updateSliderInput(session, "curated_nstudies", value = c(1, 59))
    updateSliderInput(session, "curated_nspecies", value = c(1, 29))

    updateTabItems(session, "sidebar", "filters")
  })


  # ----- Species-Cluster Tab -----
  {
    # ---- Process data ----

    # cluster selection check
    no_clusters_selected <- function(x) {
      is.null(x) || !length(x) || all(is.na(x)) || all(x == "")
    }

    # Clean and normalize plotting df
    cluster_plotting_df <- reactive({
      df0 <- combined_rows()
      if (is.null(df0) || !nrow(df0)) {
        return(NULL)
      }

      df <- df0 |>
        transmute(
          Species      = na_if(trimws(as.character(species)), ""),
          AA_original  = na_if(trimws(as.character(AA_original)), ""),
          AA_mutation  = na_if(trimws(as.character(AA_mutation)), ""),
          AA_pos       = suppressWarnings(as.integer(AA_pos)),
          AA_pos_Ecoli = suppressWarnings(as.integer(AA_pos_Ecoli))
        ) |>
        filter(
          !is.na(Species),
          !is.na(AA_pos_Ecoli),
          !is.na(AA_original),
          !is.na(AA_mutation)
        ) |>
        distinct(Species, AA_original, AA_mutation, AA_pos, AA_pos_Ecoli)

      if (!nrow(df)) {
        return(NULL)
      }
      df
    })

    # Precompute which species have all reference assets
    ref_assets_ok <- local({
      # Named logical: species -> TRUE/FALSE
      cache <- NULL

      function(refs, seqs, coordinates) {
        if (!is.null(cache)) {
          return(cache)
        }

        # only keep rows with a valid FASTA_name
        rr <- refs |>
          transmute(
            Species = as.character(.data$Species),
            FASTA   = na_if(trimws(as.character(.data$FASTA_name)), "")
          ) |>
          filter(!is.na(.data$Species), !is.na(.data$FASTA)) |>
          distinct(.data$Species, .keep_all = TRUE)

        ok <- setNames(rep(FALSE, nrow(rr)), rr$Species)

        if (nrow(rr)) {
          has_seq <- rr$FASTA %in% names(seqs)
          has_coord <- vapply(rr$FASTA, function(id) !is.null(coordinates[[id]]), logical(1))
          ok[] <- has_seq & has_coord
        }

        cache <<- ok
        cache
      }
    })

    # Species to plot (drop those without assets) + warn only when list changes
    cluster_species_to_plot <- local({
      last_dropped <- character(0)

      reactive({
        df <- cluster_plotting_df()
        if (is.null(df) || !nrow(df)) {
          return(character(0))
        }

        sp <- sort(unique(df$Species))

        ok_map <- ref_assets_ok(refs = refs_seqs, seqs = seqs, coordinates = coordinates)
        keep <- !is.na(ok_map[sp]) & ok_map[sp]

        sp_keep <- sp[keep]
        dropped <- setdiff(sp, sp_keep)

        # warn only if dropped set changed
        if (length(dropped) && !identical(sort(dropped), sort(last_dropped))) {
          warning(
            "Dropping species missing reference assets (refs/seqs/coordinates):\n",
            paste0("  - ", dropped, collapse = "\n")
          )
          last_dropped <<- dropped
        }
        if (!length(dropped)) last_dropped <<- character(0)

        sp_keep
      })
    })

    # Canonical ordered clusters (compute once per selection)
    clusters_to_plot_r <- reactive({
      if (no_clusters_selected(input$clustersToPlot)) {
        return(character(0))
      }

      clusters |>
        filter(.data$Gene == "rpoB", .data$Cluster %in% input$clustersToPlot) |>
        arrange(.data$AAstart_Ecoli) |>
        pull(.data$Cluster) |>
        unique()
    })

    # Cache plot dimensions so we don’t call plotMutations twice
    cluster_dims_r <- reactiveVal(NULL)

    observeEvent(
      list(input$apply_filters, input$clustersToPlot, input$colour_by),
      {
        df <- cluster_plotting_df()
        sp <- cluster_species_to_plot()
        cl <- clusters_to_plot_r()

        if (is.null(df) || !nrow(df) || !length(sp) || !length(cl)) {
          cluster_dims_r(NULL)
          return()
        }

        cluster_dims_r(
          plotMutations(
            muts            = df,
            clusters        = clusters,
            seqs            = seqs,
            coordinates     = coordinates,
            refs            = refs_seqs,
            gene            = "rpoB",
            clustersToPlot  = cl,
            speciesToPlot   = sp,
            colour_by       = input$colour_by,
            mergeDuplicates = TRUE,
            return_dims     = TRUE
          )
        )
      },
      ignoreInit = TRUE
    )

    # Trigger for plot redraw
    cluster_trigger <- reactive({
      input$apply_filters
      input$clustersToPlot
      input$colour_by
      TRUE
    })


    # ---- Plot ----

    # Store hover bands + plt for JS
    hover_payload_r <- reactiveVal(NULL)

    output$clusterPlotUI <- renderUI({
      d <- cluster_dims_r()
      w <- if (is.null(d)) 300 else max(300, d$width_px)
      div(
        id = "clusterWrap",
        div(
          id = "clusterPlotWrap",
          spin(plotOutput("clusterPlot", width = paste0(w, "px"), height = "85vh")),
          tags$canvas(id = "hoverCanvas"),
          uiOutput("clusterPlotMsg")
        )
      )
    })

    # Data check and message
    cluster_status <- reactive({
      cluster_trigger()

      df <- cluster_plotting_df()
      cl <- clusters_to_plot_r()
      sp <- cluster_species_to_plot()

      if (is.null(df) || nrow(df) == 0) {
        return(list(ok = FALSE, msg = "No data to show."))
      }
      if (is.null(cl) || length(cl) == 0) {
        return(list(ok = FALSE, msg = "Please select at least one cluster."))
      }
      if (is.null(sp) || length(sp) == 0) {
        return(list(ok = FALSE, msg = "No selectable species (missing reference assets)."))
      }

      list(ok = TRUE, msg = NULL)
    })

    output$clusterPlotMsg <- renderUI({
      st <- cluster_status()
      if (isTRUE(st$ok)) {
        return("NULL")
      }
      div(class = "plot-empty-msg", st$msg)
    })

    output$clusterPlot <- renderPlot(
      {
        st <- cluster_status()
        if (!isTRUE(st$ok)) {
          return(invisible())
        }

        df <- cluster_plotting_df()
        cl <- clusters_to_plot_r()
        sp <- cluster_species_to_plot()

        out <- plotMutations(
          muts             = df,
          clusters         = clusters,
          seqs             = seqs,
          coordinates      = coordinates,
          refs             = refs_seqs,
          gene             = "rpoB",
          clustersToPlot   = cl,
          speciesToPlot    = sp,
          colour_by        = input$colour_by,
          mergeDuplicates  = TRUE,
          return_hover_map = TRUE
        )

        # capture geometry from this plot
        payload <- list(
          bands = out$bands,
          plt = out$plt,
          xcols = out$xcols,
          dx = out$dx
        )

        hover_payload_r(payload)

        invisible()
      },
      width = function() {
        d <- cluster_dims_r()
        if (is.null(d)) 300 else max(300, d$width_px)
      },
      height = function() {
        d <- cluster_dims_r()
        if (is.null(d)) 300 else max(150, d$height_px)
      },
      res = 96
    )

    # Send compact bands to JS
    observeEvent(hover_payload_r(),
      {
        p <- hover_payload_r()

        req(p$bands, p$plt, p$xcols)

        session$onFlushed(function() {
          session$sendCustomMessage("hoverMap", p)
        }, once = TRUE)
      },
      ignoreInit = TRUE
    )
  }

  # ---- Legend ----

  output$legend_panel <- renderUI({
    mode <- input$colour_by

    if (is.null(mode)) {
      return(NULL)
    }

    pal <- COLOURING_PALETTES[[mode]]

    # none
    if (mode == "none") {
      return(
        tags$div(
          class = "cluster-legend",
          tags$div(class = "cluster-legend-title", "Legend"),
          legend_item("All mutations", pal["none"])
        )
      )
    }

    # distribution
    if (mode == "distribution") {
      cols <- COLOURING_PALETTES$distribution
      gradient_css <- paste(cols, collapse = ", ")

      rng <- cluster_dims_r()$distribution_range

      return(
        tags$div(
          class = "cluster-legend",
          tags$div(class = "cluster-legend-title", "Distribution across species"),
          tags$div(
            style = sprintf(
              "height:14px; border-radius:6px; margin-bottom:6px;
           background: linear-gradient(to right, %s);",
              gradient_css
            )
          ),
          tags$div(
            style = "display:flex; justify-content:space-between; font-size:0.8em;",
            tags$span(paste0(rng["min"])),
            tags$span(paste0(rng["max"]))
          )
        )
      )
    }

    # aa_type
    if (mode == "aa_type") {
      labels <- c(
        hydrophobic_aliphatic = "Hydrophobic: Aliphatic",
        hydrophobic_aromatic  = "Hydrophobic: Aromatic",
        polar_uncharged       = "Hydrophilic: Polar Uncharged",
        acidic                = "Hydrophilic: Acidic",
        basic                 = "Hydrophilic: Basic",
        unknown               = "Unknown"
      )

      return(
        tags$div(
          class = "cluster-legend",
          tags$div(class = "cluster-legend-title", "Amino acid type"),
          lapply(names(labels), function(nm) {
            legend_item(labels[[nm]], pal[[nm]])
          })
        )
      )
    }


    # aa_distance
    if (mode == "aa_distance") {
      cols <- COLOURING_PALETTES$aa_distance
      gradient_css <- paste(cols, collapse = ", ")
      rng <- cluster_dims_r()$aa_distance_range
      return(
        tags$div(
          class = "cluster-legend",
          tags$div(class = "cluster-legend-title", "Amino acid distance"),
          tags$div(
            style = sprintf(
              "height:14px; border-radius:6px; margin-bottom:6px;
     background: linear-gradient(to right, %s);",
              gradient_css
            )
          ),
          tags$div(
            style = "display:flex; justify-content:space-between; font-size:0.8em;",
            tags$span(round(rng["min"], 1)),
            tags$span(round(rng["max"], 1))
          )
        )
      )
    }

    NULL
  })

  # plot note

  output$clusterPlotNote <- renderUI({
    mode_note <- switch(input$colour_by,
      none = "",
      distribution = "Colours indicate how broadly each mutation is distributed across species.",
      aa_type = "Colours indicate the biochemical class of the mutated amino acid.",
      aa_distance = "Colours indicate Grantham distance, which measures the evolutionary distance between the original and mutated amino acids."
    )

    tags$div(
      class = "text-muted",
      style = "font-size: 0.85em; margin-bottom: 6px;",
      mode_note
    )
  })

  # ---- Mutation Frequency Tab ----

  output$heatmap_pos_mut <- ggiraph::renderGirafe({
    df0 <- combined_rows()
    req(df0, nrow(df0) > 0)

    df <- df0 |>
      transmute(
        AA_pos_Ecoli = as.character(AA_pos_Ecoli),
        AA_mutation = as.character(AA_mutation),
        Species = as.character(species)
      ) |>
      filter(
        !is.na(AA_pos_Ecoli),
        !is.na(AA_mutation),
        !is.na(Species)
      ) |>
      distinct()

    all_AAs_sorted <- sort(strsplit("HAKRGILPVFWYDESTCMNQ", "")[[1]])

    heat_df <- df |>
      filter(!is.na(AA_pos_Ecoli), !is.na(AA_mutation), !is.na(Species)) |>
      group_by(AA_pos_Ecoli, AA_mutation) |>
      summarise(
        n_species = n_distinct(Species),
        .groups = "drop"
      ) |>
      tidyr::complete(
        AA_pos_Ecoli,
        AA_mutation = all_AAs_sorted,
        fill = list(n_species = 0)
      )

    pos_levels <- heat_df$AA_pos_Ecoli |>
      unique() |>
      as.numeric() |>
      sort(decreasing = TRUE) |>
      as.character()

    n_pos <- length(pos_levels)
    bar_h_cm <- max(4, min(10, 0.15 * n_pos))

    heat_df <- heat_df |>
      mutate(
        AA_mutation = factor(AA_mutation, levels = all_AAs_sorted),
        AA_pos_Ecoli = factor(AA_pos_Ecoli, levels = pos_levels),
        tooltip_txt = paste0(
          "<b>Position:</b> ", AA_pos_Ecoli,
          "<br><b>AA:</b> ", AA_mutation,
          "<br><b>No. species:</b> ", n_species,
          "<br><i>Click tile to view species list</i>"
        ),
        tile_id = paste(AA_pos_Ecoli, AA_mutation)
      )

    max_n <- max(heat_df$n_species, na.rm = TRUE)
    if (is.na(max_n) || max_n < 0) max_n <- 0

    cols <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

    # choose fill mapping
    if (max_n < 4) {
      levs <- 0:max_n
      heat_df <- heat_df |> mutate(fill_val = factor(n_species, levels = levs))

      box_cols <- c("white", cols[seq(10, 90, length.out = max(1, max_n))])
      names(box_cols) <- as.character(levs)

      p <- ggplot(
        heat_df,
        aes(
          x = AA_mutation,
          y = AA_pos_Ecoli,
          fill = fill_val,
          tooltip = tooltip_txt,
          data_id = tile_id
        )
      ) +
        ggiraph::geom_tile_interactive(color = "grey80", linewidth = 0.2) +
        scale_x_discrete(position = "top", drop = FALSE) +
        scale_fill_manual(
          values = box_cols,
          breaks = as.character(levs),
          labels = as.character(levs),
          drop = FALSE
        ) +
        guides(
          fill = guide_legend(
            title = "No. species",
            keyheight = grid::unit(0.5, "cm"),
            keywidth = grid::unit(0.5, "cm")
          )
        )
    } else {
      if (max_n <= 0) max_n <- 1
      heat_df <- heat_df |> mutate(fill_val = n_species)

      vals <- c(0, seq(1, max_n, length.out = length(cols)))
      vals_scaled <- vals / max_n

      p <- ggplot(
        heat_df,
        aes(
          x = AA_mutation,
          y = AA_pos_Ecoli,
          fill = fill_val,
          tooltip = tooltip_txt,
          data_id = tile_id
        )
      ) +
        ggiraph::geom_tile_interactive(color = "grey80", linewidth = 0.2) +
        scale_x_discrete(position = "top", drop = FALSE) +
        scale_fill_gradientn(
          colours = c("white", cols),
          values  = vals_scaled,
          limits  = c(0, max_n),
          breaks  = c(0, 1, max_n),
          labels  = as.character(c(0, 1, max_n)),
          oob     = scales::squish
        ) +
        guides(
          fill = guide_colourbar(
            barheight = grid::unit(bar_h_cm, "cm"),
            barwidth  = grid::unit(0.6, "cm")
          )
        )
    }

    # now add common styling
    p <- p +
      labs(x = "Amino acid", y = "Position (E. coli coordinate)", fill = "No. species") +
      plotting_theme_large_axes(legend = "right", x_top = TRUE)


    # height scaling
    h_in <- max(5, min(35, 0.35 * n_pos))

    ggiraph::girafe(
      ggobj = p,
      width_svg = 12,
      height_svg = h_in,
      options = list(
        ggiraph::opts_hover(css = "fill-opacity:1; stroke:black; stroke-width:0.6px;"),
        ggiraph::opts_tooltip(css = "
        background-color: #2C3E50;
        color: #ECF0F1;
        padding: 6px;
        border-radius: 3px;
        font-family: Helvetica, sans;
        font-size: 11px;
        box-shadow: 0 0 10px rgba(0,0,0,0.5);
      "),
        ggiraph::opts_zoom(max = 3),
        ggiraph::opts_selection(type = "single", only_shiny = TRUE, css = "stroke:#ff3b30; stroke-width:2px;"),
        ggiraph::opts_sizing(rescale = TRUE)
      )
    ) |>
      ggiraph::girafe_options(
        opts_toolbar(position = "topright", saveaspng = FALSE)
      )
  })

  # In the UI, my output id is "heatmap_pos_mut" => Shiny will create an input input$heatmap_pos_mut_selected

  # Clicked tile -> species list

  selected_tile_species <- reactive({
    df <- cluster_plotting_df()
    req(df, nrow(df) > 0)

    sel_id <- input$heatmap_pos_mut_selected
    req(!is.null(sel_id), nzchar(sel_id))

    parts <- strsplit(sel_id, " ")[[1]]
    req(length(parts) >= 2)

    pos <- suppressWarnings(as.numeric(parts[1]))
    aa <- parts[2]
    req(!is.na(pos), nzchar(aa))

    spp <- df |>
      filter(AA_pos_Ecoli == pos, AA_mutation == aa) |>
      distinct(Species) |>
      arrange(Species)

    list(pos = pos, aa = aa, species = spp)
  })


  output$selected_tile_info <- renderUI({
    x <- selected_tile_species()
    tags$div(
      style = "margin-bottom:8px; font-weight:600;",
      sprintf(
        "AA_Pos_Ecoli = %s; Mutation = %s (%d species)",
        x$pos, x$aa, nrow(x$species)
      )
    )
  })


  output$selected_tile_species_table <- DT::renderDT({
    x <- selected_tile_species()
    spp <- x$species

    DT::datatable(
      spp,
      rownames = FALSE,
      options = list(
        pageLength = 10,
        scrollY = "260px",
        dom = "l tip"
      )
    )
  })


  # ---- Source Overlap Tab ----
  ## ---- Upset Plot ----
  {
    n_sources <- reactive({
      df <- combined_rows()
      if (is.null(df) || nrow(df) == 0) {
        return(0)
      }
      length(unique(na.omit(df$source)))
    })

    output$mutationUpset <- renderPlot(
      {
        req(input$apply_filters > 0)

        df0 <- combined_rows()

        if (is.null(df0) || nrow(df0) == 0) {
          return(ui_message("No data to show."))
        }

        df <- df0 |>
          transmute(
            AA_pos_Ecoli = as.character(AA_pos_Ecoli),
            AA_mutation  = as.character(AA_mutation),
            source       = as.character(source)
          ) |>
          filter(
            !is.na(AA_pos_Ecoli),
            !is.na(AA_mutation),
            !is.na(source)
          ) |>
          distinct()

        if (is.null(df) || nrow(df) == 0) {
          return(ui_message("No data to show."))
        }

        upset_df <- df |>
          distinct(AA_mutation, AA_pos_Ecoli, source) |>
          group_by(AA_mutation, AA_pos_Ecoli) |>
          summarise(
            sources = list(sort(unique(source))),
            .groups = "drop"
          )

        if (is.null(upset_df) || nrow(upset_df) == 0) {
          return(ui_message("No data to show."))
        }

        ggplot(upset_df, aes(x = sources)) +
          geom_bar() +
          ggupset::scale_x_upset() +
          labs(y = "Number of mutations", x = "Source(s)") +
          plotting_theme(legend = "none")
      },
      res = 120
    )


    output$upset_container <- renderUI({
      df <- combined_rows()

      if (is.null(df) || nrow(df) == 0) {
        return(ui_message("No data to show."))
      } else {
        ns <- n_sources()
        w <- ns * 300
        div(
          style = "
            display: flex;
            justify-content: center;
            overflow-x: auto;
            padding: 10px 0;
          ",
          div(
            style = sprintf("width: %dpx;", w),
            spin(plotOutput("mutationUpset", height = "450px", width = "100%"))
          )
        )
      }
    })
  }

  ## ---- Mutation Heatmap ----

  {
    selected_sources <- reactiveVal(NULL)

    ### ---- Helpers ----
    safe_id <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

    notify_reselect <- function(where = "heatmap") {
      shiny::showNotification("Please select the legend again", type = "message", duration = 6)
    }

    get_clickable_levels <- function(df) {
      present <- unique(as.character(df$source))
      present <- present[!is.na(present)]
      present_levels <- source_levels[source_levels %in% present]
      if (length(present_levels) == 0) present_levels <- sort(present)
      setdiff(present_levels, "None")
    }

    sync_selected_sources <- function(clickable_levels) {
      sel <- selected_sources()
      if (is.null(sel)) sel <- setNames(logical(0), character(0))

      sel2 <- setNames(rep(FALSE, length(clickable_levels)), clickable_levels)
      keep <- intersect(names(sel), clickable_levels)
      sel2[keep] <- sel[keep]

      selected_sources(sel2)
      sel2
    }

    ### ---- Positions ----
    selected_pos <- reactive({
      df0 <- combined_rows()
      if (is.null(df0) || nrow(df0) == 0) {
        return(character(0))
      }

      df0$AA_pos_Ecoli |>
        suppressWarnings(as.numeric()) |>
        na.omit() |>
        unique() |>
        sort(decreasing = TRUE) |>
        as.character()
    })

    ### ---- Plot data ----
    hm_plotting_df <- reactive({
      df0 <- combined_rows()
      if (is.null(df0) || nrow(df0) == 0) {
        return(NULL)
      }

      pos_levels <- selected_pos()
      if (!length(pos_levels)) {
        return(NULL)
      }

      df <- df0 |>
        transmute(
          AA_pos_Ecoli = as.character(AA_pos_Ecoli),
          AA_mutation  = as.character(AA_mutation),
          source       = as.character(source)
        ) |>
        filter(!is.na(AA_pos_Ecoli), !is.na(AA_mutation), !is.na(source)) |>
        distinct()

      df
    })

    clickable_levels_r <- reactive({
      df <- hm_plotting_df()
      if (is.null(df) || !nrow(df)) {
        return(character(0))
      }
      get_clickable_levels(df)
    })

    ### ---- Dynamic height ----
    heatmap_height_px <- reactive({
      pos_levels <- selected_pos()
      if (!length(pos_levels)) {
        return(300)
      }
      max(300, length(pos_levels) * 30)
    })

    ### ---- Legend UI ----
    output$heatmap_legend <- renderUI({
      req(input$apply_filters > 0)
      df <- hm_plotting_df()
      req(!is.null(df) && nrow(df) > 0)

      clickable_levels <- get_clickable_levels(df)
      sel <- sync_selected_sources(clickable_levels)

      tags$div(
        class = "hm-legend-col",
        lapply(clickable_levels, function(lvl) {
          col <- source_palette[[lvl]]
          if (is.null(col) || is.na(col) || !nzchar(col)) col <- "#D9D9D9"

          is_on <- isTRUE(sel[[lvl]])

          actionButton(
            inputId = paste0("legend__", safe_id(lvl)),
            label = tags$span(
              class = "hm-legend-btn-label",
              tags$span(class = "hm-legend-swatch", style = sprintf("--swatch:%s;", col)),
              tags$span(class = "hm-legend-text", lvl)
            ),
            class = paste("hm-legend-btn", if (is_on) "is-on" else "is-off"),
            style = sprintf("--btncol:%s;", col)
          )
        })
      )
    })

    ### ---- Legend resets ----
    observeEvent(input$apply_filters,
      {
        selected_sources(NULL)
      },
      ignoreInit = TRUE
    )

    observeEvent(input$reset_heatmap_legend,
      {
        selected_sources(NULL)
      },
      ignoreInit = TRUE
    )

    ### ---- Legend clicks ----
    legend_observers_created <- reactiveVal(character(0))

    observe({
      req(input$apply_filters > 0)
      df <- hm_plotting_df()
      if (is.null(df) || nrow(df) == 0) {
        return()
      }

      clickable_levels <- get_clickable_levels(df)
      sync_selected_sources(clickable_levels)

      ids <- paste0("legend__", safe_id(clickable_levels))
      created <- legend_observers_created()
      new_ids <- setdiff(ids, created)
      if (!length(new_ids)) {
        return()
      }

      legend_observers_created(c(created, new_ids))

      for (j in seq_along(clickable_levels)) {
        id <- ids[j]
        lvl <- clickable_levels[j]
        if (!(id %in% new_ids)) next

        local({
          lvl_local <- lvl
          id_local <- id
          observeEvent(input[[id_local]],
            {
              sel <- selected_sources()
              if (is.null(sel) || !(lvl_local %in% names(sel))) {
                return()
              }
              sel[[lvl_local]] <- !isTRUE(sel[[lvl_local]])
              selected_sources(sel)
            },
            ignoreInit = TRUE
          )
        })
      }
    })


    ### ---- Heatmap ----
    output$mutationHeatmap <- ggiraph::renderGirafe({
      req(input$apply_filters > 0)

      tryCatch(
        {
          df <- hm_plotting_df()
          validate(need(!is.null(df) && nrow(df) > 0, "No data to show."))

          pos_levels <- selected_pos()

          clickable_levels <- get_clickable_levels(df)
          sel <- sync_selected_sources(clickable_levels)

          keep_levels <- names(sel)[sel]

          if (length(keep_levels) == 0) {
            keep_levels <- clickable_levels
          }

          df <- df %>%
            filter(source %in% keep_levels)

          base_grid <- tidyr::expand_grid(
            AA_pos_Ecoli = pos_levels,
            AA_mutation  = all_AAs_sorted
          )

          df_fill <- base_grid |>
            left_join(df, by = c("AA_pos_Ecoli", "AA_mutation")) |>
            group_by(AA_pos_Ecoli, AA_mutation) |>
            mutate(source = if (all(is.na(source))) "None" else source) |>
            ungroup() |>
            filter(!is.na(source)) |>
            mutate(
              AA_mutation  = factor(AA_mutation, levels = all_AAs_sorted),
              AA_pos_Ecoli = factor(AA_pos_Ecoli, levels = pos_levels),
              source       = factor(source, levels = source_levels)
            )

          # vertical stripes inside each tile
          df_strip <- df_fill |>
            mutate(
              x_fac = factor(as.character(AA_mutation), levels = all_AAs_sorted),
              y_fac = factor(as.character(AA_pos_Ecoli), levels = pos_levels),
              x_num = as.numeric(x_fac),
              y_num = as.numeric(y_fac)
            ) |>
            group_by(x_fac, y_fac) |>
            arrange(source, .by_group = TRUE) |>
            mutate(
              all_sources = paste(sort(unique(as.character(source))), collapse = ", "),
              k = n(),
              i = row_number(),
              ymin = y_num - 0.5,
              ymax = y_num + 0.5,
              xmin = (x_num - 0.5) + (i - 1) / k,
              xmax = (x_num - 0.5) + i / k
            ) |>
            ungroup()

          # strong borders between tiles only
          df_tile_border <- df_strip |>
            distinct(x_num, y_num) |>
            mutate(
              xmin = x_num - 0.5, xmax = x_num + 0.5,
              ymin = y_num - 0.5, ymax = y_num + 0.5
            )

          p <- ggplot(df_strip) +
            ggiraph::geom_rect_interactive(
              aes(
                xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax,
                fill = source,
                tooltip = paste0(
                  "Position: ", AA_pos_Ecoli,
                  "<br>AA: ", AA_mutation,
                  "<br>Sources: ", all_sources
                ),
                data_id = paste(AA_pos_Ecoli, AA_mutation, sep = "|")
              ),
              color = "grey85", linewidth = 0.2
            ) +
            geom_rect(
              data = df_tile_border,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE,
              fill = NA, color = "grey60", linewidth = 0.3
            ) +
            scale_fill_manual(values = source_palette) +
            scale_x_continuous(
              breaks = seq_along(all_AAs_sorted),
              labels = all_AAs_sorted,
              position = "top",
              expand = c(0, 0)
            ) +
            scale_y_continuous(
              breaks = seq_along(pos_levels),
              labels = pos_levels,
              expand = c(0, 0)
            ) +
            labs(x = "Amino acid", y = "Position (E. coli coordinate)") +
            plotting_theme(legend = "none")

          ggiraph::girafe(
            ggobj = p,
            width_svg = 12,
            height_svg = heatmap_height_px() / 96,
            options = list(
              ggiraph::opts_hover(css = "fill-opacity:1;"),
              ggiraph::opts_tooltip(css = "
            background-color: #2C3E50;
            color: #ECF0F1;
            padding: 5px;
            border-radius: 2px;
            font-family: Helvetica, sans;
            font-size: 10px;
            box-shadow: 0 0 10px rgba(0,0,0,0.5);
          "),
              ggiraph::opts_zoom(max = 3),
              ggiraph::opts_selection(type = "none")
            )
          )
        },
        error = function(e) {
          notify_reselect("heatmap")
          ggiraph::girafe(ggobj = ggplot() +
            theme_void())
        }
      )
    })

    ### ---- Container ----
    output$heatmap_container <- renderUI({
      df <- hm_plotting_df()
      if (is.null(df) || nrow(df) == 0) {
        return(ui_message("No data to show."))
      }

      tags$div(
        class = "hm-layout",
        tags$div(
          class = "hm-plot",
          girafeOutput(
            outputId = "mutationHeatmap",
            width = "60vw",
            height = paste0(heatmap_height_px(), "px")
          )
        ),
        tags$div(
          class = "hm-legend",
          uiOutput("heatmap_legend")
        )
      )
    })
  }


  # ---- Data ----

  mutation_data <- reactive({
    df <- combined_rows()
    req(!is.null(df), nrow(df) > 0)

    if (!"Ref_codes" %in% names(df)) df$Ref_codes <- NA_character_
    if (!"PMID" %in% names(df)) df$PMID <- NA_character_

    pmid_lookup <- refs_paper |>
      transmute(
        PMID = as.character(PMID),
        Ref_code_from_pmid = as.character(Ref_code)
      ) |>
      filter(!is.na(PMID), PMID != "") |>
      distinct(PMID, .keep_all = TRUE)

    df_clean <- df |>
      mutate(
        Ref_codes = na_if(trimws(as.character(Ref_codes)), ""),
        PMID      = na_if(trimws(as.character(PMID)), "")
      )

    df_joined <- df_clean |>
      mutate(.row_id = dplyr::row_number()) |>
      tidyr::separate_rows(PMID, sep = "\\s*;\\s*") |>
      mutate(PMID = dplyr::na_if(trimws(PMID), "")) |>
      left_join(pmid_lookup, by = "PMID")

    pmid_refcodes <- df_joined |>
      dplyr::group_by(.row_id) |>
      dplyr::summarise(
        Ref_code_from_pmid = paste(unique(stats::na.omit(Ref_code_from_pmid)), collapse = "; "),
        PMID = paste(unique(stats::na.omit(PMID)), collapse = "; "),
        .groups = "drop"
      ) |>
      mutate(
        Ref_code_from_pmid = na_if(Ref_code_from_pmid, ""),
        PMID = na_if(PMID, "")
      )

    df_clean |>
      mutate(.row_id = dplyr::row_number()) |>
      left_join(pmid_refcodes, by = ".row_id", suffix = c("", "_collapsed")) |>
      mutate(
        Ref_codes = dplyr::coalesce(Ref_codes, Ref_code_from_pmid),
        PMID = dplyr::coalesce(PMID_collapsed, PMID),
        Gene = "rpoB",
        AA_mut_name = paste0(AA_original, AA_pos, AA_mutation),
        AA_mut_name_Ecoli = paste0(AA_original, AA_pos_Ecoli, AA_mutation)
      ) |>
      dplyr::rename(
        Species      = species,
        AA_Pos       = AA_pos,
        AA_Pos_Ecoli = AA_pos_Ecoli,
        Original     = AA_original,
        Mutation     = AA_mutation,
        Source       = source
      ) |>
      select(
        Species,
        Gene,
        AA_Pos_Ecoli,
        AA_Pos,
        Original,
        Mutation,
        AA_mut_name_Ecoli,
        AA_mut_name,
        Source,
        Ref_codes,
        PMID
      )
  })

  refcode_links_tooltip <- function(ref_codes, ref_table,
                                    code_col = "Ref_code",
                                    url_col = "DOI",
                                    title_col = "Full_ref") {
    # map: code -> url / tooltip text
    url_map <- setNames(
      as.character(ref_table[[url_col]]),
      as.character(ref_table[[code_col]])
    )

    tip_map <- setNames(
      as.character(ref_table[[title_col]]),
      as.character(ref_table[[code_col]])
    )

    vapply(as.character(ref_codes), function(s) {
      s <- trimws(s)
      if (!nzchar(s) || is.na(s)) {
        return("")
      }

      # split codes (comma/semicolon separators only)
      codes <- unlist(strsplit(s, "\\s*[,;]\\s*"))
      codes <- trimws(codes)
      codes <- codes[nzchar(codes)]

      out <- vapply(codes, function(code) {
        url <- url_map[code]
        tip <- tip_map[code]

        code_esc <- htmltools::htmlEscape(code)
        tip_esc <- htmltools::htmlEscape(ifelse(is.null(tip) || is.na(tip) || !nzchar(tip),
          code, tip
        ))

        if (!is.null(url) && !is.na(url) && nzchar(url)) {
          sprintf(
            '<a href="%s" target="_blank" rel="noopener" title="%s">%s</a>',
            htmltools::htmlEscape(url),
            tip_esc,
            code_esc
          )
        } else {
          sprintf('<span title="%s">%s</span>', tip_esc, code_esc)
        }
      }, character(1))

      paste(out, collapse = "; ")
    }, character(1))
  }

  output$mutation_table <- DT::renderDataTable({
    data <- mutation_data()

    data <- data |>
      mutate(Ref_codes = refcode_links_tooltip(
        Ref_codes,
        refs_paper,
        code_col  = "Ref_code",
        url_col   = "DOI",
        title_col = "Full_ref"
      )) |>
      select(-PMID)

    DT::datatable(
      data,
      escape = FALSE,
      rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 25)
    )
  })

  output$download_mutations <- downloadHandler(
    filename = function() paste0("RIFRMap_mutation_table_", Sys.Date(), ".csv"),
    content = function(file) {
      readr::write_csv(mutation_data(), file)
    }
  )

  filtered_references <- reactive({
    df <- mutation_data()
    req(!is.null(df), nrow(df) > 0)

    # extract all codes
    ref_codes <- df$Ref_codes |>
      as.character() |>
      trimws()

    # split multiple codes per cell ("," or ";")
    ref_codes <- unlist(strsplit(ref_codes, "\\s*[,;]\\s*"))

    ref_codes <- ref_codes[nzchar(ref_codes)]
    ref_codes <- unique(ref_codes)

    refs_paper |>
      filter(Ref_code %in% ref_codes)
  })

  output$ref_table <- DT::renderDataTable({
    data <- filtered_references()

    DT::datatable(
      data,
      escape = FALSE,
      rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 15)
    )
  })

  output$download_refs <- downloadHandler(
    filename = function() paste0("RIFRMap_reference_table_", Sys.Date(), ".csv"),
    content = function(file) {
      readr::write_csv(filtered_references(), file)
    }
  )

  # ---- Summary ----
  output$summary_table <- render_gt({
    df <- mutation_data()

    df |>
      select(c(Species, Gene, Source)) |>
      tbl_summary() |>
      as_gt()
  })


  output$bar_aa_mut_name <- renderPlot({
    df <- mutation_data()
    req(df, nrow(df) > 0, all(c("AA_mut_name_Ecoli") %in% names(df)))

    top_muts <- df |>
      mutate(AA_mut_name_Ecoli = as.character(AA_mut_name_Ecoli)) |>
      count(AA_mut_name_Ecoli, sort = TRUE) |>
      slice_head(n = 15) |>
      pull(AA_mut_name_Ecoli)

    plot_df <- df |>
      mutate(AA_mut_name_Ecoli = as.character(AA_mut_name_Ecoli)) |>
      filter(AA_mut_name_Ecoli %in% top_muts) |>
      count(AA_mut_name_Ecoli)

    ggplot(plot_df, aes(x = reorder(AA_mut_name_Ecoli, n, FUN = sum), y = n)) +
      geom_col() +
      coord_flip() +
      scale_y_continuous(
        limits = c(0, NA),
        breaks = scales::breaks_width(1),
        labels = scales::label_number(accuracy = 1),
        expand = c(0, 0)
      ) +
      plotting_theme(legend = "top") +
      theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text  = element_text(size = 12)
      ) +
      labs(
        x = "AA mutation (E. coli coordinate)",
        y = "Number of observations"
      )
  })

  output$hist_aa_pos_ecoli <- renderPlot({
    df <- mutation_data()
    req(df, nrow(df) > 0, "AA_Pos_Ecoli" %in% names(df))

    ggplot(df, aes(x = AA_Pos_Ecoli)) +
      geom_histogram(bins = 30) +
      plotting_theme() +
      theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)
      ) +
      labs(
        x = "Amino acid position (E. coli coordinate)",
        y = "Number of observations"
      )
  })
}
