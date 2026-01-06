server <- function(input, output, session) {
  
  # ---- Load all datasets (once) ----
  card    <- readRDS("data/card.rds")
  cryptic <- readRDS("data/cryptic.rds")
  dragdb  <- readRDS("data/dragdb.rds")
  yang    <- readRDS("data/yang.rds")
  curated   <- readRDS("data/curated.rds")
  
  
  # ----  Load sequencing data ----
  
  refs <- read_csv("data/rpoB_references.csv", show_col_types = FALSE)
  clusters <- read_csv("data/clusters.csv", show_col_types = FALSE)
  seqs <- readDNAStringSet("data/rpoB_references.fasta")
  load("data/coordinates.RData")
  
  # ---- Filter data -----
  ## ---- Helper: species filter ("All species" = no filter) ----
  filter_by_species <- function(df, species_selected, species_col = "species") {
    if (is.null(species_selected) || "All species" %in% species_selected) {
      return(df)
    }
    df |> dplyr::filter(.data[[species_col]] %in% species_selected)
  }
  
  ## ---- Populate species dropdowns ----
  observeEvent(TRUE, {
    
    updateSelectInput(
      session, "card_species",
      choices  = c("All species", sort(unique(card$species))),
      selected = "All species"
    )
    
    updateSelectInput(
      session, "cryptic_species",
      choices  = c("All species", sort(unique(cryptic$species))),
      selected = "All species"
    )
    
    updateSelectInput(
      session, "dragdb_species",
      choices  = c("All species", sort(unique(dragdb$species))),
      selected = "All species"
    )
    
    updateSelectInput(
      session, "yang_species",
      choices  = c("All species", sort(unique(yang$species))),
      selected = "All species"
    )
    
    updateSelectInput(
      session, "curated_species",
      choices  = c("All species", sort(unique(curated$species))),
      selected = "All species"
    )
    
  }, once = TRUE)
  
  # AA position order from ALL data 
  all_pos <- c(
    card$AA_pos_Ecoli,
    cryptic$AA_pos_Ecoli,
    dragdb$AA_pos_Ecoli,
    yang$AA_pos_Ecoli,
    curated$AA_pos_Ecoli
  ) |>
    as.numeric() |>          
    na.omit() |>
    unique() |>
    sort(decreasing = TRUE) |>
    as.character() 
  
  ## ---- Combine selected datasets ----
  combined_rows <- reactive({
    dfs <- list()
    
    # helper to apply numeric range only if column exists
    filter_range_if <- function(df, col, rng) {
      if (!col %in% names(df) || is.null(rng) || length(rng) != 2) return(df)
      df |> dplyr::filter(dplyr::between(.data[[col]], rng[1], rng[2]))
    }
    
    # helper to apply multi-select categorical only if column exists
    filter_in_if <- function(df, col, selected) {
      if (!col %in% names(df) || is.null(selected) || length(selected) == 0 || selected == "") return(df)
      df |> dplyr::filter(.data[[col]] %in% selected)
    }
    
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
      
      # numeric filters (apply only to CRyPTIC rows)
      d <- filter_range_if(d, "coef", input$cryptic_coef)
      d <- filter_range_if(d, "counts", input$cryptic_counts)
      d <- filter_range_if(d, "bh_corrected_p", input$cryptic_bhp)
      
      # categorical filter (if you use it)
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
    
    if (length(dfs) == 0) return(NULL)
    
    dplyr::bind_rows(dfs, .id = "source")
  })
  
  
  # ---- Redirect to View tab on Apply ----
  observeEvent(input$apply_filters, {
    updateTabItems(session, "sidebar", "view")
  })
  
  # ---- No data ui ----
  no_data_ui <- function(msg = "No data to show") {
    div(
      style = "padding: 40px; text-align: center; color: #777;",
      h6(msg)
    )
  }
  
  # ---- Reset button: clear inputs + go back Home
  observeEvent(input$reset_filters, {
    updateCheckboxInput(session, "card_enabled", value = FALSE)
    updateCheckboxInput(session, "cryptic_enabled", value = FALSE)
    updateCheckboxInput(session, "dragdb_enabled", value = FALSE)
    updateCheckboxInput(session, "yang_enabled", value = FALSE)
    updateCheckboxInput(session, "curated_enabled", value = FALSE)
    
    updateSelectInput(session, "card_species", selected = "All species")
    updateSelectInput(session, "cryptic_species", selected = "All species")
    updateSelectInput(session, "dragdb_species", selected = "All species")
    updateSelectInput(session, "yang_species", selected = "All species")
    
    updateSliderInput(session, "cryptic_coef", value = c(-4, 9.6))
    updateSliderInput(session, "cryptic_counts", value = c(0, 2326))
    updateSliderInput(session, "cryptic_bhp", value = c(0, 1))
    updateSelectInput(session, "cryptic_homoplasy", selected = "")
    
    updateSliderInput(session, "dragdb_provean", value = c(-10.7, -2.6))
    updateSliderInput(session, "yang_log2fc", value = c(0, 2.7))
    updateSliderInput(session, "curated_nstudies", value = c(1, 59))
    updateSliderInput(session, "curated_nspecies", value = c(1, 29))
    
    updateTabItems(session, "sidebar", "filters")
  })
  
  
  # ----- Cluster plot -----
  ## ----- Helper functions for cluster plot ----
  
  getBlockData <- function(muts, 
                           clusters, 
                           seqs, 
                           refSeq_ID,
                           coordinatesFocal,
                           species, 
                           clusterName, 
                           mergeDuplicates = TRUE) {
    
    clusterStart <- clusters$AAstart_Ecoli[clusters$Cluster == clusterName]
    clusterEnd <- clusters$AAend_Ecoli[clusters$Cluster == clusterName]
    clusterLength <- clusterEnd - clusterStart + 1
    clusterShift <- ALJEbinf::translateCoordinate(pos = clusterStart,
                                                  coordinates = coordinatesFocal,
                                                  direction = "RefToFocal",
                                                  AAinput = TRUE,
                                                  AAoutput = TRUE) - clusterStart
    block <- muts |>
      filter(Species == species) |>
      mutate(AA_mut_name_Ecoli = paste0("rpoB", "_", AA_original, AA_pos_Ecoli, AA_mutation))
    block <- block |>
      filter(!is.na(AA_pos),
             !is.na(AA_pos_Ecoli),
             !is.na(AA_original),
             !is.na(AA_mutation)) |>
      filter(AA_pos_Ecoli >= clusterStart & AA_pos_Ecoli <= clusterEnd) |>
      arrange(AA_pos_Ecoli, AA_mutation)
    
    if(mergeDuplicates) {
      block <- distinct(block, AA_pos_Ecoli, AA_mutation, .keep_all = TRUE)
    }
    
    block$RelX <- rep(NA_integer_, nrow(block))
    block$RelY <- rep(NA_integer_, nrow(block))
    y <- 1L
    for(i in seq_len(nrow(block))) {
      # the following line might need adjustment if there are indels within a cluster:
      block$RelX[i] <- block$AA_pos_Ecoli[i] - clusterStart + 1L
      if(i > 1) {
        if(block$AA_pos_Ecoli[i] == block$AA_pos_Ecoli[i - 1]) {
          y <- y + 1L
        } else {
          y <- 1L
        }
      }
      block$RelY[i] <- y
    }
    
    # add positions and original bases:
    addInfo <- data.frame(matrix(NA_integer_,
                                 ncol=ncol(block), 
                                 nrow=2 * clusterLength, 
                                 dimnames=list(NULL, colnames(block))))
    # again, the following line might need adjustment if there are indels within a cluster:
    addInfo$AA_pos <- (clusterStart:clusterEnd) + clusterShift
    addInfo$AA_pos_Ecoli <- clusterStart:clusterEnd
    addInfo$RelX <- 1:clusterLength
    addInfo$RelY <- c(rep(-1L, clusterLength), rep(0L, clusterLength))
    addInfo$AA_original <- seqs[[refSeq_ID]] |>
      Biostrings::translate(genetic.code = getGeneticCode("Bacterial", full.search = TRUE)) |>
      as.character() |>
      substr(clusterStart + clusterShift, 
             clusterEnd + clusterShift) |>
      str_split("", simplify = TRUE) |>
      as.vector()
    addInfo <- as_tibble(addInfo)
    
    block <- rbind(addInfo, block)
    block <- select(block, AA_pos, AA_pos_Ecoli, AA_original, AA_mutation, AA_mut_name_Ecoli, RelX, RelY)
    return(block)
  }
  
  renderBlock <- function(block, x0, y0, dx, dy, colcat = NULL, pars) {
    # positions:
    text(x0 + (block$RelX[block$RelY == -1L] - 1L) * dx,
         1 - y0,
         block$AA_pos[block$RelY == -1L],
         cex = , adj = c(0.5, 0.5), srt = 90, family = pars$familyPos)
    # original sequence:
    text(x0 + (block$RelX[block$RelY == 0L] - 1L) * dx,
         1 - y0 - dy * pars$dyPos,
         block$AA_original[block$RelY == 0L],
         cex = pars$cexRef, adj = c(0.5, 0.5), family = pars$familyRef)
    # mutations:
    sel <- block$RelY > 0L
    if (any(sel)) {
      text(x0 + (block$RelX[block$RelY > 0L] - 1L) * dx,
           #scaled down (×0.5) to reduce excessive vertical spacing mutation rows
           1 - y0 - dy * (pars$dyPos + pars$dyRef) - (block$RelY[block$RelY > 0L] - 1L) * dy * 0.5,
           block$AA_mutation[block$RelY > 0L],
           col = block$Colour[block$RelY > 0L],
           cex = pars$cexPos, adj = c(0.5, 0.5), family = pars$familyMuts)
    }
  }
  
  getPars <- function(plot_width_px = NULL) {
    
    # Responsive scaling
    if (is.null(plot_width_px) || is.na(plot_width_px)) {
      scale <- 1
    } else {
      scale <- plot_width_px / 1200
      scale <- max(0.65, min(1.0, scale))
    }
    
    list(
      # layout
      dxBlocks = 0.40 * scale, # relative space between blockss
      dyPos    = 0.5 * scale, # relative space for position numbers
      dyRef    = 0.50 * scale, # relative space for reference sequence
      dyBlocks = 0.5 * scale,

      # text
      cexPos     = 0.95 * scale,
      cexRef     = 1.00 * scale,
      cexSpecies = 1.00 * scale,
      cexMuts    = 0.95 * scale,
      
      # colors
      col = c(
        default     = "#000000",
        ubiquitous  = "#3e3bff",
        unique      = "#ff3e3b",
        polymorphic = "#12947F"
      ),
      
      # base graphics
      mar = pmax(0.2, c(1, 0.5, 0.5, 0.5) * scale),
      familyPos  = "sans",
      familyMuts = "Courier",
      familyRef  = "Courier",
      
      # export
      xInchPerUnit = 0.20,
      yInchPerUnit = 0.20,
      imageRes     = 300,
      
      scale = scale
    )
  }
  
  plotMutations <- function(muts, 
                            clusters, 
                            seqs,
                            coordinates,
                            refs,
                            gene = "rpoB", 
                            clustersToPlot = NULL,
                            speciesToPlot = NULL,
                            file = NULL,
                            colouring = "none",
                            mergeDuplicates = TRUE,
                            plot_width_px = getOption("cluster.plot.width", NA_real_)) {
    
    # initial checks:
    if (is.null(speciesToPlot)) {
      stop("No species to be plotted need to be provided.")
    }
    
    if (!(colouring %in% c("none", "ubi", "uni", "ubiuni"))) {
      stop("Selected option for colouring argument not supported.") 
    }
    
    
    # determine which clusters to plot:
    if (is.null(clustersToPlot) || length(clustersToPlot) == 0) {
      clustersToPlot <- clusters |>
        dplyr::filter(.data$Gene == gene) |>
        dplyr::arrange(.data$AAstart_Ecoli) |>
        dplyr::pull(.data$Cluster) |>
        unique()
    } else {
      clustersToPlot <- clustersToPlot[clustersToPlot %in% clusters$Cluster]
      if (length(clustersToPlot) == 0) stop("clustersToPlot has no valid cluster names.")
      clustersToPlot <- clusters |>
        dplyr::filter(.data$Gene == gene, .data$Cluster %in% clustersToPlot) |>
        dplyr::arrange(.data$AAstart_Ecoli) |>
        dplyr::pull(.data$Cluster) |>
        unique()
    }
    
    
    # delete species that don't have a reference
    dropSpecies <- rep(TRUE, length(speciesToPlot))

    for (i in seq_along(speciesToPlot)) {
      refstrain <- refs$RefStrain[refs$Species == speciesToPlot[[i]]]
      
      if (length(refstrain) > 0 && !is.na(refstrain[1]) && nzchar(refstrain[1])) {
        dropSpecies[i] <- FALSE
      }
    }
    
    if (any(dropSpecies)) {
      warning(paste0("No reference found for the following strains: \n",
                     paste0("  ", unique(unlist(speciesToPlot[dropSpecies])), collapse="\n"), 
                     "\nThese strains will not be included in the plot."))
      speciesToPlot <- speciesToPlot[!dropSpecies]
    }
    
    # obtain all block data:
    allBlocks <- vector(mode = "list", length = length(speciesToPlot))
    names(allBlocks) <- speciesToPlot
    
    blockDimsX <- matrix(NA, nrow = length(speciesToPlot),
                         ncol = length(clustersToPlot))
    blockDimsY <- blockDimsX
    
    for(i in 1:length(speciesToPlot)) {
      refSeq_ID <- refs$FASTA_name[refs$Species == speciesToPlot[[i]]]
      
      allBlocks[[i]] <- vector(mode = "list", length = length(clustersToPlot))
      names(allBlocks[[i]]) <- clustersToPlot
      for(j in 1:length(clustersToPlot)) {
        allBlocks[[i]][[j]] <- getBlockData(muts = muts, 
                                            clusters = clusters, 
                                            seqs = seqs,
                                            refSeq_ID = refSeq_ID,
                                            coordinatesFocal = coordinates[[refSeq_ID]],
                                            species = speciesToPlot[[i]],
                                            clusterName = clustersToPlot[j], 
                                            mergeDuplicates = mergeDuplicates)
        blockDimsX[i, j] <- max(allBlocks[[i]][[j]]$RelX)
        blockDimsY[i, j] <- max(allBlocks[[i]][[j]]$RelY)
      }
    }
    
    
    # determine whether mutations are unique or ubiquitous (for colouring them):
    if (colouring %in% c("ubi", "uni", "ubiuni")) {
      filteredMuts <- muts |>
        filter(Species %in% speciesToPlot) |>
        mutate(AA_mut_name_Ecoli = paste0(gene, "_", AA_original, AA_pos_Ecoli, AA_mutation)) |>
        distinct(Species, AA_mut_name_Ecoli, keep_all = TRUE) 
      
      nSpecies <- filteredMuts |>
        pull(Species) |>
        unique() |>
        length()
      
      mutCounts <- filteredMuts |>
        group_by(AA_mut_name_Ecoli) |>
        summarise(n = n(), .groups = "drop") |>
        select(AA_mut_name_Ecoli, n) |>
        arrange(AA_mut_name_Ecoli) |>
        mutate(category = case_when(n == 1 ~ "unique",
                                    n == nSpecies ~ "ubiquitous",
                                    TRUE ~ "polymorphic"))
      
    }
    
    # determine all plotting parameters:
    pars <- getPars(plot_width_px = plot_width_px)
    # calculate dx (dx = 1 / (nAA + (nCols - 1) * pars$dxBlocks))
    xTotal <- (sum(apply(blockDimsX, 2, max)) + ncol(blockDimsX) * pars$dxBlocks)
    dx <- 1 / xTotal
    # calculate dy (dy = 1 / (nAA + (nRows - 1) * pars$dxBlocks))
  
    yTotal <- (sum(apply(blockDimsY, 1, max)) + 
                 nrow(blockDimsY) * (pars$dyPos + pars$dyRef) +
                 nrow(blockDimsY) * pars$dyBlocks)
    dy <- 1 / yTotal
    
    # determine colours for each amino acid:
    for(i in 1:length(speciesToPlot)) {
      for(j in 1:length(clustersToPlot)) {
        allBlocks[[i]][[j]]$Colour <- pars$col["default"]
        if (colouring == "ubiuni") {
          for(k in 1:nrow(allBlocks[[i]][[j]])) {
            if (!is.na(allBlocks[[i]][[j]]$AA_mut_name_Ecoli[k])) {
              category <- mutCounts[mutCounts$AA_mut_name_Ecoli == allBlocks[[i]][[j]]$AA_mut_name_Ecoli[k], 
                                    "category"] |> pull(category)
              allBlocks[[i]][[j]]$Colour[k] <- pars$col[category]
            }
          }
        } else if (colouring == "ubi") {
          for(k in 1:nrow(allBlocks[[i]][[j]])) {
            if (!is.na(allBlocks[[i]][[j]]$AA_mut_name_Ecoli[k])) {
              category <- mutCounts[mutCounts$AA_mut_name_Ecoli == allBlocks[[i]][[j]]$AA_mut_name_Ecoli[k], 
                                    "category"] |> pull(category)
              category[category == "unique"] <- "default"
              allBlocks[[i]][[j]]$Colour[k] <- pars$col[category]
            }
          }
        } else if (colouring == "uni") {
          for(k in 1:nrow(allBlocks[[i]][[j]])) {
            if (!is.na(allBlocks[[i]][[j]]$AA_mut_name_Ecoli[k])) {
              category <- mutCounts[mutCounts$AA_mut_name_Ecoli == allBlocks[[i]][[j]]$AA_mut_name_Ecoli[k], 
                                    "category"] |> pull(category)
              category[category == "ubiquitous"] <- "default"
              allBlocks[[i]][[j]]$Colour[k] <- pars$col[category]
            }
          }
        }
      }
    }
    
    # actual plotting:
    if (!is.null(file)) {
      fileExtension <- tools::file_ext(file)
      if (fileExtension == "pdf") {
        pdf(file = file,
            width = xTotal * pars$xInchPerUnit,
            height = yTotal * pars$yInchPerUnit)
      } else if (fileExtension == "png") {
        png(file = file,
            width = xTotal * pars$xInchPerUnit,
            height = yTotal * pars$yInchPerUnit,
            units = "in",
            res = pars$imageRes)
      } else if (fileExtension %in% c("jpg", "jpeg")) {
        jpeg(file = file,
             units = "in",
             width = xTotal * pars$xInchPerUnit,
             height = yTotal * pars$yInchPerUnit,
             res = pars$imageRes)
      } else if (fileExtension == "tiff") {
        tiff(file = file,
             units = "in",
             width = xTotal * pars$xInchPerUnit,
             height = yTotal * pars$yInchPerUnit,
             res = pars$imageRes)
      } else {
        stop(paste0("Graphics file type (",fileExtension, ") not supported."))
      }
    }
    
    par(mar = pars$mar)
    plot.new()
    plot.window(
      xlim = c(0, 1),
      ylim = c(0, 1),
      xaxs = "i",
      yaxs = "i"
    )
    
    y0 <- 0
    
    n_species  <- length(speciesToPlot)
    n_clusters <- length(clustersToPlot)
 
    for (i in seq_len(n_species)) {
      # scaled down (×0.75) to reduce excessive vertical spacing between species rows
      if (i > 1) {
        y0 <- y0 + dy * (max(blockDimsY[i - 1, ]) + pars$dyPos + pars$dyRef) * 0.75
      }
      rowLabel <- speciesToPlot[i]
      par(xpd = NA)
      text(0.01, 1 - y0, rowLabel, cex = pars$cexSpecies, adj = c(0, 1), font = 3)
      par(xpd = FALSE)
      
      y0 <- y0 + pars$dyBlocks * dy
      x0 <- pars$dxBlocks * dx
      
      for (j in seq_len(n_clusters)) {
        if (j > 1) {
          x0 <- x0 + dx * (max(blockDimsX[, j - 1]) + pars$dxBlocks)
        }
        
        renderBlock(allBlocks[[i]][[j]], x0, y0, dx, dy, pars = pars)
      }
    }
    
    
    if (!is.null(file)) {
      dev.off() 
    }
  }
  
  
  ## ---- Process data ----
  
  ### ---- Cluster plot data ----
  
  legend_colour <- c(
    Reference   = "#000000",
    Ubiquitous  = "#3e3bff", 
    Unique      = "#FF0000",
    Polymorphic = "#12947F"
  )
  
  output$legend_panel <- renderUI({
    tags$div(
      tags$div(
        "Mutation type",
        style = "font-weight: 600; margin-bottom: 10px;"
      ),
      
      lapply(names(legend_colour), function(name) {
        tags$div(
          style = "
          display: flex;
          align-items: center;
          gap: 10px;
          margin: 6px 0;
        ",
          tags$span(
            style = sprintf(
              "
            width: 14px;
            height: 14px;
            border-radius: 3px;
            display: inline-block;
            background-color: %s;
            border: 1px solid rgba(0,0,0,0.3);
            ",
              legend_colour[[name]]
            )
          ),
          tags$span(
            style = "font-size: 13px;",
            tools::toTitleCase(name)
          )
        )
      })
    )
  })
  
  
  cluster_plotting_df <- reactive({
    df0 <- combined_rows()
    if (is.null(df0) || nrow(df0) == 0) return(NULL)
    
    df <- df0 |>
      dplyr::transmute(
        Species      = as.character(species),
        AA_original  = as.character(AA_original),
        AA_pos       = suppressWarnings(as.integer(AA_pos)),
        AA_pos_Ecoli = suppressWarnings(as.integer(AA_pos_Ecoli)),
        AA_mutation  = as.character(AA_mutation),
        source       = as.character(source)
      ) |>
      dplyr::filter(
        !is.na(Species),
        !is.na(AA_pos_Ecoli),
        !is.na(AA_original),
        !is.na(AA_mutation),
        AA_original != "",
        AA_mutation != ""
      ) |>
      dplyr::distinct()
    
    if (nrow(df) == 0) return(NULL)
    df
  })
  
  ### ---- helper: species has refs + seq + coords ----
  has_ref_assets <- function(species, refs, seqs, coordinates) {
    ref_row <- refs[refs$Species == species, , drop = FALSE]
    if (nrow(ref_row) == 0) return(FALSE)
    
    ref_id <- ref_row$FASTA_name[1]
    if (is.na(ref_id) || ref_id == "") return(FALSE)
    
    if (!(ref_id %in% names(seqs))) return(FALSE)
    if (is.null(coordinates[[ref_id]])) return(FALSE)
    
    TRUE
  }
  
  ### ---- Species to plot (drop those without assets) ----
  cluster_species_to_plot <- reactive({
    df <- cluster_plotting_df()
    if (is.null(df) || nrow(df) == 0) return(character(0))
    
    sp <- sort(unique(df$Species))
    
    keep <- vapply(
      sp,
      has_ref_assets,
      logical(1),
      refs = refs,
      seqs = seqs,
      coordinates = coordinates
    )
    
    sp2 <- sp[keep]
    
    # warn once per change in species selection
    dropped <- setdiff(sp, sp2)
    if (length(dropped) > 0) {
      warning(
        "Dropping species missing reference assets (refs/seqs/coordinates):\n",
        paste0("  - ", dropped, collapse = "\n")
      )
    }
    sp2
  })
  
  ### ---- Dynamic height ----
  cluster_height_px <- reactive({
    n <- length(cluster_species_to_plot())
    if (n == 0) return(250)
    
    # container width in px (falls back if not available yet)
    w <- session$clientData$output_clusterPlot_width %||% 900
    
    # smaller screens -> less height per species
    px_per_species <- dplyr::case_when(
      w < 600  ~ 220,   # phone
      w < 1000 ~ 300,   # tablet / small laptop
      TRUE     ~ 400    # desktop
    )
    
    base_px <- dplyr::case_when(
      w < 600  ~ 80,
      w < 1000 ~ 100,
      TRUE     ~ 120
    )
    
    h <- base_px + n * px_per_species
    max(250, min(h, 20000))
  })
  
  ### ---- Plot ----
  
  output$clusterPlot <- renderPlot({
    df <- cluster_plotting_df()
    validate(need(!is.null(df) && nrow(df) > 0, "No data to show"))
    
    clustersToPlot <- input$clustersToPlot 
    if (is.null(clustersToPlot) || length(clustersToPlot) == 0 || all(clustersToPlot == "")) { 
      clustersToPlot <- NULL } 
    else { 
      clustersToPlot <- clusters |> 
        dplyr::filter(Gene == "rpoB", Cluster %in% clustersToPlot) |> 
        dplyr::arrange(AAstart_Ecoli) |> dplyr::pull(Cluster) 
      }
    
    w <- session$clientData$output_clusterPlot_width %||% 900
    options(cluster.plot.width = w)
    
    plotMutations(
      muts           = df,
      clusters       = clusters,
      seqs           = seqs,
      coordinates    = coordinates,
      refs           = refs,
      gene           = "rpoB",
      clustersToPlot = clustersToPlot,
      speciesToPlot  = cluster_species_to_plot(),
      colouring      = "ubiuni",
      mergeDuplicates= TRUE
    )
  },
  res = 96,
  width  = function() session$clientData$output_clusterPlot_width %||% 900,
  height = function() cluster_height_px()
  )
  
  output$clusterPlotUI <- renderUI({
    plotOutput("clusterPlot", width = "100%")
  })
  
  
  # ---- Upset PLot ----
  n_sources <- reactive({
    df <- combined_rows()
    if (is.null(df) || nrow(df) == 0) return(0)
    length(unique(na.omit(df$source)))
  })
  
  output$mutationUpset <- renderPlot({
    req(input$apply_filters > 0)
    
    df0 <- combined_rows()
    validate(need(!is.null(df0) && nrow(df0) > 0, "No data to show"))
    
    df <- df0 |>
      dplyr::transmute(
        AA_pos_Ecoli = as.character(AA_pos_Ecoli),
        AA_mutation  = as.character(AA_mutation),
        source       = as.character(source)
      ) |>
      dplyr::filter(
        !is.na(AA_pos_Ecoli),
        !is.na(AA_mutation),
        !is.na(source),
        AA_mutation != "",
        source != ""
      ) |>
      dplyr::distinct()
    
    validate(need(nrow(df) > 0, "No data to show"))
    
    upset_df <- df |>
      dplyr::distinct(AA_mutation, AA_pos_Ecoli, source) |>
      dplyr::group_by(AA_mutation, AA_pos_Ecoli) |>
      dplyr::summarise(
        sources = list(sort(unique(source))),
        .groups = "drop"
      )
    
    validate(need(nrow(upset_df) > 0, "No data to show"))
    
    ggplot(upset_df, aes(x = sources)) +
      geom_bar() +
      ggupset::scale_x_upset() +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 8))
      ) +
      labs(y = "Number of mutation sites")
  }, res = 120)
  
  
  
  output$upset_container <- renderUI({
    df <- combined_rows()
    
    if (is.null(df) || nrow(df) == 0) {
      no_data_ui()
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
          plotOutput("mutationUpset", height = "450px", width = "100%")
        )
      )
    }
  })
  
  # ---- Mutation Heatmap ----
  
  ## ---- Process data -----
  
  all_AAs_sorted <- sort(strsplit("HAKRGILPVFWYDESTCMNQ", "")[[1]])
  
  hm_plotting_df <- reactive({
    df0 <- combined_rows()
    
    if (is.null(df0) || nrow(df0) == 0) return(NULL)
    
    df <- df0 |>
      dplyr::transmute(
        AA_pos_Ecoli = as.character(AA_pos_Ecoli),
        AA_mutation  = as.character(AA_mutation),
        source       = as.character(source)
      ) |>
      dplyr::distinct()
    
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
  

  ## ---- colors ----
  source_palette <- c(
    
    
    "CARD"    = "#ead6fd",  
    "CRyPTIC" = "#ffdfba", 
    "DRAGdb"  = "#bae1ff", 
    "Curated"   = "#fff8f0",  
    "Yang"    = "#dff2e1",  
    
   
    "CARD+CRyPTIC"   = "#D1D8BE",  
    "CARD+DRAGdb"    = "#A7C1A8", 
    "CARD+Curated"     = "#819A91",  
    "CARD+Yang"      = "#5FB3C4", 
    
    "CRyPTIC+DRAGdb" = "#B0D19A",  
    "CRyPTIC+Curated"  = "#34955e",  
    "CRyPTIC+Yang"   = "#27b481", 
    
    "Curated+DRAGdb"   = "#3cb8b2",  
    "DRAGdb+Yang"    = "#8aad47",  
    "Curated+Yang"     = "#a2cd5a", 
    
    
    "CARD+CRyPTIC+DRAGdb"  = "#ffde1a", 
    "CARD+CRyPTIC+Curated"   = "#ffce00",  
    "CARD+CRyPTIC+Yang"    = "#ffa700",  
    
    "CARD+Curated+DRAGdb"    = "#eec900",  
    "CARD+DRAGdb+Yang"     = "#e7ca52",  
    "CARD+Curated+Yang"      = "#fecc4b",  
    
    "CRyPTIC+Curated+DRAGdb" = "#fdb414",  
    "CRyPTIC+DRAGdb+Yang"  = "#cba034",  
    "CRyPTIC+Curated+Yang"   = "#F0AE6C",  
    "Curated+DRAGdb+Yang"    = "#E4C07A",  
    
    
    "CARD+CRyPTIC+Curated+DRAGdb" = "#ed453d",
    "CARD+CRyPTIC+Curated+Yang"   = "#ff7f50",  
    "CARD+Curated+DRAGdb+Yang"    = "#f40d30",  
    "CRyPTIC+Curated+DRAGdb+Yang" = "#b31b1b",  
    
   
    "CARD+CRyPTIC+Curated+DRAGdb+Yang" = "#ff0800",
    
    "None" = "#D9D9D9"
  )
  
  legend_order <- c(
    # Singles
    "CARD", "CRyPTIC", "DRAGdb", "Curated", "Yang",
    
    # Pairs
    "CARD+CRyPTIC", "CARD+DRAGdb", "CARD+Curated", "CARD+Yang",
    "CRyPTIC+DRAGdb", "CRyPTIC+Curated", "CRyPTIC+Yang",
    "Curated+DRAGdb", "DRAGdb+Yang",
    "Curated+Yang",
    
    # Triples
    "CARD+CRyPTIC+DRAGdb", "CARD+CRyPTIC+Curated", "CARD+CRyPTIC+Yang",
    "CARD+Curated+DRAGdb", "CARD+DRAGdb+Yang", "CARD+Curated+Yang",
    "CRyPTIC+Curated+DRAGdb", "CRyPTIC+DRAGdb+Yang", "CRyPTIC+Curated+Yang",
    "Curated+DRAGdb+Yang",
    
    # Quadruples
    "CARD+Curated+CRyPTIC+DRAGdbd",
    "CARD+Curated+CRyPTIC+Yang",
    "CARD+Curated+DRAGdb+Yang",
    "CRyPTIC+Curated+DRAGdb+Yang",
    
    # All
    "CARD+CRyPTIC+Curated+DRAGdb+Yang",
    
    # None 
    "None"
  )
  
  ## ---- plot ----
  output$mutationHeatmap <- ggiraph::renderGirafe({
    req(input$apply_filters > 0)
    
    df <- hm_plotting_df()
    
    req(df, nrow(df) > 0)
    
    # normalize combo names so they match palette keys
    df$source_combo <- vapply(
      strsplit(trimws(df$source_combo), "\\+"),
      function(x) paste(sort(x), collapse = "+"),
      character(1)
    )
    
    present_levels <- legend_order[legend_order %in% unique(df$source_combo)]
    df$source_combo <- factor(df$source_combo, levels = present_levels)
    
    n_pos <- length(unique(df$AA_pos_Ecoli))

    p <- ggplot(
      df,
      aes(
        x = AA_mutation,
        y = factor(AA_pos_Ecoli, levels = all_pos),
        fill = source_combo,
        tooltip = paste0(
          "Position: ", AA_pos_Ecoli,
          "<br>AA: ", AA_mutation,
          "<br>Source: ", source_combo
        )
      )
    ) +
      ggiraph::geom_tile_interactive(
        color = "grey80",
        linewidth = 0.2
      ) +
      scale_fill_manual(
        values = source_palette,
        breaks = present_levels,   
        limits = present_levels,  
        drop = TRUE
      ) + 
      scale_x_discrete(position = "top") +
      guides(
        fill = guide_legend(
          ncol = 1,
          byrow = TRUE
        )
      ) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.title.x = element_text(size = 12, margin = margin(b = 30)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text  = element_text(size = 10)
      ) +
      labs(
        x = "Amino acid",
        y = "Position",
        fill = "Source(s)"
      )
    
    ggiraph::girafe(
      ggobj = p,
      width_svg = 12,   
      height_svg = max(400, n_pos * 22)/96,
      options = list(
        ggiraph::opts_hover(css = "fill-opacity:0.85;"),
        ggiraph::opts_tooltip(css = "background-color:white; color:black;"),
        ggiraph::opts_zoom(max = 3),
        ggiraph::opts_selection(type = "none")
      )
    )
  })
  
  ## ---- heatmap container ----
  
  output$heatmap_container <- renderUI({
    df <- hm_plotting_df()
    
    if (is.null(df) || nrow(df) == 0) {
      no_data_ui()
    } else {
      girafeOutput("mutationHeatmap")
    }
  })
  
  
  
  # ---- Summary ----
  
  
  ## ---- helper: ggplot -> inline <img> (base64 png) ----
  plot_to_img <- function(p, width_px = 160, height_px = 50, dpi = 96) {
    tf <- tempfile(fileext = ".png")
    on.exit(unlink(tf), add = TRUE)
    
    ggplot2::ggsave(
      filename = tf,
      plot     = p,
      width    = width_px / dpi,
      height   = height_px / dpi,
      dpi      = dpi,
      units    = "in",
      bg       = "transparent"
    )
    
    b64 <- base64enc::base64encode(tf)
    
    htmltools::tags$img(
      src = paste0("data:image/png;base64,", b64),
      style = paste(
        "display:block;",
        "max-width:100%;",
        "height:auto;",
        "border:none;",
        "background:transparent;",
        sep = " "
      )
    )
  }
  
  ## ---- helper: build one summary row (categorical vs numeric) ----
  summarise_col <- function(df, v, hist_bins = 30) {
    x <- df[[v]]
    type <- class(x)[1]
    
    miss  <- sum(is.na(x) | (is.character(x) & x == ""))
    valid <- length(x) - miss
    
    is_cat <- v %in% c("AA_original", "AA_mutation") || is.character(x) || is.factor(x)
    
    if (is_cat) {
      x0 <- as.character(x)
      x0 <- x0[!is.na(x0) & x0 != ""]
      
      tab <- table(x0)
      tab <- sort(tab, decreasing = TRUE)
      
      freqs_html <- if (length(tab) == 0) {
        ""
      } else {
        paste(
          sprintf(
            "%d. %s&nbsp;&nbsp;%d (%.1f%%)",
            seq_along(tab),
            names(tab),
            as.integer(tab),
            100 * as.integer(tab) / sum(tab)
          ),
          collapse = "<br>"
        )
      }
      
      p <- ggplot2::ggplot(
        data.frame(cat = factor(names(tab), levels = rev(names(tab))),
                   n   = as.integer(tab)),
        ggplot2::aes(x = cat, y = n)
      ) +
        ggplot2::geom_col(fill = "grey70") +
        ggplot2::coord_flip() +
        ggplot2::theme_void() +
        ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
      
      list(
        Variable = v,
        Type     = type,
        Stats    = sprintf("%d distinct values", length(tab)),
        Freqs    = freqs_html,
        Graph    = plot_to_img(p, width_px = 300, height_px = 200),
        Valid    = valid,
        Missing  = miss
      )
      
    } else {
      x0 <- suppressWarnings(as.numeric(x))
      x0 <- x0[!is.na(x0)]
      
      stats_txt <- if (!length(x0)) {
        "No valid numeric data"
      } else {
        htmltools::HTML(sprintf(
          "Mean (sd): %.1f (%.1f)<br>
           min &lt; med &lt; max: %.1f &lt; %.1f &lt; %.1f<br>
           IQR: %.1f",
          mean(x0),
          stats::sd(x0),
          min(x0),
          stats::median(x0),
          max(x0),
          stats::IQR(x0)
        ))
      }
      
      p <- ggplot2::ggplot(data.frame(x = x0), ggplot2::aes(x)) +
        ggplot2::geom_histogram(bins = hist_bins, fill = "grey70", color = "white", linewidth = 0.1) +
        ggplot2::theme_void() +
        ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
      
      list(
        Variable = v,
        Type     = type,
        Stats    = stats_txt,
        Freqs    = sprintf("%d distinct values", length(unique(x0))),
        Graph    = plot_to_img(p, width_px = 300, height_px = 200),
        Valid    = valid,
        Missing  = miss
      )
    }
  }
  
  ## ---- UI output: summary table ----
  output$summary_table_ui <- renderUI({
    req(input$apply_filters > 0)
    
    df0 <- combined_rows() 
    validate(need(!is.null(df0) && nrow(df0) > 0, "No data to show"))
    
    df <- df0 |>
      dplyr::transmute(
        Species      = as.character(species),
        AA_original  = as.character(AA_original),
        AA_pos       = suppressWarnings(as.integer(AA_pos)),
        AA_pos_Ecoli = suppressWarnings(as.integer(AA_pos_Ecoli)),
        AA_mutation  = as.character(AA_mutation),
        Source       = as.character(source))
    
    validate(need(!is.null(df) && nrow(df) > 0, "No data to show"))
    
    cols <- intersect(
      c("Species", "AA_original", "AA_mutation", "AA_pos", "AA_pos_Ecoli", "Source"),
      names(df)
    )
    
    rows <- lapply(cols, function(v) summarise_col(df, v))
    
    htmltools::tags$table(
      class = "table table-sm table-striped",
      htmltools::tags$thead(
        htmltools::tags$tr(
          htmltools::tags$th("Variable"),
          htmltools::tags$th("Type"),
          htmltools::tags$th("Stats / Values"),
          htmltools::tags$th("Freqs (% of Valid)"),
          htmltools::tags$th("Graph"),
          htmltools::tags$th("Valid"),
          htmltools::tags$th("Missing")
        )
      ),
      htmltools::tags$tbody(
        lapply(rows, function(r) {
          htmltools::tags$tr(
            htmltools::tags$td(r$Variable),
            htmltools::tags$td(r$Type),
            htmltools::tags$td(r$Stats),
            htmltools::tags$td(style = "white-space:normal; line-height:1.25;",
                               htmltools::HTML(r$Freqs)),
            htmltools::tags$td(r$Graph),
            htmltools::tags$td(r$Valid),
            htmltools::tags$td(r$Missing)
          )
        })
      )
    )
  })
  
  
}
