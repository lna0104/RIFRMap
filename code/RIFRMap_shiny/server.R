server <- function(input, output, session) {
  
  # ---- Global plot style ----
  PLOT_STYLE <- list(
    # font family: use same name everywhere
    family = "Helvetica",
    
    # base sizes (ggplot)
    base_size   = 9,
    title_size  = 12,
    axis_title  = 10,
    axis_text   = 10,
    caption     = 8,
    
    # colours
    text_col    = "grey20",
    axis_col    = "grey20",
    tick_col    = "grey60",
    panel_bg    = "white",
    
    # margins
    plot_margin = ggplot2::margin(6, 8, 6, 8)
  )
  
  theme_app <- function(
    legend = "none",
    x_top  = FALSE
  ) {
    theme_bw(base_family = PLOT_STYLE$family, base_size = PLOT_STYLE$base_size) +
      theme(
        plot.title      = element_text(
          size = PLOT_STYLE$title_size,
          colour = PLOT_STYLE$text_col
        ),
        plot.subtitle   = element_text(
          size = PLOT_STYLE$base_size,
          colour = PLOT_STYLE$text_col
        ),
        axis.title.x    = element_text(
          size = PLOT_STYLE$axis_title,
          colour = PLOT_STYLE$axis_col,
          margin = ggplot2::margin(t = if (x_top) 8 else 0, b = if (x_top) 0 else 8)
        ),
        axis.title.y    = element_text(
          size = PLOT_STYLE$axis_title,
          colour = PLOT_STYLE$axis_col,
          margin = ggplot2::margin(r = 8)
        ),
        axis.text.x     = element_text(
          size = PLOT_STYLE$axis_text,
          colour = PLOT_STYLE$axis_col
        ),
        axis.text.y     = element_text(
          size = PLOT_STYLE$axis_text,
          colour = PLOT_STYLE$axis_col
        ),
        axis.ticks      = element_line(colour = PLOT_STYLE$tick_col),
        panel.grid      = element_blank(),
        plot.margin     = PLOT_STYLE$plot_margin,
        legend.position = legend
      )
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # ----  Plot message setting ----
  ui_message <- function(msg) {
    div(
      style = sprintf(
        paste(
          "padding:40px;",
          "font-family:%s;",
          "font-size:%dpx;",
          "font-weight:normal;",
          "color:%s;",
          "line-height:1.4;",
          sep = " "
        ),
        PLOT_STYLE$family,
        round(PLOT_STYLE$base_size * 2),
        PLOT_STYLE$text_col
      ),
      msg
    )
  }
  

  
  spin <- function(x) shinycssloaders::withSpinner(x, color = "#0dc5c1")
  
  
  
  
  # ---- Load all datasets ----
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
  # species filter ("All species" = no filter)
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
      choices  = c(sort(unique(cryptic$species))),
      selected = "Mycobacterium tuberculosis"
    )
    
    updateSelectInput(
      session, "dragdb_species",
      choices  = c("All species", sort(unique(dragdb$species))),
      selected = "All species"
    )
    
    updateSelectInput(
      session, "yang_species",
      choices  = c(sort(unique(yang$species))),
      selected = "Escherichia coli"
    )
    
    updateSelectInput(
      session, "curated_species",
      choices  = c("All species", sort(unique(curated$species))),
      selected = "All species"
    )
    
  }, once = TRUE)
  
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
      
      # numeric filters 
      d <- filter_range_if(d, "coef", input$cryptic_coef)
      d <- filter_range_if(d, "counts", input$cryptic_counts)
      d <- d |> 
        dplyr::filter(bh_corrected_p <= as.numeric(input$cryptic_bhp))
      
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
    
    if (length(dfs) == 0) return(NULL)
    
    dplyr::bind_rows(dfs, .id = "source")
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
    updateSliderInput(session, "cryptic_bhp", value = c(0, 1))
    updateSelectInput(session, "cryptic_homoplasy", selected = "")
    
    updateSliderInput(session, "dragdb_provean", value = c(-10.7, -2.6))
    updateSliderInput(session, "yang_log2fc", value = c(0, 2.7))
    updateSliderInput(session, "curated_nstudies", value = c(1, 59))
    updateSliderInput(session, "curated_nspecies", value = c(1, 29))
    
    updateTabItems(session, "sidebar", "filters")
  })
  
  
  # ----- Species-Cluster Tab -----
  {
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
           cex = pars$cexPos, adj = c(0.5, 0.5), srt = 90, family = pars$familyPos)
      # original sequence:
      text(x0 + (block$RelX[block$RelY == 0L] - 1L) * dx,
           1 - y0 - dy * pars$dyPos,
           block$AA_original[block$RelY == 0L],
           cex = pars$cexRef, adj = c(0.5, 0.5), family = pars$familyRef, font = 2)
      # mutations:
      if (length(block$AA_mutation[block$RelY > 0L]) > 0) {
        text(x0 + (block$RelX[block$RelY > 0L] - 1L) * dx,
             1 - y0 - dy * (pars$dyPos + pars$dyRef) - (block$RelY[block$RelY > 0L] - 1L) * dy,
             block$AA_mutation[block$RelY > 0L],
             col = block$Colour[block$RelY > 0L],
             cex = pars$cexMut, adj = c(0.5, 0.5), family = pars$familyMuts)
      }
    }
    
    
    getPars <- function() {
      pars <- vector(mode = "list")
      pars$dxBlocks <- 0.5  # relative space between blocks
      
      pars$dyPos <- 1.2 # relative space for position numbers
      pars$dyRef <- 1.3 # relative space for reference sequence
      pars$dyBlocks <- 2.2
      
      pars$cexPos <- 0.9 # relative font size of position numbers
      pars$cexRef <- 0.9 # relative font size of reference sequence
      pars$cexMut <- 0.9
      pars$cexSpecies <- 1 # relative font size of species labels
      
      pars$col <- c(default     = "#000000",
                    ubiquitous  = "#3e3bff", 
                    unique      = "#ff3e3b",
                    polymorphic = "#12947F")
      
      pars$mar <- c(1, 1, 1, 1) + 0.1
      pars$familyPos  <- PLOT_STYLE$family
      pars$familyMuts <- PLOT_STYLE$family
      pars$familyRef  <- PLOT_STYLE$family
      
      pars$xInchPerUnit <- 0.2
      pars$yInchPerUnit <- 0.2
      return(pars)
    }
    
    
    plotMutations <- function(
      muts,
      clusters,
      seqs,
      coordinates,
      refs,
      gene = "rpoB",
      clustersToPlot = NULL,
      speciesToPlot  = NULL,
      colouring = "none",
      mergeDuplicates = TRUE,
      plot_width_px = getOption("cluster.plot.width", NA_real_),
      return_dims = FALSE,
      return_hover_map = FALSE,
      res = 96,
      mut_lines_cap = 12L 
    ) {
      
      # initial checks 
      if (is.null(speciesToPlot) || length(speciesToPlot) == 0) {
        stop("No species to be plotted need to be provided.")
      }
      
      if (!colouring %in% c("none", "ubi", "uni", "ubiuni")) {
        stop("Selected option for colouring argument not supported.")
      }
      
      # determine clusters to plot 
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
      
      if (length(clustersToPlot) == 0) {
        stop("No clusters available to plot for this gene.")
      }
      
      # drop species without a reference 
      dropSpecies <- rep(TRUE, length(speciesToPlot))
      
      for (i in seq_along(speciesToPlot)) {
        refstrain <- refs$RefStrain[refs$Species == speciesToPlot[[i]]]
        if (length(refstrain) > 0 && !is.na(refstrain[1]) && nzchar(refstrain[1])) {
          dropSpecies[i] <- FALSE
        }
      }
      
      if (any(dropSpecies)) {
        warning(
          paste0(
            "No reference found for the following strains:\n",
            paste0("  ", unique(unlist(speciesToPlot[dropSpecies])), collapse = "\n"),
            "\nThese strains will not be included in the plot."
          )
        )
        speciesToPlot <- speciesToPlot[!dropSpecies]
      }
      
      if (length(speciesToPlot) == 0) {
        stop("All selected species were dropped because reference assets are missing.")
      }
      
      # obtain all block data 
      n_species  <- length(speciesToPlot)
      n_clusters <- length(clustersToPlot)
      
      allBlocks <- vector("list", length = n_species)
      names(allBlocks) <- speciesToPlot
      
      blockDimsX <- matrix(NA_real_, nrow = n_species, ncol = n_clusters)
      blockDimsY <- matrix(NA_real_, nrow = n_species, ncol = n_clusters)
      
      for (i in seq_along(speciesToPlot)) {
        refSeq_ID <- refs$FASTA_name[refs$Species == speciesToPlot[[i]]][1]
        
        allBlocks[[i]] <- vector("list", length = n_clusters)
        names(allBlocks[[i]]) <- clustersToPlot
        
        for (j in seq_along(clustersToPlot)) {
          blk <- getBlockData(
            muts            = muts,
            clusters        = clusters,
            seqs            = seqs,
            refSeq_ID       = refSeq_ID,
            coordinatesFocal= coordinates[[refSeq_ID]],
            species         = speciesToPlot[[i]],
            clusterName     = clustersToPlot[[j]],
            mergeDuplicates = mergeDuplicates
          )
          
          allBlocks[[i]][[j]] <- blk
          blockDimsX[i, j] <- suppressWarnings(max(blk$RelX, na.rm = TRUE))
          blockDimsY[i, j] <- suppressWarnings(max(blk$RelY, na.rm = TRUE))
        }
      }
      
      # determine mutation categories for colouring 
      mutCounts <- NULL
      if (colouring %in% c("ubi", "uni", "ubiuni")) {
        filteredMuts <- muts |>
          dplyr::filter(.data$Species %in% speciesToPlot) |>
          dplyr::mutate(AA_mut_name_Ecoli = paste0(gene, "_", AA_original, AA_pos_Ecoli, AA_mutation)) |>
          dplyr::distinct(Species, AA_mut_name_Ecoli, .keep_all = TRUE)
        
        nSpecies <- filteredMuts |>
          dplyr::pull(Species) |>
          unique() |>
          length()
        
        mutCounts <- filteredMuts |>
          dplyr::count(AA_mut_name_Ecoli, name = "n") |>
          dplyr::mutate(
            category = dplyr::case_when(
              n == 1 ~ "unique",
              n == nSpecies ~ "ubiquitous",
              TRUE ~ "polymorphic"
            )
          )
      }
      
      # plotting parameters 
      pars <- getPars()
      
      # safe max per column/row
      colMax <- apply(blockDimsX, 2, max, na.rm = TRUE)
      rowMax <- apply(blockDimsY, 1, max, na.rm = TRUE)
      
      xTotal <- sum(colMax) + n_clusters * pars$dxBlocks
      dx <- 1 / xTotal
      
      yTotal <- sum(rowMax) +
        n_species * (pars$dyPos + pars$dyRef) +
        n_species * pars$dyBlocks
      dy <- 1 / yTotal
      
      # Calculate width and height for UI sizing
      width_px  <- xTotal * pars$xInchPerUnit * res
      height_px <- yTotal * pars$yInchPerUnit * res
      
      # Build x column map  
      {
        xcols <- vector("list", n_clusters)
        
        x0c <- pars$dxBlocks * dx
        
        for (j in seq_len(n_clusters)) {
          if (j > 1) {
            x0c <- x0c + dx * (colMax[j - 1] + pars$dxBlocks)
          }
          
          n_units <- as.integer(colMax[j])
          xmin <- x0c - dx / 2
          xmax <- x0c + (max(1L, n_units) - 1L) * dx + dx / 2
          
          xcols[[j]] <- data.frame(
            cluster = clustersToPlot[[j]],
            x0      = x0c,     # center of RelX==1
            n       = n_units, # number of AA columns in this cluster
            xmin    = xmin,
            xmax    = xmax,
            stringsAsFactors = FALSE
          )
        }
        xcols <- dplyr::bind_rows(xcols)
      }

      if (isTRUE(return_dims)) {
        return(list(
          width_px  = width_px,
          height_px = height_px,
          xTotal    = xTotal,
          yTotal    = yTotal,
          pars      = pars,
          dx        = dx,
          xcols     = xcols
        ))
      }
      
      
      # assign colours 
      for (i in seq_along(speciesToPlot)) {
        for (j in seq_along(clustersToPlot)) {
          blk <- allBlocks[[i]][[j]]
          blk$Colour <- pars$col[["default"]]
          
          if (!is.null(mutCounts) && "AA_mut_name_Ecoli" %in% names(blk)) {
            # join category onto block by mutation name
            idx <- match(blk$AA_mut_name_Ecoli, mutCounts$AA_mut_name_Ecoli)
            cat <- mutCounts$category[idx]
            
            if (colouring == "ubi")    cat[cat == "unique"]      <- "default"
            if (colouring == "uni")    cat[cat == "ubiquitous"]  <- "default"
            if (colouring == "none")   cat[] <- "default"
            
            # only update where category is not NA
            ok <- !is.na(cat)
            blk$Colour[ok] <- pars$col[cat[ok]]
          }
          
          allBlocks[[i]][[j]] <- blk
        }
      }
      
      species_bands <- vector("list", n_species)
      # draw base plot 
      graphics::par(mar = pars$mar)
      graphics::plot.new()
      graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
      plt <- graphics::par("plt")
      
      y0 <- 0
      for (i in seq_len(n_species)) {
        if (i > 1) {
          y0 <- y0 + dy * (rowMax[i - 1] + pars$dyPos + pars$dyRef)
        }
        
        rowLabel <- speciesToPlot[[i]]
        
        graphics::par(xpd = NA)
        graphics::text(0, 1 - y0, rowLabel, cex = pars$cexSpecies, adj = c(0, 1), font = 2)
        
        y_line <- 1 - y0 - dy * 1.1
        graphics::segments(0, y_line, 1, y_line, col = "grey70", lwd = 1)
        graphics::par(xpd = FALSE)
        
        y0 <- y0 + pars$dyBlocks * dy
        
        # band y positions computed at the SAME y0 used by renderBlock()
        y_pos <- 1 - y0
        y_ref <- 1 - y0 - dy * pars$dyPos
        mut_top <- 1 - y0 - dy * (pars$dyPos + pars$dyRef)
        
        mut_lines <- rowMax[i]
        mut_lines_cap_i <- min(mut_lines, mut_lines_cap)
        mut_bottom <- mut_top - dy * max(0, (mut_lines_cap_i - 1))
        
        # thickness/padding tuned for visual stability across devices
        band_pos_h <- dy * 1.2
        band_ref_h <- dy * 0.90
        band_mut_pad <- dy * 0.55
        
        species_bands[[i]] <- data.frame(
          species = rowLabel,
          pos_ymin = y_pos - band_pos_h / 2, pos_ymax = y_pos + band_pos_h / 2,
          ref_ymin = y_ref - band_ref_h / 2, ref_ymax = y_ref + band_ref_h / 2,
          mut_ymin = mut_bottom - band_mut_pad, mut_ymax = mut_top + band_mut_pad,
          stringsAsFactors = FALSE
        )
        
        x0 <- pars$dxBlocks * dx
        
        for (j in seq_len(n_clusters)) {
          if (j > 1) {
            x0 <- x0 + dx * (colMax[j - 1] + pars$dxBlocks)
          }
          renderBlock(allBlocks[[i]][[j]], x0, y0, dx, dy, pars = pars)
        }
      }
      
      if (isTRUE(return_hover_map)) {
        return(list(
          bands = dplyr::bind_rows(species_bands),
          xcols = xcols,
          dx    = dx,
          plt   = plt
        ))
      }
      
      return(invisible(NULL))
    }
    
    
    # ---- Process data ----
    
    # Helper: cluster selection check
    no_clusters_selected <- function(x) {
      is.null(x) || !length(x) || all(is.na(x)) || all(x == "")
    }
    
    # 1) Clean + normalize plotting df 
    cluster_plotting_df <- reactive({
      df0 <- combined_rows()
      if (is.null(df0) || !nrow(df0)) return(NULL)
      
      df <- df0 |>
        dplyr::transmute(
          Species      = dplyr::na_if(trimws(as.character(species)), ""),
          AA_original  = dplyr::na_if(trimws(as.character(AA_original)), ""),
          AA_mutation  = dplyr::na_if(trimws(as.character(AA_mutation)), ""),
          AA_pos       = suppressWarnings(as.integer(AA_pos)),
          AA_pos_Ecoli = suppressWarnings(as.integer(AA_pos_Ecoli)),
          source       = dplyr::na_if(trimws(as.character(source)), "")
        ) |>
        dplyr::filter(
          !is.na(Species),
          !is.na(AA_pos_Ecoli),
          !is.na(AA_original),
          !is.na(AA_mutation)
        ) |>
        dplyr::distinct(Species, AA_original, AA_mutation, AA_pos, AA_pos_Ecoli, source)
      
      if (!nrow(df)) return(NULL)
      df
    })
    
    # 2) Precompute which species have all reference assets 
    ref_assets_ok <- local({
      # Named logical: species -> TRUE/FALSE
      cache <- NULL
      
      function(refs, seqs, coordinates) {
        if (!is.null(cache)) return(cache)
        
        # only keep rows with a valid FASTA_name
        rr <- refs |>
          dplyr::transmute(
            Species = as.character(.data$Species),
            FASTA   = dplyr::na_if(trimws(as.character(.data$FASTA_name)), "")
          ) |>
          dplyr::filter(!is.na(.data$Species), !is.na(.data$FASTA)) |>
          dplyr::distinct(.data$Species, .keep_all = TRUE)
        
        ok <- setNames(rep(FALSE, nrow(rr)), rr$Species)
        
        if (nrow(rr)) {
          has_seq  <- rr$FASTA %in% names(seqs)
          has_coord <- vapply(rr$FASTA, function(id) !is.null(coordinates[[id]]), logical(1))
          ok[] <- has_seq & has_coord
        }
        
        cache <<- ok
        cache
      }
    })
    
    # 3) Species to plot (drop those without assets) + warn only when list changes
    cluster_species_to_plot <- local({
      last_dropped <- character(0)
      
      reactive({
        df <- cluster_plotting_df()
        if (is.null(df) || !nrow(df)) return(character(0))
        
        sp <- sort(unique(df$Species))
        
        ok_map <- ref_assets_ok(refs = refs, seqs = seqs, coordinates = coordinates)
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
    
    # 4) Canonical ordered clusters (compute once per selection)
    clusters_to_plot_r <- reactive({
      if (no_clusters_selected(input$clustersToPlot)) return(character(0))
      
      clusters |>
        dplyr::filter(.data$Gene == "rpoB", .data$Cluster %in% input$clustersToPlot) |>
        dplyr::arrange(.data$AAstart_Ecoli) |>
        dplyr::pull(.data$Cluster) |>
        unique()
    })
    
    # 5) Cache plot dimensions so we don’t call plotMutations twice
    cluster_dims_r <- reactiveVal(NULL)
    
    observeEvent(
      list(input$apply_filters, input$clustersToPlot),
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
            refs            = refs,
            gene            = "rpoB",
            clustersToPlot  = cl,
            speciesToPlot   = sp,
            colouring       = "ubiuni",
            mergeDuplicates = TRUE,
            plot_width_px   = NA,
            return_dims     = TRUE
          )
        )
      },
      ignoreInit = TRUE
    )
    
    # 6) Trigger for plot redraw 
    cluster_trigger <- reactive({
      input$apply_filters
      input$clustersToPlot
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
      cl <- input$clustersToPlot
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
      if (isTRUE(st$ok)) return("NULL")
      div(class = "plot-empty-msg", st$msg)
    })
    
    output$clusterPlot <- renderPlot(
      {
        st <- cluster_status()
        if (!isTRUE(st$ok)) return(invisible())
        
        df <- cluster_plotting_df()
        cl <- input$clustersToPlot
        sp <- cluster_species_to_plot()
        
        out <- plotMutations(
          muts             = df,
          clusters         = clusters,
          seqs             = seqs,
          coordinates      = coordinates,
          refs             = refs,
          gene             = "rpoB",
          clustersToPlot   = cl,
          speciesToPlot    = sp,
          colouring        = "ubiuni",
          mergeDuplicates  = TRUE,
          return_hover_map = TRUE # now expected to return $bands, not huge per-letter map
        )
        
        # capture geometry from this plot
        payload <- list(
          bands = out$bands,
          plt   = out$plt,
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
    observeEvent(hover_payload_r(), {
      p <- hover_payload_r()
      
      req(p$bands, p$plt, p$xcols)
      
      session$onFlushed(function() {
        session$sendCustomMessage("hoverMap", p)
      }, once = TRUE)
    }, ignoreInit = TRUE)
    
    
    ## ---- Legend panel ----
    
    legend_reference <- c(
      Reference = "#000000"
    )
    
    legend_distribution <- c(
      Ubiquitous  = "#3e3bff",
      Polymorphic = "#12947F",
      Unique      = "#FF0000"
    )
    
    
    output$legend_panel <- renderUI({
      tags$div(
        class = "cluster-legend",
        
        ## ---- Reference section ----
        tags$div(
          class = "cluster-legend-section",
          tags$div(
            class = "cluster-legend-title",
            "Reference sequence"
          ),
          lapply(names(legend_reference), function(name) {
            tags$div(
              class = "cluster-legend-item",
              tags$span(
                class = "cluster-legend-swatch",
                style = sprintf("background-color:%s;", legend_reference[[name]])
              ),
              tags$span(
                class = "cluster-legend-label",
                "Reference"
              )
            )
          })
        ),
        
        ## ---- Distribution section ----
        tags$div(
          class = "cluster-legend-section",
          tags$div(
            class = "cluster-legend-title",
            "Mutation distribution across species"
          ),
          lapply(names(legend_distribution), function(name) {
            tags$div(
              class = "cluster-legend-item",
              tags$span(
                class = "cluster-legend-swatch",
                style = sprintf("background-color:%s;", legend_distribution[[name]])
              ),
              tags$span(
                class = "cluster-legend-label",
                name
              )
            )
          })
        )
      )
    })
    
    
  }
  
  # ---- Mutation Frequency Tab ----

  output$heatmap_pos_mut <- ggiraph::renderGirafe({
    df <- cluster_plotting_df()
    req(df, nrow(df) > 0)
    
    all_AAs_sorted <- sort(strsplit("HAKRGILPVFWYDESTCMNQ", "")[[1]])
    
    heat_df <- df |>
      dplyr::filter(!is.na(AA_pos_Ecoli), !is.na(AA_mutation), !is.na(Species)) |>
      dplyr::group_by(AA_pos_Ecoli, AA_mutation) |>
      dplyr::summarise(
        n_species = dplyr::n_distinct(Species),
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
      sort() |>
      as.character()
    
    n_pos <- length(pos_levels)
    bar_h_cm <- max(4, min(10, 0.15 * n_pos))
    
    heat_df <- heat_df |>
      dplyr::mutate(
        AA_mutation  = factor(AA_mutation, levels = all_AAs_sorted),
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
      heat_df <- heat_df |> dplyr::mutate(fill_val = factor(n_species, levels = levs))
      
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
            keywidth  = grid::unit(0.5, "cm")
          )
        )
      
    } else {
      if (max_n <= 0) max_n <- 1
      heat_df <- heat_df |> dplyr::mutate(fill_val = n_species)
      
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
      labs(x = "Amino acid", y = "Position", fill = "No. species") +
      theme_app(legend = "right", x_top = TRUE) +
      theme(
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text  = element_text(size = 12)
      )
    

    
    # height scaling: prevent insane svg height
    h_in <- max(5, min(35, 0.35 * n_pos))
    
    ggiraph::girafe(
      ggobj = p,
      width_svg  = 12,
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
    aa  <- parts[2]
    req(!is.na(pos), nzchar(aa))
    
    spp <- df |>
      dplyr::filter(AA_pos_Ecoli == pos, AA_mutation == aa) |>
      dplyr::distinct(Species) |>
      dplyr::arrange(Species)
    
    list(pos = pos, aa = aa, species = spp)
  })
  
  
  output$selected_tile_info <- renderUI({
    x <- selected_tile_species()
    tags$div(
      style = "margin-bottom:8px; font-weight:600;",
      sprintf("AA_Pos_Ecoli = %s; Mutation = %s (%d species)",
              x$pos, x$aa, nrow(x$species))
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
      if (is.null(df) || nrow(df) == 0) return(0)
      length(unique(na.omit(df$source)))
    })
    
    output$mutationUpset <- renderPlot({
      req(input$apply_filters > 0)
      
      df0 <- combined_rows()
      
      if (is.null(df0) || nrow(df0) == 0) {
        return(ui_message("No data to show."))
      }
      
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
      
      if (is.null(df) || nrow(df) == 0) {
        return(ui_message("No data to show."))
      }
      
      upset_df <- df |>
        dplyr::distinct(AA_mutation, AA_pos_Ecoli, source) |>
        dplyr::group_by(AA_mutation, AA_pos_Ecoli) |>
        dplyr::summarise(
          sources = list(sort(unique(source))),
          .groups = "drop"
        )
      
      if (is.null(upset_df) || nrow(upset_df) == 0) {
        return(ui_message("No data to show."))
      }
      
      ggplot(upset_df, aes(x = sources)) +
        geom_bar() +
        ggupset::scale_x_upset() +
        labs(y = "Number of mutation sites", x = "Source(s)") +
        theme_app(legend = "none")
    }, res = 120)
    
    
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
    ### ---- Process data -----
  
  all_AAs_sorted <- sort(strsplit("HAKRGILPVFWYDESTCMNQ", "")[[1]])
  
  selected_pos <- reactive({
    df0 <- combined_rows()
    if (is.null(df0) || nrow(df0) == 0) return(character(0))
    
    pos <- df0$AA_pos_Ecoli |>
      suppressWarnings(as.numeric()) |>
      na.omit() |>
      unique() |>
      sort(decreasing = TRUE) |>
      as.character()
    
    pos
  })
  
  
  hm_plotting_df <- reactive({
    df0 <- combined_rows()
    if (is.null(df0) || nrow(df0) == 0) return(NULL)
    
    pos_levels <- selected_pos()
    if (length(pos_levels) == 0) return(NULL)
    
    df <- df0 |>
      dplyr::transmute(
        AA_pos_Ecoli = as.character(AA_pos_Ecoli),
        AA_mutation  = as.character(AA_mutation),
        source       = as.character(source)
      ) |>
      dplyr::distinct()
    
    req(nrow(df) > 0)
    
    base_grid <- tidyr::expand_grid(
      AA_pos_Ecoli = pos_levels,
      AA_mutation = all_AAs_sorted
    )
    
    # obs <- df |>
    #   dplyr::distinct(AA_pos_Ecoli, AA_mutation, source) |>
    #   dplyr::group_by(AA_pos_Ecoli, AA_mutation) |>
    #   dplyr::summarise(
    #     source_set = list(sort(unique(source))),   
    #     .groups = "drop"
    #   ) |>
    #   dplyr::mutate(
    #     source_combo = vapply(source_set, paste, collapse = "+", FUN.VALUE = character(1))
    #   )
    obs <- df |>
      distinct(AA_pos_Ecoli, AA_mutation, source) |>
      group_by(AA_pos_Ecoli, AA_mutation) |>
      summarise(
        source_set = list(sort(unique(source))),
        .groups = "drop"
      ) |>
      mutate(
        source_combo_all = vapply(source_set, \(x) paste(x, collapse = "+"), character(1)),
        source_combo = Map(\(x, all) c(x, all), source_set, source_combo_all)
      ) |>
      select(-source_combo_all) |>
      tidyr::unnest(source_combo)
    
    base_grid |>
      dplyr::left_join(obs, by = c("AA_pos_Ecoli", "AA_mutation")) |>
      dplyr::mutate(
        source_combo = tidyr::replace_na(source_combo, "None"),
        source_set   = ifelse(is.na(source_combo), list(character(0)), source_set),
        AA_mutation  = factor(AA_mutation, levels = all_AAs_sorted),
        AA_pos_Ecoli = factor(AA_pos_Ecoli, levels = pos_levels)
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
    "CARD+Curated+DRAGdb+Yang"    = "#D34E4E",  
    "CRyPTIC+Curated+DRAGdb+Yang" = "#FF5555",  
    
   
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
    "CARD+Curated+CRyPTIC+DRAGdb",
    "CARD+Curated+CRyPTIC+Yang",
    "CARD+Curated+DRAGdb+Yang",
    "CRyPTIC+Curated+DRAGdb+Yang",
    
    # All
    "CARD+CRyPTIC+Curated+DRAGdb+Yang",
    
    # None 
    "None"
  )
  
  # Normalize "A+B" combo strings: trim, split, sort, re-join
  normalize_combo <- function(x) {
    vapply(
      strsplit(trimws(x), "\\+"),
      function(parts) paste(sort(parts), collapse = "+"),
      character(1)
    )
  }
  
  ## ---- Dynamic height -----
  
  heatmap_height_px <- reactive({
    pos_levels <- selected_pos()

    if (is.null(pos_levels) || length(pos_levels) == 0) {
      return(300)
    }
    
    n_pos <- length(pos_levels)
    # row height
    row_px <- 30 
    
    # minimum height safeguard
    max(300, n_pos * row_px)
  })
  
  ## ---- Legend ----
  selected_sources <- reactiveVal(NULL)
  
  output$heatmap_legend <- renderUI({
    req(input$apply_filters > 0)
    
    df <- hm_plotting_df()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    df$source_combo <- normalize_combo(df$source_combo)
    
    present_levels <- legend_order[legend_order %in% unique(df$source_combo)]
    if (length(present_levels) == 0) {
      present_levels <- sort(unique(df$source_combo))
    }
    
    # clickable legend levels (exclude "None")
    clickable_levels <- setdiff(present_levels, "None")
    
    # initialise selection state once
    sel <- selected_sources()
    if (is.null(sel) || !identical(sort(names(sel)), sort(clickable_levels))) {
      sel <- setNames(rep(FALSE, length(clickable_levels)), clickable_levels)
      selected_sources(sel)
    }
    
    cols <- source_palette[clickable_levels]
    
    tags$div(
      class = "hm-legend-col",
      lapply(seq_along(clickable_levels), function(i) {
        lvl <- clickable_levels[i]
        col <- cols[[lvl]]
        is_on <- isTRUE(sel[[lvl]])

        actionButton(
          inputId = paste0("legend_", i),
          label = tags$span(
            class = "hm-legend-btn-label",
            tags$span(
              class = "hm-legend-swatch",
              style = sprintf("--swatch:%s;", col) 
            ),
            lvl
          ),
          class = paste("hm-legend-btn", if (is_on) "is-on" else "is-off"),
          style = sprintf("--btncol:%s;", col)  
        )
      })
    )
  })
  
  # Select legend 
  observeEvent(input$apply_filters, {
    df <- hm_plotting_df()
    if (is.null(df) || nrow(df) == 0) return()
    
    df$source_combo <- normalize_combo(df$source_combo)
    
    present_levels <- legend_order[legend_order %in% unique(df$source_combo)]
    if (length(present_levels) == 0) {
      present_levels <- sort(unique(df$source_combo))
    }
    
    clickable_levels <- setdiff(present_levels, "None")
    
    # reset selection when filters re-apply
    selected_sources(setNames(rep(FALSE, length(clickable_levels)), clickable_levels))
    
    lapply(seq_along(clickable_levels), function(i) {
      observeEvent(input[[paste0("legend_", i)]], {
        lvl <- clickable_levels[i]
        sel <- selected_sources()
        sel[[lvl]] <- !isTRUE(sel[[lvl]])
        selected_sources(sel)
      }, ignoreInit = TRUE)
    })
  }, ignoreInit = TRUE)
  
  observeEvent(input$reset_heatmap_legend, {
    df <- hm_plotting_df()
    req(!is.null(df), nrow(df) > 0)
    
    df$source_combo <- normalize_combo(df$source_combo)
    
    present_levels <- legend_order[legend_order %in% unique(df$source_combo)]
    if (length(present_levels) == 0) present_levels <- sort(unique(df$source_combo))
    
    lvls <- setdiff(present_levels, "None")
    
    selected_sources(setNames(rep(FALSE, length(lvls)), lvls))
  }, ignoreInit = TRUE)
  
  output$mutationHeatmap <- ggiraph::renderGirafe({
    req(input$apply_filters > 0)
    
    df <- hm_plotting_df()
    req(nrow(df) > 0)
    
    df$source_combo <- normalize_combo(df$source_combo)
    
    present_levels <- legend_order[legend_order %in% unique(df$source_combo)]
    if (length(present_levels) == 0) {
      present_levels <- sort(unique(df$source_combo))
    }
    
    # selection state
    sel <- selected_sources()
    clickable_levels <- setdiff(present_levels, "None")
    
    if (is.null(sel) || !identical(sort(names(sel)), sort(clickable_levels))) {
      sel <- setNames(rep(FALSE, length(clickable_levels)), clickable_levels)
      selected_sources(sel)
    }
    
    label_to_set <- function(lbl) sort(strsplit(lbl, "\\+")[[1]])
    
    has_required <- function(source_set, required) {
      vapply(source_set, function(s) all(required %in% s), logical(1))
    }
    
    keep_levels <- names(sel)[sel]
    
    # print(sel)
    # print(keep_levels)
    
    if (length(keep_levels) > 0) {
      # union of selections: keep tiles that satisfy any selected required-set
      keep_mask <- Reduce(`|`, lapply(keep_levels, function(lbl) {
        req_set <- label_to_set(lbl)
        has_required(df$source_set, req_set)
      }))
      df$fill_group <- ifelse(keep_mask, df$source_combo, "None")
    } else {
      df$fill_group <- df$source_combo
    }
    
    # use ALL levels including None
    display_levels <- unique(c(present_levels, "None"))
    df$fill_group <- factor(df$fill_group, levels = display_levels)
    
    pos_levels <- selected_pos()
    n_pos <- length(unique(df$AA_pos_Ecoli))
    
    p <- ggplot(
      df,
      aes(
        x = AA_mutation,
        y = factor(AA_pos_Ecoli, levels = pos_levels),
        fill = fill_group,
        tooltip = paste0(
          "Position: ", AA_pos_Ecoli,
          "<br>AA: ", AA_mutation,
          "<br>Source: ", source_combo
        ),
        data_id = paste(AA_pos_Ecoli, AA_mutation)
      )
    ) +
      ggiraph::geom_tile_interactive(color = "grey80", linewidth = 0.2) +
      scale_fill_manual(
        values = source_palette[display_levels],
        breaks = display_levels,
        limits = display_levels,
        drop = FALSE
      ) +
      scale_x_discrete(position = "top", drop = FALSE) +
      labs(x = "Amino acid", y = "Position") +
      theme_app(legend = "none", x_top = TRUE) +
      theme(
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
      )
    
    h_in <- heatmap_height_px() / 96
    
    ggiraph::girafe(
      ggobj = p,
      width_svg  = 12,
      height_svg = h_in,
      options = list(
        ggiraph::opts_hover(css = "fill-opacity:1; stroke:none;"),
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
  })
  
  
  ## ---- Heatmap container ----
  
  output$heatmap_container <- renderUI({
    df <- hm_plotting_df()
    
    if (is.null(df) || nrow(df) == 0) {
      return(ui_message("No data to show."))
    }
    
    tags$div(
      class = "hm-layout",
      
      # LEFT: plot
      tags$div(
        class = "hm-plot",
        girafeOutput(
          outputId = "mutationHeatmap",
          width  = "60vw",
          height = paste0(heatmap_height_px(), "px")
        )
        
      ),
      
      # RIGHT: legend
      tags$div(
        class = "hm-legend",
        uiOutput("heatmap_legend")
      )
    )
  })
  
  }
  
  
  
  # ---- Data ----
  
  table_df <- reactive({
    df <- cluster_plotting_df()
    req(!is.null(df), nrow(df) > 0)
    
    df |>
      dplyr::mutate(
        Gene = "rpoB",
        AA_mut_name = paste0(AA_original, AA_pos, AA_mutation),
        AA_mut_name_Ecoli = paste0(AA_original, AA_pos_Ecoli, AA_mutation)
        ) |>
      dplyr::rename(
        Species       = Species,
        AA_Pos        = AA_pos,
        AA_Pos_Ecoli  = AA_pos_Ecoli,
        Original      = AA_original,
        Mutation      = AA_mutation,
        Source        = source
      ) 
  })
  
  
  output$table_data <- DT::renderDataTable({
    DT::datatable({
      data <- table_df()
      data
    })
  })
  
  output$download_table <- downloadHandler(
    filename = function() paste0("RIFRMap_table_", Sys.Date(), ".csv"),
    content  = function(file) {
      readr::write_csv(table_df(), file)
    }
  )
  
  # ---- Summary ----
  output$summary_table <- render_gt({
    data <- table_df()
    
    data |>
      select(-c(AA_Pos, AA_Pos_Ecoli, AA_mut_name, AA_mut_name_Ecoli, Original, Mutation)) |>
      tbl_summary() |>
      as_gt()
  })

  
  output$bar_aa_mut_name <- renderPlot({
    df <- table_df()
    req(df, nrow(df) > 0, all(c("AA_mut_name_Ecoli") %in% names(df)))
    
    top_muts <- df |>
      dplyr::mutate(AA_mut_name_Ecoli = as.character(AA_mut_name_Ecoli)) |>
      dplyr::count(AA_mut_name_Ecoli, sort = TRUE) |>
      dplyr::slice_head(n = 15) |>
      dplyr::pull(AA_mut_name_Ecoli)
    
    plot_df <- df |>
      dplyr::mutate(AA_mut_name_Ecoli = as.character(AA_mut_name_Ecoli)) |>
      dplyr::filter(AA_mut_name_Ecoli %in% top_muts) |>
      dplyr::count(AA_mut_name_Ecoli)
    
    ggplot(plot_df, aes(x = reorder(AA_mut_name_Ecoli, n, FUN = sum), y = n)) +
      geom_col() +
      coord_flip() +
      theme_app(legend = "top") +
      theme(
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text  = element_text(size = 12)
      ) +
      labs(
        x = "AA mutation (E. coli coordinate)",
        y = "Number of observations"
      ) 
  })
  
  
  
  output$hist_aa_pos <- renderPlot({
    df <- table_df()
    req(df, nrow(df) > 0, "AA_Pos" %in% names(df))
    
    ggplot(df, aes(x = AA_Pos)) +
      geom_histogram(bins = 30) +
      theme_app() +
      theme(
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13)
      ) +
      labs(
        x = "Amino acid position",
        y = "Number of observations"
      )
  })
  
  output$hist_aa_pos_ecoli <- renderPlot({
    df <- table_df()
    req(df, nrow(df) > 0, "AA_Pos_Ecoli" %in% names(df))
    
    ggplot(df, aes(x = AA_Pos_Ecoli)) +
      geom_histogram(bins = 30) +
      theme_app() +
      theme(
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13)
      ) +
      labs(
        x = "Amino acid position (E. coli coordinate)",
        y = "Number of observations"
      )
  })
}

