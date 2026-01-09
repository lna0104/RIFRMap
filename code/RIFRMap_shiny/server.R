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
  
  
  plot_message <- function(
    msg,
    x = 0.02,
    y = 0.98,
    size_mult = 1.1,
    col = PLOT_STYLE$text_col,
    fontface = "plain"
  ) {
    grid::grid.newpage()
    grid::grid.text(
      msg,
      x = grid::unit(x, "npc"),
      y = grid::unit(y, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(
        col        = col,
        fontsize   = PLOT_STYLE$base_size * 12/9 * size_mult,
        fontfamily = PLOT_STYLE$family,
        fontface   = fontface
      )
    )
  }
  
  
  
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
    updateTabItems(session, "sidebar", "view")
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
      pars$dxBlocks <- 0.8  # relative space between blocks
      
      pars$dyPos <- 1.1 # relative space for position numbers
      pars$dyRef <- 1.3 # relative space for reference sequence
      pars$dyBlocks <- 2.2
      
      pars$cexPos <- 0.9 # relative font size of position numbers
      pars$cexRef <- 0.9 # relative font size of reference sequence
      pars$cexMut <- 0.9
      pars$cexSpecies <- 1.1 # relative font size of species labels
      
      pars$col <- c(default     = "#000000",
                    ubiquitous  = "#3e3bff", 
                    unique      = "#ff3e3b",
                    polymorphic = "#12947F")
      
      pars$mar <- c(1, 1, 1, 1) + 0.1
      pars$familyPos  <- PLOT_STYLE$family
      pars$familyMuts <- PLOT_STYLE$family
      pars$familyRef  <- PLOT_STYLE$family
      
      pars$xInchPerUnit <- 0.3
      pars$yInchPerUnit <- 0.2
      return(pars)
    }
    
    #Create hover map for mouse movement
    block_hover_map <- function(block, x0, y0, dx, dy, pars, species, cluster) {
      out <- list()
      scale <- 1.3
      
      # AA positions (RelY == -1)
      idx <- which(block$RelY == -1L)
      if (length(idx)) {
        x <- x0 + (block$RelX[idx] - 1L) * dx
        y <- 1 - y0
        for (k in seq_along(idx)) {
          out[[length(out)+1]] <- data.frame(
            species = species, cluster = cluster,
            kind = "pos",
            label = as.character(block$AA_pos[idx[k]]),
            xmin = x[k] - dx/3, xmax = x[k] + dx/3,
            ymin = y - dy/2 * scale, ymax = y + dy/2 * scale,
            stringsAsFactors = FALSE
          )
        }
      }
      
      # reference letters (RelY == 0)
      idx <- which(block$RelY == 0L)
      if (length(idx)) {
        x <- x0 + (block$RelX[idx] - 1L) * dx
        y <- 1 - y0 - dy * pars$dyPos
        for (k in seq_along(idx)) {
          out[[length(out)+1]] <- data.frame(
            species = species, cluster = cluster,
            kind = "ref",
            label = as.character(block$AA_original[idx[k]]),
            xmin = x[k] - dx/3, xmax = x[k] + dx/3,
            ymin = y - dy/2, ymax = y + dy/2,
            stringsAsFactors = FALSE
          )
        }
      }
      
      # mutations (RelY > 0)
      idx <- which(block$RelY > 0L)
      if (length(idx)) {
        x <- x0 + (block$RelX[idx] - 1L) * dx
        y <- 1 - y0 - dy * (pars$dyPos + pars$dyRef) - (block$RelY[idx] - 1L) * dy
        for (k in seq_along(idx)) {
          out[[length(out)+1]] <- data.frame(
            species = species, cluster = cluster,
            kind = "mut",
            label = as.character(block$AA_mutation[idx[k]]),
            colour = as.character(block$Colour[idx[k]]),
            xmin = x[k] - dx/3, xmax = x[k] + dx/3,
            ymin = y - dy/2, ymax = y + dy/2,
            stringsAsFactors = FALSE
          )
        }
      }
      
      if (!length(out)) return(NULL)
      dplyr::bind_rows(out)
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
      res = 96
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
      
      # calculate width and height for UI sizing
      width_px  <- xTotal * pars$xInchPerUnit * res
      height_px <- yTotal * pars$yInchPerUnit * res
      
      if (isTRUE(return_dims)) {
        return(list(width_px = width_px, 
                    height_px = height_px, 
                    xTotal = xTotal, 
                    yTotal = yTotal, 
                    pars = pars))
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
      
      hover_map <- tibble::tibble(
        species = character(),
        cluster = character(),
        kind    = character(),
        xmin    = numeric(),
        xmax    = numeric(),
        ymin    = numeric(),
        ymax    = numeric()
      )
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
        x0 <- pars$dxBlocks * dx
        
        for (j in seq_len(n_clusters)) {
          if (j > 1) {
            x0 <- x0 + dx * (colMax[j - 1] + pars$dxBlocks)
          }
          renderBlock(allBlocks[[i]][[j]], x0, y0, dx, dy, pars = pars)
          hm <- block_hover_map(
            block   = allBlocks[[i]][[j]],
            x0      = x0, 
            y0 = y0, 
            dx = dx, 
            dy = dy,
            pars    = pars,
            species = speciesToPlot[[i]],
            cluster = clustersToPlot[[j]]
          )
          # if (!is.null(hm)) {
          #   with(hm, rect(xmin, ymin, xmax, ymax, border = "red"))
          # }
          # 
          if (!is.null(hm)) {
            hover_map <- dplyr::bind_rows(hover_map, hm)
          }
        }
      }
      
      #print(hover_map)
      
      if (isTRUE(return_hover_map)) {
        return(list(
          hover_map = hover_map,
          plt = plt
        ))
      }
      
      return(invisible(NULL))
      
    }
    
    
    # ---- Process data ----

    
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
    
    # select species that has refs + seq + coords 
    has_ref_assets <- function(species, refs, seqs, coordinates) {
      ref_row <- refs[refs$Species == species, , drop = FALSE]
      if (nrow(ref_row) == 0) return(FALSE)
      
      ref_id <- ref_row$FASTA_name[1]
      if (is.na(ref_id) || ref_id == "") return(FALSE)
      
      if (!(ref_id %in% names(seqs))) return(FALSE)
      if (is.null(coordinates[[ref_id]])) return(FALSE)
      
      TRUE
    }
    
    # Species to plot (drop those without assets) 
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
    
    
    
    # Return height and width 
    cluster_dims <- reactive({
      df <- cluster_plotting_df()
      if (is.null(df) || nrow(df) == 0) return(NULL)
      if (no_clusters_selected(input$clustersToPlot)) return(NULL)
      
      plotMutations(
        muts           = df,
        clusters       = clusters,
        seqs           = seqs,
        coordinates    = coordinates,
        refs           = refs,
        gene           = "rpoB",
        clustersToPlot = input$clustersToPlot,
        speciesToPlot  = cluster_species_to_plot(),
        colouring      = "ubiuni",
        mergeDuplicates= TRUE,
        plot_width_px  = NA,     
        return_dims    = TRUE
      )
    })
    
    
    no_clusters_selected <- function(x) {
      is.null(x) || length(x) == 0 || all(x == "")
    }
    
    cluster_trigger <- reactive({
      input$apply_filters
      input$view_tabs
      TRUE
    })
    
    # ---- Plot ----
    
    ## ---- Legend panel ----
    
    legend_colour <- c(
      Reference   = "#000000",
      Ubiquitous  = "#3e3bff", 
      Unique      = "#FF0000",
      Polymorphic = "#12947F"
    )
    
    output$legend_panel <- renderUI({
      tags$div(
        class = "cluster-legend",
        
        tags$div(
          class = "cluster-legend-title",
          "Mutation distribution across species"
        ),
        
        lapply(names(legend_colour), function(name) {
          tags$div(
            class = "cluster-legend-item",
            
            tags$span(
              class = "cluster-legend-swatch",
              style = sprintf("background-color: %s;", legend_colour[[name]])
            ),
            
            tags$span(
              class = "cluster-legend-label",
              tools::toTitleCase(name)
            )
          )
        })
      )
    })
    
    
    ## ---- Cluster plot -----
    
    hover_map_r <- reactiveVal(NULL)
    plt_r       <- reactiveVal(NULL)
    
    
    output$clusterPlotUI <- renderUI({
      div(
        id = "clusterWrap",
        div(
          id = "clusterPlotWrap",
          plotOutput("clusterPlot", width = "100%", height = "85vh"),
          tags$canvas(id = "hoverCanvas")
        )
      )
    })
    

    output$clusterPlot <- renderPlot({
      cluster_trigger() 
      
      df <- cluster_plotting_df()
      
      if (is.null(df) || nrow(df) == 0) {
        hover_map_r(NULL)
        plt_r(NULL)
        
        session$sendCustomMessage("hoverMap", "NULL")
        
        plot_message("No data to show.")
        return(invisible())
      }
      
      if (no_clusters_selected(input$clustersToPlot)) {
        hover_map_r(NULL)
        plt_r(NULL)
        
        session$sendCustomMessage("hoverMap", "NULL")
        
        plot_message("Please select at least one cluster.")
        return(invisible())
        
      }
      
      clustersToPlot <- clusters |> 
        dplyr::filter(Gene == "rpoB", Cluster %in% input$clustersToPlot) |> 
        dplyr::arrange(AAstart_Ecoli) |> dplyr::pull(Cluster) 
      
      
      out <- plotMutations(
        muts = df,
        clusters = clusters,
        seqs = seqs,
        coordinates = coordinates,
        refs = refs,
        gene = "rpoB",
        clustersToPlot = clustersToPlot,
        speciesToPlot  = cluster_species_to_plot(),
        colouring = "ubiuni",
        mergeDuplicates = TRUE,
        return_hover_map = TRUE
      )
      
      # capture geometry from this plot
      hover_map_r(out$hover_map)
      plt_r(out$plt)
      
      },
      width  = function() {
        d <- cluster_dims()
        if (is.null(d)) 300 else max(300, d$width_px)
      },
      height = function() {
        d <- cluster_dims()
        if (is.null(d)) 300 else max(150, d$height_px)
      },
      res = 96
    )
    
    # send to js 
    observeEvent(hover_map_r(), {
      
      req(hover_map_r(), plt_r())
      
      session$sendCustomMessage("hoverMap", list(
        rows = hover_map_r(),
        plt  = plt_r()
      ))
      
    }, ignoreInit = TRUE)
  }
  
  
  # ---- Source Tab ----
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
            plotOutput("mutationUpset", height = "450px", width = "100%")
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
    
    keep_levels <- names(sel)[sel]
    
    # if something selected → others become "None"
    if (length(keep_levels) > 0) {
      df$fill_group <- ifelse(df$source_combo %in% keep_levels,
                              df$source_combo, "None")
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
    if (is.null(df0) || nrow(df0) == 0) {
      return(ui_message("No data to show."))
    }
    
    df <- df0 |>
      dplyr::transmute(
        Species      = as.character(species),
        AA_original  = as.character(AA_original),
        AA_pos       = suppressWarnings(as.integer(AA_pos)),
        AA_pos_Ecoli = suppressWarnings(as.integer(AA_pos_Ecoli)),
        AA_mutation  = as.character(AA_mutation),
        Source       = as.character(source))
    
    if (is.null(df) || nrow(df) == 0) {
      return(ui_message("No data to show."))
    }
    
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
