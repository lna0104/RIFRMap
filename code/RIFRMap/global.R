library(shiny)
library(dplyr)
library(readr)
library(ggplot2)
library(Biostrings)
library(DT)
library(ggiraph)

# ---- Load data ----

card <- readRDS("data/card.rds")
cryptic <- readRDS("data/cryptic.rds")
dragdb <- readRDS("data/dragdb.rds")
yang <- readRDS("data/yang.rds")
curated <- readRDS("data/curated.rds")

refs_seqs <- read_csv("data/rpoB_references.csv", show_col_types = FALSE)
refs_paper <- read_csv("data/all_references.csv", show_col_types = FALSE)

clusters <- read_csv("data/clusters.csv", show_col_types = FALSE)
seqs <- readDNAStringSet("data/rpoB_references.fasta")
load("data/coordinates.RData")

# ----  Theme plot ----

PLOT_STYLE <- list(
  family = "Helvetica",
  base_size = 9,
  title_size = 14,
  axis_title = 11,
  axis_text = 10,
  legend_title = 12,
  legend_text = 10,
  message = 18,
  text_col = "grey20",
  axis_col = "grey20",
  tick_col = "grey60",
  panel_bg = "white",
  plot_margin = margin(6, 8, 6, 8)
)

plotting_theme <- function(legend = "none", x_top = FALSE) {
  theme_bw(
    base_family = PLOT_STYLE$family,
    base_size = PLOT_STYLE$base_size
  ) +
    theme(
      plot.title = element_text(
        size = PLOT_STYLE$title_size,
        colour = PLOT_STYLE$text_col
      ),
      axis.title.x = element_text(
        size = PLOT_STYLE$axis_title,
        colour = PLOT_STYLE$axis_col,
        margin = margin(t = if (x_top) 8 else 0, b = if (x_top) 0 else 8),
        face = "bold"
      ),
      axis.title.y = element_text(
        size = PLOT_STYLE$axis_title,
        colour = PLOT_STYLE$axis_col,
        margin = margin(r = 8),
        face = "bold"
      ),
      axis.text.x = element_text(
        size = PLOT_STYLE$axis_text,
        colour = PLOT_STYLE$axis_col
      ),
      axis.text.y = element_text(
        size = PLOT_STYLE$axis_text,
        colour = PLOT_STYLE$axis_col
      ),
      axis.ticks = element_line(colour = PLOT_STYLE$tick_col),
      panel.grid = element_blank(),
      plot.margin = PLOT_STYLE$plot_margin,
      legend.position = legend,
      legend.title = element_text(
        size = PLOT_STYLE$legend_title,
        face = "bold",
        margin = margin(b = 12)
      ),
      legend.text = element_text(
        size = PLOT_STYLE$legend_text
      )
    )
}

plotting_theme_large_axes <- function(legend = "none", x_top = FALSE) {
  plotting_theme(legend = legend, x_top = x_top) +
    theme(
      axis.title.x = element_text(
        size = 18,
        face = "bold",
        colour = PLOT_STYLE$axis_col,
        margin = margin(t = if (x_top) 8 else 0, b = if (x_top) 0 else 8)
      ),
      axis.title.y = element_text(
        size = 18,
        face = "bold",
        colour = PLOT_STYLE$axis_col,
        margin = margin(r = 8)
      ),
      legend.title = element_text(
        size = 16,
        face = "bold",
        margin = margin(b = 12)
      ),
      axis.text.x = element_text(
        size = 13,
        colour = PLOT_STYLE$axis_col
      ),
      axis.text.y = element_text(
        size = 13,
        colour = PLOT_STYLE$axis_col
      ),
      legend.text = element_text(
        size = 12
      )
    )
}

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
      round(PLOT_STYLE$message),
      PLOT_STYLE$text_col
    ),
    msg
  )
}

spin <- function(x) shinycssloaders::withSpinner(x, color = "#0dc5c1")

# ---- Filter data -----
# species filter ("All species" = no filter)
filter_by_species <- function(df, species_selected, species_col = "species") {
  if (is.null(species_selected) || "All species" %in% species_selected) {
    return(df)
  }
  df |> dplyr::filter(.data[[species_col]] %in% species_selected)
}

# helper to apply numeric range only if column exists
filter_range_if <- function(df, col, rng) {
  if (!col %in% names(df) || is.null(rng) || length(rng) != 2) {
    return(df)
  }
  df |> dplyr::filter(dplyr::between(.data[[col]], rng[1], rng[2]))
}

# helper to apply multi-select categorical only if column exists
filter_in_if <- function(df, col, selected) {
  if (!col %in% names(df) || is.null(selected) || length(selected) == 0 || selected == "") {
    return(df)
  }
  df |> dplyr::filter(.data[[col]] %in% selected)
}

# ----  Cluster plot ----

## ---- Legend panel ----

legend_reference <- c(
  Reference = "#000000"
)

legend_distribution <- c(
  Ubiquitous  = "#3e3bff",
  Polymorphic = "#12947F",
  Unique      = "#FF0000"
)


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
  clusterShift <- ALJEbinf::translateCoordinate(
    pos = clusterStart,
    coordinates = coordinatesFocal,
    direction = "RefToFocal",
    AAinput = TRUE,
    AAoutput = TRUE
  ) - clusterStart
  block <- muts |>
    filter(Species == species) |>
    mutate(AA_mut_name_Ecoli = paste0("rpoB", "_", AA_original, AA_pos_Ecoli, AA_mutation))
  block <- block |>
    filter(
      !is.na(AA_pos),
      !is.na(AA_pos_Ecoli),
      !is.na(AA_original),
      !is.na(AA_mutation)
    ) |>
    filter(AA_pos_Ecoli >= clusterStart & AA_pos_Ecoli <= clusterEnd) |>
    arrange(AA_pos_Ecoli, AA_mutation)

  if (mergeDuplicates) {
    block <- distinct(block, AA_pos_Ecoli, AA_mutation, .keep_all = TRUE)
  }

  block$RelX <- rep(NA_integer_, nrow(block))
  block$RelY <- rep(NA_integer_, nrow(block))
  y <- 1L
  for (i in seq_len(nrow(block))) {
    # the following line might need adjustment if there are indels within a cluster:
    block$RelX[i] <- block$AA_pos_Ecoli[i] - clusterStart + 1L
    if (i > 1) {
      if (block$AA_pos_Ecoli[i] == block$AA_pos_Ecoli[i - 1]) {
        y <- y + 1L
      } else {
        y <- 1L
      }
    }
    block$RelY[i] <- y
  }

  # add positions and original bases:
  addInfo <- data.frame(matrix(NA_integer_,
    ncol = ncol(block),
    nrow = 2 * clusterLength,
    dimnames = list(NULL, colnames(block))
  ))
  # again, the following line might need adjustment if there are indels within a cluster:
  addInfo$AA_pos <- (clusterStart:clusterEnd) + clusterShift
  addInfo$AA_pos_Ecoli <- clusterStart:clusterEnd
  addInfo$RelX <- 1:clusterLength
  addInfo$RelY <- c(rep(-1L, clusterLength), rep(0L, clusterLength))
  addInfo$AA_original <- seqs[[refSeq_ID]] |>
    Biostrings::translate(genetic.code = getGeneticCode("Bacterial", full.search = TRUE)) |>
    as.character() |>
    substr(
      clusterStart + clusterShift,
      clusterEnd + clusterShift
    ) |>
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
    cex = pars$cexPos, adj = c(0.5, 0.5), srt = 90, family = pars$familyPos
  )
  # original sequence:
  text(x0 + (block$RelX[block$RelY == 0L] - 1L) * dx,
    1 - y0 - dy * pars$dyPos,
    block$AA_original[block$RelY == 0L],
    cex = pars$cexRef, adj = c(0.5, 0.5), family = pars$familyRef, font = 2
  )
  # mutations:
  if (length(block$AA_mutation[block$RelY > 0L]) > 0) {
    text(x0 + (block$RelX[block$RelY > 0L] - 1L) * dx,
      1 - y0 - dy * (pars$dyPos + pars$dyRef) - (block$RelY[block$RelY > 0L] - 1L) * dy,
      block$AA_mutation[block$RelY > 0L],
      col = block$Colour[block$RelY > 0L],
      cex = pars$cexMut, adj = c(0.5, 0.5), family = pars$familyMuts
    )
  }
}


getPars <- function() {
  pars <- vector(mode = "list")
  pars$dxBlocks <- 0.5 # relative space between blocks

  pars$dyPos <- 1.2 # relative space for position numbers
  pars$dyRef <- 1.3 # relative space for reference sequence
  pars$dyBlocks <- 2.2

  pars$cexPos <- 0.9 # relative font size of position numbers
  pars$cexRef <- 0.9 # relative font size of reference sequence
  pars$cexMut <- 0.9
  pars$cexSpecies <- 1 # relative font size of species labels

  pars$col <- c(
    default = "#000000",
    ubiquitous = "#3e3bff",
    unique = "#ff3e3b",
    polymorphic = "#12947F"
  )

  pars$mar <- c(1, 1, 1, 1) + 0.1
  pars$familyPos <- PLOT_STYLE$family
  pars$familyMuts <- PLOT_STYLE$family
  pars$familyRef <- PLOT_STYLE$family

  pars$xInchPerUnit <- 0.25
  pars$yInchPerUnit <- 0.2
  return(pars)
}


plotMutations <- function(muts,
                          clusters,
                          seqs,
                          coordinates,
                          refs,
                          gene = "rpoB",
                          clustersToPlot = NULL,
                          speciesToPlot = NULL,
                          colouring = "none",
                          mergeDuplicates = TRUE,
                          return_dims = FALSE,
                          return_hover_map = FALSE,
                          res = 96,
                          mut_lines_cap = 12L) {
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
  n_species <- length(speciesToPlot)
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
        muts = muts,
        clusters = clusters,
        seqs = seqs,
        refSeq_ID = refSeq_ID,
        coordinatesFocal = coordinates[[refSeq_ID]],
        species = speciesToPlot[[i]],
        clusterName = clustersToPlot[[j]],
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
  width_px <- xTotal * pars$xInchPerUnit * res
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
        x0 = x0c, # center of RelX==1
        n = n_units, # number of AA columns in this cluster
        xmin = xmin,
        xmax = xmax,
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

        if (colouring == "ubi") cat[cat == "unique"] <- "default"
        if (colouring == "uni") cat[cat == "ubiquitous"] <- "default"
        if (colouring == "none") cat[] <- "default"

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

# ---- Mutation Heatmap ----
source_levels <- c("CARD", "CRyPTIC", "DRAGdb", "Curated", "Yang", "None")
all_AAs_sorted <- sort(strsplit("HAKRGILPVFWYDESTCMNQ", "")[[1]])

source_palette <- c(
  CARD    = "#ead6fd",
  CRyPTIC = "#ffdfba",
  DRAGdb  = "#bae1ff",
  Curated = "#fff8f0",
  Yang    = "#dff2e1",
  None    = "grey85"
)
