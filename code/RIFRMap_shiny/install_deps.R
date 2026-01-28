## ---- CRAN ----
cran_pkgs <- c(
  "bs4Dash",
  "shiny",
  "dplyr",
  "tidyr",
  "tidyverse",
  "ggiraph",
  "bslib",
  "ggplot2",
  "ggupset",
  "htmltools",
  "remotes",
  "shinyWidgets",
  "shinycssloaders",
  "gtsummary", 
  "gt",
  "DT",
  "wesanderson",
  "RColorBrewer",
  "fresh"

)

install.packages(
  setdiff(cran_pkgs, rownames(installed.packages())),
  repos = "https://cloud.r-project.org"
)

## ---- Bioconductor ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

BiocManager::install(
  setdiff(c("Biostrings", "DECIPHER"), rownames(installed.packages())),
  ask = FALSE,
  update = FALSE
)

## ---- GitHub ----
if (!requireNamespace("ALJEbinf", quietly = TRUE)) {
  remotes::install_github("JanEngelstaedter/ALJEbinf", build_vignettes = TRUE)
}

