## ---- CRAN ----
cran_pkgs <- c(
  "shiny",
  "dplyr",
  "tidyr",
  "tidyverse",
  "ggiraph",
  "ggplot2",
  "ggupset",
  "htmltools",
  "remotes",
  "shinyWidgets",
  "shinycssloaders",
  "gtsummary", 
  "gt",
  "wesanderson",
  "RColorBrewer",
  "grantham"
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
  setdiff(c("Biostrings", "DECIPHER", "BiocGenerics", "IRanges", "GenomeInfoDb", "S4Vectors", "UCSC.utils"), rownames(installed.packages())),
  ask = FALSE,
  force = TRUE,
  update = FALSE
)


## ---- GitHub ----
if (!requireNamespace("ALJEbinf", quietly = TRUE)) {
  remotes::install_github("JanEngelstaedter/ALJEbinf", build_vignettes = TRUE)
}

