## ---- CRAN ----
cran_pkgs <- c(
  "shiny",
  "shinyjs",
  "dplyr",
  "tidyr",
  "tidyverse",
  "ggiraph",
  "ggplot2",
  "ggupset",
  "plotly",
  "DT",
  "scales",
  "viridisLite",
  "htmltools",
  "remotes",
  "shinyWidgets",
  "shinycssloaders",
  "gtsummary",
  "gt",
  "wesanderson",
  "RColorBrewer",
  "grantham",
  "rentrez",
  "httr2",
  "rvest",
  "gitcreds",
  "fs",
  "janitor",
  "xml2",
  "openssl"
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
  setdiff(c("Biostrings", "DECIPHER", "BiocGenerics", "IRanges", "GenomeInfoDb", "S4Vectors", "UCSC.utils", "XVector", "pwalign"), rownames(installed.packages())),
  ask = FALSE,
  force = TRUE,
  update = FALSE
)


## ---- GitHub ----
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages(
    "remotes",
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("ALJEbinf", quietly = TRUE)) {
  remotes::install_github(
    "JanEngelstaedter/ALJEbinf"
  )
}
