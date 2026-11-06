## ---- CRAN ----
cran_pkgs <- c(
  "shiny",
  "shinyjs",
  "tidyverse", # dplyr, tidyr, readr, ggplot2, stringr, purrr
  "ggiraph",
  "ggupset",
  "ggpubr",
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
  "grantham",
  "rentrez",
  "httr2",
  "gitcreds",
  "fs",
  "janitor",
  "openssl",
  "R.utils"
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
# remotes is installed with the CRAN packages above.
if (!requireNamespace("ALJEbinf", quietly = TRUE)) {
  remotes::install_github(
    "JanEngelstaedter/ALJEbinf"
  )
}
