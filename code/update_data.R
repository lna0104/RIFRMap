library(tidyverse)
library(fs)
library(httr2)
library(rvest)
library(gitcreds)

# -------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------

token <- gitcreds_get()$password

dir_create("data/raw_data")
dir_create("data/raw_data/card_data")

# -------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------

download_github_csv <- function(repo, path, output_path, token,
                                filter_species = FALSE,
                                modify_fun = NULL) {
  repo_api <- paste0(
    "https://api.github.com/repos/",
    repo,
    "/contents/",
    path
  )
  
  resp <- request(repo_api) |>
    req_headers(
      Authorization = paste("Bearer", token),
      Accept = "application/vnd.github+json"
    ) |>
    req_perform()
  
  info <- resp_body_json(resp)
  
  dat <- read_csv(info$download_url, show_col_types = FALSE)
  
  if (filter_species && "Species" %in% names(dat)) {
    dat <- dat |> filter(!is.na(Species))
  }
  
  dat <- dat |> select(where(~ !all(is.na(.))))
  
  if (!is.null(modify_fun)) {
    dat <- modify_fun(dat)
  }
  
  dir_create(path_dir(output_path))
  write_csv(dat, output_path)
  
  message("Saved: ", output_path)
  invisible(dat)
}

download_and_extract_card <- function(
    url = "https://card.mcmaster.ca/latest/data",
    dest_file = "data/card-data.tar.bz2",
    extract_dir = "data/raw_data/card_data"
) {
  dir_create(path_dir(dest_file))
  dir_create(extract_dir)
  
  download.file(
    url = url,
    destfile = dest_file,
    mode = "wb",
    method = "libcurl"
  )
  
  untar(dest_file, exdir = extract_dir)
  
  if (file_exists(dest_file)) {
    file_delete(dest_file)
  }
  
  message("CARD data downloaded and extracted to: ", extract_dir)
}

get_drag_csv <- function(gene, save_file = NULL) {
  url <- paste0(
    "https://bicresources.jcbose.ac.in/ssaha4/drag/browse_gene.php?query=",
    utils::URLencode(gene, reserved = TRUE)
  )
  
  resp <- request(url) |>
    req_user_agent("Mozilla/5.0") |>
    req_perform()
  
  txt <- rawToChar(resp_body_raw(resp))
  href <- stringr::str_extract(txt, 'Temp/Blast/[^"]+\\.csv')
  
  if (is.na(href)) {
    stop("Could not find CSV download link in the DRAGdb page source.")
  }
  
  file_url <- paste0("https://bicresources.jcbose.ac.in/ssaha4/drag/", href)
  
  if (!is.null(save_file)) {
    dir_create(path_dir(save_file))
    download.file(file_url, destfile = save_file, mode = "wb", method = "libcurl")
    return(read_csv(save_file, show_col_types = FALSE))
  }
  
  read_csv(file_url, show_col_types = FALSE)
}

rename_drag_columns <- function(dat) {
  cols <- c(
    "gene",
    "ensembl_id",
    "species",
    "drug",
    "Nt_pos",
    "Nt_mutation",
    "AA_pos",
    "codon_mutation",
    "mutation_type",
    "provean",
    "sequencing_method",
    "strain_or_isolates",
    "sample_location",
    "pubmed_id"
  )
  
  if (ncol(dat) != length(cols)) {
    stop(
      "Column number mismatch in DRAGdb data. Expected ",
      length(cols),
      " columns but found ",
      ncol(dat),
      "."
    )
  }
  
  colnames(dat) <- cols
  dat
}

# -------------------------------------------------------------------
# Curated data
# -------------------------------------------------------------------

curated_mutations <- download_github_csv(
  repo = "JanEngelstaedter/ResistanceSurveys",
  path = "RIF/RIF_mutations_data.csv",
  output_path = "data/raw_data/RIF_curated_mutations_data.csv",
  token = token,
  filter_species = TRUE
)

curated_references <- download_github_csv(
  repo = "JanEngelstaedter/ResistanceSurveys",
  path = "RIF/RIF_mutations_references.csv",
  output_path = "data/raw_data/RIF_curated_mutations_references.csv",
  token = token,
  modify_fun = function(dat) {
    dat |>
      rename(DOI = Link) |>
      select(-Comments)
  }
)

# -------------------------------------------------------------------
# CARD
# -------------------------------------------------------------------

download_and_extract_card()
snps <- read_tsv("data/raw_data/card_data/snps.txt", show_col_types = FALSE)

snps_clean <- snps |> 
  janitor::clean_names()

write_csv(snps_clean, "data/raw_data/CARD_mutations_data.csv")

# -------------------------------------------------------------------
# DRAGdb
# -------------------------------------------------------------------

dragdb_raw_data <- get_drag_csv(
  gene = "rpoB",
  save_file = "data/raw_data/RIF_DRAGdb_mutations_data_raw.csv"
)

dragdb_clean <- dragdb_raw_data |>
  rename_drag_columns()

write_csv(
  dragdb_clean,
  "data/raw_data/RIF_DRAGdb_mutations_data.csv"
)

