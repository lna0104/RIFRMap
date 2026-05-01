#library(usethis)
#create_github_token()
#gitcreds_set()
library(tidyverse)
library(fs)
library(httr2)
library(rvest)
library(gitcreds)


# ----- Setup -----

token <- gitcreds_get()$password

# ----- Helper functions -----

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
  
  dat
}

get_pmid_from_doi <- function(doi) {
  
  doi <- stringr::str_extract(doi, "10\\..*")
  
  if (is.na(doi) || doi == "") {
    message("No valid DOI: ", raw_doi)
    return(NA_character_)
  }
  
  query <- paste0(doi, "[doi]")
  
  res <- tryCatch(
    rentrez::entrez_search(
      db = "pubmed",
      term = paste0(doi, "[doi]")
    ),
    error = function(e) NULL
  )

  
  if (is.null(res) || length(res$ids) == 0) {
    message("PMID not found for DOI: ", doi)
    return(NA_character_)
  }
  
  as.character(res$ids[1])
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
    "PMID"
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


# ----- Curated data -----

refs <- read_csv("./data/rpoB_references.csv", show_col_types = FALSE)
seqs <- readDNAStringSet("./data/rpoB_references.fasta")
load("data/coordinates.RData")

curated_mutations <- download_github_csv(
  repo = "JanEngelstaedter/ResistanceSurveys",
  path = "RIF/RIF_mutations_data.csv",
  token = token,
  filter_species = TRUE
)

curated_mutations_clean <- curated_mutations |>
  fillMutationsTable(refs, seqs, coordinates) |>
  filter(!is.na(AA_mut_name_Ecoli))

write_csv(curated_mutations_clean, "data/raw_data/RIF_curated_mutations_data.csv")

curated_references <- download_github_csv(
  repo = "JanEngelstaedter/ResistanceSurveys",
  path = "RIF/RIF_mutations_references.csv",
  token = token
)

curated_references_clean <- curated_references |>
  dplyr::mutate(
    DOI = stringr::str_trim(as.character(Link)),
    DOI = stringr::str_remove(DOI, "\\.$"),
    DOI = dplyr::na_if(DOI, ""),
    PMID = vapply(DOI, get_pmid_from_doi, FUN.VALUE = character(1)),
    PMID = dplyr::na_if(trimws(PMID), "")
  ) |>
  dplyr::select(-dplyr::any_of(c("Comments", "Link")))


write_csv(curated_references_clean, "data/raw_data/RIF_curated_mutations_references.csv")


# ----- CARD -----

download_and_extract_card()
snps <- read_tsv("data/raw_data/card_data/snps.txt", show_col_types = FALSE)

snps_clean <- snps |> 
  janitor::clean_names()

write_csv(snps_clean, "data/raw_data/CARD_mutations_data.csv")

# ----- DRAGdb -----

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

# ----- PointFinder -----

species_list <- c(
  "campylobacter",
  "enterococcus_faecalis",
  "escherichia_coli",
  "staphylococcus_faecium",
  "helicobacter_pylori",
  "klebsiella",
  "mycobacterium_tuberculosis",
  "neisseria_gonorrhoeae",
  "plasmodium_falciparum",
  "salmonella",
  "staphylococcus_aureus"
)

base_url <- "https://bitbucket.org/genomicepidemiology/pointfinder_db/raw/master/"

#dir.create("data/raw_data/pointfinder", showWarnings = FALSE)
input_dir <- "data/raw_data/pointfinder"


for (sp in species_list) {
  
  url <- paste0(base_url, sp, "/resistens-overview.txt")
  dest <- file.path(paste0(input_dir, "/", sp, "_resistens-overview.txt"))
  
  tryCatch({
    download.file(url, destfile = dest, mode = "wb")
    message("Downloaded: ", sp)
  }, error = function(e) {
    message("Failed: ", sp)
  })
}

files <- list.files(
  input_dir,
  pattern = "resistens-overview\\.txt$",
  full.names = TRUE
)

read_pointfinder <- function(file) {
  
  species <- basename(file) |>
    str_remove("_resistens-overview\\.txt$")
  
  df <- read_tsv(
    file,
    comment = "#",
    col_types = cols(.default = "c"),
    show_col_types = FALSE,
    trim_ws = TRUE
  )
  
  # Fix header if removed by comment="#"
  if (!"Gene_ID" %in% names(df)) {
    df <- read_tsv(
      file,
      col_types = cols(.default = "c"),
      show_col_types = FALSE,
      trim_ws = TRUE
    )
    
    names(df) <- str_replace(names(df), "^#", "")
  }
  
  df |>
    mutate(species = species)
}

# Merge all files
pointfinder_data <- map_dfr(files, read_pointfinder)
