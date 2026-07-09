# library(usethis)
# create_github_token()
# gitcreds_set()
library(tidyverse)
library(fs)
library(httr2)
library(rvest)
library(gitcreds)
library(Biostrings)

source("helpers.R")

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
    message("No valid DOI: ", doi)
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
  extract_dir = "raw_data/card_data"
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

list_pointfinder_species <- function(exclude = c("docs")) {
  api_url <- "https://api.bitbucket.org/2.0/repositories/genomicepidemiology/pointfinder_db/src/master/"
  species <- character(0)
  url <- api_url

  while (!is.null(url)) {
    resp <- request(url) |> req_perform()
    info <- resp_body_json(resp)
    # Keep only directory and ignore files
    dirs <- Filter(function(v) identical(v$type, "commit_directory"), info$values)
    species <- c(species, vapply(dirs, function(v) v$path, character(1)))

    url <- info[["next"]]
  }

  setdiff(species, exclude)
}

read_pointfinder <- function(file) {
  species <- basename(file) |>
    str_remove("_resistens-overview\\.txt$")

  # Read raw lines
  lines <- readLines(file)

  # Find "# Indels" line
  indel_idx <- grep("^# Indels", lines)

  # Keep only lines BEFORE it
  if (length(indel_idx) > 0) {
    lines <- lines[1:(indel_idx[1] - 1)]
  }

  # Read cleaned text
  df <- read_tsv(
    I(lines),
    col_types = cols(.default = "c"),
    show_col_types = FALSE,
    trim_ws = TRUE
  )

  # Clean header
  names(df) <- str_replace(names(df), "^#", "")

  df |>
    mutate(species = species)
}

# ----- Per-source update functions -----
# Each function performs its own retrieval, validates the result, and returns errors
# list(dataset=, status=c("ok","error"), message=, path=, rows=), so the
# orchestrator (update_all.R) can continue to the next source regardless of
# one source's failure.

update_card <- function() {
  result <- list(
    dataset = "card", status = "error", message = NA_character_,
    path = "raw_data/card_data.csv", rows = NA_integer_
  )
  tryCatch(
    {
      retry(download_and_extract_card(), max_tries = 3)
      snps <- read_tsv("raw_data/card_data/snps.txt", show_col_types = FALSE)
      snps_clean <- snps |> janitor::clean_names()

      if (nrow(snps_clean) == 0) stop("CARD download produced zero rows")

      write_csv(snps_clean, "raw_data/card_data.csv")

      result$status <- "ok"
      result$message <- "Downloaded and extracted CARD data"
      result$rows <- nrow(snps_clean)
    },
    error = function(e) {
      result$message <<- conditionMessage(e)
    }
  )
  result
}

update_dragdb <- function() {
  result <- list(
    dataset = "dragdb", status = "error", message = NA_character_,
    path = "raw_data/dragdb_data.csv", rows = NA_integer_
  )
  tryCatch(
    {
      dragdb_raw_data <- retry(
        get_drag_csv(gene = "rpoB", save_file = "raw_data/dragdb_raw.csv"),
        max_tries = 3
      )
      dragdb_clean <- dragdb_raw_data |> rename_drag_columns()

      if (nrow(dragdb_clean) == 0) stop("DRAGdb download produced zero rows")

      write_csv(dragdb_clean, "raw_data/dragdb_data.csv")

      result$status <- "ok"
      result$message <- "Downloaded DRAGdb data"
      result$rows <- nrow(dragdb_clean)
    },
    error = function(e) {
      result$message <<- conditionMessage(e)
    }
  )
  result
}

update_pointfinder <- function() {
  result <- list(
    dataset = "pointfinder", status = "error", message = NA_character_,
    path = "raw_data/pointfinder_data.csv", rows = NA_integer_
  )
  tryCatch(
    {
      species_list <- retry(list_pointfinder_species(), max_tries = 3)

      if (length(species_list) == 0) {
        stop("Could not list any species folders from PointFinder's Bitbucket repo")
      }

      base_url <- "https://bitbucket.org/genomicepidemiology/pointfinder_db/raw/master/"
      input_dir <- "./raw_data/pointfinder_data"
      dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)

      failed_species <- character(0)
      for (sp in species_list) {
        url <- paste0(base_url, sp, "/resistens-overview.txt")
        dest <- file.path(paste0(input_dir, "/", sp, "_resistens-overview.txt"))

        tryCatch(
          {
            retry(download.file(url, destfile = dest, mode = "wb"), max_tries = 3)
            message("Downloaded: ", sp)
          },
          error = function(e) {
            message("Failed: ", sp, " (", conditionMessage(e), ")")
            failed_species <<- c(failed_species, sp)
          }
        )
      }

      files <- list.files(
        input_dir,
        pattern = "resistens-overview\\.txt$",
        full.names = TRUE
      )

      if (length(files) == 0) stop("No PointFinder species downloaded successfully")

      pointfinder_data <- map_dfr(files, read_pointfinder)
      pointfinder_clean_data <- pointfinder_data |>
        filter(
          Gene_ID == "rpoB",
          str_detect(Resistance, regex("rifampicin", ignore_case = TRUE)),
          !Ref_codon %in% c("del", "ins")
        ) |>
        separate_rows(Res_codon, sep = ",") |>
        mutate(
          Res_codon = str_trim(Res_codon),
          species = species |>
            str_replace_all("_", " ") |>
            str_squish(),
          species = paste(
            str_to_title(word(species, 1)),
            str_to_lower(word(species, 2))
          ),
          mutation = paste0(Ref_codon, Codon_pos, Res_codon)
        )

      if (nrow(pointfinder_clean_data) == 0) stop("PointFinder produced zero rifampicin/rpoB rows")

      write_csv(
        pointfinder_clean_data,
        "raw_data/pointfinder_data.csv"
      )

      result$status <- "ok"
      result$message <- if (length(failed_species) > 0) {
        paste0(
          "Downloaded with ", length(failed_species), " species failing: ",
          paste(failed_species, collapse = ", ")
        )
      } else {
        "Downloaded all species successfully"
      }
      result$rows <- nrow(pointfinder_clean_data)
    },
    error = function(e) {
      result$message <<- conditionMessage(e)
    }
  )
  result
}

update_refgenes <- function() {
  result <- list(
    dataset = "refgenes", status = "error", message = NA_character_,
    path = "raw_data/refgenes_catalog.csv", rows = NA_integer_
  )
  tryCatch(
    {
      catalog_url <- "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt"

      raw_catalog <- retry(
        read_tsv(catalog_url, show_col_types = FALSE, col_types = cols(.default = "c")),
        max_tries = 3
      )

      expected_cols <- c(
        "allele", "gene_family", "product_name", "scope", "type", "subtype",
        "class", "subclass", "refseq_protein_accession", "refseq_nucleotide_accession",
        "genbank_protein_accession", "genbank_nucleotide_accession",
        "curated_refseq_start", "pubmed_reference"
      )
      missing_cols <- setdiff(expected_cols, names(raw_catalog))
      if (length(missing_cols) > 0) {
        stop(
          "ReferenceGeneCatalog.txt missing expected column(s): ",
          paste(missing_cols, collapse = ", ")
        )
      }

      # Renamed to match the column names process_data.R's parse_refgenes()
      # already expects (the pre-existing raw_data/refgenes_catalog.csv was a
      # manually-exported subset/rename of this same NCBI AMRFinderPlus file).
      refgenes_clean <- raw_catalog |>
        dplyr::select(
          `#Allele`              = allele,
          `Gene family`          = gene_family,
          `Product name`         = product_name,
          `Scope`                = scope,
          `Type`                 = type,
          `Subtype`              = subtype,
          `Class`                = class,
          `Subclass`             = subclass,
          `RefSeq protein`       = refseq_protein_accession,
          `RefSeq nucleotide`    = refseq_nucleotide_accession,
          `GenBank protein`      = genbank_protein_accession,
          `GenBank nucleotide`   = genbank_nucleotide_accession,
          `Curated RefSeq start` = curated_refseq_start,
          `Pubmed reference`     = pubmed_reference
        )

      if (nrow(refgenes_clean) == 0) stop("Reference Gene Catalog download produced zero rows")

      write_csv(refgenes_clean, "raw_data/refgenes_catalog.csv")

      result$status <- "ok"
      result$message <- "Downloaded NCBI AMRFinderPlus Reference Gene Catalog"
      result$rows <- nrow(refgenes_clean)
    },
    error = function(e) {
      result$message <<- conditionMessage(e)
    }
  )
  result
}

# RifReport's Species column sometimes names more than one species in a
# single cell - a parenthetical synonym ("Rhodococcus equi (Prescottella
# equi)"), a slash-separated pair ("Bacillus anthracis/Bacillus cereus"), or
# a comma-separated list ("Chlamydia pneumoniae, Chlamydia trachomatis").
# Left as-is, none of these match any single entry in
# raw_data/rpoB_references.csv, so fillMutationsTable() silently drops the
# row. Expand each such cell into one row per named species (duplicating the
# rest of the row) so every species gets its own chance to map.
expand_multi_species <- function(df, species_col = "Species") {
  split_species <- function(s) {
    if (is.na(s)) {
      return(NA_character_)
    }

    paren_match <- stringr::str_match(s, "^(.*?)\\s*\\(([^)]+)\\)\\s*$")
    if (!is.na(paren_match[1, 1])) {
      main <- paren_match[1, 2]
      alt <- paren_match[1, 3]
    } else {
      main <- s
      alt <- NULL
    }

    main_parts <- main |>
      stringr::str_split("\\s*[/,]\\s*") |>
      unlist() |>
      stringr::str_trim()

    all_parts <- c(main_parts, alt)
    all_parts[nzchar(all_parts)]
  }

  df[[".species_list"]] <- lapply(df[[species_col]], split_species)

  df |>
    tidyr::unnest_longer(.species_list) |>
    dplyr::mutate("{species_col}" := .species_list) |>
    dplyr::select(-.species_list)
}

# RifReport is split into two steps because mapping to E. coli coordinates
# (fillMutationsTable) needs `coordinates`, which is only produced by
# process_refseqs.R's reference-discovery/coordinate logic - and that logic
# in turn needs to know which species are referenced across all 7 datasets
# (including the RifReport species discovered here) before it runs. See
# update_all.R for how the two steps are sequenced around that dependency.

update_rifreport_download <- function(token) {
  result <- list(
    dataset = "rifreport", status = "error", message = NA_character_,
    path = "raw_data/rifreport_data.csv", rows = NA_integer_,
    raw_mutations = NULL, species = character(0)
  )
  tryCatch(
    {
      rifreport_mutations <- retry(
        download_github_csv(
          repo = "JanEngelstaedter/ResistanceSurveys",
          path = "RIF/RIF_mutations_data.csv",
          token = token,
          filter_species = TRUE
        ),
        max_tries = 3
      )

      if (is.null(rifreport_mutations) || nrow(rifreport_mutations) == 0) {
        stop("No RifReport mutation rows downloaded")
      }
      if (!"Species" %in% names(rifreport_mutations)) {
        stop("RifReport download missing expected 'Species' column")
      }

      rifreport_mutations_clean <- rifreport_mutations |>
        mutate(Gene = coalesce(Gene, "rpoB")) |>
        expand_multi_species()

      rifreport_references <- retry(
        download_github_csv(
          repo = "JanEngelstaedter/ResistanceSurveys",
          path = "RIF/RIF_mutations_references.csv",
          token = token
        ),
        max_tries = 3
      )

      rifreport_references_clean <- rifreport_references |>
        dplyr::mutate(
          DOI = stringr::str_trim(as.character(Link)),
          DOI = stringr::str_remove(DOI, "\\.$"),
          DOI = dplyr::na_if(DOI, ""),
          PMID = vapply(DOI, get_pmid_from_doi, FUN.VALUE = character(1)),
          PMID = dplyr::na_if(trimws(PMID), "")
        ) |>
        dplyr::select(-dplyr::any_of(c("Comments", "Link")))

      write_csv(rifreport_references_clean, "raw_data/rifreport_references.csv")

      result$status <- "ok"
      result$message <- "Downloaded RifReport mutations and references"
      result$rows <- nrow(rifreport_mutations_clean)
      result$raw_mutations <- rifreport_mutations_clean
      result$species <- rifreport_mutations_clean$Species |>
        unique() |>
        (\(x) x[!is.na(x) & nzchar(x)])()
    },
    error = function(e) {
      result$message <<- conditionMessage(e)
    }
  )
  result
}

update_rifreport_finalize <- function(raw_mutations, refs, seqs, coordinates) {
  result <- list(
    dataset = "rifreport", status = "error", message = NA_character_,
    path = "raw_data/rifreport_data.csv", rows = NA_integer_
  )
  tryCatch(
    {
      if (is.null(coordinates)) {
        stop("coordinates not available - run process_refseqs.R's coordinate logic first")
      }

      rifreport_mutations_clean <- raw_mutations |>
        fillMutationsTable(refs, seqs, coordinates) |>
        filter(!is.na(AA_mut_name_Ecoli))

      if (nrow(rifreport_mutations_clean) == 0) {
        stop("No RifReport rows could be mapped to E. coli coordinates")
      }

      write_csv(rifreport_mutations_clean, "raw_data/rifreport_data.csv")

      result$status <- "ok"
      result$message <- "Mapped RifReport mutations to E. coli coordinates"
      result$rows <- nrow(rifreport_mutations_clean)
    },
    error = function(e) {
      result$message <<- conditionMessage(e)
    }
  )
  result
}

# ---- Standalone execution ----
# When this file is sourced directly (not via update_all.R), replicate
# today's full refresh: run process_refseqs.R's coordinate logic first (so
# `coordinates` exists - this also fixes the previously-hidden ordering bug
# for standalone use), then run every update_*() in sequence, continuing
# past individual failures and printing a result per dataset.
# update_all.R sets `.pipeline_orchestrating` before sourcing this file, so
# this block is skipped there and only the functions above get defined.
if (!exists(".pipeline_orchestrating", inherits = FALSE)) {
  source("process_refseqs.R")

  refs <- read_csv("raw_data/rpoB_references.csv", show_col_types = FALSE)
  seqs <- readDNAStringSet("raw_data/rpoB_references.fasta")

  results <- list()
  results$card <- update_card()
  results$dragdb <- update_dragdb()
  results$pointfinder <- update_pointfinder()
  results$refgenes <- update_refgenes()

  token <- tryCatch(gitcreds::gitcreds_get()$password, error = function(e) NA_character_)
  rifreport_dl <- update_rifreport_download(token)
  results$rifreport <- if (identical(rifreport_dl$status, "ok")) {
    update_rifreport_finalize(rifreport_dl$raw_mutations, refs, seqs, coordinates)
  } else {
    rifreport_dl
  }

  for (r in results) {
    message(sprintf("[%s] %s: %s", r$dataset, r$status, r$message))
  }
}
