library(tidyverse)
library(tidyr)
library(dplyr)
library(ALJEbinf)
library(Biostrings)
library(rentrez)
library(xml2)

source("code/helpers.R")

# Global settings and parameters:
globsets <- list(
  drug_name = "Rifampicin",
  drug_code = "RIF",
  gene_name = "rpoB",
  variant_type = "single resistance variant", # filter variant type for CARD data
  cryptic_min_coef = 0, # minimum log2(MIC) for CRyPTIC data
  cryptic_max_bh_p = 0.05, # maximum BH-corrected p-value for CRyPTIC data
  provean_score = -2.5, # maximum provean score for DRAGdb data
  Yang_min_log2_fc = 0 # minimum log2(foldchange) in frequencies to include a mutation as "selected for"
)

# Loads the shared rpoB reference sequences/species table fresh from disk.
# Called explicitly (not at module level) so callers can run it *after*
# process_refseqs.R's reference-discovery step has had a chance to append
# newly-found species - otherwise the reference-gated map_*() functions
# below would use a stale species list.
load_reference_context <- function() {
  seqs <- readDNAStringSet("raw_data/rpoB_references.fasta")
  refs <- read_csv("raw_data/rpoB_references.csv", show_col_types = FALSE) |>
    filter(RefStrain == TRUE)

  list(
    seqs = seqs,
    refs = refs,
    rpoB_Ecoli = seqs[["rpoB_Escherichia_coli_MG1655"]],
    rpoB_TB = seqs[["rpoB_Mycobacterium_tuberculosis_H37Rv"]]
  )
}

# ----- Per-source processing functions -----
#
# Each function returns
# list(dataset=, status=c("ok","error"), message=, path=, rows=, data=),
# the orchestrator can continue past one
# source's failure and still combine whatever succeeded.
#
# None of these functions writes `data` to `path` themselves - update_all.R
# persists it centrally, and only if the *entire* run passes validation
# (see run_validation_checks() in update_report.R and
# run_pipeline_validation()'s `is_clean` in pipeline_run.R) - a single
# dataset's error blocks every dataset's commit this run, not just its own.
#
# The 4 sources that filter mutations to species with a known rpoB reference
# (CARD, DRAGdb, PointFinder, Reference Gene Catalog) are split into a
# parse_<name>() step (raw parsing through to a `species` column, no
# reference filtering) and a map_<name>(parsed, ctx) step (the existing
# reference-gated mapping). This exposes each source's *unfiltered* species
# list for process_refseqs.R's reference-discovery step without duplicating
# any parsing logic - see update_all.R.

## 1. CRyPTIC data (only for MTB, fixed species, no discovery needed) -----
process_cryptic <- function(ctx) {
  result <- list(
    dataset = "cryptic", status = "error", message = NA_character_,
    path = "processed_data/cryptic.rds", rows = NA_integer_, data = NULL
  )
  tryCatch(
    {
      cryptic_data <- read_csv("raw_data/cryptic_data.csv", show_col_types = FALSE)

      cryptic_processed_data <- cryptic_data |>
        janitor::clean_names() |>
        dplyr::filter(
          drug == globsets$drug_code,
          gene == globsets$gene_name,
          !stringr::str_detect(muts, "indel"),
          coef > globsets$cryptic_min_coef,
          bh_corrected_p <= globsets$cryptic_max_bh_p
        ) |>
        tidyr::separate(
          muts,
          into = c("gene_symbol", "mutations"),
          sep = "_",
          extra = "merge",
          fill = "right"
        ) |>
        tidyr::extract(
          mutations,
          into = c("AA_original", "AA_pos", "AA_mutation"),
          regex = "^(.)(-?[0-9]+)(.)$",
          remove = FALSE
        ) |>
        dplyr::filter(
          !is.na(mutations),
          mutations != "",
          !is.na(AA_original),
          !is.na(AA_pos),
          !is.na(AA_mutation),
          !is.na(gene_symbol),
          gene_symbol != ""
        ) |>
        dplyr::transmute(
          species = "Mycobacterium tuberculosis",
          AA_original,
          AA_pos = as.integer(AA_pos),
          AA_mutation,
          source = "CRyPTIC",
          DOI = "https://doi.org/10.1038/s41467-023-44325-5",
          PMID = "38216576",
          dplyr::across(dplyr::everything())
        ) |>
        dplyr::distinct()

      cryptic_mapped_ecoli <- map_to_ecoli(
        data = cryptic_processed_data,
        rpoB_target = ctx$rpoB_TB,
        rpoB_ecoli = ctx$rpoB_Ecoli
      )

      if (nrow(cryptic_mapped_ecoli) == 0) stop("CRyPTIC produced zero mapped rows")

      result$status <- "ok"
      result$message <- "Processed CRyPTIC data"
      result$rows <- nrow(cryptic_mapped_ecoli)
      result$data <- cryptic_mapped_ecoli
    },
    error = function(e) {
      result$message <<- conditionMessage(e)
    }
  )
  result
}

## 2. CARD data -----
parse_card <- function() {
  card_data <- read_csv("raw_data/card_data.csv", show_col_types = FALSE)

  card_processed_data <- card_data |>
    janitor::clean_names() |>
    tidyr::separate(
      card_short_name,
      into = c("species_abbrev", "gene", "drug"),
      sep = "_",
      fill = "right",
      extra = "merge",
      remove = TRUE
    ) |>
    dplyr::filter(
      drug == globsets$drug_code,
      parameter_type == globsets$variant_type,
      gene == globsets$gene_name
    ) |>
    dplyr::mutate(
      species = stringr::str_extract(name, "^\\S+\\s+\\S+")
    ) |>
    tidyr::extract(
      mutations,
      into = c("AA_original", "AA_pos", "AA_mutation"),
      regex = "^([A-Za-z]+)[^0-9]*([0-9]+)[^A-Za-z0-9]*([A-Za-z]+)$",
      remove = FALSE
    ) |>
    dplyr::mutate(
      AA_pos = as.integer(AA_pos),
      PMID = citation
    ) |>
    dplyr::filter(
      !is.na(species),
      !is.na(AA_original),
      !is.na(AA_pos),
      !is.na(AA_mutation),
      AA_mutation != "Var"
    ) |>
    dplyr::select(
      species,
      AA_original,
      AA_mutation,
      AA_pos,
      PMID
    ) |>
    dplyr::distinct()

  card_processed_data |>
    dplyr::group_by(
      species,
      AA_pos,
      AA_original,
      AA_mutation
    ) |>
    dplyr::summarise(
      PMID = paste(
        unique(na.omit(str_trim(PMID))),
        collapse = "; "
      ),
      .groups = "drop"
    )
}

map_card <- function(parsed, ctx) {
  result <- list(
    dataset = "card", status = "error", message = NA_character_,
    path = "processed_data/card.rds", rows = NA_integer_, data = NULL
  )
  tryCatch(
    {
      card_mapped_ecoli <- parsed |>
        dplyr::filter(species %in% ctx$refs$Species) |>
        dplyr::mutate(
          species_sq = stringr::str_replace_all(species, " ", "_")
        ) |>
        dplyr::group_by(species_sq) |>
        dplyr::group_modify(\(df, key) {
          sp <- key$species_sq[[1]]
          rpoB_target <- get_rpoB_for_species(sp, ctx$seqs)

          if (is.null(rpoB_target)) {
            return(dplyr::tibble())
          }

          map_to_ecoli(df, rpoB_target, ctx$rpoB_Ecoli)
        }) |>
        dplyr::ungroup() |>
        dplyr::select(-species_sq) |>
        dplyr::mutate(source = "CARD")

      if (nrow(card_mapped_ecoli) == 0) stop("CARD produced zero mapped rows")

      result$status <- "ok"
      result$message <- "Processed CARD data"
      result$rows <- nrow(card_mapped_ecoli)
      result$data <- card_mapped_ecoli
    },
    error = function(e) {
      result$message <<- conditionMessage(e)
    }
  )
  result
}

## 3. DRAGdb -----
parse_dragdb <- function() {
  dragdb_data <- read_csv("raw_data/dragdb_data.csv", show_col_types = FALSE)

  aa_map <- c(
    Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
    Glu = "E", Gln = "Q", Gly = "G", His = "H", Ile = "I",
    Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
    Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V"
  )

  dragdb_processed_data <- dragdb_data |>
    filter(
      str_detect(drug, globsets$drug_name),
      provean <= globsets$provean_score
    ) |>
    separate(
      codon_mutation,
      into = c("codon_original", "codon_mutation"),
      sep = "\\s*-\\s*",
      remove = FALSE
    ) |>
    mutate(
      AA_original = aa_map[codon_original],
      AA_mutation = aa_map[codon_mutation],
      AA_pos = as.integer(AA_pos),
      species = stringr::str_trim(species),
      PMID = as.character(PMID)
    ) |>
    select(species, AA_original, AA_pos, AA_mutation, provean, PMID) |>
    distinct()

  dragdb_processed_data |>
    group_by(
      species,
      AA_pos,
      AA_original,
      AA_mutation,
      provean
    ) |>
    summarise(
      PMID = paste(
        unique(na.omit(str_trim(PMID))),
        collapse = "; "
      ),
      .groups = "drop"
    )
}

map_dragdb <- function(parsed, ctx) {
  result <- list(
    dataset = "dragdb", status = "error", message = NA_character_,
    path = "processed_data/dragdb.rds", rows = NA_integer_, data = NULL
  )
  tryCatch(
    {
      dragdb_mapped_ecoli <- parsed |>
        filter(species %in% ctx$refs$Species) |>
        mutate(species_sq = str_replace(species, " ", "_")) |>
        group_by(species_sq) |>
        group_modify(\(df, key) {
          sp <- key$species_sq[1]
          rpoB_target <- get_rpoB_for_species(sp, ctx$seqs)
          map_to_ecoli(df, rpoB_target, ctx$rpoB_Ecoli)
        }) |>
        ungroup() |>
        select(!species_sq) |>
        mutate(source = "DRAGdb")

      if (nrow(dragdb_mapped_ecoli) == 0) stop("DRAGdb produced zero mapped rows")

      result$status <- "ok"
      result$message <- "Processed DRAGdb data"
      result$rows <- nrow(dragdb_mapped_ecoli)
      result$data <- dragdb_mapped_ecoli
    },
    error = function(e) {
      result$message <<- conditionMessage(e)
    }
  )
  result
}

## 4. Yang's data (fixed species, no discovery needed) -----
process_yang <- function(ctx) {
  result <- list(
    dataset = "yang", status = "error", message = NA_character_,
    path = "processed_data/yang.rds", rows = NA_integer_, data = NULL
  )
  tryCatch(
    {
      yang_data <- read_csv("raw_data/Yang2023_Fig3a_data_RIF_log2FC.csv", show_col_types = FALSE)
      rpoB_aa <- as.character(translate(ctx$rpoB_Ecoli))

      yang_mut_list <- yang_data |>
        pivot_longer(`_`:Q, names_to = "AA_mutation", values_to = "log2_fc") |>
        filter(log2_fc >= globsets$Yang_min_log2_fc) |>
        mutate(
          AA_pos_Ecoli = as.integer(paste0("5", pos)),
          AA_pos = AA_pos_Ecoli,
          species = "Escherichia coli",
          source = "Yang",
          AA_original = NA,
          PMID = "37648864"
        ) |>
        select(species, AA_original, AA_pos, AA_pos_Ecoli, AA_mutation, log2_fc, source, PMID)

      for (i in 1:nrow(yang_mut_list)) {
        pos <- as.integer(yang_mut_list$AA_pos_Ecoli[i])
        yang_mut_list$AA_original[i] <- substr(rpoB_aa, pos, pos)
      }

      if (nrow(yang_mut_list) == 0) stop("Yang produced zero rows")

      result$status <- "ok"
      result$message <- "Processed Yang data"
      result$rows <- nrow(yang_mut_list)
      result$data <- yang_mut_list
    },
    error = function(e) {
      result$message <<- conditionMessage(e)
    }
  )
  result
}

## 5. RifReport (already E. coli-native coordinates, no ctx needed) -----
process_rifreport <- function() {
  result <- list(
    dataset = "rifreport", status = "error", message = NA_character_,
    path = "processed_data/rifreport.rds", rows = NA_integer_, data = NULL
  )
  tryCatch(
    {
      rifreport_data <- read_csv("raw_data/rifreport_data.csv", show_col_types = FALSE)

      # The upstream ResistanceSurveys repo spells these two Ref_codes
      # differently between RIF_mutations_data.csv and
      # RIF_mutations_references.csv (accent dropped for Honore2001,
      # surname truncated for AlvarezNarvaez2019) - normalize to the
      # references file's spelling so build_references()'s join doesn't
      # silently drop the citation for these rows.
      rifreport_ref_code_fixes <- c(
        "Honore2001"         = "Honoré2001",
        "AlvarezNarvaez2019" = "Alvarez2019"
      )
      rifreport_data <- rifreport_data |>
        dplyr::mutate(
          Ref_code = dplyr::coalesce(
            unname(rifreport_ref_code_fixes[Ref_code]),
            Ref_code
          )
        )

      mutation_raw <- rifreport_data |>
        filter(!is.na(AA_pos_Ecoli), !is.na(AA_original), !is.na(AA_mutation), !is.na(Ref_code)) |>
        dplyr::mutate(
          AA_original = toupper(AA_original),
          AA_mutation = toupper(AA_mutation),
          AA_pos = as.integer(AA_pos),
          AA_pos_Ecoli = as.integer(AA_pos_Ecoli)
        )

      counts_species <- mutation_raw |>
        select(Species, AA_pos_Ecoli, AA_mutation) |>
        distinct() |>
        group_by(AA_pos_Ecoli, AA_mutation) |>
        summarise(n_species = n(), .groups = "drop")

      counts_studies <- mutation_raw |>
        select(Ref_code, AA_pos_Ecoli, AA_mutation) |>
        distinct() |>
        group_by(AA_pos_Ecoli, AA_mutation) |>
        summarise(n_studies = n(), .groups = "drop")

      mutation_counts <- full_join(
        counts_species,
        counts_studies,
        by = c("AA_pos_Ecoli", "AA_mutation")
      ) |>
        mutate(
          source = "RifReport",
          n_species = coalesce(n_species, 0L),
          n_studies = coalesce(n_studies, 0L)
        )

      mutation_species_rows <- mutation_raw |>
        distinct(Species, AA_pos, AA_original, AA_pos_Ecoli, AA_mutation, Ref_code) |>
        dplyr::group_by(Species, AA_pos, AA_original, AA_pos_Ecoli, AA_mutation) |>
        dplyr::summarise(
          Ref_codes = paste(sort(unique(Ref_code)), collapse = ", "),
          .groups   = "drop"
        )

      rifreport_processed <- mutation_species_rows |>
        left_join(mutation_counts, by = c("AA_pos_Ecoli", "AA_mutation")) |>
        dplyr::rename(
          species = Species
        )

      if (nrow(rifreport_processed) == 0) stop("RifReport produced zero rows")

      result$status <- "ok"
      result$message <- "Processed RifReport data"
      result$rows <- nrow(rifreport_processed)
      result$data <- rifreport_processed
    },
    error = function(e) {
      result$message <<- conditionMessage(e)
    }
  )
  result
}

## 6. PointFinder -----
parse_pointfinder <- function() {
  pointfinder_data <- read_csv("raw_data/pointfinder_data.csv", show_col_types = FALSE)

  pointfinder_data |>
    select(species, Codon_pos, Ref_codon, Res_codon, PMID) |>
    dplyr::rename(
      AA_pos = Codon_pos,
      AA_original = Ref_codon,
      AA_mutation = Res_codon
    ) |>
    mutate(
      source = "PointFinder",
      PMID = stringr::str_replace_all(PMID, ",", ";")
    )
}

map_pointfinder <- function(parsed, ctx) {
  result <- list(
    dataset = "pointfinder", status = "error", message = NA_character_,
    path = "processed_data/pointfinder.rds", rows = NA_integer_, data = NULL
  )
  tryCatch(
    {
      pointfinder_mapped_ecoli <- parsed |>
        dplyr::filter(species %in% ctx$refs$Species) |>
        dplyr::mutate(
          species_sq = stringr::str_replace_all(species, " ", "_")
        ) |>
        dplyr::group_by(species_sq) |>
        dplyr::group_modify(\(df, key) {
          sp <- key$species_sq[[1]]
          rpoB_target <- get_rpoB_for_species(sp, ctx$seqs)

          if (is.null(rpoB_target)) {
            return(dplyr::tibble())
          }

          map_to_ecoli(df, rpoB_target, ctx$rpoB_Ecoli)
        }) |>
        dplyr::ungroup() |>
        dplyr::filter(!stringr::str_detect(AA_mutation, "\\*")) |>
        dplyr::select(-species_sq)

      if (nrow(pointfinder_mapped_ecoli) == 0) stop("PointFinder produced zero mapped rows")

      result$status <- "ok"
      result$message <- "Processed PointFinder data"
      result$rows <- nrow(pointfinder_mapped_ecoli)
      result$data <- pointfinder_mapped_ecoli
    },
    error = function(e) {
      result$message <<- conditionMessage(e)
    }
  )
  result
}

## 7. Reference Gene Catalog -----
get_ncbi_organism <- function(accession) {
  s <- entrez_summary(db = "nuccore", id = accession)
  s$organism
}

parse_refgenes <- function() {
  refgenes_data <- read_csv("raw_data/refgenes_catalog.csv", col_types = cols(.default = "c"), show_col_types = FALSE)

  refgenes_data |>
    filter(
      str_detect(Subclass, regex("RIFAMPIN", ignore_case = TRUE)),
      `Gene family` == "rpoB"
    ) |>
    extract(
      `#Allele`,
      into = c("Gene", "AA_original", "AA_pos", "AA_mutation"),
      regex = "^(\\w+)_([A-Z])(\\d+)([A-Z])$",
      remove = FALSE
    ) |>
    tidyr::drop_na(AA_original, AA_pos, AA_mutation) |>
    mutate(
      organism = map_chr(`GenBank nucleotide`, get_ncbi_organism)
    ) |>
    mutate(
      species = organism |>
        stringr::word(1, 2),
      PMID = as.character(`Pubmed reference`),
      PMID = stringr::str_replace_all(PMID, ",", ";")
    ) |>
    select(species, AA_pos, AA_original, AA_mutation, PMID)
}

map_refgenes <- function(parsed, ctx) {
  result <- list(
    dataset = "refgenes", status = "error", message = NA_character_,
    path = "processed_data/refgenes.rds", rows = NA_integer_, data = NULL
  )
  tryCatch(
    {
      refgenes_mapped_ecoli <- parsed |>
        mutate(AA_pos = as.numeric(AA_pos)) |>
        dplyr::filter(species %in% ctx$refs$Species) |>
        dplyr::mutate(
          species_sq = stringr::str_replace_all(species, " ", "_")
        ) |>
        dplyr::group_by(species_sq) |>
        dplyr::group_modify(\(df, key) {
          sp <- key$species_sq[[1]]
          rpoB_target <- get_rpoB_for_species(sp, ctx$seqs)

          if (is.null(rpoB_target)) {
            return(dplyr::tibble())
          }

          map_to_ecoli(df, rpoB_target, ctx$rpoB_Ecoli)
        }) |>
        dplyr::ungroup() |>
        dplyr::select(-species_sq) |>
        mutate(source = "Reference Gene Catalog")

      if (nrow(refgenes_mapped_ecoli) == 0) stop("Reference Gene Catalog produced zero mapped rows")

      result$status <- "ok"
      result$message <- "Processed Reference Gene Catalog data"
      result$rows <- nrow(refgenes_mapped_ecoli)
      result$data <- refgenes_mapped_ecoli
    },
    error = function(e) {
      result$message <<- conditionMessage(e)
    }
  )
  result
}

# ----- Combine -----
combine_all_data <- function(frames) {
  frames <- Filter(Negate(is.null), frames)

  all_data <- frames |>
    map(~ {
      .x |> mutate(PMID = if ("PMID" %in% names(.x)) as.character(PMID) else NA_character_)
    }) |>
    bind_rows()

  saveRDS(all_data, file = "processed_data/allData.rds")
  all_data
}

# ----- Process references -----

extract_medline_field <- function(tag, lines) {
  hit <- grep(paste0("^", tag), lines)

  if (length(hit) == 0) {
    return(NA_character_)
  }

  value <- sub(paste0("^", tag, "\\s*"), "", lines[hit[1]])

  i <- hit[1] + 1
  while (i <= length(lines) && grepl("^\\s{6}", lines[i])) {
    value <- paste(value, trimws(lines[i]))
    i <- i + 1
  }

  stringr::str_squish(value)
}


extract_doi_medline <- function(lines) {
  aid_lines <- grep("^AID - ", lines, value = TRUE)

  if (length(aid_lines) == 0) {
    return(NA_character_)
  }

  doi_line <- aid_lines[grepl("\\[doi\\]$", aid_lines)]

  if (length(doi_line) == 0) {
    return(NA_character_)
  }

  doi <- doi_line[1] |>
    sub("^AID -\\s*", "", x = _) |>
    sub("\\s*\\[doi\\]$", "", x = _)

  doi <- paste0("https://doi.org/", doi)
}


extract_author_key_medline <- function(lines) {
  cn <- grep("^CN  - ", lines, value = TRUE)
  if (length(cn) > 0) {
    return(sub("^CN  -\\s*", "", cn[1]))
  }

  ir <- grep("^IR  - ", lines, value = TRUE)
  if (length(ir) > 0) {
    return(sub("^IR  -\\s*", "", ir[1]))
  }

  au <- grep("^AU  - ", lines, value = TRUE)
  if (length(au) > 0) {
    return(sub("^AU  -\\s*", "", au[1]))
  }

  NA_character_
}


extract_author_list_medline <- function(lines, max_authors = 6) {
  cn <- grep("^CN  - ", lines, value = TRUE)
  if (length(cn) > 0) {
    return(sub("^CN  -\\s*", "", cn[1]))
  }

  ir_lines <- sub("^IR  -\\s*", "", grep("^IR  - ", lines, value = TRUE))
  au_lines <- sub("^AU  -\\s*", "", grep("^AU  - ", lines, value = TRUE))

  author_vec <- if (length(ir_lines) > 0) ir_lines else au_lines

  if (length(author_vec) == 0) {
    return(NA_character_)
  }

  if (length(author_vec) <= max_authors) {
    paste(author_vec, collapse = ", ")
  } else {
    paste0(paste(author_vec[1:max_authors], collapse = ", "), ", et al")
  }
}


make_ref_code <- function(author_key, year) {
  if (is.na(author_key) || author_key == "" || is.na(year) || year == "") {
    return(NA_character_)
  }

  first_part <- stringr::str_extract(author_key, "^[A-Za-z']+")
  first_part <- stringr::str_remove_all(first_part, "[^A-Za-z]")

  paste0(first_part, year)
}


build_full_ref_medline <- function(lines, pmid_val) {
  title <- extract_medline_field("TI  - ", lines)
  so <- extract_medline_field("SO  - ", lines)
  pmcid <- extract_medline_field("PMC - ", lines)

  author_text <- extract_author_list_medline(lines)

  full_ref <- paste0(
    if (!is.na(author_text) && author_text != "") paste0(author_text, ". ") else "",
    if (!is.na(title) && title != "") paste0(title, ". ") else "",
    if (!is.na(so) && so != "") so else "",
    if (!is.na(pmid_val) && pmid_val != "") paste0(" PMID: ", pmid_val) else "",
    if (!is.na(pmcid) && pmcid != "") paste0("; PMCID: ", pmcid, ".") else "."
  )

  full_ref <- stringr::str_replace_all(full_ref, "\\.\\.", ".")
  stringr::str_squish(full_ref)
}


parse_one_pubmed_ref <- function(pmid) {
  txt <- tryCatch(
    rentrez::entrez_fetch(
      db = "pubmed",
      id = as.character(pmid),
      rettype = "medline",
      retmode = "text"
    ),
    error = function(e) NA_character_
  )

  if (is.na(txt) || txt == "") {
    message("Failed PMID: ", pmid)
    return(tibble::tibble(
      Ref_code = NA_character_,
      Year = NA_character_,
      Full_ref = NA_character_,
      DOI = NA_character_,
      PMID = as.character(pmid)
    ))
  }

  lines <- unlist(strsplit(txt, "\n"))

  dp <- extract_medline_field("DP  - ", lines)
  year <- stringr::str_extract(dp, "\\b\\d{4}\\b")
  doi <- extract_doi_medline(lines)
  author_key <- extract_author_key_medline(lines)
  ref_code <- make_ref_code(author_key, year)

  if (is.na(doi) || !nzchar(trimws(doi))) {
    doi <- paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/")
  }

  full_ref <- build_full_ref_medline(lines, pmid)

  tibble::tibble(
    Ref_code = ref_code,
    Year = year,
    Full_ref = full_ref,
    DOI = doi,
    PMID = pmid
  )
}

get_pubmed_ref_table <- function(pmids) {
  pmids |>
    as.character() |>
    stringr::str_split(";") |>
    unlist() |>
    stringr::str_trim() |>
    discard(~ .x == "" | is.na(.x)) |>
    unique() |>
    purrr::map_dfr(parse_one_pubmed_ref)
}

# Ref_code is <first author surname><year> (make_ref_code()), generated one
# PMID at a time inside parse_one_pubmed_ref(), so it carries no uniqueness
# guarantee - the generator cannot see a clash even in principle. Two failure
# modes occur in real data and need OPPOSITE fixes:
#
#   A. Two DIFFERENT papers share a code (Heep2000, Taniguchi1996,
#      Lisitsyn1984, Bland1976). Left alone, the app's reference table -
#      which filters `Ref_code %in% ref_codes` - shows every paper sharing the
#      code, so citing one mis-attributes the other alongside it. Must SPLIT.
#   B. The SAME paper appears twice with differently-worded Full_ref
#      (Ohno1996, Lavania2016). Must MERGE - suffixing would invent a second
#      paper that does not exist.
#
# Papers are identified by PMID. Rows with no PMID fold into the group's paper
# when that is unambiguous (exactly one PMID present), which is how a curated
# entry lacking a PMID picks one up from its resolved twin.
#
# Codes already in raw_data/rifreport_references.csv are join keys for
# raw_data/rifreport_data.csv's own Ref_code column, so they are NEVER
# renamed; only auto-generated entries take a suffix. Ordering is by PMID, not
# row order, so the same input always yields the same codes - otherwise
# published citations would churn between runs.
disambiguate_ref_codes <- function(refs) {
  ref_cols <- c("Ref_code", "Year", "Full_ref", "DOI", "PMID")

  suffix_at <- function(i) if (i <= 26) letters[i] else as.character(i)

  has_value <- function(x) !is.na(x) & x != ""

  # Curated first, so each field takes the curated value and falls back to the
  # generated entry only to fill a blank (e.g. a missing PMID).
  merge_paper <- function(rows) {
    rows <- rows[order(!rows$.curated), , drop = FALSE]
    out <- rows[1, , drop = FALSE]
    for (cl in intersect(ref_cols, names(rows))) {
      vals <- rows[[cl]][has_value(rows[[cl]])]
      if (length(vals) > 0) out[[cl]] <- vals[1]
    }
    out$.curated <- any(rows$.curated)
    out
  }

  resolve_group <- function(g) {
    pmids <- sort(unique(g$PMID[has_value(g$PMID)]))

    if (length(pmids) <= 1) {
      return(merge_paper(g))
    }

    papers <- lapply(pmids, function(p) merge_paper(g[has_value(g$PMID) & g$PMID == p, , drop = FALSE]))
    orphans <- g[!has_value(g$PMID), , drop = FALSE]
    if (nrow(orphans) > 0) {
      papers <- c(papers, lapply(seq_len(nrow(orphans)), function(i) orphans[i, , drop = FALSE]))
    }
    papers <- dplyr::bind_rows(papers)

    code <- g$Ref_code[1]
    curated_i <- which(papers$.curated)

    if (length(curated_i) > 0) {
      others <- setdiff(seq_len(nrow(papers)), curated_i[1])
      papers$Ref_code[curated_i[1]] <- code
      papers$Ref_code[others] <- paste0(code, vapply(seq_along(others) + 1L, suffix_at, character(1)))
    } else {
      papers$Ref_code <- paste0(code, vapply(seq_len(nrow(papers)), suffix_at, character(1)))
    }
    papers
  }

  # Rows with no Ref_code cannot be grouped or suffixed - pass them through.
  unnamed <- refs[!has_value(refs$Ref_code), , drop = FALSE]
  named <- refs[has_value(refs$Ref_code), , drop = FALSE]

  if (nrow(named) == 0) {
    return(refs)
  }

  resolved <- named |>
    dplyr::group_by(Ref_code) |>
    dplyr::group_split() |>
    lapply(resolve_group) |>
    dplyr::bind_rows()

  dplyr::bind_rows(resolved, unnamed)
}

# Builds processed_data/all_references.csv. Diffs against the existing file
# so only PMIDs not already resolved there get a fresh MEDLINE fetch -
# on a typical run where only 1-2 sources changed, this avoids re-fetching
# every reference across all 7 datasets every time. Falls back to a full
# rebuild if the existing file is missing/malformed.
build_references <- function(combined_data,
                             rifreport_refs_path = "raw_data/rifreport_references.csv",
                             output_path = "processed_data/all_references.csv") {
  rifreport_refs <- read_csv(rifreport_refs_path, show_col_types = FALSE)
  rifreport_refs_clean <- rifreport_refs |>
    dplyr::distinct() |>
    dplyr::mutate(
      Ref_code = as.character(Ref_code),
      Year     = as.character(Year),
      Full_ref = as.character(Full_ref),
      DOI      = as.character(DOI),
      PMID     = as.character(PMID)
    ) |>
    dplyr::select(dplyr::any_of(c("Ref_code", "Year", "Full_ref", "DOI", "PMID")))

  all_pmids <- combined_data |>
    dplyr::distinct(PMID) |>
    dplyr::pull(PMID) |>
    as.character() |>
    stringr::str_split(";") |>
    unlist() |>
    stringr::str_trim() |>
    purrr::discard(~ .x == "" | is.na(.x)) |>
    unique()

  existing_resolved <- tryCatch(
    {
      if (file.exists(output_path)) {
        readr::read_csv(output_path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
      } else {
        NULL
      }
    },
    error = function(e) NULL
  )

  already_resolved_pmids <- if (!is.null(existing_resolved) && "PMID" %in% names(existing_resolved)) {
    existing_resolved |>
      dplyr::filter(!is.na(Full_ref) & Full_ref != "") |>
      dplyr::pull(PMID) |>
      unique()
  } else {
    character(0)
  }

  new_pmids <- setdiff(all_pmids, already_resolved_pmids)
  reused_pmids <- intersect(all_pmids, already_resolved_pmids)

  message(
    "References: resolving ", length(new_pmids), " new PMID(s), reusing ",
    length(reused_pmids), " already-resolved PMID(s)."
  )

  newly_resolved_clean <- if (length(new_pmids) > 0) {
    get_pubmed_ref_table(new_pmids) |>
      dplyr::mutate(
        Ref_code = as.character(Ref_code),
        Year     = as.character(Year),
        Full_ref = as.character(Full_ref),
        DOI      = as.character(DOI),
        PMID     = as.character(PMID)
      )
  } else {
    tibble::tibble()
  }

  reused_resolved <- if (!is.null(existing_resolved) && length(reused_pmids) > 0) {
    existing_resolved |> dplyr::filter(PMID %in% reused_pmids)
  } else {
    tibble::tibble()
  }

  ref_cols <- c("Ref_code", "Year", "Full_ref", "DOI", "PMID")

  ref_table_clean <- dplyr::bind_rows(reused_resolved, newly_resolved_clean) |>
    dplyr::select(dplyr::any_of(ref_cols))

  combined_refs <- dplyr::bind_rows(
    rifreport_refs_clean |> dplyr::mutate(.curated = TRUE),
    ref_table_clean |> dplyr::mutate(.curated = FALSE)
  ) |>
    dplyr::arrange(dplyr::desc(.curated)) |>
    dplyr::distinct(dplyr::across(dplyr::any_of(ref_cols)), .keep_all = TRUE) |>
    disambiguate_ref_codes() |>
    dplyr::select(dplyr::any_of(ref_cols)) |>
    dplyr::arrange(Ref_code)

  write_csv(combined_refs, output_path)
  combined_refs
}

# ---- Filtering report ----
#
# Reporting only. Nothing here is part of the pipeline: neither
# run_pipeline_validation() nor the standalone driver below calls it, and it
# writes nothing. It lives in this file so the stage expressions sit beside the
# dplyr::filter() calls they mirror and read the same `globsets` thresholds -
# keep the two in sync when a filter changes.
#
# report_filtering() re-reads each raw file and applies that source's inclusion
# criteria one at a time, recording how many records survive each. It covers
# only the pre-mapping filters; mutation-code parsing and reference-gated
# mapping are not replicated and instead surface as the residual on each
# source's final line. Called by code/reports.R, and usable on its own:
#
#   .pipeline_orchestrating <- TRUE
#   source("code/process_data.R")
#   cat(format_filtering_funnels(report_filtering()))

# Applies `...` cumulatively as dplyr::filter() stages, recording how many
# records survive each one. Name every stage with the criterion it applies.
filter_funnel <- function(df, ...) {
  stages <- rlang::quos(...)
  current <- df
  out <- vector("list", length(stages))
  for (i in seq_along(stages)) {
    before <- nrow(current)
    current <- dplyr::filter(current, !!stages[[i]])
    out[[i]] <- tibble::tibble(
      stage = names(stages)[i],
      kept = nrow(current),
      dropped = before - nrow(current)
    )
  }
  dplyr::bind_rows(out)
}

# One source's rows: the raw count, one row per inclusion criterion, and the
# entry count RIFRMap serves. The gap between the last criterion and `n_final`
# is reported as that final row's loss, so each source reconciles back to its
# raw file.
funnel_rows <- function(source, path, n_raw, stages, n_final,
                        raw_label = "Records in raw file") {
  n_last <- if (nrow(stages) > 0) stages$kept[nrow(stages)] else n_raw
  dplyr::bind_rows(
    tibble::tibble(
      source = source, path = path, stage = raw_label,
      kept = n_raw, dropped = NA_integer_, kind = "raw"
    ),
    if (nrow(stages) > 0) {
      tibble::tibble(
        source = source, path = path, stage = stages$stage,
        kept = stages$kept, dropped = stages$dropped, kind = "filter"
      )
    },
    tibble::tibble(
      source = source, path = path, stage = "Entries served by RIFRMap",
      kept = n_final, dropped = n_last - n_final, kind = "final"
    )
  )
}

# Returns one row per funnel stage across all seven sources, with columns
# source / path / stage / kept / dropped / kind ("raw", "filter", "final").
# `final_counts` is a named count of the entries the app serves per source
# label; it defaults to processed_data/allData.rds, and any source missing from
# it reports NA rather than failing.
report_filtering <- function(final_counts = NULL) {
  if (is.null(final_counts)) {
    final_counts <- if (file.exists("processed_data/allData.rds")) {
      table(readRDS("processed_data/allData.rds")$source)
    } else {
      integer(0)
    }
  }
  n_final <- function(label) {
    if (label %in% names(final_counts)) as.integer(final_counts[[label]]) else NA_integer_
  }

  ## CARD ----------------------------------------------------------------
  card_raw <- read_csv("raw_data/card_data.csv", show_col_types = FALSE) |>
    janitor::clean_names() |>
    tidyr::separate(
      card_short_name,
      into = c("species_abbrev", "gene", "drug"),
      sep = "_", fill = "right", extra = "merge", remove = FALSE
    )

  card_stages <- filter_funnel(
    card_raw,
    "Drug is RIF"               = drug == globsets$drug_code,
    "Single resistance variant" = parameter_type == globsets$variant_type,
    "Gene is rpoB"              = gene == globsets$gene_name
  )

  ## CRyPTIC -------------------------------------------------------------
  cryptic_raw <- read_csv("raw_data/cryptic_data.csv", show_col_types = FALSE) |>
    janitor::clean_names()

  cryptic_stages <- filter_funnel(
    cryptic_raw,
    "Drug is RIF"      = drug == globsets$drug_code,
    "Gene is rpoB"     = gene == globsets$gene_name,
    "Excluding indels" = !stringr::str_detect(muts, "indel"),
    !!paste0("Effect size > ", globsets$cryptic_min_coef, " (log2 MIC)") :=
      coef > globsets$cryptic_min_coef,
    !!paste0("BH-corrected p <= ", globsets$cryptic_max_bh_p) :=
      bh_corrected_p <= globsets$cryptic_max_bh_p
  )

  ## DRAGdb --------------------------------------------------------------
  # The PROVEAN criterion is split across two stages because the single
  # dplyr::filter() in parse_dragdb() drops NA alongside scores above the
  # threshold. Most of DRAGdb's loss here is entries carrying no PROVEAN score
  # at all, not entries predicted to be neutral - a distinction worth seeing.
  dragdb_raw <- read_csv("raw_data/dragdb_data.csv", show_col_types = FALSE)

  dragdb_stages <- filter_funnel(
    dragdb_raw,
    "Drug is RIF"            = stringr::str_detect(drug, globsets$drug_name),
    "PROVEAN score reported" = !is.na(provean),
    !!paste0("PROVEAN score <= ", globsets$provean_score, " (deleterious)") :=
      provean <= globsets$provean_score
  )

  ## PointFinder ---------------------------------------------------------
  # Restricted to rifampicin/rpoB at download time by update_pointfinder(), so
  # the raw file carries no further inclusion filter.
  pointfinder_raw <- read_csv("raw_data/pointfinder_data.csv", show_col_types = FALSE)
  pointfinder_stages <- filter_funnel(pointfinder_raw)

  ## Reference Gene Catalog ----------------------------------------------
  refgenes_raw <- read_csv(
    "raw_data/refgenes_catalog.csv",
    col_types = cols(.default = "c"), show_col_types = FALSE
  )

  refgenes_stages <- filter_funnel(
    refgenes_raw,
    "Subclass includes rifampin" = stringr::str_detect(
      Subclass, stringr::regex("RIFAMPIN", ignore_case = TRUE)
    ),
    "Gene family is rpoB"   = `Gene family` == globsets$gene_name,
    "Parseable allele code" = stringr::str_detect(`#Allele`, "^\\w+_[A-Z]\\d+[A-Z]$")
  )

  ## RifReport -----------------------------------------------------------
  # The raw file is already rpoB/RIF-only and already expanded to one row per
  # named species by expand_multi_species() (code/update_data.R). Its large drop
  # is de-duplication, not exclusion: the survey records one row per reporting
  # study, which process_rifreport() collapses to one row per species/position/
  # substitution while counting the studies into n_studies. The dedup stage
  # mirrors that grouping key, uppercasing as process_rifreport() does.
  rifreport_raw <- read_csv("raw_data/rifreport_data.csv", show_col_types = FALSE)

  rifreport_stages <- filter_funnel(
    rifreport_raw,
    "Complete, citable E. coli-mapped record" = !is.na(AA_pos_Ecoli) &
      !is.na(AA_original) & !is.na(AA_mutation) & !is.na(Ref_code),
    "Distinct species/position/substitution"  = !duplicated(paste(
      Species, AA_pos, toupper(AA_original), AA_pos_Ecoli, toupper(AA_mutation)
    ))
  )

  ## Yang-RIF ------------------------------------------------------------
  # The raw file is a position x residue grid, so it is pivoted to one row per
  # substitution before counting; the raw figure is that pivoted total, not the
  # grid rows on disk. The stop-codon column is reported separately to reconcile
  # with the 760 standard substitutions cited in the source publication -
  # process_yang() does not exclude it explicitly, but none of its cells clear
  # the log2FC filter, so the final count is the same either way.
  yang_raw <- read_csv(
    "raw_data/Yang2023_Fig3a_data_RIF_log2FC.csv", show_col_types = FALSE
  ) |>
    tidyr::pivot_longer(`_`:Q, names_to = "AA_mutation", values_to = "log2_fc")

  yang_stages <- filter_funnel(
    yang_raw,
    "Standard amino acid substitutions" = AA_mutation != "_",
    !!paste0("log2 fold-change >= ", globsets$Yang_min_log2_fc, " (RIF-selected)") :=
      log2_fc >= globsets$Yang_min_log2_fc
  )

  dplyr::bind_rows(
    funnel_rows("CARD", "raw_data/card_data.csv",
                nrow(card_raw), card_stages, n_final("CARD")),
    funnel_rows("CRyPTIC", "raw_data/cryptic_data.csv",
                nrow(cryptic_raw), cryptic_stages, n_final("CRyPTIC")),
    funnel_rows("DRAGdb", "raw_data/dragdb_data.csv",
                nrow(dragdb_raw), dragdb_stages, n_final("DRAGdb")),
    funnel_rows("PointFinder", "raw_data/pointfinder_data.csv",
                nrow(pointfinder_raw), pointfinder_stages, n_final("PointFinder")),
    funnel_rows("Reference Gene Catalog", "raw_data/refgenes_catalog.csv",
                nrow(refgenes_raw), refgenes_stages,
                n_final("Reference Gene Catalog")),
    funnel_rows("RifReport", "raw_data/rifreport_data.csv",
                nrow(rifreport_raw), rifreport_stages, n_final("RifReport")),
    funnel_rows("Yang-RIF", "raw_data/Yang2023_Fig3a_data_RIF_log2FC.csv",
                nrow(yang_raw), yang_stages, n_final("Yang"),
                raw_label = "Substitutions in mutagenesis grid")
  )
}

# Renders report_filtering()'s tibble as one indented block per source.
format_filtering_funnels <- function(funnels) {
  blocks <- vapply(unique(funnels$source), function(src) {
    rows <- funnels[funnels$source == src, ]
    lines <- ifelse(
      is.na(rows$dropped),
      sprintf("        %-46s : %6d", rows$stage, rows$kept),
      sprintf("        %-46s : %6d  (-%d)", rows$stage, rows$kept, rows$dropped)
    )
    paste(c(paste0("    [", src, "]  ", rows$path[1]), lines), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)

  paste(blocks, collapse = "\n\n")
}

# ---- Standalone execution ----
# When this file is sourced directly (not via update_all.R), process and
# combine every source exactly as before, continuing past individual
# failures and printing a result per dataset. update_all.R sets
# `.pipeline_orchestrating` before sourcing this file, so this block is
# skipped there and only the functions above get defined.
if (!exists(".pipeline_orchestrating", inherits = FALSE)) {
  ctx <- load_reference_context()

  results <- list()
  results$cryptic <- process_cryptic(ctx)
  results$card <- map_card(parse_card(), ctx)
  results$dragdb <- map_dragdb(parse_dragdb(), ctx)
  results$yang <- process_yang(ctx)
  results$pointfinder <- map_pointfinder(parse_pointfinder(), ctx)
  results$refgenes <- map_refgenes(parse_refgenes(), ctx)
  results$rifreport <- process_rifreport()

  all_data <- combine_all_data(lapply(results, function(r) r$data))
  build_references(all_data)

  for (r in results) {
    message(sprintf("[%s] %s: %s", r$dataset, r$status, r$message))
  }
}
