library(tidyverse)
library(tidyr)
library(dplyr)
library(ALJEbinf)
library(Biostrings)
library(rentrez)
library(xml2)

source("code/translate_coordinates.R")

# Global settings and parameters:
globsets <- list(
  drug_name = "Rifampicin",
  drug_code = "RIF",
  gene_name = "rpoB",
  variant_type = "single resistance variant", # filter variant type for CARD data
  #p_value = 0.05, # maxmimum p_value for CRyPTIC data
  provean_score = -2.5, # maximum provean score for DRAGdb data
  Yang_min_log2_fc = 0 # minimum log2(foldchange) in frequencies to include a mutation as "selected for"
 )

# ----- Load database -----
cryptic_data <- read_csv("./data/raw_data/CRyPTIC_mutations_data.csv", show_col_types = FALSE)
card_data <- read_csv("./data/raw_data/CARD_mutations_data.csv", show_col_types = FALSE)
dragdb_data <- read_csv("./data/raw_data/RIF_DRAGdb_mutations_data.csv", show_col_types = FALSE)
yang_data <- read_csv("./data/raw_data/Yang2023_Fig3a_data_RIF_log2FC.csv", show_col_types = FALSE)
curated_data <- read_csv("./data/raw_data/RIF_curated_mutations_data.csv", show_col_types = FALSE) 
curated_refs <- read_csv("./data/raw_data/RIF_curated_mutations_references.csv", show_col_types = FALSE) 

# ----- Load sequencing data -----
seqs <- readDNAStringSet("./data/rpoB_references.fasta")
rpoB_Ecoli <- seqs[["rpoB_Escherichia_coli_MG1655"]]
rpoB_TB <- seqs[["rpoB_Mycobacterium_tuberculosis_H37Rv"]]
refs <- read_csv("data/rpoB_references.csv") |>
  filter(RefStrain == TRUE)
  
# ----- Process database -----
## 1. CRyPTIC data (only for MTB)  -----
cryptic_processed_data <- cryptic_data |>
  janitor::clean_names() |>
  dplyr::filter(
    drug == globsets$drug_code,
    gene == globsets$gene_name,
    !stringr::str_detect(muts, "indel")
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
  rpoB_target = rpoB_TB,
  rpoB_ecoli = rpoB_Ecoli
)

write_csv(
  cryptic_mapped_ecoli,
  "data/processed_data/CRyPTIC_processed.csv"
)

saveRDS(
  cryptic_mapped_ecoli,
  file = "code/RIFRMap/data/cryptic.rds"
)

## 2. CARD data  -----
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

card_group_pmid <- card_processed_data |> 
  group_by(
    species,
    AA_pos,
    AA_original,
    AA_mutation
  ) |>
  summarise(
    PMID = paste(
      unique(na.omit(str_trim(PMID))),
      collapse = "; "
    ),
    .groups = "drop"
  )

card_mapped_ecoli <- card_group_pmid |>
  dplyr::filter(species %in% refs$Species) |>
  dplyr::mutate(
    species_sq = stringr::str_replace_all(species, " ", "_")
  ) |>
  dplyr::group_by(species_sq) |>
  dplyr::group_modify(\(df, key) {
    sp <- key$species_sq[[1]]
    rpoB_target <- get_rpoB_for_species(sp, seqs)
    
    if (is.null(rpoB_target)) {
      return(dplyr::tibble())
    }
    
    map_to_ecoli(df, rpoB_target, rpoB_Ecoli)
  }) |>
  dplyr::ungroup() |>
  dplyr::select(-species_sq) |>
  dplyr::mutate(source = "CARD")

write_csv(
  card_mapped_ecoli,
  "data/processed_data/CARD_processed.csv"
)

saveRDS(
  card_mapped_ecoli,
  file = "code/RIFRMap/data/card.rds"
)

## 3. DRAGdb -----

aa_map <- c(
  Ala="A", Arg="R", Asn="N", Asp="D", Cys="C",
  Glu="E", Gln="Q", Gly="G", His="H", Ile="I",
  Leu="L", Lys="K", Met="M", Phe="F", Pro="P",
  Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V"
)

dragdb_processed_data <- dragdb_data |>
  filter(
    str_detect(drug, globsets$drug_name),
    provean < globsets$provean_score) |>
  separate(
    codon_mutation,
    into = c("codon_original", "codon_mutation"),
    sep = "\\s*-\\s*",
    remove = FALSE
  ) |>
  mutate(
    AA_original   = aa_map[codon_original],
    AA_mutation = aa_map[codon_mutation],
    AA_pos = as.integer(AA_pos), 
    species =  stringr::str_trim(species),
    PMID = as.character(PMID)
  ) |>
  select(species, AA_original, AA_pos, AA_mutation, provean, PMID) |>
  distinct()

dragdb_group_pmid <- dragdb_processed_data |> 
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

dragdb_mapped_ecoli <- dragdb_group_pmid |>
  filter(species %in% refs$Species) |>
  mutate(species_sq = str_replace(species, " ", "_")) |>
  group_by(species_sq) |>
  group_modify(\(df, key) {
    sp <- key$species_sq[1]
    rpoB_target <- get_rpoB_for_species(sp, seqs)
    map_to_ecoli(df, rpoB_target, rpoB_Ecoli)
  }) |>
  ungroup() |>
  select(!species_sq) |>
  mutate(source = "DRAGdb")

write.csv(dragdb_mapped_ecoli, "data/processed_data/DRAGdb_processed.csv", row.names = FALSE)
saveRDS(dragdb_mapped_ecoli, file = "code/RIFRMap/data/dragdb.rds")


## 4. Yang's data -----
aa_seq <- translate(rpoB_Ecoli)
rpoB_aa <- as.character(translate(rpoB_Ecoli))

yang_mut_list <- yang_data |>
  pivot_longer(`_`:Q, names_to = "AA_mutation", values_to = "log2_fc") |>
  filter(log2_fc >= globsets$Yang_min_log2_fc) |>
  mutate(AA_pos_Ecoli = as.integer(paste0("5", pos)),
         AA_pos = AA_pos_Ecoli,
         species = "Escherichia coli",
         source = "Yang",
         AA_original = NA, 
         PMID = "37648864") |>
  select(species, AA_original, AA_pos, AA_pos_Ecoli, AA_mutation, log2_fc, source, PMID)

for (i in 1:nrow(mut_list_Yang)) {
  pos <- as.integer(mut_list_Yang$AA_pos_Ecoli[i])
  yang_mut_list$AA_original[i] <- substr(rpoB_aa, pos, pos)
}

write.csv(yang_mut_list, "data/processed_data/Yang_processed.csv", row.names = FALSE)
saveRDS(yang_mut_list, file = "code/RIFRMap/data/yang.rds")


## 5. Curated data (Negin's data) -----

mutation_raw <- curated_data |>
  filter(!is.na(AA_pos_Ecoli), !is.na(AA_original), !is.na(AA_mutation), !is.na(Ref_code)) |>
  dplyr::mutate(AA_original = toupper(AA_original),
                AA_mutation = toupper(AA_mutation),
                AA_pos = as.integer(AA_pos),
                AA_pos_Ecoli = as.integer(AA_pos_Ecoli))
  
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
    source = "Curated",
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
  
curated_processed <- mutation_species_rows |>
  left_join(mutation_counts, by = c("AA_pos_Ecoli", "AA_mutation")) |>
  dplyr::rename(
    species = Species
  )

write.csv(curated_processed, "data/processed_data/Curated_processed.csv", row.names = FALSE)
saveRDS(curated_processed, file = "code/RIFRMap/data/curated.rds")

## PointFinder -----




## Combined data -----

all_data <- dplyr::bind_rows(
  cryptic_mapped_ecoli,
  card_mapped_ecoli,
  dragdb_mapped_ecoli,
  yang_mut_list,
  curated_processed
)

# Process references

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
  so    <- extract_medline_field("SO  - ", lines)
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
  
  dp         <- extract_medline_field("DP  - ", lines)
  year       <- stringr::str_extract(dp, "\\b\\d{4}\\b")
  doi        <- extract_doi_medline(lines)
  author_key <- extract_author_key_medline(lines)
  ref_code   <- make_ref_code(author_key, year)
  
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

ref_table <- all_data |>
  dplyr::distinct(PMID) |>
  dplyr::pull(PMID) |>
  get_pubmed_ref_table()

ref_table_clean <- ref_table |>
  dplyr::mutate(
    Ref_code = as.character(Ref_code),
    Year     = as.character(Year),
    Full_ref = as.character(Full_ref),
    DOI      = as.character(DOI),
    PMID     = as.character(PMID)
  ) |>
  dplyr::select(Ref_code, Year, Full_ref, DOI, PMID)

curated_refs_clean <- curated_refs |>
  dplyr::mutate(
    Ref_code = as.character(Ref_code),
    Year     = as.character(Year),
    Full_ref = as.character(Full_ref),
    DOI      = as.character(DOI),
    PMID     = as.character(PMID)
  ) |>
  dplyr::select(dplyr::any_of(c("Ref_code", "Year", "Full_ref", "DOI", "PMID")))

combined_refs <- dplyr::bind_rows(
  curated_refs_clean,
  ref_table_clean
)

write_csv(
  combined_refs,
  "data/processed_data/all_references.csv"
)

write_csv(
  combined_refs,
  "code/RIFRMap/data/all_references.csv"
)
