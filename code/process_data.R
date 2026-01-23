library(tidyverse)
library(tidyr)
library(dplyr)
library(ALJEbinf)
library(Biostrings)

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

# Load database 
cryptic_raw_data <- read.csv("./data/raw_data/CRyPTIC_muts.csv", stringsAsFactors = FALSE)
card_raw_data <- read.table("./data/raw_data/CARD_snps.txt", header = TRUE, sep = "\t")
dragdb_raw_data <- read.csv("./data/raw_data/DRAGdb_RIF.csv", header = FALSE,
                            stringsAsFactors = FALSE)
muts <- read.csv("./data/raw_data/muts.csv") #curated data

# Load sequencing data
seqs <- readDNAStringSet("./data/rpoB_references.fasta")
rpoB_Ecoli <- seqs[["rpoB_Escherichia_coli_MG1655"]]
rpoB_TB <- seqs[["rpoB_Mycobacterium_tuberculosis_H37Rv"]]
refs <- read_csv("data/rpoB_references.csv") |>
  filter(RefStrain == TRUE)
  
# Process database
## 1. CRyPTIC data (only for MTB)
cryptic_processed_data <- cryptic_raw_data |>
  rename_with(tolower) |>
  filter(drug == globsets$drug_code,
         #bh_corrected_p < globsets$p_value,
         gene == globsets$gene_name,
         !str_detect(muts, "indel") # get rif of indel mutation
      ) |> 
  separate(
    muts,
    into = c("gene_symbol", "mutations"),
    sep = "_"
  ) |>
  extract(
    mutations,
    into  = c("AA_original", "AA_pos", "AA_mutation"),
    regex = "^(.)(-?[0-9]+)(.)$",
    remove = FALSE
  ) |>
  # remove rows with missing or empty mutation information 
  filter(!is.na(mutations) & !is.na(AA_original) & gene_symbol != "") |>
  # remove redundant columns
  select(!c("drug", "gene_symbol", "mutations", "gene", "position")) |>
  mutate(species = "Mycobacterium tuberculosis",
         AA_pos = as.integer(AA_pos),
         source = "CRyPTIC") |>
  distinct()

cryptic_mapped_ecoli <- map_to_ecoli(cryptic_processed_data, rpoB_TB, rpoB_Ecoli) 
  
write.csv(cryptic_mapped_ecoli, "data/processed_data/CRyPTIC_rpoB_processed.csv", row.names = FALSE)
saveRDS(cryptic_mapped_ecoli, file = "code/RIFRMap_shiny/data/cryptic.rds")

## 2. CARD data
card_processed_data <- card_raw_data |>
  rename_with(tolower) |>
  separate(
    card.short.name,
    into = c("species_abbrev", "gene", "drug"),
    sep = "_",
    fill = "right",    
    remove = TRUE
  ) |>
  filter(drug == globsets$drug_code,
         parameter.type == globsets$variant_type,
         gene == globsets$gene_name) |>
  mutate(
    species = str_extract(name, "^\\S+\\s+\\S+") 
         ) |>
  extract(
    mutations,
    into = c("AA_original", "AA_pos", "AA_mutation"),
    regex = "([A-Za-z]+)[^0-9]*([0-9]+)[^A-Za-z0-9]*([A-Za-z]+)",
    remove = FALSE 
  ) |>
  mutate(AA_pos = as.integer(AA_pos)) |>
  filter(AA_mutation != "Var") |>
  select(species, AA_original, AA_mutation, AA_pos) |>
  distinct()

card_mapped_ecoli <- card_processed_data |>
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
  mutate(source = "CARD")

write.csv(card_mapped_ecoli, "data/processed_data/CARD_rpoB_processed.csv", row.names = FALSE)
saveRDS(card_mapped_ecoli, file = "code/RIFRMap_shiny/data/card.rds")

## 3. DRAGdb

colnames(dragdb_raw_data) <- c(
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

aa_map <- c(
  Ala="A", Arg="R", Asn="N", Asp="D", Cys="C",
  Glu="E", Gln="Q", Gly="G", His="H", Ile="I",
  Leu="L", Lys="K", Met="M", Phe="F", Pro="P",
  Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V"
)

dragdb_processed_data <- dragdb_raw_data |>
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
    species =  stringr::str_trim(species)
  ) |>
  select(species, AA_original, AA_pos, AA_mutation, provean) |>
  distinct()

dragdb_mapped_ecoli <- dragdb_processed_data |>
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

write.csv(dragdb_mapped_ecoli, "data/processed_data/DRAGdb_rpoB_processed.csv", row.names = FALSE)
saveRDS(dragdb_mapped_ecoli, file = "code/RIFRMap_shiny/data/dragdb.rds")


## 4. Yang's data
aa_seq <- translate(rpoB_Ecoli)
rpoB_aa <- as.character(translate(rpoB_Ecoli))

mut_list_Yang <- read_csv("./data/raw_data/Yang2023_Fig3a_data_RIF_log2FC.csv", show_col_types = FALSE) |>
  pivot_longer(`_`:Q, names_to = "AA_mutation", values_to = "log2_fc") |>
  filter(log2_fc >= globsets$Yang_min_log2_fc) |>
  mutate(AA_pos_Ecoli = as.integer(paste0("5", pos)),
         AA_pos = AA_pos_Ecoli,
         species = "Escherichia coli",
         source = "Yang",
         AA_original = NA) |>
  select(species, AA_original,  AA_pos, AA_pos_Ecoli, AA_mutation, log2_fc, source)

for (i in 1:nrow(mut_list_Yang)) {
  pos <- as.integer(mut_list_Yang$AA_pos_Ecoli[i])
  mut_list_Yang$AA_original[i] <- substr(rpoB_aa, pos, pos)
}

write.csv(mut_list_Yang, "data/processed_data/Yang_processed.csv", row.names = FALSE)
saveRDS(mut_list_Yang, file = "code/RIFRMap_shiny/data/yang.rds")


## 5. Curated data (Negin's data)

mutation_raw <- muts |>
  filter(!is.na(AA_pos_Ecoli), !is.na(AA_original), !is.na(AA_mutation)) |>
  filter(is.na(Warning)) |>
  dplyr::mutate(AA_original = toupper(AA_original),
                AA_mutation = toupper(AA_mutation))
  
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
  
final_mutation_list <- mutation_species_rows |>
  left_join(mutation_counts, by = c("AA_pos_Ecoli", "AA_mutation")) |>
  mutate(species = Species) |>
  select(!Species)

write.csv(final_mutation_list, "data/processed_data/curated_processed.csv", row.names = FALSE)
saveRDS(final_mutation_list, file = "code/RIFRMap_shiny/data/curated.rds")

# Combine data
combined_all_rows <- bind_rows(
  cryptic_mapped_ecoli,
  card_mapped_ecoli,
  dragdb_mapped_ecoli,
  mut_list_Yang, 
  final_mutation_list
) 

write_csv(combined_all_rows, "data/processed_data/RIF_all_sources.csv")


