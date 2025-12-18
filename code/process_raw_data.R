library(tidyverse)
library(tidyr)
library(dplyr)
library(ALJEbinf)
library(Biostrings)

source("code/process_mutations.R")

# Global settings and parameters:
globsets <- list(
  drug_name = "Rifampicin",
  drug_code = "RIF",
  gene_name = "rpoB",
  variant_type = "single resistance variant", # filter variant type for CARD data
  p_value = 0.05, # maxmimum p_value for CRyPTIC data
  provean_score = -2.5, # maximum provean score for DRAGdb data
  Yang_min_log2_fc = 0 # minimum log2(foldchange) in frequencies to include a mutation as "selected for"
 )

# Load database 
cryptic_raw_data <- read.csv("./data/raw_data/CRyPTIC_muts.csv", stringsAsFactors = FALSE)
card_raw_data <- read.table("./data/raw_data/CARD_snps.txt", header = TRUE, sep = "\t")
dragdb_raw_data <- read.csv("./data/raw_data/DRAGdb_RIF.csv", header = FALSE,
                            stringsAsFactors = FALSE)
muts <- read_csv("./data/raw_data/muts.csv", show_col_types = FALSE) #reported data from Negin's study

# Load sequencing data
seqs <- readDNAStringSet("./data/rpoB_references.fasta")
rpoB_Ecoli <- seqs[["rpoB_Escherichia_coli_MG1655"]]
rpoB_TB <- seqs[["rpoB_Mycobacterium_tuberculosis_H37Rv"]]

# Process database
## 1. CRyPTIC data (only for MTB)
cryptic_processed_data <- cryptic_raw_data |>
  rename_with(tolower) |>
  filter(drug == globsets$drug_code,
         bh_corrected_p < globsets$p_value,
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
  select(!c("drug", "gene_symbol", "mutations", "AA_original", "gene", "position")) |>
  mutate(species = "Mycobacterium tuberculosis",
         AA_pos = as.integer(AA_pos)) |>
  distinct()

cryptic_mapped_ecoli <- map_to_ecoli(cryptic_processed_data, rpoB_TB, rpoB_Ecoli)
  
write.csv(cryptic_mapped_ecoli, "data/processed_data/CRyPTIC_rpoB_processed.csv", row.names = FALSE)

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
  select(species, AA_mutation, AA_pos) |>
  distinct()

card_mapped_ecoli <- card_processed_data |>
  mutate(species_sq = str_replace(species, " ", "_")) |>
  group_by(species_sq) |>
  group_modify(\(df, key) {
    sp <- key$species_sq[1]
    rpoB_target <- get_rpoB_for_species(sp, seqs)
    map_to_ecoli(df, rpoB_target, rpoB_Ecoli)
  }) |>
  ungroup() |>
  select(!species_sq) 

write.csv(card_mapped_ecoli, "data/processed_data/CARD_rpoB_processed.csv", row.names = FALSE)


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
  ) |>
  mutate(
    AA_pos = as.integer(AA_pos)
  ) |>
  select(species, AA_pos, AA_mutation, provean) |>
  distinct()

dragdb_mapped_ecoli <- map_to_ecoli(dragdb_processed_data, rpoB_TB, rpoB_Ecoli)

write.csv(dragdb_mapped_ecoli, "data/processed_data/DRAGdb_rpoB_processed.csv", row.names = FALSE)


## 4. Yang's data

mut_list_Yang <- read_csv("./data/raw_data/Yang2023_Fig3a_data_RIF_log2FC.csv", show_col_types = FALSE) |>
  pivot_longer(`_`:Q, names_to = "AA_mutation", values_to = "log2_fc") |>
  filter(log2_fc >= globsets$Yang_min_log2_fc) |>
  mutate(AA_pos_Ecoli = as.integer(paste0("5", pos)),
         species = "Escherichia coli") |>
  select(species, AA_pos_Ecoli, AA_mutation, log2_fc)

write.csv(mut_list_Yang, "data/processed_data/Yang_processed.csv", row.names = FALSE)


## 5. Negin's data

mutation_list <- muts |>
  filter(!is.na(AA_pos_Ecoli), !is.na(AA_original), !is.na(AA_mutation)) 

mutation_list_studies <- mutation_list |>
  select(Ref_code, AA_pos_Ecoli, AA_mutation) |>
  distinct() |>
  group_by(AA_pos_Ecoli, AA_mutation) |>
  summarise(n_studies = n(), .groups = "drop")

mutation_list_final <- mutation_list |>
  mutate(species = Species) |>
  select(species, AA_pos_Ecoli, AA_mutation) |>
  distinct() |>
  left_join(mutation_list_studies, by = c("AA_pos_Ecoli", "AA_mutation"))

write.csv(mutation_list_final, "data/processed_data/Negin_processed.csv", row.names = FALSE)


## Create the mutation comparison dataframe
combined_data <- bind_rows(
  mutation_list |> mutate(source = "Negin"),
  cryptic_mapped_ecoli |> mutate(source = "CRyPTIC"),
  card_mapped_ecoli    |>  mutate(source = "CARD"),
  dragdb_mapped_ecoli  |> mutate(source = "DRAGdb"), 
  mut_list_Yang |> mutate(source = "Yang")) |>  # Combine sources
  distinct() 

write.csv(combined_data, "data/processed_data/combined_data.csv", row.names = FALSE)


