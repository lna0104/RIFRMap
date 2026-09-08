# RIFRMap — Data Summary Script
library(dplyr)
library(tidyr)
library(readr)

file_name <- "output/RIFRMap_data_summary.txt"

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
all_data <- readRDS("processed_data/allData.rds")

# -----------------------------------------------------------------------------
# 2. Compute summary statistics — overall
# -----------------------------------------------------------------------------
n_species <- n_distinct(all_data$species)
n_positions <- n_distinct(all_data$AA_pos_Ecoli)
n_AA_mutations <- n_distinct(paste(all_data$AA_pos_Ecoli, all_data$AA_mutation))
total_entries <- nrow(all_data)
pos_min <- min(all_data$AA_pos_Ecoli, na.rm = TRUE)
pos_max <- max(all_data$AA_pos_Ecoli, na.rm = TRUE)

# -----------------------------------------------------------------------------
# 3. Top species by unique mutations
# -----------------------------------------------------------------------------
top_species <- all_data |>
  mutate(mut_id = paste(AA_pos_Ecoli, AA_mutation)) |>
  group_by(species) |>
  summarise(n_unique_mutations = n_distinct(mut_id), .groups = "drop") |>
  arrange(desc(n_unique_mutations)) |>
  slice_head(n = 5)

top_species_sources <- all_data |>
  filter(species %in% top_species$species) |>
  group_by(species) |>
  summarise(sources = paste(sort(unique(source)), collapse = ", "), .groups = "drop")

top_species <- top_species |>
  left_join(top_species_sources, by = "species")

top_species_text <- paste(
  sprintf(
    "    %d. %-45s %d\n       Sources: %s",
    seq_len(nrow(top_species)),
    top_species$species,
    top_species$n_unique_mutations,
    top_species$sources
  ),
  collapse = "\n\n"
)

# -----------------------------------------------------------------------------
# 4. Per-source summary
# -----------------------------------------------------------------------------
per_source <- all_data |>
  group_by(source) |>
  summarise(
    n_species = n_distinct(species),
    n_positions = n_distinct(AA_pos_Ecoli),
    n_mutations = n_distinct(paste(AA_pos_Ecoli, AA_mutation)),
    n_entries = n(),
    .groups = "drop"
  ) |>
  arrange(match(source, c(
    "CARD", "CRyPTIC", "DRAGdb", "PointFinder",
    "RefGenesCatalog", "Curated", "Yang"
  )))

# -----------------------------------------------------------------------------
# 5. Cross-source overlap
# -----------------------------------------------------------------------------
overlap <- all_data |>
  mutate(mut_id = paste(AA_pos_Ecoli, AA_mutation)) |>
  group_by(mut_id) |>
  summarise(n_sources = n_distinct(source), .groups = "drop")

n_pairs_1_source <- sum(overlap$n_sources == 1)
n_pairs_2plus <- sum(overlap$n_sources >= 2)
n_pairs_3plus <- sum(overlap$n_sources >= 3)
n_pairs_all7 <- sum(overlap$n_sources == 7)

# Mutations shared across all 7 sources — compact list only
all7_mutations <- all_data |>
  mutate(mut_id = paste(AA_pos_Ecoli, AA_mutation)) |>
  group_by(mut_id, AA_pos_Ecoli, AA_original, AA_mutation) |>
  summarise(n_sources = n_distinct(source), .groups = "drop") |>
  filter(n_sources == 7) |>
  arrange(AA_pos_Ecoli) |>
  mutate(label = paste0(AA_original, AA_pos_Ecoli, AA_mutation)) |>
  pull(label)

all7_text <- if (length(all7_mutations) == 0) {
  "    None found.\n"
} else {
  paste0("    ", paste(all7_mutations, collapse = ", "), "\n")
}

# -----------------------------------------------------------------------------
# 6. References cited
# -----------------------------------------------------------------------------
refs_lookup <- read_csv("processed_data/all_references.csv", show_col_types = FALSE)

n_refs_total <- nrow(dplyr::distinct(refs_lookup))

# -----------------------------------------------------------------------------
# 7. Per-source exclusive mutations (not in any other source)
# -----------------------------------------------------------------------------
mut_source_membership <- all_data |>
  mutate(mut_id = paste(AA_pos_Ecoli, AA_mutation)) |>
  group_by(mut_id) |>
  summarise(n_sources = n_distinct(source), .groups = "drop")

unique_per_source <- all_data |>
  mutate(mut_id = paste(AA_pos_Ecoli, AA_mutation)) |>
  left_join(mut_source_membership, by = "mut_id") |>
  filter(n_sources == 1) |>
  group_by(source) |>
  summarise(n_unique_pairs = n_distinct(mut_id), .groups = "drop") |>
  arrange(match(source, c(
    "CARD", "CRyPTIC", "DRAGdb", "PointFinder",
    "RefGenesCatalog", "Curated", "Yang"
  )))

per_source <- per_source |>
  left_join(unique_per_source, by = "source") |>
  mutate(n_unique_pairs = ifelse(is.na(n_unique_pairs), 0, n_unique_pairs))

per_source_text <- paste(
  apply(per_source, 1, function(r) {
    paste0(
      "    [", r["source"], "]\n",
      "        Species                                : ", r["n_species"], "\n",
      "        Unique E. coli positions               : ", r["n_positions"], "\n",
      "        Unique mutations                       : ", r["n_mutations"], "\n",
      "        Mutations exclusive to source          : ", r["n_unique_pairs"], "\n",
      "        Total entries                          : ", r["n_entries"]
    )
  }),
  collapse = "\n\n"
)

# -----------------------------------------------------------------------------
# 8. Build and write report
# -----------------------------------------------------------------------------
report <- paste0(
  "RIFRMap - Data Summary Report\n",
  "Date: ", Sys.time(), "\n",
  "--------------------------------------------------------------------\n\n",
  "OVERALL\n",
  "    Total entries (all sources, before dedup)  : ", total_entries, "\n",
  "    Total unique species                        : ", n_species, "\n",
  "    Total unique E. coli positions covered      : ", n_positions, "\n",
  "    Total unique mutations                      : ", n_AA_mutations, "\n",
  "    Total distinct references                   : ", n_refs_total, "\n",
  "    E. coli position range                      : ", pos_min, " -- ", pos_max, "\n\n",
  "--------------------------------------------------------------------\n\n",
  "TOP 5 SPECIES BY UNIQUE MUTATIONS\n",
  top_species_text, "\n\n",
  "--------------------------------------------------------------------\n\n",
  "PER SOURCE\n\n",
  per_source_text, "\n\n",
  "--------------------------------------------------------------------\n\n",
  "CROSS-SOURCE OVERLAP\n",
  "    Mutations in exactly 1 source       : ", n_pairs_1_source, "\n",
  "    Mutations in >= 2 sources           : ", n_pairs_2plus, "\n",
  "    Mutations in >= 3 sources           : ", n_pairs_3plus, "\n",
  "    Mutations in all 7 sources          : ", n_pairs_all7, "\n\n",
  "MUTATIONS REPORTED IN ALL 7 SOURCES\n",
  all7_text, "\n",
  "--------------------------------------------------------------------\n"
)

write_file(report, file_name)
cat(report)
cat(paste0("\nThis summary has been saved in file ", file_name, ".\n"))
