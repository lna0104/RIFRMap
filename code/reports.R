# RIFRMap — Data Summary Script
library(dplyr)
library(readr)

file_name <- "output/RIFRMap_data_summary.txt"

# Display order of the per-source sections (values of allData$source);
# anything not listed sorts last.
SOURCE_ORDER <- c(
  "CARD", "CRyPTIC", "DRAGdb", "PointFinder",
  "Reference Gene Catalog", "RifReport", "Yang-RIF"
)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
all_data <- readRDS("processed_data/allData.rds")

# -----------------------------------------------------------------------------
# 2. Raw inputs and filtering funnels
# -----------------------------------------------------------------------------
# report_filtering() traces every source from its raw file to the entries the
# app serves, one line per inclusion criterion. It lives in code/process_data.R
# so it sits beside the filters it mirrors and reads the same `globsets`
# thresholds - see the "Filtering report" section at the bottom of that file.
#
# Setting `.pipeline_orchestrating` suppresses process_data.R's standalone
# driver block, so sourcing it defines the reporting functions without running
# the pipeline (see the guarded block at the bottom of code/process_data.R).
.pipeline_orchestrating <- TRUE
suppressPackageStartupMessages(source("code/process_data.R"))
rm(.pipeline_orchestrating) # only checked at source time; don't leave it set

filtering_funnels <- report_filtering(final_counts = table(all_data$source))
raw_summary_text <- format_filtering_funnels(filtering_funnels)
raw_total <- sum(filtering_funnels$kept[filtering_funnels$kind == "raw"])

# -----------------------------------------------------------------------------
# 3. Compute summary statistics — overall
# -----------------------------------------------------------------------------
n_species <- n_distinct(all_data$species)
n_positions <- n_distinct(all_data$AA_pos_Ecoli)
n_AA_mutations <- n_distinct(paste(all_data$AA_pos_Ecoli, all_data$AA_mutation))
total_entries <- nrow(all_data)
pos_min <- min(all_data$AA_pos_Ecoli, na.rm = TRUE)
pos_max <- max(all_data$AA_pos_Ecoli, na.rm = TRUE)

# -----------------------------------------------------------------------------
# 4. Top species by unique mutations
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
# 5. Per-source summary
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
  arrange(match(source, SOURCE_ORDER))

# -----------------------------------------------------------------------------
# 6. Cross-source overlap
# -----------------------------------------------------------------------------
# Mutation identity is (AA_pos_Ecoli, AA_mutation) - AA_original is not part
# of the key, since species differ in which residue sits at a given E. coli
# position. The differing originals are collapsed into the label (e.g. D/Q516K).
overlap <- all_data |>
  mutate(mut_id = paste(AA_pos_Ecoli, AA_mutation)) |>
  group_by(mut_id, AA_pos_Ecoli, AA_mutation) |>
  summarise(
    n_sources = n_distinct(source),
    label = paste0(
      paste(sort(unique(na.omit(AA_original))), collapse = "/"),
      dplyr::first(AA_pos_Ecoli), dplyr::first(AA_mutation) # S4Vectors masks first()
    ),
    .groups = "drop"
  )

n_sources_total <- n_distinct(all_data$source)
n_pairs_1_source <- sum(overlap$n_sources == 1)
n_pairs_2plus <- sum(overlap$n_sources >= 2)
n_pairs_3plus <- sum(overlap$n_sources >= 3)
n_pairs_all <- sum(overlap$n_sources == n_sources_total)

# Mutations shared across every source — compact list only
all_sources_mutations <- overlap |>
  filter(n_sources == n_sources_total) |>
  arrange(AA_pos_Ecoli, AA_mutation) |>
  pull(label)

all_sources_text <- if (length(all_sources_mutations) == 0) {
  "    None found.\n"
} else {
  paste0("    ", paste(all_sources_mutations, collapse = ", "), "\n")
}

# -----------------------------------------------------------------------------
# 7. References cited
# -----------------------------------------------------------------------------
refs_lookup <- read_csv("processed_data/all_references.csv", show_col_types = FALSE)

n_refs_total <- nrow(dplyr::distinct(refs_lookup))

# -----------------------------------------------------------------------------
# 8. Per-source exclusive mutations (not in any other source)
# -----------------------------------------------------------------------------
unique_per_source <- all_data |>
  mutate(mut_id = paste(AA_pos_Ecoli, AA_mutation)) |>
  left_join(select(overlap, mut_id, n_sources), by = "mut_id") |>
  filter(n_sources == 1) |>
  group_by(source) |>
  summarise(n_unique_pairs = n_distinct(mut_id), .groups = "drop")

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
# 9. Build and write report
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
  "RAW DATA AND FILTERING\n",
  "    Records across all raw input files          : ", raw_total, "\n",
  "    Entries served by RIFRMap                   : ", total_entries, "\n\n",
  "    Each block below traces one source from its raw file to the entries\n",
  "    the app serves. Named stages are the inclusion criteria applied in\n",
  "    code/process_data.R; the drop on the final line covers mutation-code\n",
  "    parsing, de-duplication and reference-gated mapping, which are not\n",
  "    broken out individually.\n\n",
  raw_summary_text, "\n\n",
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
  "    Mutations in all ", n_sources_total, " sources          : ", n_pairs_all, "\n\n",
  "MUTATIONS REPORTED IN ALL ", n_sources_total, " SOURCES\n",
  all_sources_text, "\n",
  "--------------------------------------------------------------------\n"
)

write_file(report, file_name)
cat(report)
cat(paste0("\nThis summary has been saved in file ", file_name, ".\n"))
