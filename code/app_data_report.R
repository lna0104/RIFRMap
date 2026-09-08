# RIFRMap — App Data Report
#
# Prints what is actually in processed_data/, organised the way the Shiny app
# presents it, so this report can be read side by side with the running app and
# any mismatch spotted by eye.
# Run from the repo root:
#   source("code/app_data_report.R")
#
# Writes output/app_data_report_<timestamp>.txt (gitignored via *.txt) and
# echoes the same content to the console. Reads only; touches nothing else.

source("global.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# ---- Output plumbing ----------------------------------------------------

report_lines <- character(0)

emit <- function(...) {
  lines <- c(...)
  if (length(lines) == 0) lines <- ""
  report_lines <<- c(report_lines, lines)
  writeLines(lines)
}

rule <- function(ch = "-", n = 78) paste(rep(ch, n), collapse = "")

section <- function(title) {
  emit("", rule("="), title, rule("="))
}

# The seven datasets, in the order the Data Filters tab lists them. `label` is
# the value carried in each .rds's own `source` column.
DATASETS <- tribble(
  ~name,         ~label,
  "card",        "CARD",
  "cryptic",     "CRyPTIC",
  "dragdb",      "DRAGdb",
  "yang",        "Yang",
  "rifreport",   "RifReport",
  "pointfinder", "PointFinder",
  "refgenes",    "Reference Gene Catalog"
)

get_dataset <- function(name) get(name, envir = globalenv())

# formatC() pads to `digits` width, which wrecks column alignment - trim it.
fmt <- function(x) trimws(formatC(x, format = "g", digits = 6))

# The app keys mutation identity on (AA_pos_Ecoli, AA_mutation) 
mut_label <- function(orig, pos, mut) {
  paste0(paste(sort(unique(orig)), collapse = "/"), pos[1], mut[1])
}

# ---- Header -------------------------------------------------------------

stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

emit(rule("="))
emit("RIFRMap - App Data Report")
emit(rule("="))
emit(sprintf("Generated:          %s", stamp))
emit(sprintf("App version:        %s", APP_VERSION))
emit(sprintf("Data last updated:  %s", DATA_LAST_UPDATED))
emit("")
emit("Read this next to a running app (shiny::runApp()). Each section names the")
emit("app surface it corresponds to. Nothing here is a pass/fail verdict.")

# ---- 1. Filter controls -------------------------------------------------

section("1. FILTER CONTROLS  (Explore -> Data Filters)")

emit("")
emit("Data range = the actual min/max in processed_data. UI range = the bounds")
emit("the control offers, from SLIDER_RANGES in global.R. If the UI range does")
emit("not cover the data range, the uncovered rows are unreachable in the app.")
emit("")

# The numeric filters applied in server.R's combined_rows(). cryptic_bhp is a
# radioButtons rather than a slider, so it has no SLIDER_RANGES entry.
FILTER_CONTROLS <- tribble(
  ~control,             ~dataset,    ~column,
  "cryptic_coef",       "cryptic",   "coef",
  "cryptic_counts",     "cryptic",   "counts",
  "cryptic_bhp",        "cryptic",   "bh_corrected_p",
  "dragdb_provean",     "dragdb",    "provean",
  "yang_log2fc",        "yang",      "log2_fc",
  "rifreport_nstudies", "rifreport", "n_studies",
  "rifreport_nspecies", "rifreport", "n_species"
)

emit(sprintf(
  "%-20s %-22s %-16s %-9s %s",
  "control", "data range", "UI range", "covers?", "NAs"
))
emit(rule("-"))

uncovered <- list()

for (i in seq_len(nrow(FILTER_CONTROLS))) {
  ctl <- FILTER_CONTROLS$control[i]
  d <- get_dataset(FILTER_CONTROLS$dataset[i])
  col <- FILTER_CONTROLS$column[i]
  x <- d[[col]]

  n_na <- sum(is.na(x))
  data_rng <- range(x, na.rm = TRUE)
  ui_rng <- SLIDER_RANGES[[ctl]]

  if (is.null(ui_rng)) {
    # radioButtons: choices are 0.05 / 0.01 / 0.005 / 0.001 / None (= 1.0),
    # defaulting to None, so the widest reachable filter is <= 1.
    ui_txt <- "<= 1 (radio)"
    covers <- data_rng[2] <= 1
  } else {
    ui_txt <- sprintf("%s - %s", fmt(ui_rng[1]), fmt(ui_rng[2]))
    covers <- ui_rng[1] <= data_rng[1] && ui_rng[2] >= data_rng[2]
    if (!covers) uncovered[[ctl]] <- list(d = d, col = col, ui = ui_rng)
  }

  emit(sprintf(
    "%-20s %-22s %-16s %-9s %d",
    ctl,
    sprintf("%s - %s", fmt(data_rng[1]), fmt(data_rng[2])),
    ui_txt,
    if (covers) "yes" else "NO",
    n_na
  ))
}

# Spell out exactly what an uncovered control costs - a row count alone would
# understate it, since the lost rows cluster on the most-reported mutations.
if (length(uncovered) > 0) {
  emit("")
  emit(rule("-"))
  emit("UNREACHABLE ROWS - the app cannot show these at any setting:")
  for (ctl in names(uncovered)) {
    u <- uncovered[[ctl]]
    lost <- u$d |>
      filter(!dplyr::between(.data[[u$col]], u$ui[1], u$ui[2]))
    emit("")
    emit(sprintf("  %s: %d of %d rows", ctl, nrow(lost), nrow(u$d)))
    muts <- lost |>
      count(AA_pos_Ecoli, AA_original, AA_mutation, sort = TRUE)
    for (j in seq_len(min(10, nrow(muts)))) {
      emit(sprintf(
        "    %s%s%s  (%d rows)",
        muts$AA_original[j], muts$AA_pos_Ecoli[j], muts$AA_mutation[j], muts$n[j]
      ))
    }
    if (nrow(muts) > 10) emit(sprintf("    ... and %d more", nrow(muts) - 10))
  }
}

emit("")
emit("cryptic_bhp is a radioButtons (0.05 / 0.01 / 0.005 / 0.001 / None), and")
emit("defaults to None = 1.0, so by default it filters nothing.")
emit("")

# ---- 2. Species per source ----------------------------------------------

section("2. SPECIES PER SOURCE  (Explore -> Data Filters, species dropdowns)")

ref_species <- unique(refs_seqs$Species)

emit("")
emit(sprintf(
  "Reference catalog (raw_data/rpoB_references.csv): %d species",
  length(ref_species)
))
emit("")
emit("A species with no reference sequence is dropped from the Mutation Map")
emit("(cluster_species_to_plot() warns about this) but still appears in the")
emit("tables and the other plots. Expected, not a fault - listed for awareness.")

for (i in seq_len(nrow(DATASETS))) {
  d <- get_dataset(DATASETS$name[i])
  sp <- sort(unique(d$species))
  missing_ref <- setdiff(sp, ref_species)

  emit("")
  emit(rule("-"))
  emit(sprintf("%s - %d species, %d rows", DATASETS$label[i], length(sp), nrow(d)))
  emit(rule("-"))
  emit(paste("  ", paste(sp, collapse = ", ")))
  if (length(missing_ref) > 0) {
    emit("")
    emit(sprintf(
      "  no reference sequence (%d): %s",
      length(missing_ref), paste(missing_ref, collapse = ", ")
    ))
  }
}

# ---- 3. Per-source totals -----------------------------------------------

section("3. PER-SOURCE TOTALS  (Explore -> Summary, and the Data tab)")

emit("")
emit("The Data tab is unfiltered (all_rows()), so its row counts should match")
emit("the 'rows' column exactly. The Explore tab is filtered, so it will match")
emit("only when every dataset is enabled at default filter settings.")
emit("")

emit(sprintf("%-24s %8s %10s %11s %11s", "source", "rows", "species", "positions", "mutations"))
emit(rule("-"))

for (i in seq_len(nrow(DATASETS))) {
  d <- get_dataset(DATASETS$name[i])
  emit(sprintf(
    "%-24s %8d %10d %11d %11d",
    DATASETS$label[i],
    nrow(d),
    n_distinct(d$species),
    n_distinct(d$AA_pos_Ecoli),
    n_distinct(paste(d$AA_pos_Ecoli, d$AA_mutation))
  ))
}

emit(rule("-"))
emit(sprintf(
  "%-24s %8d %10d %11d %11d",
  "allData.rds (combined)",
  nrow(allData <- readRDS("processed_data/allData.rds")),
  n_distinct(allData$species),
  n_distinct(allData$AA_pos_Ecoli),
  n_distinct(paste(allData$AA_pos_Ecoli, allData$AA_mutation))
))


# ---- 4. Mutation Map coverage -------------------------------------------

section("4. MUTATION MAP COVERAGE  (Explore -> Mutation Map)")

pos_in_cluster <- function(pos) {
  vapply(
    pos,
    function(p) isTRUE(any(p >= clusters$AAstart_Ecoli & p <= clusters$AAend_Ecoli)),
    logical(1)
  )
}

emit("")
emit("Clusters defined in raw_data/clusters.csv:")
for (i in seq_len(nrow(clusters))) {
  emit(sprintf(
    "  %-4s %d - %d", clusters$Cluster[i],
    clusters$AAstart_Ecoli[i], clusters$AAend_Ecoli[i]
  ))
}
emit("")
emit("Positions outside every cluster are not drawn in the default cluster")
emit("view - they appear only under the \"Show all mutations\" toggle. That is")
emit("by design, and is the usual explanation for a mutation looking absent.")
emit("")

allData$in_cluster <- pos_in_cluster(allData$AA_pos_Ecoli)

emit(sprintf("Distinct positions overall:      %d", n_distinct(allData$AA_pos_Ecoli)))
emit(sprintf(
  "  inside a cluster:              %d",
  n_distinct(allData$AA_pos_Ecoli[allData$in_cluster])
))
emit(sprintf(
  "  outside every cluster:         %d",
  n_distinct(allData$AA_pos_Ecoli[!allData$in_cluster])
))
emit(sprintf(
  "Rows outside every cluster:      %d of %d",
  sum(!allData$in_cluster), nrow(allData)
))
emit("")

by_species <- allData |>
  group_by(species) |>
  summarise(
    positions = n_distinct(AA_pos_Ecoli),
    in_cl     = n_distinct(AA_pos_Ecoli[in_cluster]),
    out_cl    = n_distinct(AA_pos_Ecoli[!in_cluster]),
    .groups   = "drop"
  ) |>
  mutate(has_ref = species %in% ref_species) |>
  arrange(desc(positions))

emit(sprintf("%-42s %10s %8s %8s %6s", "species", "positions", "in", "out", "ref?"))
emit(rule("-"))
for (i in seq_len(nrow(by_species))) {
  emit(sprintf(
    "%-42s %10d %8d %8d %6s",
    by_species$species[i], by_species$positions[i],
    by_species$in_cl[i], by_species$out_cl[i],
    if (by_species$has_ref[i]) "yes" else "NO"
  ))
}

# ---- 5. Top positions and mutations -------------------------------------

section("5. TOP POSITIONS AND MUTATIONS  (Explore -> Mutation Frequency)")

emit("")
emit("Most frequent positions (E. coli numbering):")
emit("")
emit(sprintf("%-10s %8s %10s %10s", "position", "rows", "species", "sources"))
emit(rule("-"))

top_pos <- allData |>
  group_by(AA_pos_Ecoli) |>
  summarise(
    rows = n(), species = n_distinct(species),
    sources = n_distinct(source), .groups = "drop"
  ) |>
  arrange(desc(rows)) |>
  slice_head(n = 15)

for (i in seq_len(nrow(top_pos))) {
  emit(sprintf(
    "%-10d %8d %10d %10d", top_pos$AA_pos_Ecoli[i],
    top_pos$rows[i], top_pos$species[i], top_pos$sources[i]
  ))
}

emit("")
emit("Most frequent mutations:")
emit("")
emit(sprintf("%-12s %8s %10s %10s", "mutation", "rows", "species", "sources"))
emit(rule("-"))

top_mut <- allData |>
  group_by(AA_pos_Ecoli, AA_mutation) |>
  summarise(
    label = mut_label(AA_original, AA_pos_Ecoli, AA_mutation),
    rows = n(), species = n_distinct(species),
    sources = n_distinct(source), .groups = "drop"
  ) |>
  arrange(desc(rows)) |>
  slice_head(n = 15)

for (i in seq_len(nrow(top_mut))) {
  emit(sprintf(
    "%-12s %8d %10d %10d",
    top_mut$label[i], top_mut$rows[i], top_mut$species[i], top_mut$sources[i]
  ))
}

# The RRDR hotspots are listed unconditionally, whether or not they make the
# top 15 above: these are the mutations whose disappearance matters most, so
# they get a fixed place in the report rather than a ranked one.
RRDR_POSITIONS <- c(513, 516, 522, 526, 531, 533)

emit("")
emit(rule("-"))
emit(sprintf(
  "RRDR hotspots (positions %s)",
  paste(RRDR_POSITIONS, collapse = ", ")
))
emit(rule("-"))
emit("")
emit(sprintf("%-12s %8s %10s   %s", "mutation", "rows", "species", "sources"))
emit(rule("-"))

hot <- allData |>
  filter(AA_pos_Ecoli %in% RRDR_POSITIONS) |>
  group_by(AA_pos_Ecoli, AA_mutation) |>
  summarise(
    label = mut_label(AA_original, AA_pos_Ecoli, AA_mutation),
    rows = n(), species = n_distinct(species),
    src = paste(sort(unique(source)), collapse = ", "), .groups = "drop"
  ) |>
  arrange(AA_pos_Ecoli, desc(rows))

for (i in seq_len(nrow(hot))) {
  emit(sprintf(
    "%-12s %8d %10d   %s",
    hot$label[i], hot$rows[i], hot$species[i], hot$src[i]
  ))
}

# ---- 6. Cross-source overlap --------------------------------------------

section("6. CROSS-SOURCE OVERLAP  (Explore -> Dataset Comparison)")

emit("")
emit("How many distinct mutations are shared by how many sources. Compare")
emit("against the upset plot's bar heights.")
emit("")

overlap <- allData |>
  group_by(AA_pos_Ecoli, AA_mutation) |>
  summarise(n_sources = n_distinct(source), .groups = "drop")

emit(sprintf("%-12s %12s", "in N sources", "mutations"))
emit(rule("-"))
tab <- overlap |>
  count(n_sources) |>
  arrange(n_sources)
for (i in seq_len(nrow(tab))) {
  emit(sprintf("%-12d %12d", tab$n_sources[i], tab$n[i]))
}
emit(rule("-"))
emit(sprintf("%-12s %12d", "total", nrow(overlap)))

emit("")
emit("Mutations present in 4 or more sources:")
emit("")

shared <- allData |>
  group_by(AA_pos_Ecoli, AA_mutation) |>
  summarise(
    label = mut_label(AA_original, AA_pos_Ecoli, AA_mutation),
    n_sources = n_distinct(source),
    src = paste(sort(unique(source)), collapse = ", "), .groups = "drop"
  ) |>
  filter(n_sources >= 4) |>
  arrange(desc(n_sources), AA_pos_Ecoli)

if (nrow(shared) == 0) {
  emit("  (none)")
} else {
  for (i in seq_len(nrow(shared))) {
    emit(sprintf(
      "  %-10s %d sources: %s",
      shared$label[i], shared$n_sources[i], shared$src[i]
    ))
  }
}

# ---- Write --------------------------------------------------------------

if (!dir.exists("output")) dir.create("output", recursive = TRUE)

out_path <- file.path(
  "output",
  sprintf("app_data_report_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
writeLines(report_lines, out_path)

emit("")
emit(rule("="))
emit(sprintf("Report written to %s", out_path))
emit(rule("="))


