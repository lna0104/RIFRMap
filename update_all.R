# update_all.R
#
# Single entry point to refresh every RIFRMap dataset end to end:
#   - downloads whatever has a scripted source (CARD, DRAGdb, PointFinder,
#     RifReport, RefGenes), validating each and continuing past individual failures;
#   - discovers species newly seen across all 7 datasets that are missing an
#     rpoB reference, and attempts to fetch one from NCBI for each;
#   - regenerates the E. coli coordinate map if the reference fasta changed;
#   - skips reprocessing any dataset whose raw data hasn't changed since the
#     last successful run (processed_data/pipeline_state.rds);
#   - regenerates processed_data/allData.rds and all_references.csv;
#   - prints a per-dataset summary and appends it to
#     processed_data/last_update_log.csv.
#
# Run with:
#   source("update_all.R")
#
# For debugging a single source instead of a full refresh, source
# update_data.R / process_data.R / process_refseqs.R directly - each still
# works standalone.

library(tidyverse)
library(Biostrings)
library(openssl)

.pipeline_orchestrating <- TRUE

source("helpers.R")
source("update_data.R")
source("process_data.R")
source("process_refseqs.R")

PIPELINE_STATE_PATH <- "processed_data/pipeline_state.rds"
LOG_PATH <- "processed_data/last_update_log.csv"

RAW_PATHS <- list(
  card        = "raw_data/card_data.csv",
  cryptic     = "raw_data/cryptic_data.csv",
  dragdb      = "raw_data/dragdb_data.csv",
  yang        = "raw_data/Yang2023_Fig3a_data_RIF_log2FC.csv",
  rifreport   = "raw_data/rifreport_data.csv",
  pointfinder = "raw_data/pointfinder_data.csv",
  refgenes    = "raw_data/refgenes_catalog.csv"
)

SCRIPTED_DATASETS <- c("card", "dragdb", "pointfinder", "rifreport", "refgenes")
MANUAL_DATASETS <- c("cryptic", "yang")

hash_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  as.character(openssl::sha1(file(path)))
}

load_pipeline_state <- function() {
  if (file.exists(PIPELINE_STATE_PATH)) readRDS(PIPELINE_STATE_PATH) else list()
}

save_pipeline_state <- function(state) {
  saveRDS(state, file = PIPELINE_STATE_PATH)
}

load_existing_rds <- function(path) {
  if (file.exists(path)) readRDS(path) else NULL
}

safe_parse <- function(dataset, fn) {
  tryCatch(
    list(ok = TRUE, value = fn()),
    error = function(e) {
      message("Could not parse ", dataset, " for species discovery: ", conditionMessage(e))
      list(ok = FALSE, value = NULL)
    }
  )
}

message("=== RIFRMap data pipeline refresh starting ===")

# ---- 1. Update scripted sources ----
update_results <- list()
update_results$card <- update_card()
update_results$dragdb <- update_dragdb()
update_results$pointfinder <- update_pointfinder()
update_results$refgenes <- update_refgenes()

# usethis::create_github_token()
# Generate a new token select all scopes
# gitcreds::gitcreds_set()
token <- tryCatch(gitcreds::gitcreds_get()$password, error = function(e) NA_character_)
rifreport_dl <- update_rifreport_download(token)
update_results$rifreport <- rifreport_dl

# ---- 2. Lightly parse whichever raw data is available (pre-reference-filter) ----
parsed_card <- safe_parse("card", parse_card)
parsed_dragdb <- safe_parse("dragdb", parse_dragdb)
parsed_pointfinder <- safe_parse("pointfinder", parse_pointfinder)
parsed_refgenes <- safe_parse("refgenes", parse_refgenes)

species_pool <- c(
  if (parsed_card$ok) parsed_card$value$species else character(0),
  if (parsed_dragdb$ok) parsed_dragdb$value$species else character(0),
  if (parsed_pointfinder$ok) parsed_pointfinder$value$species else character(0),
  if (parsed_refgenes$ok) parsed_refgenes$value$species else character(0),
  if (identical(rifreport_dl$status, "ok")) rifreport_dl$species else character(0)
)

# ---- 3. Discover & fetch rpoB references for any newly-seen species ----
discovery <- discover_and_fetch_new_references(species_pool)
if (length(discovery$added_species) > 0) {
  message("Added rpoB references for: ", paste(discovery$added_species, collapse = ", "))
}
if (length(discovery$failed_species) > 0) {
  message("Could not find rpoB references for: ", paste(discovery$failed_species, collapse = ", "))
}

# ---- 4. Regenerate coordinates (picks up any newly-appended sequences) ----
coordinates <- regenerate_coordinates_if_needed()
DR
# ---- 5. Finish RifReport now that coordinates is available ----
refs <- read_csv("raw_data/rpoB_references.csv", show_col_types = FALSE)
seqs <- readDNAStringSet("raw_data/rpoB_references.fasta")

if (identical(rifreport_dl$status, "ok")) {
  update_results$rifreport <- update_rifreport_finalize(rifreport_dl$raw_mutations, refs, seqs, coordinates)
}

# ---- 6. Change detection: decide which datasets need reprocessing ----
state <- load_pipeline_state()
current_hashes <- lapply(RAW_PATHS, hash_file)

needs_processing <- vapply(names(RAW_PATHS), function(ds) {
  prev <- state[[ds]]$hash
  cur <- current_hashes[[ds]]
  is.null(prev) || is.na(cur) || !identical(prev, cur)
}, logical(1))

# ---- 7. Process (or skip) each dataset ----
ctx <- load_reference_context()
process_results <- list(
  cryptic = if (needs_processing[["cryptic"]]) process_cryptic(ctx) else NULL,
  card = if (needs_processing[["card"]] && parsed_card$ok) map_card(parsed_card$value, ctx) else NULL,
  dragdb = if (needs_processing[["dragdb"]] && parsed_dragdb$ok) map_dragdb(parsed_dragdb$value, ctx) else NULL,
  yang = if (needs_processing[["yang"]]) process_yang(ctx) else NULL,
  pointfinder = if (needs_processing[["pointfinder"]] && parsed_pointfinder$ok) map_pointfinder(parsed_pointfinder$value, ctx) else NULL,
  refgenes = if (needs_processing[["refgenes"]] && parsed_refgenes$ok) map_refgenes(parsed_refgenes$value, ctx) else NULL,
  rifreport = if (needs_processing[["rifreport"]] && identical(update_results$rifreport$status, "ok")) process_rifreport() else NULL
)

# ---- 8. Combine + rebuild references ----
frames <- lapply(names(RAW_PATHS), function(ds) {
  pr <- process_results[[ds]]
  if (!is.null(pr)) pr$data else load_existing_rds(sprintf("processed_data/%s.rds", ds))
})

all_data <- combine_all_data(frames)
build_references(all_data)

# ---- 9. Persist state for every dataset that succeeded this run ----
for (ds in names(RAW_PATHS)) {
  pr <- process_results[[ds]]
  if (!is.null(pr) && identical(pr$status, "ok")) {
    state[[ds]] <- list(hash = current_hashes[[ds]], last_success_time = Sys.time())
  }
}
save_pipeline_state(state)

# ---- 10. Build and print the summary ----
summary_rows <- lapply(names(RAW_PATHS), function(ds) {
  is_scripted <- ds %in% SCRIPTED_DATASETS
  update_res <- update_results[[ds]]
  process_res <- process_results[[ds]]

  if (is_scripted && !is.null(update_res) && !identical(update_res$status, "ok")) {
    status <- "failed"
    reason <- paste0(update_res$message, " (using last-known-good processed data, if any)")
  } else if (!needs_processing[[ds]]) {
    status <- if (is_scripted) "skipped" else "manual-unchanged"
    reason <- "raw data unchanged since last successful run"
  } else if (!is.null(process_res) && identical(process_res$status, "ok")) {
    status <- if (is_scripted) "updated" else "manual-updated"
    reason <- process_res$message
  } else {
    status <- "failed"
    reason <- if (!is.null(process_res)) process_res$message else "raw data unavailable or unparseable"
  }

  tibble::tibble(
    dataset = ds,
    status = status,
    reason = reason,
    rows = if (!is.null(process_res)) process_res$rows else NA_integer_,
    timestamp = Sys.time()
  )
})

summary_df <- dplyr::bind_rows(summary_rows)

message("=== RIFRMap data pipeline refresh finished ===")
print(summary_df, width = Inf)

readr::write_csv(summary_df, LOG_PATH, append = file.exists(LOG_PATH))
