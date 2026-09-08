# pipeline_run.R
#
# Shared core for the two pipeline entry points, code/validate_all.R and
# code/update_all.R: run_pipeline_validation() does everything through
# validation - download, parse, discover/fetch missing rpoB references,
# regenerate coordinates, finish RifReport, decide what needs reprocessing,
# process each dataset entirely in memory, and run every validation check
# (run_validation_checks(), update_report.R) against the result. It commits
# nothing to processed_data/, pipeline_state.rds, last_update_log.csv, or
# changelog.csv - reference discovery and coordinate regeneration are the
# only writes here, and those are deliberately ungated.
#
# Both entry points call this, then decide separately what to do with the
# result: validate_all.R always stops after reporting; update_all.R commits
# only if the whole run came back clean (see its header for the exact
# all-or-nothing rule).
#
# Sourced by validate_all.R/update_all.R only

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

# Runs download -> parse -> reference discovery -> coordinates -> RifReport
# finalize -> change detection -> per-dataset processing -> validation, and
# returns everything the caller needs to decide what to do next. See the
# file header for exactly what this does and doesn't touch on disk.
run_pipeline_validation <- function() {
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
  # Deliberately NOT gated - append-only. Problems here surface via
  # check_refstrain_gaps() in the validation report instead.
  discovery <- discover_and_fetch_new_references(species_pool)
  if (length(discovery$added_species) > 0) {
    message("Added rpoB references for: ", paste(discovery$added_species, collapse = ", "))
  }
  if (length(discovery$failed_species) > 0) {
    message("Could not find rpoB references for: ", paste(discovery$failed_species, collapse = ", "))
  }

  # ---- 4. Regenerate coordinates (picks up any newly-appended sequences) ----
  coordinates <- regenerate_coordinates_if_needed()
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

  # Snapshot each about-to-be-reprocessed dataset's current .rds *before*
  # processing (nothing overwrites it here - only the caller's commit phase,
  # if any, does) - this is the only point old data is available for
  # update_report.R's diff_dataset() to compare against.
  old_data_snapshot <- setNames(
    lapply(names(RAW_PATHS), function(ds) {
      if (needs_processing[[ds]]) load_existing_rds(sprintf("processed_data/%s.rds", ds)) else NULL
    }),
    names(RAW_PATHS)
  )

  # ---- 7. Process (or skip) each dataset - held only in memory ----
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

  # ---- 8. Validate everything processed this run ----
  pre_species <- list(
    card = if (parsed_card$ok) parsed_card$value$species else NULL,
    dragdb = if (parsed_dragdb$ok) parsed_dragdb$value$species else NULL,
    pointfinder = if (parsed_pointfinder$ok) parsed_pointfinder$value$species else NULL,
    refgenes = if (parsed_refgenes$ok) parsed_refgenes$value$species else NULL,
    rifreport = if (identical(rifreport_dl$status, "ok")) rifreport_dl$species else NULL
  )

  validation_findings <- run_validation_checks(process_results, pre_species, refs, coordinates, needs_processing)

  # Diff every dataset that processed successfully against its pre-processing
  # snapshot (see step 6) - used by build_release_report()/append_changelog()
  # if and when the caller actually commits.
  diffs <- setNames(
    lapply(names(RAW_PATHS), function(ds) {
      pr <- process_results[[ds]]
      if (!is.null(pr) && identical(pr$status, "ok")) {
        diff_dataset(ds, old_data_snapshot[[ds]], pr$data)
      } else {
        NULL
      }
    }),
    names(RAW_PATHS)
  )

  # ---- 9. Build the summary ----
  # Status here describes what *would* happen if this dataset's processing
  # is considered in isolation - it's the caller's job (update_all.R) to
  # decide whether the whole run is clean enough to actually apply it.
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
  run_time <- Sys.time()

  # Validation is report-only: every check emits "warning", and findings never
  # gate a release. The single gating condition is a source that CRASHED and so
  # produced no data at all - that is not a validation opinion, it is an absent
  # dataset, and publishing around it would silently serve stale data for that
  # source. See build_validation_report() for what the report shows.
  n_failures <- sum(summary_df$status == "failed")

  message("=== RIFRMap pipeline run finished (validation only - nothing committed yet) ===")
  print(summary_df, width = Inf)

  list(
    update_results = update_results,
    process_results = process_results,
    validation_findings = validation_findings,
    summary_df = summary_df,
    diffs = diffs,
    state = state,
    current_hashes = current_hashes,
    needs_processing = needs_processing,
    run_time = run_time,
    n_failures = n_failures,
    is_clean = (n_failures == 0)
  )
}
