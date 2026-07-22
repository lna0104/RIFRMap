# update_all.R
#
# Stage 2 (release) of the two-command update workflow (see
# code/validate_all.R for Stage 1). Re-runs the exact same pipeline +
# validation as validate_all.R (run_pipeline_validation(),
# code/pipeline_run.R) - so this is self-contained and safe to run even if
# you never ran validate_all.R first, since it always re-validates before
# ever touching disk - and then:
#
#   - if the run is fully clean (zero processing failures AND zero
#     validation errors, anywhere - warnings never block): commits every
#     successfully-processed dataset to processed_data/, rebuilds
#     allData.rds/all_references.csv, advances
#     processed_data/pipeline_state.rds, and publishes the release report
#     (output/release_reports/) + processed_data/changelog.csv (read by the
#     Shiny app's Data tab's "What's new" panel);
#
#   - otherwise: stops immediately. This is whole-run, all-or-nothing - if
#     *anything* is wrong *anywhere*, nothing commits for *any* dataset,
#     even ones with zero problems of their own. Every dataset that would
#     have updated is relabelled "blocked" and logged to
#     processed_data/last_update_log.csv for the audit trail, but nothing
#     else happens - no processed_data/ writes, no pipeline_state.rds
#     update, no release report, no changelog entry. See
#     output/validation_reports/ for what's wrong, fix it, and re-run
#     (source("code/validate_all.R") first then re-run this script).
#
# Run with (working directory = repo root):
#   source("code/update_all.R")

library(tidyverse)
library(Biostrings)
library(openssl)

.pipeline_orchestrating <- TRUE

source("code/helpers.R")
source("code/update_data.R")
source("code/process_data.R")
source("code/process_refseqs.R")
source("code/update_report.R")
source("code/pipeline_run.R")

message("=== RIFRMap release run starting ===")
result <- run_pipeline_validation()

validation_report_path <- build_validation_report(result$validation_findings, result$summary_df, result$run_time)
message("Validation report written to ", validation_report_path)

if (!result$is_clean) {
  blocked_summary <- result$summary_df |>
    dplyr::mutate(
      reason = ifelse(
        status %in% c("updated", "manual-updated"),
        "Blocked: this run's validation found unresolved error(s)/failure(s) elsewhere - see the validation report.",
        reason
      ),
      status = ifelse(status %in% c("updated", "manual-updated"), "blocked", status)
    )
  print(blocked_summary, width = Inf)
  readr::write_csv(blocked_summary, LOG_PATH, append = file.exists(LOG_PATH))

  message(
    "=== RIFRMap release run STOPPED: ", result$n_failures, " processing failure(s) and ",
    result$n_errors, " validation error(s) found - nothing published this run. See ",
    validation_report_path, ". Fix the underlying issue(s) and re-run to publish. ==="
  )
} else {
  # ---- Commit phase: only ever reached when the whole run is clean ----
  for (ds in names(RAW_PATHS)) {
    pr <- result$process_results[[ds]]
    if (!is.null(pr) && identical(pr$status, "ok")) {
      saveRDS(pr$data, file = pr$path)
    }
  }

  frames <- lapply(names(RAW_PATHS), function(ds) {
    pr <- result$process_results[[ds]]
    if (!is.null(pr) && identical(pr$status, "ok")) pr$data else load_existing_rds(sprintf("processed_data/%s.rds", ds))
  })

  all_data <- combine_all_data(frames)
  build_references(all_data)

  state <- result$state
  for (ds in names(RAW_PATHS)) {
    pr <- result$process_results[[ds]]
    if (!is.null(pr) && identical(pr$status, "ok")) {
      state[[ds]] <- list(hash = result$current_hashes[[ds]], last_success_time = Sys.time())
    }
  }
  save_pipeline_state(state)

  print(result$summary_df, width = Inf)
  readr::write_csv(result$summary_df, LOG_PATH, append = file.exists(LOG_PATH))

  release_report_path <- build_release_report(result$summary_df, result$diffs, result$run_time)
  message("Release report written to ", release_report_path)
  append_changelog(result$summary_df, result$diffs, result$run_time)

  message("=== RIFRMap release run finished: published successfully ===")
}
