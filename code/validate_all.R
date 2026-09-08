# validate_all.R
#
# Stage 1 of the two-command update workflow: runs the full pipeline
# (download through processing, see run_pipeline_validation() in
# code/pipeline_run.R) in memory and validates the result, without touching
# processed_data/, processed_data/pipeline_state.rds,
# processed_data/last_update_log.csv, or processed_data/changelog.csv.
# Always safe to run repeatedly - after fixing something the validation, re-run
# this until it reports clean, then run source("code/update_all.R") to actually
# publish.
#
# Run with (working directory = repo root):
#   source("code/validate_all.R")

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

message("=== RIFRMap validation run starting ===")
result <- run_pipeline_validation()

validation_report_path <- build_validation_report(result$validation_findings, result$summary_df, result$run_time)
message("Validation report written to ", validation_report_path)

if (result$is_clean) {
  n_ready <- sum(result$summary_df$status %in% c("updated", "manual-updated"))
  message(
    "Validation clean - ", n_ready, " dataset(s) ready to publish. ",
    "Run source(\"code/update_all.R\") to release."
  )
} else {
  message(
    "Validation found ", result$n_failures, " processing failure(s) - see ",
    validation_report_path, ". Fix the underlying issue(s) and re-run ",
    "source(\"code/validate_all.R\") until clean before releasing. ",
    "Any warnings in the report are informational and do not block a release."
  )
}
