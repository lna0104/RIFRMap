# update_report.R
#
# Implements the reporting/validation half of the two-command update
# workflow split across code/validate_all.R (Stage 1) and code/update_all.R
# (Stage 2) - see code/pipeline_run.R's run_pipeline_validation(), which
# both call:
#
#   Stage 1 (validation, internal): run_validation_checks() inspects each
#   freshly-processed-but-not-yet-persisted dataset (still only in memory -
#   see process_data.R) plus the reference catalog as a whole, and returns a
#   tibble of findings. Every finding is a "warning": validation reports, it
#   does not adjudicate, and no finding has ever blocked a release since the
#   error severity was retired. build_validation_report() always writes this to
#   output/validation_reports/, even when there's nothing to report - this
#   file is never read by the Shiny app. validate_all.R always stops after
#   writing it; update_all.R writes the same report before deciding whether
#   to proceed.
#
#   Stage 2 (release, public, whole-run gate): update_all.R only proceeds
#   past its own validation_report call if the *entire* run is clean - which
#   now means exactly one thing: zero processing failures, anywhere. A crashed
#   source produced no data, so publishing around it would silently serve that
#   source's stale data. Findings themselves never block. If clean,
#   diff_dataset() compares each committed dataset's new data against a
#   pre-processing snapshot of the old data (`old_data_snapshot`, built in
#   run_pipeline_validation()), and build_release_report() +
#   append_changelog() turn that into a detailed report
#   (output/release_reports/) and a durable, git-tracked one-line-per-
#   dataset changelog (processed_data/changelog.csv) that the Shiny app's
#   Data tab reads as its "what's new" summary. If anything in the run
#   wasn't clean, update_all.R stops before ever calling either of those -
#   nothing commits for *any* dataset, even ones with zero problems of their
#   own - so only a fully-validated run is ever published.
#
# Sourced by validate_all.R/update_all.R only - diffing needs a before/after
# snapshot that only run_pipeline_validation() (code/pipeline_run.R)
# captures, so unlike update_data.R/process_data.R/process_refseqs.R this
# file has no standalone driver block.

DATASET_DISPLAY_NAMES <- c(
  card        = "CARD",
  cryptic     = "CRyPTIC",
  dragdb      = "DRAGdb",
  yang        = "Yang-RIF",
  rifreport   = "RifReport",
  pointfinder = "PointFinder",
  refgenes    = "Reference Gene Catalog"
)


DIFF_KEY_COLS <- c("species", "AA_original", "AA_pos_Ecoli", "AA_mutation")

# Row-wise "did this value change" that tolerates NA on either side and
# compares via character representation so type drift between runs doesn't
# register as a spurious change.
.values_differ <- function(a, b) {
  a_na <- is.na(a)
  b_na <- is.na(b)
  ifelse(a_na & b_na, FALSE, ifelse(a_na | b_na, TRUE, as.character(a) != as.character(b)))
}

.mutation_label <- function(df) {
  paste0(df$species, ": ", df$AA_original, df$AA_pos_Ecoli, df$AA_mutation)
}

.empty_findings <- function() {
  tibble::tibble(
    dataset = character(0), severity = character(0),
    subject = character(0), reason = character(0)
  )
}

# Diffs one dataset's new processed data frame against its pre-processing
# snapshot (NULL if the dataset never existed before). Returns a list
# describing what changed; used by both build_release_report() (verbose) and
# append_changelog() (one-liner).
diff_dataset <- function(dataset, old_df, new_df) {
  base <- list(
    dataset = dataset,
    is_new_dataset = is.null(old_df),
    new_species = character(0),
    removed_species = character(0),
    n_added = 0L,
    n_removed = 0L,
    n_modified = 0L,
    modified_cols = character(0),
    added_examples = character(0),
    removed_examples = character(0)
  )

  if (is.null(old_df)) {
    base$new_species <- sort(unique(new_df$species))
    base$n_added <- nrow(new_df)
    base$added_examples <- head(.mutation_label(new_df), 10)
    return(base)
  }

  if (!all(DIFF_KEY_COLS %in% names(old_df)) || !all(DIFF_KEY_COLS %in% names(new_df))) {
    # Defensive fallback: schema doesn't support key-level diffing this run -
    # still surface a row-count-level change rather than erroring the run.
    base$n_added <- max(nrow(new_df) - nrow(old_df), 0L)
    base$n_removed <- max(nrow(old_df) - nrow(new_df), 0L)
    return(base)
  }

  old_df <- dplyr::mutate(old_df, dplyr::across(dplyr::all_of(DIFF_KEY_COLS), as.character))
  new_df <- dplyr::mutate(new_df, dplyr::across(dplyr::all_of(DIFF_KEY_COLS), as.character))

  base$new_species <- sort(setdiff(unique(new_df$species), unique(old_df$species)))
  base$removed_species <- sort(setdiff(unique(old_df$species), unique(new_df$species)))

  added <- dplyr::anti_join(new_df, old_df, by = DIFF_KEY_COLS)
  removed <- dplyr::anti_join(old_df, new_df, by = DIFF_KEY_COLS)

  base$n_added <- nrow(added)
  base$n_removed <- nrow(removed)
  base$added_examples <- head(.mutation_label(added), 10)
  base$removed_examples <- head(.mutation_label(removed), 10)

  compare_cols <- setdiff(intersect(names(old_df), names(new_df)), DIFF_KEY_COLS)

  if (length(compare_cols) > 0) {
    matched <- dplyr::inner_join(new_df, old_df, by = DIFF_KEY_COLS, suffix = c(".new", ".old"))

    if (nrow(matched) > 0) {
      diff_matrix <- vapply(compare_cols, function(col) {
        .values_differ(matched[[paste0(col, ".new")]], matched[[paste0(col, ".old")]])
      }, logical(nrow(matched)))

      # vapply drops to a bare vector (no dim) when there's exactly one
      # matched row - reshape so rowSums()/colSums() below stay valid.
      if (is.null(dim(diff_matrix))) {
        diff_matrix <- matrix(diff_matrix, nrow = nrow(matched), dimnames = list(NULL, compare_cols))
      }

      row_changed <- rowSums(diff_matrix) > 0
      base$n_modified <- sum(row_changed)
      base$modified_cols <- compare_cols[colSums(diff_matrix[row_changed, , drop = FALSE]) > 0]
    }
  }

  base
}

# One line summarising a dataset's diff, used both in the
# detailed report and the Data-tab changelog.
.diff_summary_line <- function(diff) {
  if (isTRUE(diff$is_new_dataset)) {
    return(sprintf(
      "New dataset: %d record(s) across %d species.",
      diff$n_added, length(diff$new_species)
    ))
  }

  parts <- character(0)

  if (diff$n_added > 0) {
    species_note <- if (length(diff$new_species) > 0) {
      sprintf(" (%d new species: %s)", length(diff$new_species), paste(diff$new_species, collapse = ", "))
    } else {
      ""
    }
    parts <- c(parts, sprintf("%d new mutation record(s) added%s", diff$n_added, species_note))
  }

  if (diff$n_removed > 0) {
    species_note <- if (length(diff$removed_species) > 0) {
      sprintf(" (%d species dropped: %s)", length(diff$removed_species), paste(diff$removed_species, collapse = ", "))
    } else {
      ""
    }
    parts <- c(parts, sprintf("%d mutation record(s) removed%s", diff$n_removed, species_note))
  }

  if (diff$n_modified > 0) {
    cols_note <- if (length(diff$modified_cols) > 0) {
      sprintf(" (%s)", paste(diff$modified_cols, collapse = ", "))
    } else {
      ""
    }
    parts <- c(parts, sprintf("%d record(s) updated%s", diff$n_modified, cols_note))
  }

  if (length(parts) == 0) {
    return("Reprocessed with no net record-level changes detected.")
  }

  paste0(paste(parts, collapse = "; "), ".")
}

# ---- Stage 1: validation checks ----
#
# Each check_*() returns a tibble(dataset, severity, subject, reason).
# `severity` is always "warning" - findings are informational and never block
# a commit. The "error" severity was retired deliberately: because the release
# gate is whole-run, a single over-eager error blocked EVERY dataset from ever
# publishing (an upstream citation gap did exactly that, indefinitely).

# Flags every raw_data/rpoB_references.csv row that isn't a designated
# RefStrain. Always a warning, never gates a commit - even the case with no
# coordinates.RData entry (the failure mode that crashes
# ALJEbinf::fillMutationsTable()) is only ever informational here, so the
# reason is legible in the report instead of a raw dplyr stack trace; the
# crash itself, if it happens, still surfaces separately as a processing
# failure (see run_pipeline_validation()'s `n_failures`, code/pipeline_run.R).
check_refstrain_gaps <- function(refs, coordinates) {
  non_ref <- refs |> dplyr::filter(is.na(RefStrain) | RefStrain != TRUE)
  if (nrow(non_ref) == 0) {
    return(.empty_findings())
  }

  has_coords <- non_ref$FASTA_name %in% names(coordinates)

  tibble::tibble(
    dataset = NA_character_,
    severity = "warning",
    subject = non_ref$FASTA_name,
    reason = ifelse(
      has_coords,
      paste0(
        "Species '", non_ref$Species, "' reference sequence is not marked as a RefStrain, ",
        "but does have a coordinate mapping - usable, just not the designated best reference."
      ),
      paste0(
        "Species '", non_ref$Species, "' reference sequence is not marked as a RefStrain and ",
        "has no entry in processed_data/coordinates.RData - any mutation record resolving to ",
        "this sequence cannot be mapped to E. coli coordinates and will be excluded, or can ",
        "crash mapping outright."
      )
    )
  )
}

# Species present in a source's pre-mapping parsed input but absent from its
# final mapped output. map_to_ecoli()/fillMutationsTable() already silently
# drop rows that can't be translated to E. coli coordinates -, not a bug,
# so this stays a warning rather than gating the dataset;
# it just makes the drop visible instead of silent.
check_dropped_species <- function(dataset, pre_species, post_species) {
  dropped <- setdiff(pre_species, post_species)
  dropped <- dropped[!is.na(dropped) & nzchar(dropped)]
  if (length(dropped) == 0) {
    return(.empty_findings())
  }

  tibble::tibble(
    dataset = dataset,
    severity = "warning",
    subject = dropped,
    reason = paste0(
      "Present in the source's parsed input but no rows survived mapping to E. coli ",
      "coordinates (no usable reference sequence, or all rows failed coordinate translation)."
    )
  )
}

# Rows with incomplete citation metadata. A complete record requires at
# least one of DOI or PMID - either one alone is enough. A column entirely
# absent from this dataset counts as missing for every row.
check_missing_metadata <- function(dataset, data) {
  if (is.null(data) || nrow(data) == 0) {
    return(.empty_findings())
  }

  if (!("DOI" %in% names(data)) && !("PMID" %in% names(data))) {
    return(.empty_findings())
  }

  is_missing_col <- function(col) {
    if (!(col %in% names(data))) {
      return(rep(TRUE, nrow(data)))
    }
    vals <- data[[col]]
    is.na(vals) | !nzchar(trimws(as.character(vals)))
  }

  incomplete <- is_missing_col("DOI") & is_missing_col("PMID")

  n_missing <- sum(incomplete)
  if (n_missing == 0) {
    return(.empty_findings())
  }

  tibble::tibble(
    dataset = dataset,
    severity = "warning",
    subject = sprintf("%d row(s)", n_missing),
    reason = sprintf("%d of %d row(s) are missing both DOI and PMID.", n_missing, nrow(data))
  )
}

# Orchestrates every check above across all 7 datasets into one findings
# tibble. `pre_species` is a named list (dataset -> character vector of
# species seen before mapping/reference-filtering); datasets with no
# meaningful "before" species list (cryptic, yang - fixed single-species,
# no discovery step) are simply absent from it. Per-dataset checks only run
# for datasets reprocessed this run (`needs_processing`) that finished with
# status "ok" - a dataset that errored out has no `data` to check anyway.
run_validation_checks <- function(process_results, pre_species, refs, coordinates, needs_processing) {
  findings <- list(check_refstrain_gaps(refs, coordinates))

  for (ds in names(process_results)) {
    pr <- process_results[[ds]]
    if (is.null(pr) || !identical(pr$status, "ok") || !isTRUE(needs_processing[[ds]])) {
      next
    }

    findings[[length(findings) + 1]] <- check_missing_metadata(ds, pr$data)

    if (!is.null(pre_species[[ds]])) {
      post <- if (!is.null(pr$data) && "species" %in% names(pr$data)) unique(pr$data$species) else character(0)
      findings[[length(findings) + 1]] <- check_dropped_species(ds, pre_species[[ds]], post)
    }

    if (identical(ds, "rifreport") && !is.null(pr$excluded_rows) && nrow(pr$excluded_rows) > 0) {
      ex <- pr$excluded_rows
      findings[[length(findings) + 1]] <- tibble::tibble(
        dataset = "rifreport",
        severity = "warning",
        subject = paste0(ex$Species, ifelse(is.na(ex$Strain), "", paste0(" (", ex$Strain, ")"))),
        reason = paste0("Excluded by fillMutationsTable(): ", dplyr::coalesce(ex$Warning, "unresolved"))
      )
    }
  }

  dplyr::bind_rows(findings)
}

# Builds the Stage 1 (internal) validation report and writes it to
# output/validation_reports/, returning the file path. Always called, always
# written - even when clean, so there's a durable "nothing was wrong this
# run" record. Never feeds processed_data/changelog.csv or the Data tab.
build_validation_report <- function(validation_findings, summary_df, run_time,
                                    out_dir = "output/validation_reports") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  display <- function(ds) {
    if (is.na(ds) || !ds %in% names(DATASET_DISPLAY_NAMES)) {
      "Reference catalog (all datasets)"
    } else {
      DATASET_DISPLAY_NAMES[[ds]]
    }
  }

  # Findings are warnings only - the report informs, it does not adjudicate.
  # A crashed dataset is reported separately below, since that is a missing
  # dataset rather than a finding about one.
  failed <- summary_df |> dplyr::filter(status == "failed")
  warnings_ <- validation_findings |> dplyr::filter(severity == "warning")

  sections <- character(0)

  if (nrow(failed) > 0) {
    lines <- sprintf("    - %s: %s", vapply(failed$dataset, display, character(1)), failed$reason)
    sections <- c(sections, paste0("  PROCESSING FAILURES\n", paste(lines, collapse = "\n"), "\n"))
  }

  if (nrow(warnings_) > 0) {
    lines <- sprintf(
      "    - [%s] %s\n        %s",
      vapply(warnings_$dataset, display, character(1)), warnings_$subject, warnings_$reason
    )
    sections <- c(sections, paste0(
      "  WARNINGS (informational, does not block)\n", paste(lines, collapse = "\n"), "\n"
    ))
  }

  body <- if (length(sections) == 0) "    No issues found this run.\n" else paste(sections, collapse = "\n")

  report <- paste0(
    "RIFRMap - Validation Report (Stage 1, internal)\n",
    "Run date: ", format(run_time), "\n",
    "This report is never published to the app's Data page - see the release\n",
    "report for what actually got published.\n",
    "--------------------------------------------------------------------\n\n",
    body, "\n",
    "--------------------------------------------------------------------\n"
  )

  path <- file.path(out_dir, sprintf("validation_report_%s.txt", format(run_time, "%Y%m%d_%H%M%S")))
  writeLines(report, path)
  path
}

# ---- Stage 2: release report ----
#
# Builds the full per-run text report and writes it to output/release_reports/,
# returning the file path. Only ever called from update_all.R's clean branch
# (see run_pipeline_validation()'s `is_clean`, code/pipeline_run.R)
build_release_report <- function(summary_df, diffs, run_time,
                                 out_dir = "output/release_reports") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  display <- function(ds) DATASET_DISPLAY_NAMES[[ds]]

  unchanged <- summary_df |> dplyr::filter(status %in% c("skipped", "manual-unchanged"))
  updated <- summary_df |> dplyr::filter(status %in% c("updated", "manual-updated"))
  failed <- summary_df |> dplyr::filter(status == "failed")

  unchanged_text <- if (nrow(unchanged) == 0) {
    "    None.\n"
  } else {
    paste0("    - ", vapply(unchanged$dataset, display, character(1)), "\n", collapse = "")
  }

  updated_text <- if (nrow(updated) == 0) {
    "    None.\n"
  } else {
    paste(vapply(updated$dataset, function(ds) {
      diff <- diffs[[ds]]
      lines <- sprintf("    [%s]", display(ds))
      lines <- c(lines, sprintf("        Rows: %s", summary_df$rows[summary_df$dataset == ds]))

      if (!is.null(diff)) {
        lines <- c(lines, sprintf("        %s", .diff_summary_line(diff)))

        if (length(diff$added_examples) > 0) {
          lines <- c(lines, "        Added (examples):", paste0("            ", diff$added_examples))
        }
        if (length(diff$removed_examples) > 0) {
          lines <- c(lines, "        Removed (examples):", paste0("            ", diff$removed_examples))
        }
      }

      paste(lines, collapse = "\n")
    }, character(1)), collapse = "\n\n")
  }

  failed_text <- if (nrow(failed) == 0) {
    "    None.\n"
  } else {
    paste0("    - ", vapply(failed$dataset, display, character(1)), ": ", failed$reason, "\n", collapse = "")
  }

  report <- paste0(
    "RIFRMap - Release Report (Stage 2, public)\n",
    "Run date: ", format(run_time), "\n",
    "--------------------------------------------------------------------\n\n",
    "UNCHANGED DATASETS\n", unchanged_text, "\n",
    "--------------------------------------------------------------------\n\n",
    "UPDATED DATASETS (published to processed_data/changelog.csv)\n\n", updated_text, "\n\n",
    "--------------------------------------------------------------------\n\n",
    "FAILED DATASETS\n", failed_text, "\n",
    "--------------------------------------------------------------------\n"
  )

  path <- file.path(out_dir, sprintf("release_report_%s.txt", format(run_time, "%Y%m%d_%H%M%S")))
  writeLines(report, path)
  path
}

# Appends one row per changed dataset to processed_data/changelog.csv - the
# durable, git-tracked "what's new" history the Data tab reads. Only called
# when at least one dataset actually changed this run (see update_all.R), so
# "latest release" always means the last run with real data changes, not
# just the last time the pipeline happened to execute.
append_changelog <- function(summary_df, diffs, run_time,
                             path = "processed_data/changelog.csv") {
  changed <- summary_df |> dplyr::filter(status %in% c("updated", "manual-updated"))
  if (nrow(changed) == 0) {
    return(invisible(NULL))
  }

  entries <- tibble::tibble(
    run_date = as.character(as.Date(run_time)),
    dataset = vapply(changed$dataset, function(ds) DATASET_DISPLAY_NAMES[[ds]], character(1)),
    status = changed$status,
    summary = vapply(changed$dataset, function(ds) {
      diff <- diffs[[ds]]
      if (is.null(diff)) "Reprocessed." else .diff_summary_line(diff)
    }, character(1))
  )

  readr::write_csv(entries, path, append = file.exists(path))
  invisible(path)
}
