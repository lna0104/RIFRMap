library(tidyverse)
library(Biostrings)
library(ALJEbinf)
library(rentrez)
library(openssl)

source("helpers.R")

# ---- Reference discovery & auto-fetch ----
#
# Given a list of species observed across the 7 mutation datasets, fetches an
# rpoB reference sequence + metadata row for any species not already present
# in raw_data/rpoB_references.csv, and appends the result
# to raw_data/rpoB_references.csv/.fasta. Species already present are left
# untouched, preserving manually-curated rows
#
# For each missing species: search NCBI's assembly database for a reference
# genome (falling back to any latest RefSeq assembly if no genome is flagged
# "reference genome"), download it, and pull out the sequence whose header
# matches `gene_match_pattern`. Failures for one species (no assembly found,
# no rpoB match in the downloaded genome, network error, ...) are caught and
# reported without stopping the rest of the species list.
#
# `gene_match_pattern` is inherently approximate - RefSeq CDS header wording
# for a given gene varies by assembly - and may need tuning against real
# output rather than being correct on the first try.
discover_and_fetch_new_references <- function(species_list,
                                              gene = "rpoB",
                                              gene_match_pattern = paste0(
                                                "\\[gene=rpoB\\]",
                                                "|rpoB($|[^a-zA-Z])",
                                                "|DNA-directed RNA polymerase subunit beta($|[^'])",
                                                "|DNA-directed RNA polymerase subunit beta chain",
                                                "|RNA polymerase, beta subunit",
                                                "|DNA-directed RNA polymerase subunit beta/beta'"
                                              ),
                                              refs_path = "raw_data/rpoB_references.csv",
                                              fasta_path = "raw_data/rpoB_references.fasta",
                                              download_dir = "output/rpoB_references") {
  species_list <- species_list |>
    stringr::str_squish() |>
    unique()
  species_list <- species_list[!is.na(species_list) & nzchar(species_list)]

  existing_refs <- readr::read_csv(refs_path, show_col_types = FALSE)
  missing_species <- setdiff(species_list, existing_refs$Species)

  added_species <- character(0)
  failed_species <- character(0)

  if (length(missing_species) == 0) {
    message("No new species to fetch rpoB references for.")
    return(list(added_species = added_species, failed_species = failed_species))
  }

  message("Attempting to fetch rpoB references for ", length(missing_species), " new species.")

  new_rows <- list()
  new_sequences <- list()

  for (species in missing_species) {
    fetch_result <- tryCatch(
      {
        term_ser <- sprintf(
          '"%s"[ORGN] AND ("latest refseq"[filter] AND "reference genome"[filter])',
          species
        )

        assembly_results <- get_summaries(db = "assembly", term = term_ser)

        if (length(assembly_results) != 0) {
          assembly_final <- assembly_results[[1]]
          ref_strain <- TRUE
        } else {
          term_ser <- sprintf('"%s"[ORGN] AND "latest refseq"[filter]', species)
          assembly_results <- get_summaries(db = "assembly", term = term_ser)
          if (length(assembly_results) == 0) {
            stop("No assembly found in NCBI for this species")
          }
          assembly_final <- assembly_results[[1]]
          ref_strain <- FALSE
        }

        link_to_nuc <- entrez_link(dbfrom = "assembly", db = "nuccore", id = assembly_final$uid)
        n_link_ref <- length(link_to_nuc$links$assembly_nuccore_refseq)
        n_link <- length(link_to_nuc$links$assembly_nuccore)

        if (n_link == 0) {
          stop("No nucleotide records linked to this assembly")
        }

        nuc_ids <- if (n_link_ref >= 1) link_to_nuc$links$assembly_nuccore_refseq else link_to_nuc$links$assembly_nuccore
        nuc_summary <- entrez_summary(db = "nuccore", id = nuc_ids)
        nuc_first <- if (length(nuc_ids) > 1) nuc_summary[[1]] else nuc_summary

        strain <- nuc_first$strain
        ref_seq <- nuc_first$accessionversion

        species_us <- gsub(" ", "_", species)
        strain_us <- gsub(" ", "_", strain)
        fasta_name <- paste0(gene, "_", species_us, "_", strain_us)

        genome_path <- download_file(assembly_final, dir = download_dir)
        genome_fasta <- sub("\\.gz$", "", genome_path)
        sequences <- readDNAStringSet(genome_fasta)
        matching_seq <- sequences[grepl(gene_match_pattern, names(sequences), ignore.case = TRUE)]

        if (length(matching_seq) == 0) {
          stop(paste0("No ", gene, " sequence found in downloaded genome"))
        } else if (length(matching_seq) > 1) {
          matching_seq <- matching_seq[1]
        }
        names(matching_seq) <- fasta_name

        list(
          row = tibble::tibble(
            Species = species,
            Strain = strain,
            RefStrain = ref_strain,
            Gene = gene,
            Assembly_ID = assembly_final$assemblyname,
            FASTA_name = fasta_name,
            `NCBI Reference sequence` = ref_seq,
            Origen = "Auto-fetched via NCBI assembly search"
          ),
          seq = matching_seq
        )
      },
      error = function(e) {
        message("Failed to fetch rpoB reference for ", species, ": ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(fetch_result)) {
      failed_species <- c(failed_species, species)
    } else {
      new_rows[[species]] <- fetch_result$row
      new_sequences[[species]] <- fetch_result$seq
      added_species <- c(added_species, species)
    }
  }

  if (length(new_rows) > 0) {
    combined_new_rows <- dplyr::bind_rows(new_rows)
    updated_refs <- dplyr::bind_rows(existing_refs, combined_new_rows)
    readr::write_csv(updated_refs, refs_path)

    existing_seqs <- readDNAStringSet(fasta_path)
    combined_new_seqs <- do.call(c, unname(new_sequences))
    all_seqs <- c(existing_seqs, combined_new_seqs)
    writeXStringSet(all_seqs, fasta_path)
  }

  list(added_species = added_species, failed_species = failed_species)
}

# ---- Coordinate regeneration ----
#
# Regenerates the per-reference-sequence coordinate map used to translate
# amino-acid positions between each species' native numbering and the
# canonical E. coli numbering. Cached/skipped based on a SHA1 hash of
# raw_data/rpoB_references.fasta - call this *after*
# discover_and_fetch_new_references() so a newly-expanded fasta is picked up
# (its changed hash correctly triggers regeneration).
regenerate_coordinates_if_needed <- function(fasta_path = "raw_data/rpoB_references.fasta",
                                             hash_path = "processed_data/fastahash.Rds",
                                             coords_path = "processed_data/coordinates.RData") {
  seqs <- readDNAStringSet(fasta_path)
  current_hash <- as.character(openssl::sha1(file(fasta_path)))

  if (file.exists(hash_path) && file.exists(coords_path) &&
    current_hash == readRDS(hash_path)) {
    message("Sequences file has not changed, loading original coordinates")
    load(coords_path)
  } else {
    message("Sequences file has changed, regenerating coordinates")
    saveRDS(current_hash, file = hash_path)
    coordinates <- ALJEbinf::getAllCoordinates(seqs, "rpoB_Escherichia_coli_MG1655")
    save(coordinates, file = coords_path)
  }

  coordinates
}

# ---- Standalone execution ----
# When this file is sourced directly (not via update_all.R), just regenerate
# coordinates if the reference fasta changed - the same behavior as before
# this file gained the reference-discovery function above. Reference
# discovery itself needs a species_list assembled from all 7 datasets, which
# only update_all.R has; it is not run automatically here.
if (!exists(".pipeline_orchestrating", inherits = FALSE)) {
  coordinates <- regenerate_coordinates_if_needed()
}
