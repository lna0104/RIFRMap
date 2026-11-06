#' Retry an expression until it succeeds
#'
#' Evaluates `expr`, re-evaluating it after any error, up to `max_tries` times
#' in total. The final attempt runs without error handling, so its error
#' propagates to the caller.
#'
#' @param expr the expression to evaluate.
#' @param max_tries the maximum number of attempts.
#' @param envir the environment in which to run the expression.
#'
#' @returns the return value of the expression on its successful execution
#' @export
#'
#' @examples
#' retry(stopifnot(runif(1) > 0.95))
retry <- function(expr, max_tries = Inf, envir = parent.frame()) {
  expr <- substitute(expr)
  while(max_tries > 1) {
    result <- tryCatch(eval(expr, envir), error = function(e) {return(e)})
    if(!is(result, "error")) {
      return(result)
    } else {
      cat(sprintf("Error in %s : %s\n", deparse(result$call), result$message), file = stderr())
      cat("Retrying...\n", file = stderr())
      max_tries <- max_tries - 1
    }
  }
  #If this is our last try, just run without error handling.
  eval(expr, envir)
}

#' Download summaries from an Entrez database search
#'
#' Searches an Entrez database for a given search term, downloading the summaries for each result.
#'
#' @param db a character string giving the name of the Entrez database to search.
#' @param term a character string giving the search term to pass to the database.
#' @param retmax the maximum number of summaries to request from the database at a time.
#' @param max_tries the maximum number of times to try each web request before throwing an error.
#'
#' @details
#' `retmax` limits the number of summaries requested in each web request to the database.
#' `get_summaries()` will make multiple requests, each of size not exceeding `retmax`, until all available summaries have been downloaded.
#' Each web request will be attempted `max_tries` times to protect against intermittent internet failure.
#'
#' @returns a list of esummary records
#' @export
#'
#' @examples
#' get_summaries('assembly', 'GCA_000006945.2[Assembly Accession] OR GCA_016028495.1[Assembly Accession] AND (latest[filter] AND all[filter] NOT anomalous[filter])')
get_summaries <- function(db, term, retmax = 500, max_tries = 10) {
  #Firstly, we perform a search, but set retmax = 0 to have it return no entries.
  #This just loads entries into the history server.
  #We will then use web_history to get the summaries.
  search_result <- rentrez::entrez_search(db = db, term = term, retmax = 0, use_history = TRUE)
  #Entrez limits the number of entries it will return in a single search.
  #As such, we need to make repeated requests to get all the summaries.
  summaries <- list()
  while(length(summaries) < search_result$count) {
    #The web request may fail, and we don't want to lose all our progress so far if it does, so use retry() for it.
    these_summaries <- retry(
      rentrez::entrez_summary(db = db, web_history = search_result$web_history, retstart = length(summaries), retmax = retmax, always_return_list = TRUE),
      max_tries = max_tries
    )
    summaries <- append(summaries, these_summaries)
  }
  return(summaries)
}

# extracts the file path to a genome from a summary:
get_genome_file_path <- function(summary, suffix, dir) {
  name <- gsub(".*/([^/]+)/?", "\\1", summary$ftppath_refseq)
  file_name <- paste0(name, suffix)
  file_path <- file.path(dir, file_name)
  return(file_path)
}


#Takes an Entrez summary, and downloads and unzips the fasta file from refseq.
#Each summary object must have an `ftppath_refseq` entry.
#Downloaded files are placed in the `dir` directory, which is created if it does not exist.
#' Download Fasta File From RefSeq
#'
#' Given an esummary object, downloads (and unzips if necessary) a corresponding fasta files of a given type from RefSeq.
#'
#' @param summary an esummary object with an `ftppath_refseq` entry.
#' @param suffix a character string giving the suffix of the file to be downloaded (following the assembly name).
#' @param dir a character string giving the directory to place the downloaded file in (created if it does not exist).
#' @param max_tries the maximum number of times to try each web request before throwing an error.
#'
#' @details
#' The URL from which to download the file is constructed from `summary`'s `ftppath_refseq` entry and the given `suffix`.
#' `suffix` must include the file extension, including `.gz` if applicable.
#' Each web request will be attempted `max_tries` times to protect against intermittent internet failure.
#' The downloaded file is unzipped if it ends in `gz`, and the original `.gz` files is removed.
#'
#' @returns a character string giving the path to the downloaded (and unzipped if applicable) file
#' @export
#'
#' @examples
#' withr::with_tempfile("td", {
#'   summaries <- get_summaries('assembly', 'GCA_000006945.2[Assembly Accession] OR GCA_016028495.1[Assembly Accession] AND (latest[filter] AND all[filter] NOT anomalous[filter])')
#'   files <- list()
#'   for(i in 1:length(summaries)) {
#'     try( #Continue if any download fails
#'       files[[i]] <- download_file(summaries[[i]], dir = td)
#'     )
#'   }
#' })
download_file <- function(summary, suffix = "_cds_from_genomic.fna.gz", dir = "output/genomes", max_tries = 10, overwrite = FALSE) {
  #Create the output directory if it does not exist.
  if(!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  file_path <- get_genome_file_path(summary, suffix = suffix, dir = dir)
  file_name <- basename(file_path)

  if ((!overwrite) && (file.exists(sub("\\.gz$", "", file_path)))) {
    cat("  --> File already exists, skipped.\n")
  } else {
    url <- paste0(summary$ftppath_refseq, "/", file_name)
    url <- gsub("^ftp:", "https:", url)
    print(url)
    retry(
      download.file(url, file_path),
      max_tries = max_tries
    )
    if(grepl("gz$", suffix)) {
      file_path <- R.utils::gunzip(file_path, overwrite = TRUE)
    }
  }
  return(file_path)
}

map_to_ecoli <- function(data, rpoB_target, rpoB_ecoli) { 
  
  # get target coordinates relative to E. coli in nucleotide positions:
  target_coords <- ALJEbinf::getCoordinates(rpoB_target, rpoB_ecoli, aligOutput = FALSE)

  positions <- sort(unique(data$AA_pos))
  aa_map <- data.frame(
    AA_pos = positions,
    AA_pos_Ecoli = vapply(
      positions,
      function(pos) {
        ALJEbinf::translateCoordinate(
          pos = pos,
          coordinates = target_coords,
          direction = "FocalToRef",
          AAinput = TRUE,
          AAoutput = TRUE
        )
      },
      integer(1)
    )
  ) |>
    dplyr::filter(!is.na(AA_pos_Ecoli))
  
  mapped_muts <- data |>
    dplyr::left_join(aa_map, by = "AA_pos") |>
    dplyr::filter(!is.na(AA_pos_Ecoli)) |>
    dplyr::distinct()
  
  return(mapped_muts)
}

get_rpoB_for_species <- function(species, seqs) {
  nm <- names(seqs)
  
  hit <- grep(species, nm, ignore.case = TRUE)
  
  if (length(hit) == 0) {
    message("No rpoB sequence found for species: ", species)
    return(NULL)
  }
  
  if (length(hit) > 1) {
    message(
      "Multiple rpoB sequences found for ", species,
      "; using first match: ", nm[hit[1]]
    )
  }
  
  seqs[[hit[1]]]
}
