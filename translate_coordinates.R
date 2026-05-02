map_to_ecoli <- function(data, rpoB_target, rpoB_ecoli) { 
  
  # get target coordinates relative to E. coli in nucleotide positions:
  target_coords <- getCoordinates(rpoB_target, rpoB_ecoli, aligOutput = FALSE)
  
  aa_map <- data.frame(
    AA_pos = sort(unique(data$AA_pos)),
    AA_pos_Ecoli = vapply(
      sort(unique(data$AA_pos)),
      function(pos) {
        translateCoordinate(
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
    filter(!is.na(AA_pos_Ecoli))|>
    distinct() 
  
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





