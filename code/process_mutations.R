map_to_ecoli <- function(data, rpoB_target, rpoB_Ecoli) { 
  
  # get target coordinates relative to E. coli in nucleotide positions:
  target_coords <- getCoordinates(rpoB_target, rpoB_Ecoli, aligOutput = FALSE)
  
  # get target coordinates relative to E. coli in AA positions:
  target_coords_AA <- data.frame(
    AA_pos = sapply(target_coords$posFocal, function(pos) {
      translateCoordinate(pos, target_coords, direction = "FocalToRef", AAinput = FALSE, AAoutput = TRUE)
    }),
    AA_pos_Ecoli = sapply(target_coords$posRef, function(pos) {
      translateCoordinate(pos, target_coords, direction = "FocalToRef", AAinput = FALSE, AAoutput = TRUE)
    })
  ) |>
    distinct() |>
    filter(!is.na(AA_pos))
  
  mapped_muts <- data |>
    left_join(target_coords_AA, by = join_by(AA_pos)) |>
    filter(!is.na(AA_pos_Ecoli))
  
  return(mapped_muts)
}

get_rpoB_for_species <- function(species, seqs) {
  # flexible match: species string anywhere in the name
  hit <- which(grepl(species, names(seqs), ignore.case = TRUE))
  if (length(hit) == 0) {
    stop("No rpoB sequence found for species: ", species)
  }
  if (length(hit) > 1) {
    message("Multiple rpoB sequences for ", species, ", using first match")
  }
  seqs[[hit[1]]]
}





