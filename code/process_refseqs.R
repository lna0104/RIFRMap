library(tidyverse)
library(Biostrings)
library(ALJEbinf)
library(tidyverse)
library(rentrez)

source("code/bioinformatics.R")

# -----Generate reference sequences-----
#Load data
data <- read.csv("data/processed_data/RIF_all_sources.csv")
#Load data from rpoB_comprehensive git
refs <- read_csv("data/rpoB_references.csv") |>
  filter(RefStrain == TRUE)
seqs <- readDNAStringSet("./data/rpoB_references.fasta")
#List unique species
species_list <- data |>
  dplyr::mutate(species = stringr::str_squish(species)) |>  # clean species names
  dplyr::filter(!species %in% refs$Species) |>              # exclude reference species
  dplyr::distinct(species) |>                               # unique species only
  dplyr::pull(species)
  

# Select unique species from the mutations data and create a gene reference data frame
selected_muts <- tibble(
    Species = species_list) |>  
  mutate(
    ID = NA,  # Initialize the 'ID' column with NA
    Strain = NA,  # Initialize the 'Strain' column with NA
    RefStrain = NA,  # Initialize the 'RefStrain' column with NA
    Gene = "rpoB",  # Add the gene name as 'rpsL'
    Assembly_ID = NA,  # Initialize the 'Assembly_ID' column with NA
    NCBI_Reference_sequence = NA  # Initialize the 'NCBI_Reference_sequence' column with NA
  ) |>
  distinct()

summaries<-list()

# Loop through each species in the selected_muts data frame
for (i in 1:nrow(selected_muts)){
  species = selected_muts$Species[i]
  cat(paste0("Working on i = ",i, " - ", species, " species\n"))
  
  {
    # Search for reference genomes for the species in the NCBI assembly database
    # Check if the species is Escherichia coli and include strain in the search term
    if (grepl("Escherichia coli", species, ignore.case = TRUE)) {
      term_ser <- sprintf('"%s"[ORGN] AND "MG1655"[Strain] AND ("latest refseq"[filter] AND "reference genome"[filter])', species)
    } else {
      term_ser <- sprintf('"%s"[ORGN] AND ("latest refseq"[filter] AND "reference genome"[filter])', species)
    }
    assembly_search_results <- get_summaries(db="assembly", term=term_ser)
    
    # Check if any assembly search results were found
    if (length(assembly_search_results) == 0) {
      cat(paste0("Cannot find assembly data for ", species, ". Leaving NA.\n"))
      # record failure in the table (no crash, no next)
      selected_muts$RefStrain[i] <- FALSE
      next
    } else {
      assembly_search_final_result <- assembly_search_results[[1]]
      selected_muts$RefStrain[i] <- TRUE
    }

    # Store the assembly summary in the 'summaries' list using the assembly UID as the key
    summaries[[assembly_search_final_result$uid]]<-assembly_search_final_result
    selected_muts$ID[i]<-assembly_search_final_result$uid
  }
  # Store the assembly name in the 'Assembly_ID' column
  selected_muts$Assembly_ID[i]<-assembly_search_final_result$assemblyname
  
  {
    # Fetch links to nucleotide sequence data for the selected assembly
    link_to_nuc <- entrez_link(dbfrom = "assembly", db = "nuccore", id = assembly_search_final_result$uid)
    n_link_ref <- length(link_to_nuc$links$assembly_nuccore_refseq)
    n_link <- length(link_to_nuc$links$assembly_nuccore)
    
    # Check if any nucleotide sequences are found
    if (n_link == 0) {
      cat("No nucleotide information found! Skipping\n")
      next
    }
    
    # Function to extract strain and reference sequence
    extract_strain_info <- function(nuc_result) {
      strain <- nuc_result$strain
      ref_seq <- nuc_result$accessionversion
      
      # If the species is E. coli, force strain to be MG1655
      if (grepl("Escherichia coli", species, ignore.case = TRUE)) {
        strain <- "MG1655"
      }
      
      return(list(strain = strain, ref_seq = ref_seq))
    }
    
    # If no refseq sequences were found
    if (n_link_ref == 0) {
      cat("No refseq information found! Skipping\n")
      next
    }
    
    # If refseq sequences are found
    if (n_link_ref >= 1) {
      nuc_search_results <- entrez_summary(db = "nuccore", id = link_to_nuc$links$assembly_nuccore_refseq)
      if (n_link_ref > 1) {
        strain_info <- extract_strain_info(nuc_search_results[[1]])
      } else {
        strain_info <- extract_strain_info(nuc_search_results)
      }
    }
    
    # Assign extracted strain and reference sequence
    selected_muts$Strain[i] <- strain_info$strain
    selected_muts$NCBI_Reference_sequence[i] <- strain_info$ref_seq
  }
}

#Generate the fasta name
added_name_muts <- selected_muts |>
  mutate(Species_name = gsub(" ", "_", Species),  # Replace spaces with underscores
         Strain_name = gsub(" ", "_", Strain),
         FASTA_name = paste0(Gene, "_", Species_name, "_", Strain_name),
         `NCBI Reference sequence` = NCBI_Reference_sequence) |>
  select(-c(Species_name, Strain_name, NCBI_Reference_sequence)) 

combined <- dplyr::bind_rows(added_name_muts, refs) |>
  select(-ID)

# Save the table as a CSV file
write.csv(combined, "data/rpoB_references.csv", row.names = FALSE)

# Initialize an empty list to store the sequences from all downloads
all_sequences <- list()
# Loop through each genome summary in the 'summaries' list
for (j in 1:length(summaries)) {
  id <- names(summaries[j])
  cat(paste0("Downloading genome for species #", j, "/", length(summaries), "\n..."))
  file_path <- download_file(summaries[[j]], dir = "data/rpoB_references")
  fasta_file <- sub("\\.gz$", "", file_path)
  sequences <- readDNAStringSet(fasta_file)
  matching_seq <- sequences[grep("rpoB|DNA-directed RNA polymerase subunit beta($|[^'])|DNA-directed RNA polymerase subunit beta chain|RNA polymerase, beta subunit|DNA-directed RNA polymerase subunit beta/beta", names(sequences))]
  if (length(matching_seq) == 0){
    cat(paste0("No rpsL found for species ID ", id, "\n"))
    next
  } else if (length(matching_seq) > 1) {
    cat("There are ", length(matching_seq), "copies of rpsL. \n")
    matching_seq <- matching_seq[1]
  }
  names(matching_seq) <- added_name_muts |>
    filter(ID == id) |>
    pull(FASTA_name)
  all_sequences <- c(all_sequences, matching_seq)
}
# Combine all sequences into a single DNAStringSet object
add_seqs <- do.call(base::c, all_sequences)
seqs <- c(seqs, add_seqs)

# Write the combined sequences to a new FASTA file
writeXStringSet(seqs, "data/rpoB_references.fasta")

# Generate coordinates 
if(file.exists("./data/fastahash.Rds") && as.character(openssl::sha1(file("./data/rpoB_references.fasta"))) == readRDS("./data/fastahash.Rds")){
  print("Sequences file has not changed, loading original coordinates")
  load(file = "./data/coordinates.RData")
} else {
  old_fastahash <- as.character(openssl::sha1(file("./data/rpoB_references.fasta")))
  saveRDS(old_fastahash, 
          file = "./data/fastahash.Rds")
  print("Sequences file has changed, regenerating coordinates")
  #get coordinates
  coordinates <- ALJEbinf::getAllCoordinates(seqs, "rpoB_Escherichia_coli_MG1655")
  save(coordinates, file = "./data/coordinates.RData")
}



