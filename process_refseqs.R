library(tidyverse)
library(Biostrings)
library(ALJEbinf)
library(tidyverse)
library(rentrez)

source("helpers.R")

# Generate coordinates 
if(file.exists("processed_data/fastahash.Rds") && as.character(openssl::sha1(file("raw_data/rpoB_references.fasta"))) == readRDS("processed_data/fastahash.Rds")){
  print("Sequences file has not changed, loading original coordinates")
  load(file = "processed_data/coordinates.RData")
} else {
  old_fastahash <- as.character(openssl::sha1(file("raw_data/rpoB_references.fasta")))
  saveRDS(old_fastahash, 
          file = "processed_data/fastahash.Rds")
  print("Sequences file has changed, regenerating coordinates")
  #get coordinates
  coordinates <- ALJEbinf::getAllCoordinates(seqs, "rpoB_Escherichia_coli_MG1655")
  save(coordinates, file = "processed_data/coordinates.RData")
}



