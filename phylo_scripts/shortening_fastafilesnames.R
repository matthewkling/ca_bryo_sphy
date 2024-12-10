########################################
## Making fasta files with short names
########################################

#loading libraries
library(dplyr)
library(tidyverse)
#install.packages("phylotools")
library(phylotools)

# DO EVERYTHING BELOW FOR BOTH LIVERS AND MOSSES

# STEP 1: Constructing reference tables for all the markers,

# get a list of alignments
align      <- list.files("output/sequences_per_group/mosses", full.names = T)
gene_names <- tools::file_path_sans_ext(list.files("output/sequences_per_group/mosses"))

# ref_table
refs <- vector(mode = "list", length = length(align))

# for each alignment
for ( i in 1:length(align)){
#for ( i in 1:1){

  # get the old names
  old_names <- data.frame( get.fasta.name(align[i]) )

  # get the new names
  new_names <- old_names %>%  separate(1, into = c("name", "accesion") , sep = "_")
  new_names <- new_names %>%  separate(name, into = c("gen", "sp", "extra"), sep = " ")%>%
    unite("short_name", c(gen, sp), sep = "_")

  # make a reference table
  reftable <- data_frame(old_names, new_names$short_name, new_names$accesion)
  colnames(reftable) <- c("old_names", "short_names", "accession")

  refs[[i]] <- reftable
}


# save to files
for (i in 1:length(refs)){

  # write the ref table
  write_csv(refs[[i]], paste0("output/sequences_per_group/mosses/reftabs/", gene_names[i], ".csv"))

}




#STEP 2:
#changing the names of the fasta files

# for each alignment
for (i in 1:length(gene_names)){

  # read the ref table
  ref_table <- read.csv( paste0("output/sequences_per_group/mosses/reftabs/", gene_names[i], ".csv") )

  # rename the fasta
  rename.fasta(infile = paste0("output/sequences_per_group/mosses/", gene_names[i], ".fasta"),
               ref_table = refs[[i]],
               outfile = paste0("data/raw_short_seq/mosses/", gene_names[i], ".fasta"))

}

