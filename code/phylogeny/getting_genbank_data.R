####################################
# Getting sequences from gene bank #
####################################

#install.packages('rentrez')
#install.packages('taxize')
library(rentrez)
library(dplyr)
library(tidyverse)
library(ape)
library(rdrr)
library(taxize)

#----------------#
# Read the data: #
#----------------#

# synonyms dataset
source("scripts/format_names_for_querry.R")

# create character vector with desired genes
#gene_list <- c('atpb', 'trnK', 'trnL', 'matK', 'matR', 'ndhF', 'rbcL', 'trnL-trnF', 'internal transcribed spacer')


# TANNER: change this to your key
set_entrez_key("1374ce1f440c8e235f330bdfbf6ba623b109")
Sys.getenv("ENTREZ_KEY")

#-------------------------------------------------------------
# Source some functions

source("scripts/functions/create_term.R")

source("scripts/functions/get_ID.R")

#####################################
# Getting sequences ids in for loop #
#####################################

# NOTE: the for-loop requires very stable internet connection and
# takes a long time to run. So it might be necessary to run it
# gene by gene instead.


# gene list
gene_list <- c( "rbcL", 'matK', "trnK", 'trnL', 'ndhF', 'internal transcribed spacer')


########----------------  for liverworts!

# empty container
ids.df <- data.frame(hep_wide_names$ACCtrinom)
ids.df$rbcL <- NA
ids.df$matK <- NA
ids.df$trnK <- NA
ids.df$trnL <- NA
ids.df$ndhF <- NA
ids.df$ITS     <- NA


# progress bar
bar <- txtProgressBar(style = 3, width = 40)
nsteps <- nrow(hep_wide_names)


# for each gene
#for (g in 1:length(gene_list)){
for (g in 4:6){

  # print gene
  print( paste0("working on gene ", gene_list[g]))


  # for each species
  for(i in 1:nrow(hep_wide_names)){

    # progress bar
    setTxtProgressBar(bar, i / nsteps)

    # get a temporary vector with all the names for that species
    temp_name_vect   <- as.character(hep_wide_names[i,])[hep_wide_names[i,] != ""]

    # for each synonym
    for (s in 1:length(temp_name_vect)){

      # if there is no ID yet
      if ( is.na(ids.df[i, (g+1) ]) == T ){

        # create a search term
        term_temp <- create_term(temp_name_vect[s], gene_list[g])

        # and save results in the dataframe
        ids.df[i, g+1] <- get_max_ID_all_small(term_temp)
        Sys.sleep(0.1)

      }

    }


  }

}



# save the output
#write.csv(ids.df, "output/per_group/liverwort_genbank_accessions.csv", row.names = F)




########----------------  for mosses!

# empty container
ids.df <- data.frame(moss_wide_names$ACCtrinom)
ids.df$rbcL <- NA
ids.df$matK <- NA
ids.df$trnK <- NA
ids.df$trnL <- NA
ids.df$ndhF <- NA
ids.df$ITS     <- NA


# progress bar
bar <- txtProgressBar(style = 3, width = 40)
nsteps <- nrow(moss_wide_names)


# for each gene
#for (g in 1:length(gene_list)){
for (g in 4:6){

  # print gene
  print( paste0("working on gene ", gene_list[g]))


  # for each species
  for(i in 1:nrow(moss_wide_names)){

    # progress bar
    setTxtProgressBar(bar, i / nsteps)

    # get a temporary vector with all the names for that species
    temp_name_vect   <- as.character(moss_wide_names[i,])[moss_wide_names[i,] != ""]

    # for each synonym
    for (s in 1:length(temp_name_vect)){

      # if there is no ID yet
      if ( is.na(ids.df[i, (g+1) ]) == T ){

        # create a search term
        term_temp <- create_term(temp_name_vect[s], gene_list[g])

        # and save results in the dataframe
        ids.df[i, g+1] <- get_max_ID_all_small(term_temp)
        Sys.sleep(0.1)

      }

    }


  }

}



# save the output
#write.csv(ids.df, "output/per_group/moss_genbank_accessions.csv", row.names = F)

