##############################
# Download Genbank sequences #
##############################


# RE- RUNING THIS SCRIPT WILL OVERWRITE THE SEQUENCES DOWNLOADED

set_entrez_key("1374ce1f440c8e235f330bdfbf6ba623b109")
Sys.getenv("ENTREZ_KEY")



# read the accession numbers dataset
acc_num <- read.csv("output/per_group/liverwort_genbank_accessions.csv")

# get rid of the marker i wont use: ndhF and matK
acc_num <- subset(acc_num, select = -c(ndhF, matK))


# keep only the entries with DNA
accessions <- acc_num %>%
  filter(is.na(rbcL) == F | is.na(trnK) == F | is.na(trnL) == F | is.na(ITS) == F )

# # read list of taxa with dna
# taxa_w_dna <- read.csv("output/2023/taxa_with_dna.csv")
# taxa_w_dna <- unlist(taxa_w_dna$taxon)

# keep only the taxa with dna
accessions <- acc_num %>% filter(taxon %in% taxa_w_dna)

sequences <- vector(mode = "list",
                    length = ncol(accessions)-1)



# for each gene
for (i in 2:ncol(accessions)){

  gene_name <- colnames(accessions)[i] # get gene name

  print(paste0("working on gene ", gene_name))

  gene_accessions        <- accessions[ , c(1,i)] # get sequences

  gene_accessions_no_na  <- na.omit(gene_accessions)

  gene_accessions_vector <- gene_accessions_no_na[ ,2] # get rid of NAs

  # make names with binomial+acc number
  names <- paste0(gene_accessions_no_na[ ,1], "_", gene_accessions_no_na[ ,2])

  # download the sequences
  sequences[[i-1]] <- read.GenBank(as.character(gene_accessions_vector),
                                   chunk.size = 300, quiet = F)
  sequences[[i-1]] <- updateLabel(sequences[[i-1]],
                                  old = names(sequences[[i-1]]),
                                  new = names)
}


# name elements of the list with gene names
names(sequences) <- colnames(accessions)[2:ncol(accessions)]



# write fastas

# for (i in 1:length(sequences)){
#   name <- names(sequences)[i]
#   print(name)
#   write.FASTA(sequences[[i]],
#               file = paste0("output/sequences_per_group/liverworts/", name, ".fasta"))
# }


# DOWNLOAD MOSSES ---------------------------------

# read the accession numbers dataset
acc_num <- read.csv("output/per_group/moss_genbank_accessions.csv")

# get rid of the marker i wont use: ndhF and matK
acc_num <- subset(acc_num, select = -ndhF)


# keep only the entries with DNA
accessions <- acc_num %>%
  filter(is.na(rbcL) == F | is.na(matK) == F | is.na(trnK) == F | is.na(trnL) == F | is.na(ITS) == F )

# create an empty list
sequences <- vector(mode = "list",
                    length = ncol(accessions)-2) # -2 cause not doing ITS



# for each gene except ITS which is being problematic
#for (i in 2:ncol(accessions)){
for (i in 2:5){

  gene_name <- colnames(accessions)[i] # get gene name

  print(paste0("working on gene ", gene_name))

  gene_accessions        <- accessions[ , c(1,i)] # get sequences

  gene_accessions_no_na  <- na.omit(gene_accessions)

  gene_accessions_vector <- gene_accessions_no_na[ ,2] # get rid of NAs

  # make names with binomial+acc number
  names <- paste0(gene_accessions_no_na[ ,1], "_", gene_accessions_no_na[ ,2])

  # download the sequences
  sequences[[i-1]] <- read.GenBank(as.character(gene_accessions_vector),
                                   chunk.size = 100, quiet = F)
  sequences[[i-1]] <- updateLabel(sequences[[i-1]],
                                 old = names(sequences[[i-1]]),
                                  new = names)
}



# name elements of the list with gene names
names(sequences) <- colnames(accessions)[2:(ncol(accessions)-1)]



# write fastas

# for (i in 1:length(sequences)){
#   name <- names(sequences)[i]
#   print(name)
#   write.FASTA(sequences[[i]],
#               file = paste0("output/sequences_per_group/mosses/", name, ".fasta"))
# }









#### I need to find a different way for ITS below!

# Manual download for ITS because some sequences can't be downloaded and it affects the loop

ITS_accessions <- accessions[ , c(1,6)] # get sequences

ITS_accessions_no_na  <- na.omit(gene_accessions) # get rid of NAs

ITS_accessions_vector <- as.character(ITS_accessions_no_na[ ,2] )


# create empty list
ITS <- vector(mode = "list", length = nrow(ITS_accessions_no_na))


# find sequences on file
for (i in 1:length(ITS)){

   ITS[[i]] <- read.GenBank(as.character(ITS_accessions_vector[i]))

   name <- paste0(ITS_accessions_no_na[i ,1], "_", ITS_accessions_no_na[i ,2])

   ITS[[i]] <- updateLabel(ITS[[i]],
                           old = names(ITS[[i]]),
                           new = name)

 }



# # export file in a very hacky way!!
#
# for (i in 1:length(ITS)){
#
#   # start the file
#   if (i == 1){
#     write.FASTA(ITS[[i]],
#                 file = "output/sequences_per_group/mosses/ITS.fasta")
#   } else {
#
#     # append to the file
#     write.FASTA(ITS[[i]],
#                 file = "output/sequences_per_group/mosses/ITS.fasta", append = T)
#   }
#
# }











