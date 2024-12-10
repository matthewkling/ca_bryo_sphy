##############################
#  Concatenating alignements
################################

#install.packages("devtools")
#devtools::install_github("tardipede/concatipede")


library(concatipede)


# for liverworts:

path <- "data/alignments/liverworts/"

files <- find_fasta(path)

  # preparing matching file
concatipede_prepare(files, out = paste0(path, "seqnames"))

  # check the matching file makes sense before running line below

  # concatenate
concatipede(filename = paste0(path, "seqnames_checked.xlsx"),
            dir = path,
            out = "con_align_liverwort",
            plotimg = T)
