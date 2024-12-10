
# function to create GenBank query term based on input taxon and gene

create_term <- function(taxon, gene){
  term <- paste0(taxon, "[ORGN] AND ", gene, "[WORD]")
  term
}
