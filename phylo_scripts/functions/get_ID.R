# function to get GenBank IDs based on desired criteria.
#   This function only gets IDs for sequences that are
#   less than 10,000 base pairs long in order for alignment to work.

get_max_ID_all_small <- function(term){
  search <- entrez_search(db = "nuccore", term = term, retmax = 1)
  if (search$count == 0){
    NA
  } else if (search$count > 0){
    max_ID = NA # variable to store ID associated with entry that has the most base pairs
    max_len_bp = 0 # number of basepairs for longest sequence
    for (ID in search$ids) {
      search_sum <- entrez_summary(db = "nuccore",
                                   id = ID)
      if (search_sum$slen < 10000){
        if (max_len_bp < search_sum$slen){
          max_ID = ID
          max_len_bp = search_sum$slen
        }
      }
    }
    max_ID
  }
}
#---------------------------------------------------------
