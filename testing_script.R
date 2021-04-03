library(tidyverse)
library(Biostrings)

set.seed(1999)

data <- Biostrings::DNAStringSet(replicate(100000, paste0(sample(c("A", "C", "T", "G"), 24, replace = TRUE), collapse = ""))) %>% 
  translate()

# Removes any stop codons "*" and saves new data as n_data. Also provides seq names.
# Saves the n_data as a .fasta file
n_data <- data[vcountPattern("*", data) == 0]
names(n_data) <- paste0("seq_",seq_along(n_data))

# Finds a random sequence from n_data
# duplicates it 100 times and adds one random AA mutation at position 1,3,6 or 8
seq_e <- n_data[[sample(seq_along(n_data), 1)]]
enrich_set <- AAStringSet()
for(i in 1:100){
  enrich_set[i] <- as(seq_e, "AAStringSet")
}
places <- c(1,3,6,8)
for(i in 1:100){
  j <- sample(places,1)
  x <- AAString(sample(AA_STANDARD, 1))
  enrich_set[[i]][j] <- x
}
names(enrich_set) <- paste0("test_",seq_along(enrich_set))

# n_data is randomised to a training and testing set, both with 3^4 sequences
# testing set contains the 100 mutated sequences 
r_data <- sample(n_data)

bg <- r_data[1:30000]
test <- c(r_data[30001:59900], enrich_set)

# Files are outputed as a .fasta
writeXStringSet(bg, "bg.fasta")
writeXStringSet(test, "test.fasta")
writeXStringSet(n_data, "n_data.fasta")
writeXStringSet(enrich_set, "enrich_set.fasta")

# This could be optomised with a lapply. 
#However, I had issues with retaining argument names. 


