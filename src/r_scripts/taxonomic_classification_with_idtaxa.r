library(DECIPHER)


# Set working directory
wor_dir <- "path/to/your/working/directory"

fasta_file_path <- "/home/maurice/projects/ellie/fasta_files/representative_sequences/16s_leaf/otus_sequences_16s_leaf.fasta"


seqs <- readDNAStringSet(fasta_file_path)


load("/home/maurice/Downloads/GTDB_r226-mod_April2025.RData")

ids <- IdTaxa(seqs,
   trainingSet,
   strand="both", # or "top" if same as trainingSet
   threshold=60, # 60 (cautious) or 50 (sensible)
   processors=NULL)