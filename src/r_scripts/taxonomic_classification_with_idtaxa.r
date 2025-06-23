library(DECIPHER)


# Set working directory
wor_dir <- "path/to/your/working/directory"

fasta_file_path <- "/home/maurice/projects/ellie/fasta_files/representative_sequences/16s_leaf/otus_sequences_16s_leaf.fasta"


seqs <- readDNAStringSet(fasta_file_path)



if (amplicon == "its") {
  ref_file <- "/Users/elliegascoyne/resources/tempobiome/taxonmic_classifcation/its/idtaxa/UNITE_v2025_February2025.RData"
} else if (amplicon == "16s_leaf") {
  ref_file <- "/Users/elliegascoyne/resources/tempobiome/taxonmic_classifcation/16s/idtaxa/GTDB_r226-mod_April2025.RData"
} else if (amplicon == "18s") {
  ref_file <- "/Users/elliegascoyne/resources/tempobiome/taxonmic_classifcation/18s/idtaxa/PR2_v4.14.0.RData"
} else {
  stop("Invalid amplicon type. Please choose 'its' or  '16s_leaf'")
}


load(ref_file)

ids <- IdTaxa(seqs,
   trainingSet,
   strand="both", # or "top" if same as trainingSet
   threshold=50,
   processors=NULL)


ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest


taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))


colnames(taxid) <- ranks

taxid <- as.data.frame(taxid)

sum(!is.na(taxid$genus))
