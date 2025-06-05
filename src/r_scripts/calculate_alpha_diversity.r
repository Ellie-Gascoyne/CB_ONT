library(vegan)
 

# read in table

path_to_otu_table <- ""

# Read in OTU table were rows are features and columns are samples
otu_table <- read.table(path_to_otu_table,
                          header = TRUE, sep = "\t", row.names=1)


# Remove unwanted columns (samples)

# Convert to matrix
otu_table_mat <- as.matrix(otu_table)


# transpose the matrix so rows are samples and columns are features
transposed_table <- t(otu_table_mat)



# remove samples with depth less than 10000
sample_depths <- rowSums(transposed_table)
transposed_table <- transposed_table[sample_depths >= 10000, ]

# Finding the sample with the lowest read depth
min_depth <- min(rowSums(transposed_table))

# Perform rarefaction
rarefied_table <- rrarefy(transposed_table, sample = min_depth)

# 
simpsons_diversity<- diversity(rarefied_table, index = "simpson")

#
observed_richness <- specnumber(rarefied_table)

inverse_simpson <- diversity(rarefied_table, index = "invsimpson")


#  Faith's phylogenetic diversity
# Combine the results into a data frame
results <- data.frame(
  Sample = rownames(rarefied_table),
  Observed_Richness = observed_richness,
    Simpson_Diversity = simpsons_diversity,
    Inverse_Simpson = inverse_simpson
)



output_path <- ""

# Write the results to a CSV file
write.csv(results, output_path, row.names = FALSE)
