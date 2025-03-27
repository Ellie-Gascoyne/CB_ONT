library(tidyr)

# Define working directory
wor_dir <- "/home/maurice/projects/ellie/databases/emu_gtdb"


# Define the path to the taxonomy file
taxa_file <- file.path(wor_dir, "ssu_all_r220_taxonomy.tsv")

taxonomy_df <- read.table(taxa_file,
    check.names = FALSE,
    header = FALSE,
    sep = "\t",
    comment.char = "", # This tells R to ignore the comment character
    fill = TRUE
) # This will fill missing values with NA if any row has fewer columns


# Convert to data frame
taxonomy_df <- as.data.frame(taxonomy_df)



#
has_na <- any(is.na(taxonomy_df))
print(has_na)



# Extract the GCF_ID from the first column
taxonomy_df$GCF_ID <- sub("^[^_]*_([^~]*)~.*", "\\1", taxonomy_df$V1)


has_na <- any(is.na(taxonomy_df))
print(has_na)



# Rename columns
colnames(taxonomy_df) <- c("Full_ID", "taxonomy", "ID")


# Define the path to the taxid map file
taxid_path <- file.path(wor_dir, "gtdb-taxdump-R220/taxid.map")



taxid_map_df <- as.data.frame(read.table(taxid_path,
    check.names = FALSE,
    header = F,
    sep = "\t"
))


# Rename
colnames(taxid_map_df) <- c("ID", "tax_id")


# Merge the data frames
merged_df <- merge(taxonomy_df, taxid_map_df, by = "ID", all.x = TRUE)


# Select the first and fourth columns
seq2tax_df <- merged_df[, c(2, 4)] # Adjust column indices if necessary


# Define the path to the output file
output_file <- file.path(wor_dir, "seq2tax.map.tsv")

# Write the selected columns to a TSV file
write.table(seq2tax_df,
    file = output_file,
    sep = "\t", # Set tab as the separator
    row.names = FALSE, # Set to TRUE if you want to include row names
    col.names = FALSE, # Set to TRUE if you want to include column headers
    quote = FALSE
) # Avoid quoting text


# Select the first and fourth columns
tax_list_df <- merged_df[, c(4, 3)]


# Remove duplicate rows from the data frame
tax_list_df_unique <- unique(tax_list_df)


# Separate the taxonomy column into multiple columns
tax_list_mod_df <- separate(tax_list_df_unique, col = taxonomy, into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "merge")


# Define the path to the output file
output_file <- file.path(wor_dir, "taxonomy.tsv")

write.table(tax_list_mod_df,
    file = output_file, # Specify your path and file name
    sep = "\t", # Set tab as the separator
    row.names = FALSE, # Set to TRUE if you want to include row names
    col.names = TRUE,
    quote = FALSE
) # Avoid quoting text
