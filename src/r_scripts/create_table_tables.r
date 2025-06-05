

path_to_otu_table <- ""

# Read in OTU table were rows are features and columns are samples
otu_table <- read.table(path_to_otu_table,
                          header = TRUE, sep = "\t", row.names=1)


otu_table_100<-apply(otu_table_table,2,function(x){(x/sum(x))*100})


#Make taxonomic tables
phylum_table <- apply(otu_table,2,function(xx){tapply(xx,taxonomy_table[rownames(otu_table),"phylum"],sum)})
class_table <- apply(otu_table,2,function(xx){tapply(xx,taxonomy_table[rownames(otu_table),"class"],sum)})
order_table <- apply(otu_table,2,function(xx){tapply(xx,taxonomy_table[rownames(otu_table),"order"],sum)})
family_table <- apply(otu_table,2,function(xx){tapply(xx,taxonomy_table[rownames(otu_table),"family"],sum)})
genus_table <- apply(otu_table,2,function(xx){tapply(xx,taxonomy_table[rownames(otu_table),"genus"],sum)})
species_table <- apply(otu_table,2,function(xx){tapply(xx,taxonomy_table[rownames(otu_table),"species"],sum)