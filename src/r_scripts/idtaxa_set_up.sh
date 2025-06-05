cd /home/maurice/resources/marker_gene_databases/16s/gtdb/226.0/full_length/idtaxa

wget https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_all/ssu_all.fna.gz

#Decompress fasta file
gunzip ssu_all.fna.gz

grep -E '^>' ssu_all.fna |
    tr -d '>' |
    sed 's/ d__/\td__/' |
    sed 's/\[.*//' >ssu_all_taxonomy.tsv
