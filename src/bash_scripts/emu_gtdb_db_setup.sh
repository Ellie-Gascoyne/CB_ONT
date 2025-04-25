# define working directory
wor_dir=/home/maurice/projects/ellie

db_dir="${wor_dir}/databases"

# Define emu gtdb database path
emu_db_path="${db_dir}/emu_gtdb"

mkdir -p "$emu_db_path"

cd "$emu_db_path" || exit

wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_all/ssu_all_r220.fna.gz

gunzip ssu_all_r220.fna.gz

grep -E '^>' ssu_all_r220.fna | tr -d '>' | sed 's/ d__/\td__/' | sed 's/\[.*//' >ssu_all_r220_taxonomy.tsv

wget https://github.com/shenwei356/gtdb-taxdump/releases/download/v0.5.0/gtdb-taxdump-R220.tar.gz

tar -xvzf gtdb-taxdump-R220.tar.gz

# Run R script

emu build-database "${emu_db_path}/custom_db" \
    --sequences "${emu_db_path}/ssu_all_r220.fna" \
    --seq2tax "${emu_db_path}/seq2tax.map.tsv" \
    --taxonomy-list "${emu_db_path}/taxonomy.tsv"
