#!/bin/bash
source ~/.bashrc

pixi shell


# Define working dircetory
wor_dir=/Users/elliegascoyne/Desktop/Projects/ONT

threads=12

phred_quality_score=20

emu_db_path=/Users/elliegascoyne/Desktop/Projects/ONT/databases/custom_db

# Define fastq demultiplexed directory
fastq_dir="${wor_dir}/fastq_files/demultiplexed/16s_leaf"

# Make fastq demultiplexed directory
mkdir -p "${fastq_dir}"

# Manually copy fastq files to the demultiplexed directory

##################
#Primer removal
##################

# Make a directory for primer removal
mkdir -p "${wor_dir}/fastq_files/primers_removed/16s_leaf"

#fwd_primer for its = "TCCGTAGGTGAACCTGCGG"
# rev_primer for its = "GCATATCAATAAGCGGAGGA"
# fwd_primer="AGAGTTTGATCMTGGCTCAG"
# rev_primer="CTACCVGGGTATCTAATCCBG"

for file in "${wor_dir}/fastq_files/demultiplexed/16s_leaf"/*.fastq.gz; do
    [ -e "$file" ] || continue
    i=$(basename "$file" .fastq.gz)

    cutadapt \
        -j "${threads}" \
        -e 0.2 \
        -O 15 \
        --match-read-wildcards \
        --revcomp \
        --discard-untrimmed \
        -g AGAGTTTGATCMTGGCTCAG...VGGATTAGATACCCBGGTAG \
        -o "${wor_dir}/fastq_files/primers_removed/16s_leaf/${i}.fastq.gz" \
        "${wor_dir}/fastq_files/demultiplexed/16s_leaf/${i}.fastq.gz"

done

#######################################
#######################################
#Size filtering
#######################################
#######################################

min_length=650
max_length=950

mkdir -p "${wor_dir}/fastq_files/size_filt/16s_leaf"

for file in "${wor_dir}/fastq_files/primers_removed/16s_leaf/"*.fastq.gz; do
    [ -e "$file" ] || continue
    i=$(basename "$file" .fastq.gz)

    echo "Performing size filtering on sample: $i"

    seqkit seq \
        -m "${min_length}" \
        -M "${max_length}" \
        -g \
        -o "${wor_dir}/fastq_files/size_filt/16s_leaf/${i}.fastq.gz" \
        "${wor_dir}/fastq_files/primers_removed/16s_leaf/${i}.fastq.gz"

done

#######################################
#######################################
#Quality filtering
#######################################
#######################################

#Create directory for output of chopper
mkdir -p "${wor_dir}/fastq_files/chopper/16s_leaf"

for file in "${wor_dir}/fastq_files/size_filt/16s_leaf/"*.fastq.gz; do
    [ -e "$file" ] || continue
    i=$(basename "$file" .fastq.gz)

    echo "Performing quality filtering on sample: $i"

    chopper \
        --input "${wor_dir}/fastq_files/size_filt/16s_leaf/${i}.fastq.gz" \
        --threads "${threads}" \
        --headcrop 0 \
        --tailcrop 0 \
        -q "${phred_quality_score}" \
        -l "${min_length}" \
        --maxlength "${max_length}" |
        gzip > \
            "${wor_dir}/fastq_files/chopper/16s_leaf/${i}.fastq.gz"

done

#############################
# Classification with EMU
#############################

for file in "${wor_dir}/fastq_files/size_filt/16s_leaf/"*.fastq.gz; do
    [ -e "$file" ] || continue
    i=$(basename "$file" .fastq.gz)

    emu abundance \
        --type lr:hq \
        --keep-counts \
        --keep-read-assignments \
        --threads "${threads}" \
        --min-abundance 0.001 \
        --db "${emu_db_path}" \
        --output-basename "${i}" \
        --output-dir "${wor_dir}/taxonomic_classification/emu/" \
        "${wor_dir}/fastq_files/chopper/16s_leaf/${i}.fastq.gz"
done





#Make count and relative abundance taxonomic tables using the 'emu combine-outputs' command


#Create a directory for output of Emu
mkdir -p "${wor_dir}/tables/emu"



#Create a vector of each taxonomic rank
taxa=("species" "genus" "family" "order" "class" "phylum")


#Loop through each taxonomic rank to create the respective table
for t in ${taxa[@]}
do 
emu combine-outputs "${wor_dir}/taxonomic_classification/emu/" "$t"
done



#Create a directory for relative abundance tables produced by the 'emu combine-outputs' command
mkdir -p $wor_dir/tables/emu/rel_abun


#Move tables to relative abundance table
mv $wor_dir/taxonomic_classification/emu/emu-combined*  $wor_dir/tables/emu/rel_abun




#Loop through each taxonomic rank to create the respective table
for t in ${taxa[@]}
do 
emu combine-outputs "${wor_dir}/" "$t" --counts 
done


#Create a directory for count tables produced by the 'emu combine-outputs' command
mkdir -p $wor_dir/tables/emu/counts


#Move tables to count tables to emu/counts directory 
mv $wor_dir/emu_output/emu-combined*    $wor_dir/tables/emu/counts
