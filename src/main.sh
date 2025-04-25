#!/bin/bash

# Define working dircetory
wor_dir=/home/maurice/projects/ellie

threads=6

phred_quality_score=20

emu_db_path=/home/maurice/projects/ellie/databases/emu_gtdb/custom_db

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
