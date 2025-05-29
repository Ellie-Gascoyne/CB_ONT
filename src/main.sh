#!/bin/bash
source ~/.bashrc

pixi shell

# Define working dircetory
# Ellie's working directory
wor_dir=/Users/elliegascoyne/Desktop/Projects/ONT
# Maurice's working directory
wor_dir=/home/maurice/projects/ellie

threads=12

phred_quality_score=20

# Define
# cargo install fqkit

# Get the directory of current script
SCRIPT_DIR="$(dirname "$0")"

# Define fastq demultiplexed directory
fastq_dir="${wor_dir}/fastq_files/demultiplexed/$amplicon"

# Make fastq demultiplexed directory
mkdir -p "${fastq_dir}"

# Manually copy fastq files to the demultiplexed directory

for amplicon in "${wor_dir}/fastq_files/demultiplexed/"*; do
    amplicon=$(basename "$amplicon")

    amplicon=16s_leaf

    # QC check the demultiplexed fastq files
    echo "Performing QC check on demultiplexed fastq files for $amplicon"

    #########################################
    #########################################
    # Demultiplexed samples
    #########################################
    #########################################

    # Make a directory for QC checks
    mkdir -p "${wor_dir}/qc/demultiplexed/${amplicon}"

    # Check the size and length distribution of each fastq file
    # Loop through each fastq file in the demultiplexed directory
    for file in "${wor_dir}/fastq_files/demultiplexed/${amplicon}"/*.fastq.gz; do
        [ -e "$file" ] || continue
        i=$(basename "$file" .fastq.gz)

        echo "Performing QC check on sample: $i"
        # Check the size
        fqkit size \
            "${wor_dir}/fastq_files/demultiplexed/${amplicon}/${i}.fastq.gz" > \
            "${wor_dir}/qc/demultiplexed/${amplicon}/${i}_stats.txt"

        fqkit length \
            "${wor_dir}/fastq_files/demultiplexed/${amplicon}/${i}.fastq.gz" > \
            "${wor_dir}/qc/demultiplexed/${amplicon}/${i}_length_distribution.txt"

        Rscript "$SCRIPT_DIR"/r_scripts/plot_sequence_length_distribution.r \
            -s demultiplexed \
            -a "$amplicon" \
            -i "${wor_dir}/qc/demultiplexed/${amplicon}/${i}_length_distribution.txt" \
            -o "${wor_dir}/qc/demultiplexed/${amplicon}/${i}_length_distribution.pdf"

    done

    #######################################
    #######################################
    # Primer removal with Cutadapt
    #######################################
    #######################################
    echo "Performing primer removal with Cutadapt for $amplicon"

    # Make a directory for primer removal
    mkdir -p "${wor_dir}/fastq_files/noprimers/${amplicon}"

    # Dtermine of amplion is 16S or ITS

    if [ "$amplicon" = "16s" ]; then

        primers="AGAGTTTGATCMTGGCTCAG...AAGTCGTAACAAGGTAACCG"
        action=trim

    elif [ "$amplicon" = "16s_leaf" ]; then

        # fwd_primer="AGAGTTTGATCMTGGCTCAG"
        # rev_primer="CTACCVGGGTATCTAATCCBG"

        primers="AGAGTTTGATCMTGGCTCAG...CVGGATTAGATACCCBGGTG"
        action=trim

    elif [ "$amplicon" = "its" ]; then
        primers="TCCGTAGGTGAACCTGCGG...GCATATCAATAAGCGGAGGA"
        action=retain
    else
        echo "Unrecognized amplicon type."
        exit 1
    fi

    amplicon_upper=${amplicon^^}

    echo "Trimming $amplicon primers from the demultiplexed samples with Cutadapt"

    # Run cutadapt to remove primers
    for file in "${wor_dir}/fastq_files/demultiplexed/${amplicon}"/*.fastq.gz; do
        [ -e "$file" ] || continue
        i=$(basename "$file" .fastq.gz)

        echo "Performing primer trimming on sample: $i"
        # Run cutadapt to remove primers
        cutadapt \
            -j "$threads" \
            -e 0.05 \
            -O 12 \
            --revcomp \
            --action="$action" \
            --match-read-wildcards \
            --discard-untrimmed \
            -g "$primers" \
            -o "${wor_dir}/fastq_files/noprimers/${amplicon}/${i}.fastq.gz" \
            "${wor_dir}/fastq_files/demultiplexed/${amplicon}/${i}.fastq.gz"

        #QC check the primer removal
        echo "Performing QC check on primer removal for sample: $i"
        mkdir -p "${wor_dir}/qc/primer_removal/${amplicon}"

        # Check the size
        fqkit size \
            "${wor_dir}/fastq_files/noprimers/${amplicon}/${i}.fastq.gz" > \
            "${wor_dir}/qc/primer_removal/${amplicon}/${i}_stats.txt"

        fqkit length \
            "${wor_dir}/fastq_files/noprimers/${amplicon}/${i}.fastq.gz" > \
            "${wor_dir}/qc/primer_removal/${amplicon}/${i}_length_distribution.txt"

        Rscript "$SCRIPT_DIR"/r_scripts/plot_sequence_length_distribution.r \
            -s noprimers \
            -a "$amplicon" \
            -i "${wor_dir}/qc/primer_removal/${amplicon}/${i}_length_distribution.txt" \
            -o "${wor_dir}/qc/primer_removal/${amplicon}/${i}_length_distribution.pdf"

    done

    #######################################
    #######################################
    #Size filtering
    #######################################
    #######################################

    if [ "$amplicon" = "16s" ]; then
        min_length=1200
        max_length=1800
    elif [ "$amplicon" = "16s_leaf" ]; then
        min_length=600
        max_length=1000
    elif [ "$amplicon" = "its" ]; then
        min_length=200
        max_length=1200
    else
        echo "Unrecognized amplicon type."
        exit 1
    fi

    mkdir -p "${wor_dir}/fastq_files/size_filt/$amplicon"

    for file in "${wor_dir}/fastq_files/noprimers/$amplicon/"*.fastq.gz; do
        [ -e "$file" ] || continue
        i=$(basename "$file" .fastq.gz)

        echo "Performing size filtering on sample: $i"

        seqkit seq \
            -m "${min_length}" \
            -M "${max_length}" \
            -g \
            -o "${wor_dir}/fastq_files/size_filt/$amplicon/${i}.fastq.gz" \
            "${wor_dir}/fastq_files/noprimers/$amplicon/${i}.fastq.gz"

        #QC check the size filtering
        echo "Performing QC check on size filtering for sample: $i"
        mkdir -p "${wor_dir}/qc/size_filtering/$amplicon"
        # Check the size
        fqkit size \
            "${wor_dir}/fastq_files/size_filt/$amplicon/${i}.fastq.gz" > \
            "${wor_dir}/qc/size_filtering/$amplicon/${i}_stats.txt"

        fqkit length \
            "${wor_dir}/fastq_files/size_filt/$amplicon/${i}.fastq.gz" > \
            "${wor_dir}/qc/size_filtering/$amplicon/${i}_length_distribution.txt"

        # Plot the sequence length distribution
        Rscript "$SCRIPT_DIR"/r_scripts/plot_sequence_length_distribution.r \
            -s size_filt \
            -a "$amplicon" \
            -i "${wor_dir}/qc/size_filtering/$amplicon/${i}_length_distribution.txt" \
            -o "${wor_dir}/qc/size_filtering/$amplicon/${i}_length_distribution.pdf"

    done

    #######################################
    #######################################
    #Quality filtering
    #######################################
    #######################################

    #Create directory for output of chopper
    mkdir -p "${wor_dir}/fastq_files/chopper/$amplicon"

    for file in "${wor_dir}/fastq_files/size_filt/$amplicon/"*.fastq.gz; do
        [ -e "$file" ] || continue
        i=$(basename "$file" .fastq.gz)

        echo "Performing quality filtering on sample: $i"

        chopper \
            --input "${wor_dir}/fastq_files/size_filt/$amplicon/${i}.fastq.gz" \
            --threads "${threads}" \
            --headcrop 0 \
            --tailcrop 0 \
            -q "${phred_quality_score}" \
            -l "${min_length}" \
            --maxlength "${max_length}" |
            gzip > \
                "${wor_dir}/fastq_files/chopper/$amplicon/${i}.fastq.gz"

        # QC check the quality filtering
        echo "Performing QC check on quality filtering for sample: $i"
        mkdir -p "${wor_dir}/qc/chopper/$amplicon"
        # Check the size
        fqkit size \
            "${wor_dir}/fastq_files/chopper/$amplicon/${i}.fastq.gz" > \
            "${wor_dir}/qc/chopper/$amplicon/${i}_stats.txt"

        # Check the length distribution
        fqkit length \
            "${wor_dir}/fastq_files/chopper/$amplicon/${i}.fastq.gz" > \
            "${wor_dir}/qc/chopper/$amplicon/${i}_length_distribution.txt"

        # Plot the sequence length distribution
        Rscript "$SCRIPT_DIR"/r_scripts/plot_sequence_length_distribution.r \
            -s chopper \
            -a "$amplicon" \
            -i "${wor_dir}/qc/chopper/$amplicon/${i}_length_distribution.txt" \
            -o "${wor_dir}/qc/chopper/$amplicon/${i}_length_distribution.pdf"

    done

    #######################################
    #######################################
    #ITSxpress
    #######################################
    #######################################

    if [ "${amplicon}" == "its" ]; then
        echo "Performing ITSxpress on the chopper ${amplicon_upper} samples"
        send_slack_message "Performing ITSxpress on the chopper ${amplicon_upper} samples"

        # Create directory for output of chopper
        mkdir -p "${wor_dir}/fastq_files/itsxpress/${amplicon}"

        for file in "${wor_dir}/fastq_files/chopper/${amplicon}"/*.fastq.gz; do
            [ -e "$file" ] || continue
            i=$(basename "$file" .fastq.gz)

            if find "${wor_dir}/qc/itsxpress/${amplicon}" -maxdepth 3 -type f -name "${i}*" | grep -q .; then
                echo "QC files already exist for $i. Skipping ITSxpress step and QC step."
                continue
            fi

            echo "Performing ITSxpress on sample: $i"

            # Find ITS start and end positions with ITSxpress
            itsxpress \
                --fastq "${wor_dir}/fastq_files/chopper/${amplicon}/${i}.fastq.gz" \
                --threads "$threads" \
                --single_end \
                --region ALL \
                --taxa All \
                --outfile "${wor_dir}/fastq_files/itsxpress/${amplicon}/${i}.fastq.gz"

        done
    fi

done

##################################
# Dereplication and OTU clustering
##################################

for amplicon in "${wor_dir}/fastq_files/demultiplexed/"*; do
    amplicon=$(basename "$amplicon")

    amplicon=16s_leaf

    amplicon_upper=${amplicon^^}

    ################
    # Dereplication
    ################

    echo "Dereplicating quality filtered ${amplicon_upper} fastq files"

    if [ "$amplicon" == "16s" ] || [ "$amplicon" == "16s_leaf" ]; then
        dir=chopper
    elif [ "$amplicon" == "its" ]; then
        dir=itsxpress
    else
        echo "Unrecognized amplicon type."
        exit 1
    fi

    # Make a directory for the dereplicated sequences
    mkdir -p "${wor_dir}/fasta_files/dereplicate/${amplicon}/"

    for file in "${wor_dir}/fastq_files/${dir}/${amplicon}/"*.fastq.gz; do
        [ -e "$file" ] || continue
        i=$(basename "$file" .fastq.gz)

        output_file="${wor_dir}/fasta_files/dereplicate/${amplicon}/${i}.fasta"

        if [ -s "$output_file" ]; then
            echo "Output file for $i already exists. Skipping dereplication for sample: $i"
            continue
        fi

        vsearch \
            --fastx_uniques "${wor_dir}/fastq_files/${dir}/${amplicon}/${i}.fastq.gz" \
            --fastaout "$output_file" \
            --sizeout \
            --relabel "$i". \
            --fasta_width 0
    done

done

#Pool all dereplicated samples into one file

for amplicon in "${wor_dir}/fastq_files/demultiplexed/"*; do
    amplicon=$(basename "$amplicon")

    amplicon=16s_leaf

    #Make directory for representative sequences
    mkdir -p "${wor_dir}/fasta_files/representative_sequences/${amplicon}"

    # Check if output file already exists
    output_file="${wor_dir}/fasta_files/representative_sequences/${amplicon}/all_dereplicate_samples_pooled.fasta"
    if [ -s "$output_file" ]; then
        echo "Output file already exists. Skipping concatenation of ${amplicon} samples."
    else
        echo "Concatenating ${amplicon} samples..."
        cat "${wor_dir}/fasta_files/dereplicate/${amplicon}/"*.fasta >"$output_file"
    fi

done

# Create OTU tables and representative sequences

for amplicon in "${wor_dir}/fastq_files/demultiplexed/"*; do
    amplicon=$(basename "$amplicon")

    amplicon=16s_leaf

    vsearch \
        --cluster_size "${wor_dir}/fasta_files/representative_sequences/${amplicon}/all_dereplicate_samples_pooled.fasta" \
        --threads "$threads" \
        --id 0.98 \
        --strand plus \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --centroids "${wor_dir}/fasta_files/representative_sequences/${amplicon}/otu_representative_sequences_${amplicon}.fasta"

    # Sort the representative sequences by size and remove singletons
    vsearch \
        --sortbysize "${wor_dir}/fasta_files/representative_sequences/${amplicon}/otu_representative_sequences_${amplicon}.fasta" \
        --threads "$threads" \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --minsize 3 \
        --output "${wor_dir}/fasta_files/representative_sequences/${amplicon}/otu_representative_sequences_${amplicon}_sorted.fasta"

    # Remove chimeras from the representative sequences
    vsearch --uchime_denovo "${wor_dir}/fasta_files/representative_sequences/${amplicon}/otu_representative_sequences_${amplicon}_sorted.fasta" \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --qmask none \
        --nonchimeras "${wor_dir}/fasta_files/representative_sequences/${amplicon}/otu_representative_sequences_${amplicon}_sorted_nonchimeras.fasta" \
        --chimeras "${wor_dir}/fasta_files/representative_sequences/${amplicon}/otu_representative_sequences_${amplicon}_sorted_chimeras.fasta"

    # Remove host sequences from the representative sequences
    # to do: add host sequences to the database

    # Add label to the representative sequences
    vsearch --fastx_filter "${wor_dir}/fasta_files/representative_sequences/${amplicon}/otu_representative_sequences_${amplicon}_sorted_nonchimeras.fasta" \
        --threads "$threads" \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --relabel OTU_ \
        --fastaout "${wor_dir}/fasta_files/representative_sequences/${amplicon}/otus_sequences_${amplicon}.fasta"

    # Create OTU table from the dereplicated sequences
    vsearch --usearch_global "${wor_dir}/fasta_files/representative_sequences/${amplicon}/all_dereplicate_samples_pooled.fasta" \
        --threads "$threads" \
        --db "${wor_dir}/fasta_files/representative_sequences/${amplicon}/otus_sequences_${amplicon}.fasta" \
        --id 0.98 \
        --strand plus \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --qmask none \
        --dbmask none \
        --otutabout "${wor_dir}/tables/sequence_tables/${amplicon}/otu_table_${amplicon}.txt"

    # Use IDTAXA to assign taxonomy to the representative sequences

    # Create taxa tables

done
