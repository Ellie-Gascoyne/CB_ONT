#!/bin/bash
source ~/.bashrc

pixi shell


# Define working dircetory
wor_dir=/Users/elliegascoyne/Desktop/Projects/ONT

threads=12

phred_quality_score=20

emu_db_path=/Users/elliegascoyne/Desktop/Projects/ONT/databases/custom_db

# Define fastq demultiplexed directory
fastq_dir="${wor_dir}/fastq_files/demultiplexed/$amplicon"

# Make fastq demultiplexed directory
mkdir -p "${fastq_dir}"

# Manually copy fastq files to the demultiplexed directory


############################################
#Primer removal
############################################

amplicon=its

# Make a directory for primer removal
mkdir -p "${wor_dir}/fastq_files/noprimers/${amplicon}"

# Dtermine of amplion is 16S or ITS

if [ "$amplicon" = "16s" ]; then

    primers="AGAGTTTGATCMTGGCTCAG...AAGTCGTAACAAGGTAACCG"
    action=trim

elif [ "$amplicon" = "16s_leaf" ]; then

    primers="AGAGTTTGATCMTGGCTCAG...AAGTCGTAACAAGGTAACCG"
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
        -e 0.10 \
        -O 10 \
        --revcomp \
        --action="$action" \
        --match-read-wildcards \
        --discard-untrimmed \
        -g "$primers" \
        -o "${wor_dir}/fastq_files/noprimers/${amplicon}/${i}.fastq.gz" \
        "${wor_dir}/fastq_files/demultiplexed/${amplicon}/${i}.fastq.gz"

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
    min_length=650
    max_length=950
elif [ "$amplicon" = "its" ]; then
    min_length=200
    max_length=1200
else
    echo "Unrecognized amplicon type."
    exit 1
fi

pixi add --platform "osx-arm64" "seqkit==2.10.0"


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

done

#############################
# Classification with EMU
#############################

if [ "$amplicon" = "16s" ]; then
    emu_db_path=/Users/elliegascoyne/Desktop/Projects/ONT/databases/custom_db
elif [ "$amplicon" = "16s_leaf" ]; then
    emu_db_path=/Users/elliegascoyne/Desktop/Projects/ONT/databases/custom_db
elif [ "$amplicon" = "its" ]; then
    emu_db_path=/Users/elliegascoyne/Desktop/Projects/ONT/databases/unite_emu_db
else
    echo "Unrecognized amplicon type."
    exit 1
fi



for file in "${wor_dir}/fastq_files/size_filt/$amplicon/"*.fastq.gz; do
    [ -e "$file" ] || continue
    i=$(basename "$file" .fastq.gz)

    emu abundance \
        --type lr:hq \
        --keep-counts \
        --keep-read-assignments \
        --threads "${threads}" \
        --min-abundance 0.0001 \
        --db "${emu_db_path}" \
        --output-basename "${i}" \
        --output-dir "${wor_dir}/taxonomic_classification/emu/$amplicon" \
        "${wor_dir}/fastq_files/chopper/$amplicon/${i}.fastq.gz"
done




emu collapse-taxonomy "${wor_dir}/taxonomic_classification/emu/$amplicon/Pos_Ctrl_1_ITS_146_1_rel-abundance.tsv"  'genus'


#Make count and relative abundance taxonomic tables using the 'emu combine-outputs' command


#Create a directory for output of Emu
mkdir -p "${wor_dir}/tables/emu/${amplicon}"







#Create a vector of each taxonomic rank
taxa=("species" "genus" "family" "order" "class" "phylum")


#Loop through each taxonomic rank to create the respective table
for t in ${taxa[@]}
do 
emu combine-outputs "${wor_dir}/taxonomic_classification/emu/${amplicon}"
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
