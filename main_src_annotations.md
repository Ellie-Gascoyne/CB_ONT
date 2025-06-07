# Main src annotations
This markdown file contains in depth annotations for the code being used. The main src code contains notes for the aim of the overall step, however this markdown file contains notes on the syntax itself. 

## Working in bash
Bash (Bourne Again shell) is a shell that can be used on macOS. You can also run Z Shell (Zsh) which is the default for later mac software that can autocorrect for commands, and more powerful auto-completion.

To change the shell from zsh to bash, run the following command in the terminal:

```
chsh -s /bin/bash
```
This can be changed back to zsh by using the same command but change 'bash' to 'zsh'. 

Other common commands for bash include listing all files and directories for current directory 

```
ls
```
making a new directory:
```
mkdir directory_name
```
navigating directories:
```
cd /path/to/directory
```
displaying current directory
```
pwd
```
removing files:
```
rm file_name.txt
```
Finding files:
```
find /path-name "filename"
```

## Set up pixi shell

It has been noted that sometimes bash does not automatically understand pixi commands. Hence, to open the pixi shell at the start of every session, the following code can be used:

```
source ~/.bashrc

pixi shell
```

The pixi shell is where the code will be run from. However, a seperate terminal with a bash shell can be used to check that the outputs of each step, e.g. that the files are put into the correct directory.

## Main Script notes 

### Step 1 Labelling Objects

The following objects should be labelled as follows the rest of the script: 

```
# Ellie's working directory
wor_dir=/Users/elliegascoyne/Desktop/Projects/ONT

ref_db_dir=/Users/elliegascoyne/resources/tempobiome

threads=12

phred_quality_score=20

SCRIPT_DIR=/Users/elliegascoyne/Documents/GitHub/CB_ONT/src

fastq_dir="${wor_dir}/fastq_files/demultiplexed/$amplicon"
```
### Step 2- Run a QC on the demultiplexed files

The overall aim of this step is to check the number of sequences is reduced at each step, and that the size of the fragements being kept is within the correct thresholds. The code first checks that the files exist, then run a qc using fqkit. A length distriubtion plot (looks a little like a mass spec output) is then generated via R.

Firstly, make sure there to make the output directory 

```
 mkdir -p "${wor_dir}/qc/demultiplexed/${amplicon}"
```
Then, the qc is run for each file with the specific amplicon in its file name as a forloop using fqkit.
```
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
```
#### Step 2 in Detail

```
for file in "${wor_dir}/fastq_files/demultiplexed/${amplicon}"/*.fastq.gz; do
        [ -e "$file" ] || continue
        i=$(basename "$file" .fastq.gz)
echo "Performing QC check on sample: $i"
```
${wor_dir} and ${amplicon} are variables for the working directory and the name of the amplicon group.
This line loops over every .fastq.gz file in the specified directory (.../demultiplexed/${amplicon}/).
If the file doesn't exist (-e checks for existence), skip this iteration (continue).
Extracts the base filename (i.e., removes the .fastq.gz extension) to use it as a sample identifier (i).
Prints a message indicating that QC is starting for the sample.

```
        fqkit size \
            "${wor_dir}/fastq_files/demultiplexed/${amplicon}/${i}.fastq.gz" > \
            "${wor_dir}/qc/demultiplexed/${amplicon}/${i}_stats.txt"
```
Runs fqkit size to get size and metadata of the fastq file then saves the result in a .txt file in the qc/demultiplexed/... folder.

```
        fqkit length \
            "${wor_dir}/fastq_files/demultiplexed/${amplicon}/${i}.fastq.gz" > \
            "${wor_dir}/qc/demultiplexed/${amplicon}/${i}_length_distribution.txt"
```

Computes the distribution of sequence lengths in the fastq file then saves it to another .txt file for downstream plotting.

```
        Rscript "$SCRIPT_DIR"/r_scripts/plot_sequence_length_distribution.r \
            -s demultiplexed \
            -a "$amplicon" \
            -i "${wor_dir}/qc/demultiplexed/${amplicon}/${i}_length_distribution.txt" \
            -o "${wor_dir}/qc/demultiplexed/${amplicon}/${i}_length_distribution.pdf"
```
Calls an R script to generate a PDF plot of the read length distribution.

Inputs:

-s: the sample type (demultiplexed)
-a: the amplicon name
-i: the input stats file (output text file from the fqkit length step)
-o: the output PDF file (saved to the same destination as the fqkit files.

**Note- the QC from each step will be in its own seperate repo.**

### Step 3- Removing primer sequences with Cutadapt

The overall aim of this step determines which primer set to use and whether to cut or retain them. The primers are then removed via Cutadapt. Trimmed reads are then saved into a new output repo. The same qc step as seen in Step 2 is then used to make sure the primers have been removed, this should be seen in the R output plot as the length should be reduced.

```
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
```

#### Step 3 (primer removal) in detail:

```
echo "Performing primer removal with Cutadapt for $amplicon"
```
Prints which amplicon is being used.
```
mkdir -p "${wor_dir}/fastq_files/noprimers/${amplicon}"
```
Creates the output directory for the amplicon in 'noprimers' repo if not made before. -p ensures it creates parent directories if needed, and doesn't error if it already exists.
```
if [ "$amplicon" = "16s" ]; then
    primers="AGAGTTTGATCMTGGCTCAG...AAGTCGTAACAAGGTAACCG"
    action=trim

```
In this example, the 16s primers have been chosen so the 'action=trim' means Cutadapt will remove these primer sequence from reads. (16s_leaf primers should be used for the Arabidopsis work).

```
elif [ "$amplicon" = "its" ]; then
    primers="TCCGTAGGTGAACCTGCGG...GCATATCAATAAGCGGAGGA"
    action=retain

```
For its primers, the action=retain command is used, meaning primers are kept in the reads. This is because ITSxpress later in the pipeline for the ITS samples. This needs the primers to recognize the sites to trim off the conserved regions. For this instance, ITSexpress will trim the primers and conserved regions. 
```
else
    echo "Unrecognized amplicon type."
    exit 1

```
if an unknown amplicon type is detected, the code will exit with an error.
```
echo "Trimming $amplicon primers from the demultiplexed samples with Cutadapt"

```
prints which primers are being trimmed from the sample with Cutadapt.
```
for file in "${wor_dir}/fastq_files/demultiplexed/${amplicon}"/*.fastq.gz; do
    [ -e "$file" ] || continue
    i=$(basename "$file" .fastq.gz)

```
This loops over each fastq.gz file in the demultiplexed directory for the desired amplicon. It then extracts the sample name (without the fast.gz extension) as a variable (i)

```
cutadapt \
    -j "$threads" \                     # Number of threads to use
    -e 0.05 \                           # Max error rate (5%)
    -O 12 \                             # Minimum overlap with primer
    --revcomp \                         # Also match reverse complements
    --action="$action" \               # trim or retain based on amplicon
    --match-read-wildcards \           # Allow ambiguous bases (e.g., N)
    --discard-untrimmed \              # Discard reads that don’t match primers
    -g "$primers" \                    # Primer sequences
    -o "${wor_dir}/fastq_files/noprimers/${amplicon}/${i}.fastq.gz" \
    "${wor_dir}/fastq_files/demultiplexed/${amplicon}/${i}.fastq.gz"

```
runs cutadapt on the desired primers. Please find notes for each line in the code snippet. The trimmed files are then sent to the output repo 'noprimers' with each amplicon having a seperate folder.
```
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
```
The same qc step as before is carried out for all samples, but this time putting them into a new output qc folder for 'primer removal'. In this instance the 'sample-type' (-s) is noprimers for the desired amplicon (-a) and taken from the primer-removal text file output from fqkit length output and creates a PDF plot.

### Step 4- Size Filtering

The overall aim of this step is to set the thresholds for expected size of the fragments once the primers have been removed. Seqkit can then remove the sequences that fall outside of the range. The same qc step is used again to confirm that the unwanted sequences have been removed. 
```
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
```
#### Step 4 (size filtering) in detail:
```
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

```
Sets the minimum and maximum length of sequences for each amplicon type (depending on input) and will exit out with an error for unrecognised sequences. (elif = else if, e.g. if not 16s, then do 16s_leaf, if not then do its). 
```
mkdir -p "${wor_dir}/fastq_files/size_filt/$amplicon"
```
Makes size_filt directory for amplicon (making any necessary parent directories). 
```
for file in "${wor_dir}/fastq_files/noprimers/$amplicon/"*.fastq.gz; do
    [ -e "$file" ] || continue
    i=$(basename "$file" .fastq.gz)

```
Makes a forloop for each file in the 'noprimers' repo for the desired amplicon to extract the sample name as a variable (i). 
```
seqkit seq \
    -m "${min_length}" \   # Minimum sequence length
    -M "${max_length}" \   # Maximum sequence length
    -g \                   # Output in gzipped format
    -o "${wor_dir}/fastq_files/size_filt/$amplicon/${i}.fastq.gz" \
    "${wor_dir}/fastq_files/noprimers/$amplicon/${i}.fastq.gz"

```
performs size filtering using Seqkit. Please find notes for each line in code snippet. Overall, it discards any reads that fall outside the desired length range from the no primer files and sends the output to the size_filt directory. 
```
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

```
Carried out same qc step with fqkit as before, to make sure that the undesired reads (falling outside the length range) have been removed.

### Step 5- Quality Filtering

The overall aim of this step is to filter out undesired reads that do not meet the desired phred score of 20. The phred quality score measures the identification of the nucleobases generated by the ONT and denotes how confident we are in the assignment of each base call by the sequencer. For many purposes, a Phred Score of 20 or above is acceptable, because this means that whatever it qualifies is 99% accurate, with a 1% chance of error. This is especially important for ONT as lower accuracy (95%) than illumina (99%). 
```
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
```
#### Step 5 (Quality Filtering) in detail: 
```
mkdir -p "${wor_dir}/fastq_files/chopper/$amplicon"

```
Firstly, a directory is made for the amplicon in the new 'chopper' file. Any non-existing parent files are created (-p).
```
for file in "${wor_dir}/fastq_files/size_filt/$amplicon/"*.fastq.gz; do
    [ -e "$file" ] || continue
    i=$(basename "$file" .fastq.gz)

```
This loops through the size filt directory for the desired amplicon and extracts the sample name as a variable (i) and the continue check skips if no files match the pattern (avoids errors).
```
chopper \
    --input "${wor_dir}/fastq_files/size_filt/$amplicon/${i}.fastq.gz" \ #Specifies the input FASTQ file (from size filtering).
    --threads "${threads}" \                                             #Uses multiple CPU threads (from ${threads} variable).
    --headcrop 0 \                                           
    --tailcrop 0 \                                                       # No cropping from either end of the reads.
    -q "${phred_quality_score}" \                                        # Minimum acceptable Phred quality score.
    -l "${min_length}" \                                                 # Minimum length to retain a read after trimming.
    --maxlength "${max_length}" |                                        # Maximum length allowed for a read after trimming.
    gzip > \
    "${wor_dir}/fastq_files/chopper/$amplicon/${i}.fastq.gz"             # Compresses the output and saves it to the chopper output directory.

```
This code uses Seqkit to keep reads that have a phred score of 20 or higher (i.e have 99% basecalled accuracy) within the reads that have already been filtered for desired length. Any reads that don't reach the quality threshold are discarded.
```
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
```
The same qc step as before using 'chopper' as a sample name. the output plots and files from the fqkit and rscript are sent to the newly created chopper directory for the desired amplicon in qc directory. The number of reads should have decreased, but the length of reads should remain within the desired thresholds. 

### Step 6- ITSxpress (merging and trimming paired-end sequences for ITS reads only)
The overall aim of this step is to extract just the ITS region from the sequencing reads from the chopper directory. This is similar to the 'Cutadapt' function performed for 16s samples but ITxpress is specifically designed for ITS and can trim the ITS sequences at a faster rate.
```
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
```
#### Step 6 - ITSxpress in detail 
```
if [ "${amplicon}" == "its" ]; then

```
This is a conditional check to make sure the ITSxpress only runs on files labelled as with the 'its' amplicon.
```
    echo "Performing ITSxpress on the chopper ${amplicon_upper} samples"
    send_slack_message "Performing ITSxpress on the chopper ${amplicon_upper} samples"

```
This outputs a log message to the terminal to indicate that the ITSexpress code has started to run.
```
    mkdir -p "${wor_dir}/fastq_files/itsxpress/${amplicon}"

```
Creates a new repo in the fastq_files directory for the itsxpress output files.
```
    for file in "${wor_dir}/fastq_files/chopper/${amplicon}"/*.fastq.gz; do
        [ -e "$file" ] || continue
        i=$(basename "$file" .fastq.gz)
```
The start of the forloop which scans all the files in the chopper directory for any files that have the 'its' amplicon. If the file doesn't exist the '-e' function skips to the next one. Finally, the 'i=' step extracts the sample name from the fast.gz files whilst leaving the extension. 
```
        if find "${wor_dir}/qc/itsxpress/${amplicon}" -maxdepth 3 -type f -name "${i}*" | grep -q .; then
            echo "QC files already exist for $i. Skipping ITSxpress step and QC step."
            continue
        fi
```
This makes sure the ITSexpress qc file exists already for the sample. If the file is found it is skipped to avoid multiple of the same file. 
```
        echo "Performing ITSxpress on sample: $i"
        itsxpress \
            --fastq "${wor_dir}/fastq_files/chopper/${amplicon}/${i}.fastq.gz" \  #input file from the filtered 'chopper' file.
            --threads "$threads" \                                                # desiered number of threads to use as set out in the variable at the start of the pipeline
            --single_end \                                                        # Describes reads as single reads as opposed to paired reads
            --region ALL \                                                        # Extracts all ITS regions (e.g. ITS1, ITS2 etc.)
            --taxa All \                                                          # Works for all taxa with no restrictions.
            --outfile "${wor_dir}/fastq_files/itsxpress/${amplicon}/${i}.fastq.gz" #Describes the output path directory.
    done
fi

```
The message 'Performing ITSxpress on sample (SAMPLENAME)' is printing so we know which sample is being processed. Find details for each step in the code snippet. ITS regions have now been extracted and and merged for fungal identification. 

### Step 7- OTU clustering a
