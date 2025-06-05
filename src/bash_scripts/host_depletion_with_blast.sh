#!/bin/bash
set -euo pipefail

# Function to show usage
usage() {
    echo "Usage: $0 -s <input_fasta> -o <output_fasta> -n <host_fasta> -t <threads> -r <ref_db_path> [-i <identity_threshold>] [-c <coverage_threshold>] [-d <temp_dir>] [-k] [-h]"
    echo "  -s    Input fasta file (required)"
    echo "  -o    Host depleted fasta file (required)"
    echo "  -n    Host sequences fasta file (required)"
    echo "  -t    Number of threads (default: 12)"
    echo "  -i    Identity threshold (default: 0.97)"
    echo "  -c    Coverage threshold (default: 0.90)"
    echo "  -r    Reference databases path (required)"
    echo "  -d    Temporary directory (default: /tmp/host_depletion)"
    echo "  -k    Keep temporary directory (do not delete at end)"
    echo "  -h    Show this help message"
    exit 1
}

# Default values
input_fasta=""
output_fasta=""
host_fasta=""
threads=12
identity_threshold=0.97
coverage_threshold=0.90
ref_db_path=""
temp_dir="/tmp/host_depletion"
keep_temp=false

# Parse command-line arguments
while getopts ":s:o:n:t:i:c:r:d:kh" opt; do
    case $opt in
    s) input_fasta="$OPTARG" ;;
    o) output_fasta="$OPTARG" ;;
    n) host_fasta="$OPTARG" ;;
    t) threads="$OPTARG" ;;
    i) identity_threshold="$OPTARG" ;;
    c) coverage_threshold="$OPTARG" ;;
    r) ref_db_path="$OPTARG" ;;
    d) temp_dir="$OPTARG" ;;
    k) keep_temp=true ;;
    h) usage ;;
    :)
        echo "Option -$OPTARG requires an argument." >&2
        exit 1
        ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    esac
done

# Check required arguments
if [[ -z "$input_fasta" || -z "$output_fasta" || -z "$host_fasta" || -z "$ref_db_path" ]]; then
    echo "Missing required argument(s)." >&2
    usage
fi

mkdir -p "$temp_dir"

blast_db="${ref_db_path}/host_depletion/GCF_002870075.4_Lsat_Salinas_v11_genomic"
blast_output="${temp_dir}/blast_results.tsv"

echo "Running BLAST to identify host sequences in sample: $input_fasta"

blastn \
    -db "${blast_db}" \
    -num_threads "${threads}" \
    -query "$input_fasta" \
    -outfmt "6 qseqid sallseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" \
    -max_target_seqs 5 \
    -out "$blast_output"

# Filter out sequences which have >= identity_threshold identity and >= coverage_threshold coverage
awk -v idt="$identity_threshold" -v covt="$coverage_threshold" '{if ($3 >= idt * 100 && $4/$9 >= covt) print $0}' \
    "$blast_output" |
    awk '!x[$1]++' |
    awk '{print $1}' >"${temp_dir}/filtered_sequences.txt"

echo "Filtering representative sequences to remove host derived amplicons from sample: $input_fasta"

# Filter representative sequences to remove host derived amplicons
seqkit grep --invert-match \
    -f "${temp_dir}/filtered_sequences.txt" \
    "$input_fasta" > \
    "$output_fasta"

echo "Host depleted fasta file created: $output_fasta"

# Extract host-matching sequences for reference
seqkit grep \
    -f "${temp_dir}/filtered_sequences.txt" \
    "$input_fasta" > \
    "$host_fasta"

echo "Host sequences extracted to: $host_fasta"

# Calculate percentage of host depletion
total_sequences=$(grep -c "^>" "$input_fasta")
kept_sequences=$(grep -c "^>" "$output_fasta")
host_sequences=$(grep -c "^>" "$host_fasta")
echo "Total sequences in input: $total_sequences"
echo "Sequences kept after host depletion: $kept_sequences"
echo "Host sequences extracted: $host_sequences"
depletion_percentage=$(echo "scale=2; 100 * ($total_sequences - $kept_sequences) / $total_sequences" | bc)
echo "Host depletion percentage: $depletion_percentage%"

# Clean up temporary files
if [ "$keep_temp" = false ]; then
    rm -rf "$temp_dir"
    echo "Temporary files cleaned up from: $temp_dir"
else
    echo "Temporary directory kept at: $temp_dir"
fi
echo "Host depletion process completed successfully."
