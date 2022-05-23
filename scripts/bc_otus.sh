#!/bin/bash

# Requirements:
# seqkit >= 2.1.0
# vsearch >= 2.18.0
# minimap2 >= 2.24-r1122

# Must have a reference genome file for S. albidoflavus J1074 named
# "SalbusJ1074.fasta
# Must have a directory with the read data, set it below

###################################################################
# Variables
###################################################################
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
fq_dir="fastq"
gen_ref="../SalbusJ1074.fasta"

###################################################################
# Functions
###################################################################
function die() {
    printf >&2 '[%(%Y-%m-%d %H:%M:%S)T ERROR]: %s\n' -1 "$@"
    exit 1
}

function usage() {
    echo "Usage: $0 <directory containing FASTQ.gz files>"
}

function is_command() {
    command -v "$1" >/dev/null 2>&1
}

function logthis() {
  printf '[%(%Y-%m-%d %H:%M:%S)T INFO]: %s\n' -1 "$@" > /dev/stderr
}

###################################################################
# Main script
###################################################################
# Check directory was given and exists
[[ "$#" -eq 1 ]] || {
    usage
    exit 1
}
[ -d "$1" ] || die "Directory '$1' does not exist!"
is_command seqkit || die "Error: command 'seqkit' not in PATH!"
is_command vsearch || die "Error: command 'vsearch' not in PATH!"
is_command minimap2 || die "Error: command 'minimap2' not in PATH!"

ln -s "$1" "$fq_dir"

# Combine R1 data from two sequencing runs
logthis "Pooling R1 sequencing data..."
seqkit fq2fa \
    -o comb_L001_R1_001.fasta \
    "$fq_dir"/20TAG_S1_L001_R1_001.fastq.gz \
    "$fq_dir"/40TAG_S2_L001_R1_001.fastq.gz

# Use the fixed pattern "CTAATTGGCCGTCGA" to identify sequences with a barcode
logthis "Locating barcode sequences..."
seqkit locate -m 1 -p 'CTAATTGGCCGTCGA' \
    --bed comb_L001_R1_001.fasta \
    -o comb_L001_R1_001_fished.bed

# Extract the barcode portion of the read using the BED file created above,
# discard the unnecessary info inserted by "subseq" in the header
seqkit subseq -d 0 -u 35 -f \
    --bed comb_L001_R1_001_fished.bed \
    comb_L001_R1_001.fasta |
    seqkit replace -p '(.*)_\d+-\d+:[+-]?_usf.*' \
    -r '$1' \
    -o comb_L001_R1_001_barcodes.fasta

# Dereplicate and cluster at 94%
logthis "Dereplicating and clustering barcode sequences..."
vsearch --derep_fulllength \
    comb_L001_R1_001_barcodes.fasta \
    --strand plus \
    --output comb_L001_R1_001_barcodes_derep.fna \
    --sizeout \
    --fasta_width 0 \
    --uc comb_L001_R1_001_barcodes_derep.uc
vsearch --sortbylength \
    comb_L001_R1_001_barcodes_derep.fna \
    --output comb_L001_R1_001_barcodes_sorted.fna
vsearch  --threads 60 \
    --cluster_size comb_L001_R1_001_barcodes_sorted.fna \
    --id 0.94 --iddef 1 \
    --sizein --sizeout  \
    --centroids comb_L001_R1_001_barcodes_cluster94.fna \
    --uc comb_L001_R1_001_barcodes_cluster94.uc

# Keep only clusters with >= 10 members (min10 OTUs)
vsearch --fastx_filter comb_L001_R1_001_barcodes_cluster94.fna \
    --fastaout comb_L001_R1_001_barcodes_cluster94_min10.fna \
    --minsize 10

# Recover all original (before derep/cluster) read headers for min10 OTUs
logthis "Recovering original reads for best OTUs..."
perl "$DIR"/junction_mapping.pl comb_L001_R1_001_barcodes_cluster94.uc \
    comb_L001_R1_001_barcodes_derep.uc |
    sort -k1b,1 \
    > min10_recov_repl_reads.tsv

# Save a map from cluster ID to barcode
grep '^C' comb_L001_R1_001_barcodes_cluster94.uc |
    awk -v OFS='\t' \
    '($3>=10){sub(/;size=[[:digit:]]+/, "", $9); print $0}' |
    cut -f2,9 |
    sort -k2,2 |
    join -1 2 -2 1 -t$'\t' - \
    <(seqkit replace -p ';size=\d+' -r '' \
        comb_L001_R1_001_barcodes_cluster94_min10.fna | 
        seqkit fx2tab | sort) |
    cut -f2,3 > cluster94id_to_bc.tsv

# Trim R2 reads and convert to a table
cat fastq/*_R2_*.fastq.gz |
    seqtk trimfq -L 100 - |
    seqkit fx2tab -Qi |
    sort -k1b,1 \
    > comb_L001_R2_001_trim.tsv

# Mapping
# Plasmid/short sequence assemblies failed; let's try mapping to the reference
# genome instead
logthis "Mapping reads to reference..."
mkdir -p minimap2/ref
mkdir -p minimap2/all_trimmed
ln "$gen_ref" minimap2/ref/

# Join R2 read data table with cluster table created above. Rename headers to
# keep track of OTU and ensure ID uniqueness
join -j 1 -t $'\t' \
    comb_L001_R2_001_trim.tsv \
    min10_recov_repl_reads.tsv |
    sort -k3n |
    perl -F"\t" -anE \
    'state %c_count;
    $c_count{$F[2]}++;
    say ">OTU_$F[0]:$c_count{$F[2]}\n$F[1]"' |
    pigz -9 \
    > minimap2/all_trimmed/all_trimmed.fastq.gz
minimap2 -t 4 -ax sr minimap2/ref/SalbusJ1074.fasta \
    minimap2/all_trimmed/all_trimmed.fastq.gz \
    > minimap2/trimmed_ref_overlaps.sam

logthis "Done!"
logthis "Result SAM file is in minimap2/trimmed_ref_overlaps.sam"

