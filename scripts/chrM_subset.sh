#!/bin/bash 

########################################################################
## This script subsets BAM files to only preserve chrM_reads ###########
########################################################################

## load modules on HPC
module load samtools 
#define input and output directories
INPUT_BAM=$1
OUTPUT_DIR="../filtered_bams/"
#loop through each bam file and get only chrM reads and write out
SAMPLE_ID_TEMP="${INPUT_BAM##*/}"
SAMPLE_ID=${SAMPLE_ID_TEMP%.bam}
samtools view --threads 16 -b -h ${INPUT_BAM} chrM | samtools sort --threads 16 > ${OUTPUT_DIR}${SAMPLE_ID}.chrM.sorted.bam
samtools index -@ 16 "${OUTPUT_DIR}${SAMPLE_ID}.chrM.sorted.bam"