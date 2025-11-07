#!/bin/bash 

################################################################################
## This script detects NUMTs from whole genome sequencing BAM files
## samtools, samblaster and blat need to be installed to run the pipeline
## samtools can be downloaded at http://www.htslib.org/download/
## samblaster can be downloaded at https://github.com/GregoryFaust/samblaster
## blat can be downloaded at http://hgdownload.soe.ucsc.edu/admin/exe/
################################################################################

INPUT_BAM=$1 # Input WGS bam file
OUTPUT_DIR="../NUMT_results" # output folder path
REF_GRCh38="/PATH/TO/GRCH38/REFERENCE/FASTA" # human reference genome (downloaded from GATK)

# These python scripts are called in the pipeline
CLUSTER_SCRIPT='./searchNumtCluster_fromDiscordantReads.py'
BREAKPOINT_SCRIPT='./searchBreakpoint_fromblatoutputs.py'

# Define sampleid and file variables
SAMPLE_ID1="${INPUT_BAM##*/}" #removes everything in file path to the last slash
SAMPLE_ID2=${SAMPLE_ID1%.bam}  # '%' removes the file name
OUTPUT="${OUTPUT_DIR}/${SAMPLE_ID2}" #new output path in orginal output directory based on sample ID
OUTPUT2="${OUTPUT_DIR}/${SAMPLE_ID2}"
INPUT_DISC="${OUTPUT}.mt.disc.sam"
INPUT_SPLIT="${OUTPUT}.mt.split.sam" 

# Filter input bam for only chrM reads and then use samblaster to generate discorant and split chrM reads
samtools view -@ 28 -m 10G -h -F 2 $INPUT_BAM | grep -e @ -e MT -e chrM | samtools sort -@ 16 -n  | samtools view -h | samblaster --ignoreUnmated -e -d $INPUT_DISC -s $INPUT_SPLIT -o /dev/null

# Run python script to determine approximate genome coordinates of NUMTs based on presence discordant reads pairs (min 5 to define a cluster)
python3 $CLUSTER_SCRIPT ${SAMPLE_ID2} ${INPUT_BAM} ${INPUT_DISC}
echo "look for cluster's done, move to look for breakpoints"

# This takes the input to Breakpoints and generates a multifasta fule of all reads for each NUMT region
# Each read in the multifasta is then blatted against the GRCh38 reference and the .psl outputs of the blat 
# are read into the breakpoint python script to identify breakpoints
while read line; do
    echo $line
    INPUT_Dis="$(echo $line | cut -d ' ' -f 3)"
    INPUT_Split="$(echo $line | cut -d ' ' -f 4)"
    INPUT_WGS="$(echo $line | cut -d ' ' -f 5)"
    CHR="$(echo $line | cut -d ' ' -f 6)"
    START="$(echo $line | cut -d ' ' -f 7)"
    END="$(echo $line | cut -d ' ' -f 8)"
    sampleID="$(echo $line | cut -d ' ' -f 1)"
    REGION="${CHR}:${START}-${END}"
    OUTPUT="${OUTPUT_DIR}/${sampleID}_${CHR}.${START}.${END}"
    samtools view ${INPUT_WGS} ${REGION} >${OUTPUT}.sam

    awk '$6 !~ /150M|149M|148M|149S|148S/' ${OUTPUT}.sam | cut -f1,10 >${OUTPUT}.fasta
    perl -pi -e 's/^/>/g' ${OUTPUT}.fasta
    perl -pi -e 's/\t/\n/g' ${OUTPUT}.fasta
    pblat -threads=28 ${REF_GRCh38}  ${OUTPUT}.fasta  ${OUTPUT}.psl
    python3 ${BREAKPOINT_SCRIPT} ${OUTPUT}.psl ${sampleID} ${CHR} ${START} ${END} ${OUTPUT}
    rm ${OUTPUT}.fasta
    rm ${OUTPUT}.sam
done < ${INPUT_DISC}.breakpointINPUT.tsv

echo "job done"

# Make placeholder file for snakemake pipeline
echo "NUMT Calling Done <3" > ${OUTPUT2}_calling_finished.txt
