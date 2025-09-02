#!/bin/bash 

##################################################################################
## This script generates bed files with regions to exclude from the analysis    ##
## of mtDNA due to presence of NUMTs > length of mean paired read length        ##
##################################################################################


#################
## Version 1.0 ##
#################


# load samtools
module load samtools

# load picard
module load picard/2.25.2

# load R 
module load R

# Load original bam_file
bam_file=$1

## MODIFY ## Location of bam files
bam_dir='../samples/'
## MODIFY ## Location of NUMT directory for sample
numt_dir='../NUMT_results/'
## MODIFY ## Location of NUMT directory for sample
bed_dir='../exclusion_beds/'

## MODIFY ## Read Length
read_length=150

# Get read length of bam file
if [ -e "$bam_file" ]; then
    echo "Processing $bam_file"

    # Get sample id
    sampid_temp="${bam_file##*/}"
    sampid=${sampid_temp%.chrM.sorted.bam}

    # First make the bed files from breakpoint files
    for brkpt_file in ${numt_dir}${sampid}*Breakpoints.tsv; do
        python3 make_bed.py $brkpt_file
    done
    
    # Get mean insert size for reads for the Bam file
    picard CollectInsertSizeMetrics -I $bam_file \
    -O insert_size_metrics.txt -H insert_size_histogram.pdf -M 0.5
    mean_insert_size=$(grep -A 1 MEDIAN_INSERT_SIZE insert_size_metrics.txt | tail -n 1 | cut -f 1)

    #add read length to mean insert size to determine max acceptable NUMT size
    max_NUMT_size=$((mean_insert_size+2*read_length ))


    # concatenate bed files and filter for NUMTs that are >= max_NUMT_size
    # these are the NUMTs that should be excluded
    for bed_file in ${bed_dir}/${sampid}*.bed; do
        awk -v max_NUMT_size="$max_NUMT_size" '{if ($3-$2 >= max_NUMT_size) print $0}' $bed_file > "${bed_dir}/${sampid}.filtered.bed"
    done
else
    echo "File not found: $bam_file"
fi

# Move all of the files to ./exclusion_beds/
mv *NUMT_regions.bed ../exclusion_beds/

# Make placeholder file for snakemake pipeline
echo "NUMT bed exclusion generation successful" > ${bed_dir}${sampid}_bed_finished.txt