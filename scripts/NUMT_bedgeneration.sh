#!/bin/bash -l
#$ -P icamp
#$ -j y 
#$ -N filtered_out_discorant_and_split_reads
#$ -m n
#$ -l h_rt=48:00:00
#$ -pe omp 16
#load samtools
#match files in given folder to NUMT folder

##################################################################################
## This script generates bed files with regions to exclude from the analysis    ##
## of mtDNA due to presence of NUMTs > length of mean paired read length        ##
##################################################################################


#################
## Version 1.0 ##
#################


# load samtools
module load samtools/1.9

# load picard
module load picard/2.25.2

## MODIFY ## Location of bam files
bam_dir='BAMPATH'
## MODIFY ## Location of NUMT directory for sample
numt_dir='NUMTPATH'
## MODIFY ## Read Length
read_length=150

#loop through all bam files in the directory to get the length
for bam_file in ${bam_dir}/*.sorted.bam; do
    #check if file exists 
    if [ -e "$bam_file" ]; then
        echo "Processing $bam_file"
        #get NUMT directory name from bam file name based on subj
        sampid_temp="${bam_file##*/}"
        sampid=${sampid_temp%.sorted.bam}
        res_dir=$(echo "$sampid" | cut -d'-' -f2)
        #get mean insert size
        picard CollectInsertSizeMetrics -I $bam_file \
        -O insert_size_metrics.txt 
        mean_insert_size=$(grep -A 1 MEDIAN_INSERT_SIZE insert_size_metrics.txt | tail -n 1 | cut -f 1)
        #add read length to mean insert size
        max_NUMT_size=$((mean_insert_size+2*read_length ))

        # now loop through breakpoint files in associated NUMT directory
        # create NUMT bed files
        # concatenate and filter for NUMTs that are >= max_NUMT_size

        #create bed file for NUMTs
        for brkpt_file in ${numt_dir}/${res_dir}/*.bed; do
            #get NUMT directory name from bam file name based on subj
            sampid_temp="${brkpt_file##*/}"
            sampid=${sampid_temp%.bed}
            #filter for NUMTs that are >= max_NUMT_size
            awk -v max_NUMT_size="$max_NUMT_size" '{if ($3-$2 >= max_NUMT_size) print $0}' $brkpt_file > "${numt_dir}/${res_dir}/${sampid}.filtered.bed"
        done
done        
        
        
        
