#!/bin/bash 

##################################################################################################
## Thus script removes discordant, split and pblat associated reads from aligned WGS bam files  ##
## Before running this script you need to run the NUMT identification pipeline                  ##
##################################################################################################

##################################
## Version 2.0 ###################
##################################

## MODIFY ## Location of NUMT directory for sample
numt_dir='../NUMT_results/'               
bam_dir='../samples/'
filtered_bam_dir='../filtered_bams/'

# Get bam file from stdin
bam_file=$1

# Check if file exists 
if [ -e "$bam_file" ]; then
    echo "Processing $bam_file"
    # Get NUMT directory name from bam file name based on subj
    sampid_temp="${bam_file##*/}"
    sampid=${sampid_temp%.chrM.sorted.bam}                                                                                                       

    #make text file of all discordant, split, and pblat associated reads to remove
    #discordant reads
    samtools view -@ 16 "${numt_dir}${sampid}.mt.disc.sam" | cut -f1 | sort -u > ${numt_dir}${sampid}.disc_to_rm.txt
    #split reads
    samtools view -@ 16 "${numt_dir}/${sampid}.mt.split.sam"| cut -f1 | sort -u > ${numt_dir}${sampid}.split_to_rm.txt
    #pblat associated reads (this theoretically should not matter since pblat reads don't map initially to chrM, but somehow some do)
    cat ${numt_dir}${sampid}*.psl | grep 'chrM' | awk '{print $10}' > ${numt_dir}${sampid}.pblat_to_rm.txt
    cat ${numt_dir}${sampid}.disc_to_rm.txt ${numt_dir}${sampid}.split_to_rm.txt ${numt_dir}${sampid}.pblat_to_rm.txt \
    > ${numt_dir}${sampid}.remove_list.txt
    #remove all NUMT-associated reads
    samtools view -@ 16 ${bam_file} -N ${numt_dir}${sampid}.remove_list.txt -U "${filtered_bam_dir}/${sampid}_no_NUMT_reads.bam"

    #calculate the number of reads removed and write to a summary file
    input_len=$(wc -l < "$bam_file")
    output_len=$(wc -l < "${filtered_bam_dir}/${sampid}_no_NUMT_reads.bam")
    dif=$((input_len - output_len))
    echo "$dif NUMT associated reads removed from $sampid" > "${filtered_bam_dir}/${sampid}_numt_read_removal_summary.txt"
    echo "Number of reads in original bam: $input_len" >> "${filtered_bam_dir}/${sampid}_numt_read_removal_summary.txt"
    echo "Number of reads in filtered bam: $output_len" >> "${filtered_bam_dir}/${sampid}_numt_read_removal_summary.txt"
    echo "Number of discordant reads removed: $(wc -l < ${numt_dir}${sampid}.disc_to_rm.txt)" >> "${filtered_bam_dir}/${sampid}_numt_read_removal_summary.txt"
    echo "Number of split reads removed: $(wc -l < ${numt_dir}${sampid}.split_to_rm.txt)" >> "${filtered_bam_dir}/${sampid}_numt_read_removal_summary.txt"
    echo "Number of pblat reads removed: $(wc -l < ${numt_dir}${sampid}.pblat_to_rm.txt)" >> "${filtered_bam_dir}/${sampid}_numt_read_removal_summary.txt"
    
    # remove intermediate files
    rm -f ${numt_dir}${sampid}.disc_to_rm.txt ${numt_dir}${sampid}.split_to_rm.txt ${numt_dir}${sampid}.pblat_to_rm.txt ${numt_dir}${sampid}.remove_list.txt
else
    echo "File not found: $bam_file"
fi


