### iCamP Lab NUMT Filtering Pipeline

NUMT Calling, NUMT-associated read filtering, and NUMT exclusion region identification Pipeline

This is a pipeline to call NUMTs developed by Wei et. al 2022: https://doi.org/10.1038/s41586-022-05288-7

I then added three components:
1) Filter aligned WGS BAM files to include only reads mapping to chrM
2) Filter chrM-only BAM files to remove any reads poentially attributed to NUMTs by Identification Pipeline
3) Identify NUMT regions and generate exclusion bed files for each sample in chrM that are larger than 2*(read length)+(insert size) 
such that it is impossible to determine if reads mapping to this region are chrM or NUMT derived

Steps to run the pipeline:
1) Copy the contents to the directory in which you want to run NUMT pipeline.
```shell
cp -r ./* $DESTINATION_DIRECTORY
```
2) Create and activate base conda environment.
'''shell
conda env create -f ./envs/snakemake_base_ENV.yml
conda activate snakemake_base_ENV
```
3) Add your bam files to samples directory.
(would probably be easier to change ./samples to current location of your bam files if you have a lot of files.
In order to do this change input file paths in scripts - I can add this for the next version)
```shell
cp $BAM_SOURCE_DIR/*.bam ./samples
```
4) Run the pipeline:
```shell
# Go to scripts directory
cd scripts

# Test that all of the file names are correct
# and get number of jobs that will be run
snakemake -s NUMT_pipeline_v1.snake --dry-run

# Run the pipeline using cluster-generic
snakemake -s NUMT_pipeline_v1.snake --executor cluster-generic --sdm conda \
--cluster-generic-submit-cmd "qsub -P icamp -pe omp {threads}" --jobs ### fill this in based on dry-run output ###
```


