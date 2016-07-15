####Instructions for how to run the dada2 pipeline.

####1. Dada2 requires separate fastq files for each sample. We therefore rely on the qiime library function to label every read with the sample it came from
###split_libraries_fastq.py -i /scratch/users/mmdavid/autism_mic_backup/Argonne_Dec2015/Undetermined_S0_L001_R2_001.fastq -o /scratch/users/ctataru5/microbiome/dada2/Argonne1_Dec/R2_r -m /scratch/users/mmdavid/autism_mic_backup/Argonne_Dec2015/maudedavid-1.txt --store_qual_scores --store_demultiplexed_fastq -b /scratch/users/mmdavid/autism_mic_backup/Argonne_Dec2015/Undetermined_S0_L001_I1_001.fastq --barcode_type goley_12
###This command needs to be run FOR EVERY R1 and R2 file of raw reads

testing testing
