#Instructions for how to run the dada2 pipeline.

###1. Dada2 requires separate fastq files for each sample. We therefore rely on the qiime library function to label every read with the sample it came from
#####split_libraries_fastq.py -i /scratch/users/mmdavid/autism_mic_backup/Argonne_Dec2015/Undetermined_S0_L001_R2_001.fastq -o /scratch/users/ctataru5/microbiome/dada2/Argonne1_Dec/R2_r -m /scratch/users/mmdavid/autism_mic_backup/Argonne_Dec2015/maudedavid-1.txt --store_qual_scores --store_demultiplexed_fastq -b /scratch/users/mmdavid/autism_mic_backup/Argonne_Dec2015/Undetermined_S0_L001_I1_001.fastq --barcode_type goley_12
###This command needs to be run FOR EVERY R1 and R2 file of raw reads


###2. Use the R script split_fastq_bySample.R to extract separate fastq files for eah sample.
####ml load R/3.3.0 
#### Rscript split_fastq_bySample.R -i /scratch/users/ctataru5/microbiome/dada2/Argonne1_Dec/R2_r/seqs.fastq
### This command needs to be run FOR EVERY R1 and R2 seqs.fastq file output by the previous step
#### ** To reflect manual changes made in the mapping file during combination, sample 81 is changed to sample 181.1 and sample 180 in the Argonne2_March batch is changed to 180.1.

###3. Delete reads that are not present in both the forward and reverse read files for each sample. Dada2 can't deal with incompatible read identifiers:
#### cd {path to fastq files by sample}
#### Rscript deleteDiffs.R -p {path to fastq files by sample}
#### ex: Rscript deleteDiffs.R -p /scratch/users/ctataru5/microbiome/dada2/ 

###4. Run the dada2 algorithm passing in the path name to the folder with fastq by sample
  * You will need to install the phyloseq and dada2 packages
  
      `source("https://bioconductor.org/biocLite.R")`
      `biocLite('dada2')`
      `biocLite('phyloseq')` 
  * Rscript runDada2.R -p {path to fastq files} -o {path to output directory} -m {path to mapping file}
  
####  Outputs an otu table with readSeqs x sampleId in text format, a mapping file that matches the sample ids in the otu table
