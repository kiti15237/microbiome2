#Instructions for how to run the dada2 pipeline.

###1. Dada2 requires separate fastq files for each sample. We therefore rely on the qiime library function to label every read with the sample it came from
  * This command needs to be run FOR EVERY R1 and R2 file of raw reads  
`split_libraries_fastq.py -i /scratch/users/mmdavid/autism_mic_backup/Argonne_Dec2015/Undetermined_S0_L001_R2_001.fastq -o /scratch/users/ctataru5/microbiome/dada2/Argonne1_Dec/R2_r -m /scratch/users/mmdavid/autism_mic_backup/Argonne_Dec2015/maudedavid-1.txt --store_qual_scores --store_demultiplexed_fastq -b /scratch/users/mmdavid/autism_mic_backup/Argonne_Dec2015/Undetermined_S0_L001_I1_001.fastq --barcode_type goley_12`  

###2. Use the R script split_fastq_bySample.R to extract separate fastq files for eah sample.
  * This command needs to be run FOR EVERY R1 and R2 seqs.fastq file output by the previous step  
  * To reflect manual changes made in the mapping file during combination, sample 81 is changed to sample 181.1 and sample 180 in the Argonne2_March batch is changed to 180.1.
`ml load R/3.3.0 `  
`Rscript split_fastq_bySample.R -i /scratch/users/ctataru5/microbiome/dada2/Argonne1_Dec/R2_r/seqs.fastq`  

###3. Delete reads that are not present in both the forward and reverse read files for each sample. Dada2 can't deal with incompatible read identifiers:
`cd {path to fastq files by sample}`  
`Rscript deleteDiffs.R -p {path to fastq files by sample}`  

###4. Run the dada2 algorithm passing in the path name to the folder with fastq by sample
  * You will need to install the phyloseq and dada2 packages
  
      `source("https://bioconductor.org/biocLite.R")`  
      `biocLite('dada2')`  
      `biocLite('phyloseq')`   
  * Outputs an otu table with readSeqs x sampleId in text format, a mapping file that matches the sample ids in the otu table  
      `Rscript runDada2.R -p {path to fastq files} -o {path to output directory} -m {path to mapping file}`  
  * Note: you need to run on big mem cluster. Don't know the exact mem requirements, but regular srun login with no mem specifications will fail

###5. You now have an sequence x sample table, a phylogentic tree, a pure sequence fasta file, and a corresponding mapping file in the out directory. Now use blast to assign the reads to otu ids using green genes database 13_5 (or silva if using tax4fun)
  * download the gg_13_5.fasta.gz databse fasta file from : http://greengenes.secondgenome.com/downloads/database/13_5
  * install blast locally: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
  * `makeblastdb -in {path to gg_13_5.fasta} -parse_seqids -dbtype nucl`
  * `blastn -db {path to gg_13_5.fasta}(same as prev. step ^) -query {path to sequences.fasta}(from step 4) -num_alignments 1 -o {path to output directory}`
  * Parse the output of blast into a table of sequenceid (assigned in uniquestoFasta() function in runDada2) and otuId. Output in same folder as input.
     `clean_blastOutput.py -i {path to blastn output.txt} `
  * The order to sequences remains the same, so otuIds can be directly pasted back into the rownames of the otu_table output by runDada2.
     `assign_blast_otuid.R -b {blast file from clean_blastOutput.py} -s {sequence.fasta from runDada2.py} -o {path to outdir}`
  * Write newly formed seqtab with otu ids into biom format to use with picrust:
     `write_otu_biom -i {path to .txt file you want to convert}`
  * transfer biom files to sherlock
  * Note on normalization: normalizing by copy number should theoretically still allow us to compare between tables. It takes the predicted copy number of each OTU found in nature (using IMG but not clear which datasets from IMG) and divides your abundance by the copy number. This way, if your abundance is 'normal' you get a new abundance of 1. I feel intuitively that this is not appropriate for our analysis, but I don't know why)
  * Create KO table  
     `predict_metagenomes.py -i {path to unnormalized seqtab by otu id .biom} -o {KO_table.biom}`  


***************************************************************
  * Determine contributions of each OTU to each sample  
    *  I'm currently using metagenome_contributions rather blindly, and without any normalization. Need to read more to see what norms need to be done  
     `metagenome_contributions.py -i {seqtab.biom} -o {output file} -gg_version 13_5`   
  * transfer back locally  
  * Convert the list based metagenome_contributions output into table of ko abundance x sample  
     `write_KO_table -i {metagenome_contributions.tab} -m {mapping_dada2.txt}`  
  * Run analysis at will, ko_gene_contribution_analysis.R doesn't seem good.  
