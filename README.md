#Important Data and Location
###Otu Tables: 
######in dada2 folder. Use filtered_otu_table_normCSS.biom, read in to R with phyloseq import_biom() function. Look in scripts/getTables.R for example of importing. 
### Combined otu table and relevant mapping file info: 
######filtered_otuTable_normCSS_info.txt. contig by sample(and mapping info). In the class column, 1 means case, 0 means control
###Combind otu table and relevant mapping filer into AND tax assignments for contigs
######filtered_otuTable_normCSS_info_taxa.txt. contig by sample and mapping info and taxanomic assignment
### Tax Tables and Trees: 
######in dada2 folder. It is recommended to use tree file tree_midpointRoot.tre until a better rooting method is found. Can be read in as done in first chunk of scripts/Summary.Rmd
### /scripts/Summary.Rmd. 
######Easiest way to read in data is displayed in first chunk. Most relevant analyses displayed here

#Pipeline summary:
We received the paired end sequencing information from Argonne. I'm not sure what protocol was used. We processed the sequencing information using dada2. The goal of dada2 is to assemble reads into denoised, chimera-free contigs, and then assign taxonomy to those contigs. Their claim to fame is that they resolve fine scale variation instead of lumping contigs into OTU's and taking a representative sequence. The details of this method can be found here: http://www.nature.com/nmeth/journal/v13/n7/pdf/nmeth.3869.pdf. The output of dada2 is a table of abundances, taxa by samples. We normalize these abundances using CSS after removing the outliers. You can see details about CSS here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4010126/, and details about outlier removal further in this document. We then build a phylogenetic tree of the taxa using phangorn. We do a multiple sequence alignment for all cleaned contigs present after dada2. Phangorn then builds a tree based on the similarities between sequences, with details here: https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btq706#1989321. What we have in the end is an unrooted tree, stored in dada2/tree_unrooted.tre. See "important data and location" for locations of all the aforementioned files. Additionally, 


#Pipeline
##Receive the raw files. Files come from Argonne in fastq format with an accompanying mapping file. Find raw data in /scratch/users/ctataru5/microbiome/raw. Find samples from...
####Batch 1 : /scratch/users/ctataru5/microbiome/raw/December
####Batch 2 : /scratch/users/ctataru5/microbiome/raw/March
####Batch 3 : /scratch/users/ctataru5/microbiome/raw/April
####Batch 4 : /scratch/users/ctataru5/microbiome/raw/August
####Batch 5 : /scratch/users/ctataru5/microbiome/raw/October
Every folder should contain its mapping file. This will take the form of a .txt document, most likely named something with "maude" in it. These come straight from Argonne.

## Concatenate Mapping Files:
####1. Using text editor(excel is easiest), concatenate all the mapping files together.
####2. run qiime's `validate_mapping_file.py -m mapping.txt` to check that all is well. Might have to add an empty description column at the end, or check that all sample ids are unique, all barcode and linker sequences are present, etc. Check http://qiime.org/scripts/validate_mapping_file.html for details.

## Run Dada2 pipeline on the raw data. See README.md in dada2 folder for more details. This process will leave us with:
#### 1. otu_table.txt
#### 2. a seqs.txt document which is a fasta file detailing all the sequences used in the otu file. This is useful for constructing a phylogenetic tree if the dada2 pipeline fails to do so
#### 3. tax_table.txt which is a sequence by taxonomic classification table
#### 4. tree.tre which is the output of optim.pml from the R package "phangorn". 
*All these outputs can be found in dada2/ folder
*The major reason from using dada2 is error correction and data cleaning. This paper : http://biorxiv.org/content/early/2015/08/06/024034 : argues why this is a good method.

## Alternatively, run qiime
#### 1. Take the seqs.txt file from dada2 output, and run command pick_otus.py -i 'path to seqs.txt' -o 'path to qiime output' -m uclust. In this case, it was pick_otus.py -i $SCRATCH/microbiome/dada2/dada2_output_oct/seqs.fna -o $SCRATCH/microbiome/qiime/picked_otus_oct_open_2/ -m uclust
*If you do this, substitute path to qiime folder every time you see path to dada2


##Plot PCoA:
To view your results in a pretty easy format, run qiime command beta_diversity_through_plots.py.
####1. If using complete dada2 protocol, it will be `beta_diversity_through_plots.py -i dada2/otu_table.biom -m mapping.txt -t dada2/tree.tre -o plots_raw -f`
####2. To view your results, navigate to plots_raw/weighted_unifrac_emperor_pcoa_plot and open index.html. Most likely, the data will look spread out, and in general need of normalization


![PCoA_noNorm](https://github.com/kiti15237/microbiome2/blob/master/figures/PCoA/dada2_noNorm_colPair.png)



##Normalize:
Here you have options. Rarefaction, DeSeq, of CSS. We chose CSS as per : http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531 - Susan Holmes
####1. IF you are using dada2 solely, you will need to convert otu_table.txt to biom format. You can do this using the script scripts/write_otu_biom.R. Simply open R script, change path variable to dada2 folder, and params variables to basename of the file. If you followed above steps, path = dada2/ params = otu_table. 
####2. Call the qiime command `normalize_table.py -i dada2/otu_table.biom -o dada2/otu_table_normCSS.biom`
####3. For convenience, call `biom convert -i otu_table_normCSS.biom -o otu_table_normCSS.txt --to-tsv`


## Replot normalized PCoA:
run qiime command beta_diversity_through_plots.py.
####1.`beta_diversity_through_plots.py -i dada2/otu_table_normCSS.biom -m mapping.txt -t dada2/tree.tre -o plots_raw -f`
####2. To view your results, navigate to plots_raw/weighted_unifrac_emperor_pcoa_plot and open index.html. Most likely, the data will look spread out, and in general need of normalization
Most likely you will see some outliers. These are mostly samples that did not have a lot of mergable reads in the dada2 process, and therefore look very funky on the PCoA plot. 

![PCoA_noNorm](https://github.com/kiti15237/microbiome2/blob/master/figures/PCoA/dada2_normCSS_colPair.png)

##Remove the outliers
####1. Find out what sample ids correspond to the outliers on PCoA plot. Simplest way is to click on samples in the Key tab of emperor plot until the little white arrow that appears points to one of the outliers
####2. create a 'samples_that_suck.txt' document, with one sample id per line
####3. run qiime's `filter_samples_from_otu_table.py -i otu_table.biom -o filtered_otu_table.biom --negate_sample_id_fp --sample_id_fp samplesThatSuck.txt`
####4. Re-normalize, and re-plot using the above 2 steps. Check out your new, filtered, normalized plots. Look good, right? :)

![PCoA_noNorm](https://github.com/kiti15237/microbiome2/blob/master/figures/PCoA/dada2_normCSS_filtered_colPair.png)



`


#####open_oct
Otu tables constructed using qiime: original, normCSS (normalized with cumulative sum scoring) , pca_filtered, etc.


#####Qiime
Lab stuff, from otu picking, etc. TODO: Need to put instructions of how to use qiime for my own notes

#####Dada2
Instructions on how to run
http://biorxiv.org/content/early/2015/08/06/024034











####Results for parents####
1. Run the script resultsForParents.R to produce individual sibling pair excel files
2. Run the script emperor_mapping.R to produce individual sibling pair mapping files that will make personalizing with emperor much easier
3. Copy the new mapping files to sherlock into resultsForParents folder
4. Run the bash script makeEmperorPlots after replacing the appropriate folders. Use a filtered otu table that doesn't have 168.1 and 367 because these were sequenced badly and screw up the PCoA plot
5. Copy the plot folders to a local machine to view them


#Create powerpoint
1. Open old pp to use as a template and delete all graphs and charts
2. copy the first three columns on the excel sheet into pp and resize correctly
3. Select the top and bottom borders on tool bar, change font to 12pt, and bold and right align each of the column headers
4. Return to excel, select all data, and insert chart. I have template. Without template, select stacked bar, switch rows and columns, and copy into pp and adjust size
5. Open appropriate emperor_mapping_sid1_sid2.txt.plot
6. Select color by Description and swap colors so most of the points are blue and two are red. 
7. Select size by age_group. This is an extra column intended to make this step easier. You only select 15 groups instead of 91 samples
8. change background and axis colors, screen shot, copy to paint, grab just the graph, and paste into pp
9. Arrange labels accordingly
10. SAVE
