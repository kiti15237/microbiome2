#####
##Receive the raw files. Files come from Argonne in fastq format with an accompanying mapping file. Find raw data in /scratch/users/ctataru5/microbiome/raw. Find samples from...
####Batch 1 : /scratch/users/ctataru5/microbiome/raw/December
####Batch 2 : /scratch/users/ctataru5/microbiome/raw/March
####Batch 3 : /scratch/users/ctataru5/microbiome/raw/April
####Batch 4 : /scratch/users/ctataru5/microbiome/raw/August
####Batch 5 : /scratch/users/ctataru5/microbiome/raw/October
Every folder should contain its mapping file. This will take the form of a .txt document, most likely named something with "maude" in it. These come straight from Argonne.

## Run Dada2 pipeline on the raw data. See README.md in dada2 folder for more details. This process will leave us with:
#### 1. otu_table.txt
#### 2. a seqs.txt document which is a fasta file detailing all the sequences used in the otu file. This is useful for constructing a phylogenetic tree if the dada2 pipeline fails to do so
#### 3. tax_table.txt which is a sequence by taxonomic classification table
#### 4. tree.rds which is the output of optim.pml from the R package "phangorn". To access to tree itself, use tree = readRDS("outdir/tree.rds")$tree
*All these outputs can be found in dada2/ folder
*The major reason from using dada2 is error correction and data cleaning. This paper : http://biorxiv.org/content/early/2015/08/06/024034 : argues why this is a good method.

## Alternatively, run qiime
#### 1. Take the seqs.txt file from dada2 output, and run command pick_otus.py -i 'path to seqs.txt' -o 'path to qiime output' -m uclust. In this case, it was pick_otus.py -i $SCRATCH/microbiome/dada2/dada2_output_oct/seqs.fna -o $SCRATCH/microbiome/qiime/picked_otus_oct_open_2/ -m uclust
*If you do this, substitute path to qiime folder every time you see path to dada2

##Normalize:
Here you have options. Rarefaction, DeSeq, of CSS. We chose CSS as per : http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531 - Susan Holmes
####1. IF you are using dada2 solely, you will need to convert otu_table.txt to biom format. You can do this using the script scripts/write_otu_biom.R. Simply open R script, change path variable to dada2 folder, and params variables to basename of the file. If you followed above steps, path = dada2/ params = otu_table. 
####2. Call the qiime command `normalize_table.py -i dada2/otu_table.biom -o dada2/otu_table_normCSS.biom`
####3. For convenience, call `biom convert -i otu_table_normCSS.biom -o otu_table_normCSS.txt --to-tsv`



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
