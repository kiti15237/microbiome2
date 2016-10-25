#####

#####open_oct
Otu tables of all sorts: original, normCSS (normalized with cumulative sum scoring) , pca_filtered, etc.


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
