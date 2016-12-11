# Install fpc which contains the method clusterboot()
#install.packages("fpc", repos = "http://cran.cnr.berkeley.edu/")
library(fpc)
library(phyloseq)
library(ape)
library(phytools)

source("~/Lab/microbiome2/scripts/getTables.R")
mapping <- getMapping()
table <- getPcaFilteredOtuTable()
tree = read.tree("~/Lab/16S/otuTables/open_oct/rep_set.tre")

temp <- splitAutControl(table, mapping)
aut <- temp[[1]]
control <- temp[[2]]
autMap <- temp[[3]]
cMap <- temp[[4]]


metadata <- sample_data(data.frame( SampleID = mapping$SampleID,
                                    pair=mapping$Pair, 
                                    classifier=mapping$classifier,
                                    age=mapping$age_month_ok, 
                                    batch=factor(mapping$batch),
                                    treatment= factor(mapping$Treatment, levels = c("Aut", "Control"))))
rownames(metadata) <- mapping$SampleID
metadata_aut <- sample_data(data.frame(pair = autMap$Pair,
                                       classifier = autMap$classifier,
                                       age = autMap$age_month_ok, 
                                       batch = autMap$batch, 
                                       treatment = autMap$Treatment))
rownames(metadata_aut) <- autMap$SampleID

ps <- phyloseq(otu_table(table, taxa_are_rows = T), 
               sample_data(metadata),
               phy_tree(tree))

ps_aut <- phyloseq(otu_table(aut, taxa_are_rows=T), 
                   sample_data(metadata_aut),
                   phy_tree(tree))



distance_measure <- "wunifrac"
dists <- as.matrix(distance(ps_aut, method = distance_measure))
k = 3

# Clustering algorithm; to be edited such that the clustering algorithm can be input from cmd line.
clusterboot_kmeans <- clusterboot(dists, distances = T, bootmethod="boot", clustermethod=kmeansCBI, krange= k)
#clusterboot_kmeans <- clusterboot(, bootmethod="boot", clustermethod=kmeansCBI, krange= 4)

# Print out the clusters
clusters <- clusterboot_kmeans$result$partition

print (clusters)

clusters_treatment = c(0,0,0)
#Iterate through the cluster assignments; if the patient's diagnosis contains ASD-related terms, add to ASD count for the cluster.
for (i in 1:length(clusters)){
  if(mapping$Treatment[i] == 'Aut'){
    clusters_treatment[clusters[i]] = clusters_treatment[clusters[i]] + 1
  }

}
#cluster_and_pair = rbind(clusters, mapping$Pair)
#siblings_together <- sum(duplicated(cluster_and_pair, MARGIN=2 )) / length(clusters)

autMap$classifier[clusters == 1] 
autMap$classifier[clusters == 2]


# Print the stability result for each cluster (a score closer to 1 represents high stability, a score closer to 0 represents low stability)
clusterboot_kmeans$bootmean

# Print out the number of times each cluster was dissolved
clusterboot_kmeans$bootbrd

# % of people in group who were clustered with their sibling
#interesting that we don't generally have siblings clustering together
#for (id in 1:k){
#  print(sum(duplicated(mapping$Pair[mapping$SampleID %in% names(clusters[clusters==id])])) / length(clusters[clusters==id]))
#}
#print("\n")
#print("\n")
#Are phenotypes enriched in any 1 cluster...no...
#for (id in 1:k){
#  print(sum(mapping$Treatment[mapping$SampleID %in% names(clusters[clusters==id])] == "Aut") / length(clusters[clusters==id]))
#}

#See some very vague enrichment of more severe cases in one cluster, less severe in the other two. It's a strech though...
for (id in 1:k){
  print(sum(autMap$classifier[autMap$SampleID %in% names(clusters[clusters==id])] < -5) / length(clusters[clusters==id]))
}
