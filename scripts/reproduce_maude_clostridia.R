library(phyloseq)
source("~/Lab/microbiome2/scripts/getTables.R")
mapping <- getMapping()
table <- getOtuTableL5()

temp <- splitAutControl(table, mapping)
aut <- temp[[1]]
control <- temp[[2]]
autMap <- temp[[3]]
cMap <- temp[[4]]
tree <- getTree()

metadata_aut <- sample_data(data.frame(pair = autMap$Pair,
                                       classifier = autMap$classifier,
                                       age = autMap$age_month_ok, 
                                       batch = autMap$batch, 
                                       treatment = autMap$Treatment))
rownames(metadata_aut) <- autMap$SampleID
ps_aut <- phyloseq(aut, 
                   sample)



distance_measure = "wunifrac"
clust <- hclust(phyloseq::distance(ps_use, method = distance_measure))
plot(clust)

library(dendextend)
dend <- as.dendrogram(clust)
classifier <- autMap$classifier[match(labels(dend), autMap$SampleID)]
treatments <- mapping$Treatment[match(labels(dend), mapping$SampleID)]
labels_colors(dend) <- ifelse(classifier < -6, "red", "blue")
dend <- hang.dendrogram(dend, hang = 0.3)
plot(dend)


#The difference between pairs should be much closer than the general distances aut-control distances
diffall = abs(aut - control)
mag_differences <- (apply(diffall, 2, sum))
max(mag_differences)
min(mag_differences)
plot(density(mag_differences))
cids <- cMap$SampleID[cMap$Pair %in% autMap$Pair[autMap$SampleID %in% colnames(sorted_by_differences)]]
sorted_by_differences <- diffall[ , order(mag_differences)]
cids_sorted <- cids[order(mag_differences)]

#let's take the 20 most variance:
sd(diffall[1,])
sfdiffall<-apply(diffall, MARGIN = 1, sd)
select<-diffall[sfdiffall > .01,]



#let's cluster
#install.packages("pheatmap")
library(pheatmap)
#heatmap 
colfunc <- colorRampPalette(c("white", "black"))
#ok we're doing the clustering on the euclidean distances between abundance vectors
heatmap_test<-pheatmap(t(diffall), clustering_method ="ward.D2")
heatmap_test<-pheatmap(select, clustering_method ="ward.D", )

#diffall.clust <- cutree(heatmap_test, h=2)
diffall.clust <- data.frame( cbind(t(diffall), cluster = cutree(heatmap_test$tree_col, k = 2)))

cluster_1.6 = rownames(diffall.clust)[diffall.clust$cluster == 1]
cluster_2.6 = rownames(diffall.clust)[diffall.clust$cluster == 2]


clust1_con <- intersect(intersect(cluster_1.3, cluster_1.4), cluster_1.5)
clust2_con <- intersect(intersect(cluster_2.3, cluster_2.4), cluster_2.5)

