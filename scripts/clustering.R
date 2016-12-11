library(phyloseq)
source("~/Lab/microbiome/scripts/getTables.R")
table <- getOtuTable()
mapping <- getMapping()

#awkward third sibling
table <- table[ , mapping$SampleID != "33"]
mapping <- mapping[ mapping$SampleID != "33", ]
temp <- splitAutControl(table, mapping)
aut <- temp[[1]]
control <- temp[[2]]
autMap <- temp[[3]]
cMap <- temp[[4]]
tree <- getTree()


metadata <- sample_data(data.frame( pair=mapping$Pair, 
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
ps_aut <- phyloseq(otu_table(aut, taxa_are_rows=T), 
                   sample_data(metadata_aut),
                   phy_tree(tree))

ps <- phyloseq(table, 
               sample_data(metadata),
               phy_tree(tree))
ps_aut <- phyloseq(aut, 
                   sample, 
                   phy_tree(tree))


#########Kmeans#############
library(ggplot2)
library(ggfortify)
km <- kmeans(t(aut), centers = 2)
pca <- prcomp(t(aut))
autoplot(km, data = autMap)
autoplot(pca, data = km, colour = "cluster", size = 10, loadings = T)
autoplot(pca, data = autMap, colour = "classifier", size = 10, loadings = T)

km$cluster
for(i in seq(1,max(km$cluster))){
  sids_in_cluster <- names(km$cluster)[ km$cluster == i]
  treatments <- mapping$Treatment[mapping$SampleID %in% sids_in_cluster]
  print(treatments)
  print(sum(treatments == "Aut") - sum(treatments == "Control"))
  print('\n')
}

write.table(aut[,colnames(aut) %in% names(km$cluster)[km$cluster == 1]], "~/Lab/data/autCluster1.txt")
write.table(aut[,colnames(aut) %in% names(km$cluster)[km$cluster == 2]], "~/Lab/data/autCluster2.txt")


##########Hierarchical###############
library(GUniFrac)

ps_use <- ps_aut
dists <- GUniFrac(t(aut), tree)$unifracs
a <- dists[lower.tri(dists)]
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

groups <- cutree(clust, k=2)
write.table(aut[,colnames(aut) %in% names(groups)[groups == 1]], "~/Lab/data/autCluster_hier1.txt")
write.table(aut[,colnames(aut) %in% names(groups)[groups == 2]], "~/Lab/data/autCluster_hier2.txt")

