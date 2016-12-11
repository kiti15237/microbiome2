library(ade4)
library(ggfortify)

source("~/Lab/microbiome2/scripts/getTables.R")
mapping <- getMapping()
table <- getOtuTable()


temp <- splitAutControl(table, mapping)
aut <- temp[[1]]
control <- temp[[2]]
autMap <- temp[[3]]
cMap <- temp[[4]]


differences <- aut - control
#pca where each otu is a vector of samples
pca <- prcomp((differences))
sdev <- pca$sdev
#Take a look at the otu distribution
plot(pca$x[,1], pca$x[,2])

#The idea is we keep everything that is not super congregated in the middle
pc12 = cbind(pca$x[,1], pca$x[,2])
toKeep = pca$x[,1] > 2.2 | pca$x[,1] < -2 | pca$x[,2] > 2 | pca$x[,2] < -2
prunedTable_byDiffs <- table[toKeep, ]

#look at the pca after we've removed the 'extraneous' otus
plot(pca$x[,1][toKeep], pca$x[,2][toKeep])

pca_treatment = prcomp(t(prunedTable_byDiffs))
autoplot(pca_treatment, data = mapping, colour = "Treatment", size = 10, loadings = T)
pca_orig = prcomp(t(table))
autoplot(pca_orig, dat = mapping, colour = "Treatment", size = 10, loadings = T)

write.table(prunedTable_byDiffs , "~/Lab/16S/otuTables/open_oct/otu_table_pcaFiltered_byDiffs.txt", quote = F)



#prunedDifferences <- differences[rownames(differences) %in% rownames(otus_toKeep),]

#write.table(prunedTable_byDiffs , "~/Lab/16S/otuTables/open_oct/otu_table_pcaFiltered_byDiffs.txt", quote = F)

pca_treatment <- prcomp(t(prunedTable_byDiffs))
autoplot(pca_treatment, data = mapping, colour = "Treatment", size = 10, loadings = T)



#pca where each otu is a vector of samples
pca <- prcomp((table))
sdev <- pca$sdev
autoplot(pca, loadings = T)

#looks like the projection onto pc1 should be above 0 to capture the majority of the variability
otus_toKeep <- rownames(pca$x)[pca$x[,1] > 0 ]
pruned_table1 <- table[rownames(table) %in% otus_toKeep , ]
#write.table(pruned_table1 , "~/Lab/data/otuTables_Sept_13.5/otu_table_pcaFiltered.txt", quote = F)

####Round two#####
pca <- prcomp(pruned_table1)
autoplot(pca, loadings = T)

otus_toKeep <- rownames(pca$x)[pca$x[,1] < 0 ]
pruned_table2 <- table[rownames(table) %in% otus_toKeep , ]
#write.table(pruned_table2 , "~/Lab/data/otuTables_sept_13.5/ko_table_pcaFiltered_round2.txt", quote = F)

#pca 
pca3 <- prcomp(t(table))
autoplot(pca3, data = mapping, colour = "Treatment", size = 10, loadings = T)

pca4 <- (prcomp((pruned_table2)))
autoplot(pca4, loadings = T, data = mapping, colour = "Treatment", size = 10)

otu_table_top12Taxa <- read.csv("~/Lab/data/otuTables_Sept_13.8/otu_table_top20Taxa.txt", row.names=1, sep="")
colnames(otu_table_top12Taxa) <- gsub("X", "", colnames(otu_table_top12Taxa))
pca5 <- prcomp(t(otu_table_top12Taxa[-12,]))
autoplot(pca5, data = mapping, colour = "Treatment", size = 10, loadings = T)
