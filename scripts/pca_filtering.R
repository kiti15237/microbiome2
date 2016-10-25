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

toKeep = pca$x[,1] > 2.2 | pca$x[,1] < -2 | pca$x[,2] > 2 | pca$x[,2] < -2
prunedTable_byDiffs <- table[toKeep, ]

#look at the pca after we've removed the 'extraneous' otus
plot(pca$x[,1][toKeep], pca$x[,2][toKeep])

pca_pruned= prcomp(t(prunedTable_byDiffs))
autoplot(pca_pruned, data = mapping, colour = "Treatment", size = 10, loadings = T)

pca_orig = prcomp(t(table))
autoplot(pca_orig, data = mapping, colour = "Treatment", size = 10, loadings = T)

write.table(prunedTable_byDiffs , "~/Lab/16S/otuTables/open_oct/otu_table_pcaFiltered_byDiffs.txt", quote = F)