source("~/Lab/scripts/getTables.R")
table <- getOtuTablePCAFiltered()
#table <- data.matrix(table[1:12,])
mapping <- getMapping()

#awkward third sibling
table <- table[ , mapping$SampleID != "33"]
mapping <- mapping[ mapping$SampleID != "33", ]
temp <- splitAutControl(table, mapping)
aut <- temp[[1]]
control <- temp[[2]]
autMap <- temp[[3]]
cMap <- temp[[4]]


aut1 <- getAutClust1Hier()
aut1 <- table[,colnames(table) %in% colnames(aut1)]
autMap1 <- mapping[mapping$SampleID %in% colnames(aut1),]
aut2 <- getAutClust2Hier()
aut2 <- table[,colnames(table) %in% colnames(aut2)]
autMap2 <- mapping[mapping$SampleID %in% colnames(aut2),]

cMap1 <- cMap[cMap$Pair %in% autMap1$Pair,]
con1 <- control[, colnames(control) %in% cMap1$SampleID]
cMap2 <- cMap[cMap$Pair %in% autMap2$Pair,]
con2 <- control[, colnames(control) %in% cMap2$SampleID]

aut1 <- aut1[,order(autMap1$Pair)]
con1 <- con1[,order(cMap1$Pair)]
aut2 <- aut2[,order(autMap2$Pair)]
con2 <- con2[,order(cMap2$Pair)]

aut1 <- aut1[rownames(aut1) %in% rownames(table),]
con1 <- con1[rownames(con1) %in% rownames(table),]
aut2 <- aut2[rownames(aut2) %in% rownames(table),]
con2 <- con2[rownames(con2) %in% rownames(table),]

write.table(aut1, "~/Lab/data/maude/genus/aut_cluster1.txt")
write.table(aut2, "~/Lab/data/maude/genus/aut_cluster2.txt")
write.table(con1, "~/Lab/data/maude/genus/con_cluster1.txt")
write.table(con2, "~/Lab/data/maude/genus/con_cluster2.txt")
write.table(aut, "~/Lab/data/maude/genus/aut.txt")
write.table(control, "~/Lab/data/maude/genus/control.txt")

differences <- t(aut) - t(control)
tests <- apply(differences, 2, t.test)
pvalues <- sapply(tests, function(test) return(test$p.value))
adj <- p.adjust(as.numeric(pvalues), method = "fdr")
adj <- adj[!is.nan(adj)]
print(sum(adj < 0.15))
adj

pvalues <- c()
for(i in seq(1,nrow(aut1))){
  pvalues <- c(pvalues, (wilcox.test(as.numeric(aut2[i,]), as.numeric(con2[i,]))$p.value))
}
adj <- p.adjust(as.numeric(pvalues), method = "fdr")
adj <- adj[!is.nan(adj)]
print(sum(adj < 0.05))

