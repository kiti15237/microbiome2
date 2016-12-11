source("~/Lab/microbiome2/scripts/getTables.R")
mapping <- getMapping()
temp <- getPcaFilteredOtuTable()
table = temp[[1]]
mapping = temp[[2]]
#awkward third sibling
table <- table[ , mapping$SampleID != "33"]
mapping <- mapping[ mapping$SampleID != "33", ]
temp <- splitAutControl(table, mapping)
aut <- temp[[1]]
control <- temp[[2]]
autMap <- temp[[3]]
cMap <- temp[[4]]



numTaxaAut <- apply(aut, 2, function(sample) return(sum(sample != 0)))
numTaxaControl <- apply(control, 2, function(sample) return(sum(sample!=0)))
numTaxa_pvalue <-wilcox.test(numTaxaAut, numTaxaControl, paired=TRUE)$p.value
print(paste("Number of taxa present p value: ", numTaxa_pvalue))
hist(numTaxaAut, col = rgb(1,0,0,0.5), breaks=8, main="Number of taxa in ASD samples", xlab = "# of Taxa")
hist(numTaxaControl, col = rgb(0,0,1,0.5), breaks = 8, main="Number of taxa in Control samples", xlab = "# of Taxa")

#Rare if they appear is fewer than threshold samples
threshold <- 12
rare_taxa <- apply(table, 1, function(taxa) return(sum(taxa > 0) < threshold )) # counts are greater than 0 in fewer than threshold samples
abundance_rare_perAutSample <- apply(aut[rare_taxa, ], 2, function(sample) return(sum(sample)))#for every sample, go through every taxa we've identified as rare and sum their abundances 
abundance_rare_perConSample <- apply(control[rare_taxa, ], 2, function(sample) return(sum(sample)))
rare_pvalue <- t.test(abundance_rare_perAutSample, abundance_rare_perConSample, paired=T, alternative="less")$p.value
t.test(abundance_rare_perAutSample - abundance_rare_perConSample)$p.value
#apparently a paired t test is the sample as checking if their difference is sig. different from 0, who knew?
print(paste("Abundance of rare taxa p value: ", rare_pvalue))
#There are generally fewer rare taxa in autism samples as compared to their sibling controls


library(vegan)
diversity_measure <- "shannon"
autDiversity <- apply(aut, 2, function(sample) return(diversity(sample, index = diversity_measure)))
controlDiversity <- apply(control, 2, function(sample) return(diversity(sample, index = diversity_measure)))
div_p <- t.test(autDiversity, controlDiversity, paired = T)$p.value
print(paste("Diversity p value: ", div_p))
hist(autDiversity, col = rgb(1,0,0,0.5), breaks=8, main="Shannon Diversity in ASD samples", xlab = "diversity", xlim=c(3.4, 5.5))
hist(controlDiversity, col = rgb(0,0,1,0.5), breaks = 8, main="Shannon Diversity in ASD samples", xlab = "diversity", xlim=c(3.4,5.5))




#By pairi
#check commented out filtering above. Took out the groups that were appearing in just 1 or 2 sample, but nothing
differences <- aut- control
tests <- apply(differences, 1, t.test)
p_diffs <- sapply(tests, function(test) return(test$p.value))
p_adj <- p.adjust(p_diffs, method="fdr")
hist(p_adj)
print(as.numeric(p_adj))


#print out rare_taxa seqtab
#write.table(table[rare_taxa, ], "~/Lab/data/dada2/output/0mismatch_25overlap/seqtab_otuId_rareOnly_thresh3.txt")




############Specific Genera Tests and general diverisy test based on Kang et. al
tax <- read.delim("~/Lab/16S/otuTables/open_oct/rep_set_tax_assignments.txt", header=FALSE)
rownames(table) <- tax$V2[match(rownames(table), tax$V1)]

prev_control <- control[grepl("g__Prevotella", rownames(table)), ]
prev_aut <- aut[grepl("g__Prevotella", rownames(table)), ]
ks.test(colSums(prev_aut), colSums(prev_control), paired = T)$p.value
plot(density(prev_aut), col="red")
lines(density(prev_control), col = "blue")


cop_control <- control[grepl("g__Coproco", rownames(table)), ]
cop_aut <- aut[grepl("g__Coproco", rownames(table)), ]
shapiro.test(cop_aut)
t.test(colSums(cop_aut), colSums(cop_control), paired = T)$p.value
wilcox.test(cop_aut, cop_control, paired= F)$p.value
plot(density(cop_aut), col="red")
lines(density(cop_control), col = "blue")


veil_control <- colSums(control[grepl("f__Veillonellaceae", rownames(table)), ])
veil_aut <- colSums(aut[grepl("f__Veillonellaceae", rownames(table)), ])
shapiro.test(veil_aut)
t.test(veil_aut, veil_control, paired = T)$p.value
wilcox.test(veil_aut, veil_control, paired= F)$p.value
plot(density(veil_aut), col="red")
lines(density(veil_control), col = "blue")

bac_control <- control[grepl("g__Bacteroides", rownames(table)), ]
bac_aut <- aut[grepl("g__Bacteroides", rownames(table)), ]
shapiro.test(bac_aut)
t.test(bac_aut, bac_control, paired = F, alternative = "less")$p.value
wilcox.test(bac_aut, bac_control, paired= F)$p.value
plot(density(bac_aut), col="red")
lines(density(bac_control), col = "blue")

lac_control <- control[grepl("g__Lactobacillus", rownames(table)), ]
lac_aut <- aut[grepl("g__Lactobacillus", rownames(table)), ]
shapiro.test(lac_aut)
t.test(lac_aut, lac_control, paired = F, alternative = "less")$p.value
wilcox.test(lac_aut, lac_control, paired= F)$p.value
plot(density(lac_aut), col="red")
lines(density(lac_control), col = "blue")

prev_to_bac_aut <- prev_aut / bac_aut
prev_to_bac_control <- prev_control / bac_control
t.test(prev_to_bac_aut, prev_to_bac_control, paired = T)$p.value
plot(density(prev_to_bac_aut), col="red")
lines(density(prev_to_bac_control), col = "blue")

boxplot(log(prev_control, 10), log(prev_aut, 10), horizontal = T, col=c("blue", "red"), range=0)
boxplot(log(lac_control, 10), log(lac_aut, 10), horizontal = T, col=c("blue", "red"))
boxplot(log(cop_control, 10), log(cop_aut, 10), horizontal = T, col=c("blue", "red"))
boxplot(log(cop_control, 10), log(cop_aut, 10), horizontal = T, col=c("blue", "red"))
boxplot(log(cop_control, 10), log(cop_aut, 10), horizontal = T, col=c("blue", "red"))

library(asbio)
index = "shan"
alpha_div_aut <- apply(aut, 2, function(col) return( alpha.div(col, index)))
alpha_div_control <- apply(control, 2, function(col) return( alpha.div(col, index)))
t.test(alpha_div_aut, alpha_div_control, paired = T)$p.value
#diversity_p <- t.test(diversity(t(aut)), diversity(t(control)), paired=T, alternative = "less")$p.value
#print(diversity_p)

#pairwise pearson correlation tests on only neurotypical samples
library(igraph)
for(i in 1:nrow(table)){
  corrs_control <- sapply(seq(1,nrow(control)), function(j) return(cor.test(control[i,], control[j,], type="pearson")$statistic))
  corrs_aut <- sapply(seq(1,nrow(aut)), function(j) return(cor.test(aut[i,], aut[j,], type="pearson")$statistic))
}
corrs_control[is.na(corrs_control)] <- 0
corrs_aut[is.na(corrs_aut)] <- 0
adjacency_control <- matrix(corrs_control, nrow = nrow(table), ncol = nrow(table))
#adjacency_control[abs(adjacency_control) < 0.8] <- 0
#adjacency_control [is.na(adjacency_control)] <- 0
adjacency_aut <- matrix(corrs_aut, nrow = nrow(table), ncol = nrow(table))
#adjacency_aut[abs(adjacency_aut) < 0.8] <- 0
#adjacency_aut[is.na(adjacency_aut)] <- 0

sum(corrs_control > 0.8 & corrs_aut < 0.3)
sum(corrs_aut > 0.8 & corrs_control < 0.3)

g <- graph.adjacency(adjacency_control, weighted=T)


edgeList <- get.data.frame(g)
colnames(edgeList) <- c("Source", "Target", "Weight")
#edgeList_pruned <- edgeList[abs(edgeList$Weight) > 0.75,]
write.csv(edgeList, "~/Lab/data/dada2/output/0mismatch_25overlap/corr_edgeList_control_0.8.csv", row.names=F)
nodeList<- data.frame(Id = rownames(table))
write.csv(nodeList, "~/Lab/data/dada2/output/0mismatch_25overlap/corr_nodeList_control_0.8.csv", row.names=F)
