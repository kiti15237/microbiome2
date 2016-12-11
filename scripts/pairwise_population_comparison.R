
library(phyloseq)
tree = read_tree_greengenes("~/Lab/16S/dada2/tree.tre")
tax <- read.csv("~/Lab/16S/dada2/taxa_table.txt", row.names=1, sep="")
#tree = read_tree_greengenes("~/Lab/16S/otuTables/open_oct/rep_set.tre")
source("~/Lab/microbiome2/scripts/getTables.R")
mapping <- getMapping()
#table <- getOtuTable()
table<- getOtuTable()
#table <- getPcaFilteredOtuTable()

table <- otu_table(table, taxa_are_rows = T)
#table <- table[,colnames(table) != "210"]
#mapping <- mapping[mapping$SampleID != "210", ]
#otus <- getOtusGreatestMI(table, mapping, 300)
#table <- table[rownames(table) %in% otus, ]


mapping$classifier[is.na(mapping$classifier)] <-    4


temp <- splitAutControl(table, mapping)
aut <- temp[[1]]
control <- temp[[2]]
autMap <- temp[[3]]
cMap <- temp[[4]]



autLow <- table[, mapping$classifier < -6]
autLowMapping <- mapping[mapping$classifier < -6, ]
autHigh <- table[, mapping$classifier > -6 & mapping$classifier < 0]
autHighMapping <- mapping[ mapping$classifier > -6 & mapping$classifier < 0,]

#Pair 2 has two autism samples and messes everything up
autHigh <- autHigh[,autHighMapping$SampleID!= 33.0]
autHighMapping <- autHighMapping[autHighMapping$SampleID!=33.0,]

cLow <- table[, (mapping$Pair %in% autLowMapping$Pair) & mapping$Treatment == "Control"]
cLowMapping <- mapping[(mapping$Pair %in% autLowMapping$Pair) & mapping$Treatment == "Control",]
cHigh <- table[,(mapping$Pair %in% autHighMapping$Pair) & mapping$Treatment == "Control"]
cHighMapping <- mapping[ (mapping$Pair %in% autHighMapping$Pair) & mapping$Treatment == "Control",]



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

metadata_control <- sample_data(data.frame(pair = cMap$Pair,
                                       classifier = cMap$classifier,
                                       age = cMap$age_month_ok, 
                                       batch = cMap$batch, 
                                       treatment = cMap$Treatment))
rownames(metadata_control) <- cMap$SampleID

ps <- phyloseq(otu_table(otu_table, taxa_are_rows=T), 
               sample_data(metadata),
               phy_tree(tree), 
               tax_table(tax))

ps_aut <- phyloseq(otu_table(aut, taxa_are_rows=T), 
                   sample_data(metadata_aut),
                   phy_tree(tree))

ps_control <- phyloseq(otu_table(control, taxa_are_rows = T),
                       sample_data(metadata_control))




##################### AGP STUFF ##################

getAGPObjects <- function(){
  mapping_agp <- read.delim("~/Lab/Secure__file_for_results/ag_under_8_mapping.txt", row.names=1)
  table_agp <-  read.delim("~/Lab/Secure__file_for_results/ag_under_8_otu_table_even_10k.txt", row.names=1, comment.char="#")
  colnames(table_agp)<-gsub("X","",colnames(table_agp))
  table_agp <- apply(table_agp, 2, function(col) return(col / sum(col)))
  table_agp <- otu_table(table_agp, taxa_are_rows=T)
  metadata_agp <- sample_data(data.frame(pair = rep(NA, nrow(mapping_agp)), 
                                         classifier = rep(NA, nrow(mapping_agp)),
                                         age = mapping_agp$AGE,
                                         batch = rep(NA, nrow(mapping_agp)), 
                                         treatment = rep("AGP", nrow(mapping_agp))))
  rownames(metadata_agp) <- rownames(mapping_agp)
  ps_agp <- phyloseq(table_agp,
                     sample_data(metadata_agp))
  
  ps_aut_agp <- merge_phyloseq(ps_aut, ps_agp, phy_tree(tree))
  ps_control_agp <- merge_phyloseq(ps_control, ps_agp, phy_tree(tree))
  return(list(ps_aut_agp, ps_control_agp))
  
}
plotAGPComps <- function(ps_aut_agp, ps_control_agp, distance_measure){
    
    dists_aut_agp <- as.matrix(distance(ps_aut_agp, method= distance_measure))
    treatment_pairs1 <- expand.grid(sample_data(ps_aut_agp)$treatment, sample_data(ps_aut_agp)$treatment)
    dists_aut_agp_list <- as.vector(dists_aut_agp)
    treatment_pairs1 <- treatment_pairs1[order(dists_aut_agp_list, decreasing = T),]
    treatment_pairs1 <- paste(treatment_pairs1[,1], treatment_pairs1[,2])
    dists_aut_agp_list <- dists_aut_agp_list[order(dists_aut_agp_list, decreasing = T)]
    dists_aut_agp_plot <- dists_aut_agp_list[treatment_pairs1 == "AGP Aut" ]
    
    dists_control_agp <- as.matrix(distance(ps_control_agp, method= distance_measure))
    treatment_pairs2 <- expand.grid(sample_data(ps_control_agp)$treatment, sample_data(ps_control_agp)$treatme)
    dists_control_agp_list <- as.vector(dists_control_agp)
    treatment_pairs2 <- treatment_pairs2[order(dists_control_agp_list, decreasing = T),]
    treatment_pairs2 <- paste(treatment_pairs2[,1], treatment_pairs2[,2])
    dists_control_agp_list <- dists_control_agp_list[order(dists_control_agp_list, decreasing = T)]
    dists_control_agp_plot <- dists_control_agp_list[treatment_pairs2 == "AGP Control" ]
    return(list(dists_aut_agp_plot, dists_control_agp_plot))
    ###Plot agp stuff
#     par(xpd=TRUE, mar = c(6,4,6,3))
#     plot(x = seq(1,length(dists_aut_agp_plot)), y=dists_aut_agp_plot,
#        col = "red", ylab="Distance", xlab="pair ids")
#     points(x = seq(1,length(dists_control_agp_plot)), y = dists_control_agp_plot, col="blue")
#     p <- t.test(dists_aut_agp_plot, dists_control_agp_plot)$p.value
#     title(main = "Distance btw. aut/agp vs. control/agp samples", cex.main=0.8,
#         sub = paste("P value aut-agp vs. control-agp: ", p), cex.sub = 0.75, line=4.5)
#     legend(x=-275, y=0.87, 
#          c("Distances btw Aut/AGP samples", "Distances btw Control/AGP samples"),
#          col = c('red', 'blue'), pch=c(1,1,1), cex=0.75)
}
#tempPss <- getAGPObjects()
#ps_aut_agp <- tempPss[[1]]
#ps_control_agp <- tempPss[[2]]
#temp <- plotAGPComps(ps_aut_agp, ps_control_agp, distance_measure = distance_measure )
#dists_aut_agp_plot <- temp[[1]]
#dists_control_agp_plot <- temp[[2]]
#######Relevant phyloseq analysis#####

ps_use <- readRDS("~/Lab/16S/dada2/phyloseqObj.rds")
distance_measure <- "wunifrac"
dists <- as.matrix(phyloseq::distance(ps_use, method = distance_measure))
pairIds <- expand.grid(sample_data(ps_use)$pair, sample_data(ps_use)$pair)
pairs <- expand.grid(rownames(sample_data(ps_use)), rownames(sample_data(ps_use)))
treatment_pairs3 <- expand.grid(sample_data(ps_use)$treatment, sample_data(ps_use)$treatment)
dists_list <- as.vector(dists)
pairIds <- pairIds[order(dists_list, decreasing=T),]
pairs <- pairs[order(dists_list, decreasing = T), ]
treatment_pairs3 <- treatment_pairs3[order(dists_list, decreasing = T),]
dists_list <- dists_list[order(dists_list, decreasing = T)]

treatment_pairs3 <- paste(treatment_pairs3[,1], treatment_pairs3[,2])
pairs <- paste(pairs[,1], pairs[,2])
#pairIds <- paste(pairIds[,1], pairIrds[,2])

dists_aut <- dists_list[treatment_pairs3=="Aut Aut"]
pairs_aut <- pairs[treatment_pairs3=="Aut Aut"]

dists_mixed <- dists_list[treatment_pairs3=="Aut Control" ]
pairs_mixed <- pairs[treatment_pairs3=="Aut Control" ]

dists_control <- dists_list[treatment_pairs3=="Control Control"]
pairs_control <- pairs[treatment_pairs3=="Control Control"]

dists_pairs <- dists_list[pairIds[,1] == pairIds[,2] & treatment_pairs3 == "Aut Control"]

#Here, we want to verify that the distances between pairs is significantly closer than we would expect from the average population (aut-control comparisons)
#hist(sample(dists_mixed, 100000, replace = T) - mean(dists_mixed), col = rgb(1,0,0,0.5), breaks = 20)
#hist(sample(dists_pairs, 100000, replace = T) - mean(dists_pairs), col = rgb(0,0,1,0.5), add = T, breaks = 20)

##First we want to test that the shape of the distributions is similar: we get a very high p value meaning that we cannot reject the null hypothesis
#and the shapes are decently similar
#ks.test(dists_mixed - mean(dists_mixed), dists_pairs - mean(dists_pairs))
#then, we want to check if the means are different
#ks.test(dists_mixed, dists_pairs) 


aut_top <- dists_aut[dists_aut != 0]
control_top <- dists_control[dists_control != 0]
mixed_top <- dists_mixed
p1 <- wilcox.test(aut_top, control_top, paired = F)$p.value
p2 <- ks.test(aut_top, mixed_top, paired = F)$p.value


hist(sample(aut_top, length(dists_control), replace = F), col = rgb(1,0,0,0.5), breaks = 100,
     main="Unifrac Distances within same phenotype groups \n PCA filtered data", xlab="Weighted Unifrac Distance")
hist(control_top, col = rgb(0,0,1,0.5), add=T, breaks = 100,  xlab="Weighted Unifrac Distance")





par(xpd=TRUE, mar = c(5,6,7,3))
plot(density(control_top), col = "blue", main= "Unifrac Distances within same phenotype groups \n PCA filtered data", xlab = "Weighted Unifrac Distance")
lines(density(sample(aut_top, length(dists_control))), col = "red")



################
distance_measure = "wunifrac"
ps_aut <- subset_samples(ps, treatment == "Aut")
ps_control <- subset_samples(ps, treatment == "Control")
dists_aut <- phyloseq::distance(ps_aut, method = distance_measure)
dists_control <- phyloseq::distance(ps_control, method = distance_measure)
dists <-as.matrix(phyloseq::distance(ps, method = distance_measure))
dists_aut_control <- as.vector(dists[sample_data(ps)$treatment=="Aut", sample_data(ps)$treatment=="Control"])

par(xpd=TRUE, mar = c(5,6,7,3))
plot(density(as.vector(dists_control)), col = "blue", 
     main= "Wunifrac Distances within same phenotype groups \n PCA filtered data", xlab = "Weighted Unifrac Distance",
     xlim=c(0, 0.8))
lines(density(as.vector(dists_aut)), col = "red")
lines(density(dists_aut_control), col = "purple")

hist(dists_control, col = rgb(0,0,1,0.5), breaks = 100, 
     main = "Weighted Unifrac Distances Between Control Samples and Autism Samples",
    xlab = "Weighted Unifrac Distances", xlim = c(0, 0.6))
hist(dists_aut, col=rgb(1,0,0,0.5), add=T, breaks=100)


#SET ROOT IN PHYLOGENETIC TREE
