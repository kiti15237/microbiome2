#table <- read.csv("~/Lab/data/dada2/output/0mismatch_25overlap/otu_kpcof.txt", row.names=1, sep="")
library(phyloseq)
data <- import_biom("~/Lab/data/dada2/output/0mismatch_25overlap/qiime/closed_13.8/otu_table.biom")
table_agp <-  read.delim("~/Lab/Secure__file_for_results/ag_under_8_otu_table_even_10k.txt", row.names=1, comment.char="#")
mapping <- read.table("~/Lab/data/dada2/mapping_dada2.txt", comment.char="", header=T)
mapping_agp <- read.delim("~/Lab/Secure__file_for_results/ag_under_8_mapping.txt", row.names=1)

#Necessary formatting
table <- otu_table(data)
colnames(table)<-gsub("X","",colnames(table))
mapping$SampleID <- as.character(mapping$SampleID)
table <- table[, colnames(table)!="168"]
mapping <- mapping[mapping$SampleID != "168",]
colnames(table_agp)<-gsub("X","",colnames(table_agp))

#order
mapping <- mapping[order(mapping$SampleID),]
table <- table[,order(colnames(table))]

table <- apply(table, 2, function(col) return(col / sum(col)))
table <- otu_table(table, taxa_are_rows=T)


table_agp <- apply(table_agp, 2, function(col) return(col / sum(col)))
table_agp <- otu_table(table_agp, taxa_are_rows=T)

mapping$classifier[is.na(mapping$classifier)] <-    4


aut <- table[ , mapping$Treatment == "Aut"]
control <- table[ , mapping$Treatment == "Control"]
autMap <- mapping[mapping$Treatment == "Aut",]
cMap <- mapping[mapping$Treatment == "Control",]


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

tree = read_tree_greengenes("~/Lab/data/dada2/output/0mismatch_25overlap/qiime/closed_13.8/97_otus_tree.tree")

metadata <- sample_data(data.frame( pair=mapping$Pair, 
                                    classifier=mapping$classifier,
                                    age=mapping$age_month_ok, 
                                    batch=factor(mapping$batch),
                                    treatment= factor(mapping$Treatment, levels = c("Aut", "Control"))))
rownames(metadata) <- mapping$SampleID

metadata_agp <- sample_data(data.frame(pair = rep(NA, nrow(mapping_agp)), 
                                       classifier = rep(NA, nrow(mapping_agp)),
                                       age = mapping_agp$AGE,
                                       batch = rep(NA, nrow(mapping_agp)), 
                                       treatment = rep("AGP", nrow(mapping_agp))))
rownames(metadata_agp) <- rownames(mapping_agp)

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

ps <- phyloseq(table, 
               sample_data(metadata))

ps_aut <- phyloseq(otu_table(aut, taxa_are_rows=T), 
                   sample_data(metadata_aut))

ps_control <- phyloseq(otu_table(control, taxa_are_rows = T),
                       sample_data(metadata_control))

ps_agp <- phyloseq(table_agp,
                   sample_data(metadata_agp))

ps_aut_agp <- merge_phyloseq(ps_aut, ps_agp, phy_tree(tree))
ps_control_agp <- merge_phyloseq(ps_control, ps_agp, phy_tree(tree))

#ig = make_network(ps_aut, type = "samples", distance = "bray", max.dist = 0.75)
#plot_network(ig, ps_aut, color = "treatment", line_weight = 0.4, 
#             label = NULL)


#aut_nt_distances <- sapply(seq(1,ncol(control)), function(j) return(distance(aut[,1], control[,j], method = "jaccard")))
#for(i in seq(2,ncol(aut))){
#  print(i)
#  distances <- sapply(seq(1,ncol(control)), function(j) return(jaccard(aut[,i]* control[,j])))
#  print(distances)
#  aut_nt_distances <- rbind(aut_nt_distances, distances)
#}
#colnames(aut_nt_distances) <- colnames(control)
#rownames(aut_nt_distances) <- colnames(aut)

dists_aut_agp <- as.matrix(distance(ps_aut_agp, method="dpcoa"))
treatment_pairs <- expand.grid(sample_data(ps_aut_agp)$treatment, sample_data(ps_aut_agp)$treatme)
dists_aut_agp_list <- as.vector(dists_aut_agp)
treatment_pairs <- treatment_pairs[order(dists_aut_agp_list, decreasing = T),]
treatment_pairs <- paste(treatment_pairs[,1], treatment_pairs[,2])
dists_aut_agp_list <- dists_aut_agp_list[order(dists_aut_agp_list, decreasing = T)]
dists_aut_agp_plot <- dists_aut_agp_list[treatment_pairs == "AGP Aut" | treatment_pairs == "Aut AGP"]

dists_control_agp <- as.matrix(distance(ps_control_agp, method="dpcoa"))
treatment_pairs <- expand.grid(sample_data(ps_control_agp)$treatment, sample_data(ps_control_agp)$treatme)
dists_control_agp_list <- as.vector(dists_control_agp)
treatment_pairs <- treatment_pairs[order(dists_control_agp_list, decreasing = T),]
treatment_pairs <- paste(treatment_pairs[,1], treatment_pairs[,2])
dists_control_agp_list <- dists_control_agp_list[order(dists_control_agp_list, decreasing = T)]
dists_control_agp_plot <- dists_control_agp_list[treatment_pairs == "AGP Control" | treatment_pairs == "Control AGP"]



par(xpd=TRUE, mar = c(6,4,6,3))
plot(x = seq(1,length(dists_aut_agp_plot)), y=dists_aut_agp_plot,
     col = "red", ylab="Distance", xlab="pair ids")
points(x = seq(1,length(dists_control_agp_plot)), y = dists_control_agp_plot, col="blue")
p <- t.test(dists_aut_agp_plot, dists_control_agp_plot)$p.value
title(main = "Distance btw. aut/agp vs. control/agp samples", cex.main=0.8,
      sub = paste("P value aut-agp vs. control-agp: ", p), cex.sub = 0.75, line=4.5)
legend(x=-275, y=0.87, 
       c("Distances btw Aut/AGP samples", "Distances btw Control/AGP samples"),
       col = c('red', 'blue'), pch=c(1,1,1), cex=0.75)


#######Relevant phyloseq analysis#####
dists <- as.matrix(distance(ps_merged, method = "dpcoa"))
pairs <- expand.grid(rownames(sample_data(ps_merged)), rownames(sample_data(ps_merged)))
treatment_pairs <- expand.grid(sample_data(ps_merged)$treatment, sample_data(ps_merged)$treatment)
dists_list <- as.vector(dists)
pairs <- pairs[order(dists_list, decreasing = T), ]
treatment_pairs <- treatment_pairs[order(dists_list, decreasing = T),]
dists_list <- dists_list[order(dists_list, decreasing = T)]

treatment_pairs <- paste(treatment_pairs[,1], treatment_pairs[,2])
pairs <- paste(pairs[,1], pairs[,2])

dists_aut <- dists_list[treatment_pairs=="Aut Aut"]
pairs_aut <- pairs[treatment_pairs=="Aut Aut"]

dists_mixed <- dists_list[treatment_pairs=="Aut Control" | treatment_pairs == "Control Aut"]
pairs_mixed <- pairs[treatment_pairs=="Aut Control" | treatment_pairs == "Control Aut"]

dists_control <- dists_list[treatment_pairs=="Control Control"]
pairs_control <- pairs[treatment_pairs=="Control Control"]


p1 <- t.test(dists_aut, dists_control, paired = F)$p.value
p2 <- t.test(dists_aut, dists_mixed, paired = F)$p.value
par(xpd=TRUE, mar = c(10,4,7,3))
plot(x=seq(1,1000.5, by=1000/length(pairs_mixed)), y=dists_mixed, ylab="Distances", xlab="pair ids")
title(main = "Pairwise Distances by condition", line = 5.5,
      sub = c(paste("Pvalue aut-aut vs. control-control: ", p1, '\n', "Pvalue aut-aut vs. aut-control: ", p2)))
points(x = seq(1,1000, by = 1000/length(pairs_aut)), y = dists_aut, col='red')
points(x= seq(1,1000, by = 1000/length(pairs_control)), y = dists_control, col='blue')
legend(x=300, y=0.87, c("Distances btw Aut/Control samples", "Distances btw Aut samples", "Distances btw Control samples"),col = c('black','red', 'blue'), pch=c(1,1,1))






#par(xpd=TRUE, mar = c(3,4,7,3))
#plot(x= rep(1, length(dists_aut)), y = dists_aut, col= 'red', xlim=c(0.5,2.5), ylab = "distances", xlab = "")
#title("Distances between same condition pairwise", line = 4)
#points(x = rep(2, length(dists_control)), y=dists_control, col = 'blue')
#legend(x=1, y=1.047, c("Distances btw Aut samples", "Distances btw Control samples"),col = c('red', 'blue'), pch=c(1,1))
