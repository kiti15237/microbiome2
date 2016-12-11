source("~/Lab/scripts/getTables.R")
#let's work with differences
#aut<-read.table("/Users/mmdavid/Documents/argonne_data/microbiome_analysis1_19samples/differences/genus/aut.txt", header = T)
#controls<-read.table("/Users/mmdavid/Documents/argonne_data/microbiome_analysis1_19samples/differences/genus/control.txt", header = T)
aut<-read.table("~/Lab/data/maude/class/aut.txt", header = T)
controls<-read.table("~/Lab/data/maude/class/control.txt", header = T)

diffall<-aut-controls
library(ade4)


#let's take the 20 most variance:
sd(diffall[1,])
sfdiffall<-apply(diffall, MARGIN = 1, sd)
select<-diffall[sfdiffall> 0.02,]

#let's do a t-test
pval=c()
for (i in 1:dim(select)[1]){
  a<-t.test(select[i,])$p.value
  pval<-c(pval,a)
}
adj<-p.adjust(pval,method="fdr")
adj[adj<0.15]

#let's cluster
install.packages("pheatmap")
library(pheatmap)
#heatmap 
colfunc <- colorRampPalette(c("white", "black"))
#ok we're doing the clustering on the distance unifrac distance 
heatmap_test<-pheatmap(t(diffall), clustering_method ="ward.D2")
heatmap_test<-pheatmap(select, clustering_method ="ward.D")

diffall.clust <- cutree(heatmap_test, h=2)
diffall.clust <- cbind(t(diffall), 
                       cluster = cutree(heatmap_test$tree_row, 
                                        k = 2))
#so here I have the two clusters
diffall.clust[,36] #give the clusters

cluster_1_diff<-names(diffall.clust[diffall.clust[,36]==1,36])
cluster_2_diff<-names(diffall.clust[diffall.clust[,36]==2,36])


#ok let's do the matrix diff on these samples only
diffallcluster1<-diffall[,c(cluster_1_diff)]
diffallcluster2<-diffall[,c(cluster_2_diff)]
pval=c()
for (i in 1:dim(diffallcluster)[1]){
  a<-t.test(diffallcluster1[i,])$p.value
  pval<-c(pval,a)
}
adj<-p.adjust(pval,method="fdr")
adj[adj<0.15]
diffallcluster1[15,] #significant 0.004 after correction
boxplot(as.numeric(diffallcluster1[15,]))

pval=c()
for (i in 1:dim(diffallcluster2)[1]){
  a<-t.test(diffallcluster2[i,])$p.value
  pval<-c(pval,a)
}
adj<-p.adjust(pval,method="fdr")
adj[adj<0.15]
diffallcluster2[15,] #significant 0.004 after correction
boxplot(as.numeric(diffallcluster2[15,]))

#are they differences in the classifier
colnamcluster2<-gsub("X","",colnames(diffallcluster2))
colnamcluster1<-gsub("X","",colnames(diffallcluster1))
wilcox.test(ASD_otu_CSS_mapping[colnamcluster2,]$classifier,ASD_otu_CSS_mapping[colnamcluster1,]$classifier)
boxplot(ASD_otu_CSS_mapping[colnamcluster2,]$classifier,ASD_otu_CSS_mapping[colnamcluster1,]$classifier)
boxplot(ASD_otu_CSS_mapping[colnamcluster2,]$age_month_ok,ASD_otu_CSS_mapping[colnamcluster1,]$age_month_ok)
wilcox.test(ASD_otu_CSS_mapping[colnamcluster2,]$age_month_ok,ASD_otu_CSS_mapping[colnamcluster1,]$age_month_ok[-11])


(ASD_otu_CSS_mapping[colnamcluster2,]$Aut_status,ASD_otu_CSS_mapping[colnamcluster1,]$Aut_status)

library(pvclust)
result1000 <- pvclust(diffall, method.hclust="ward", nboot=1000)
plot(result1000, )


##Are the difference normals?
tdiffall<-t(diffall)
pval=c()
for (i in 1:35){
  a<-ks.test(tdiffall [,i],"pnorm",mean=mean(tdiffall [,i]),sd=sd(tdiffall [,i]))$p.value
  pval<-c(a,pval)
}

library(ggplot2)
require (reshape)
tthresh_10pourcent<-t(thresh_10pourcent)
long = melt(tdiffall)
long<- long[,-1]
colnames(long)<-c("variable","value")
ggplot(long, aes(x=value, color=variable)) + geom_density()

rownames(tdiffall)<-gsub("X","",rownames(tdiffall))
tdiffall <- tdiffall[ order(row.names(tdiffall)), ]
diffall.pca<-dudi.pca(diffall)
Cond<-factor(clean_mapping_file_91[rownames(clean_mapping_file_91) %in% rownames(tdiffall),]$Aut_status)

plot(diffall.pca$co, col = Cond)

#looks good let's try to do a interclass pca
diffall.pca.bet<-bca(diffall.pca,Cond, scan = FALSE, nf = 2)
colText=
  s.arrow(diffall.pca.bet$co)
Pol.pca.bet<-bca(Pol.pca,Cond)

pval=C()
for (i in 1:dim(ASD_otu_CSS_mapping)[2]){
  cor.test(ASD_otu_CSS_mapping$age_month_ok,ASD_otu_CSS[,i], method = "spearman")
  pval<-c(pval,a)}

pval[p.adjust(pval, method = "fdr") < 0.15]