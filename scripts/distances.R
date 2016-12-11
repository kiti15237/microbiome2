ps <- readRDS("~/Lab/16S/dada2/phyloseqObj_midRootTree.rds")
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
     xlim=c(0, 0.5))
lines(density(as.vector(dists_aut)), col = "red")
lines(density(dists_aut_control), col = "purple")

hist(dists_control, col = rgb(0,0,1,0.5), breaks = 100, 
     main = "Weighted Unifrac Distances Between Control Samples and Autism Samples",
     xlab = "Weighted Unifrac Distances", xlim = c(0, 0.6))
hist(dists_aut, col=rgb(1,0,0,0.5), add=T, breaks=100)

boxplot(dists_aut, dists_control, dists_aut_control, 
        names=c("Aut-Aut distances", "Control-Control distances", "Aut-Control distances"),
        main = "Unifrac Distances Between Samples", 
        sub=paste(ks.test(dists_aut, dists_control)$p.value, ks.test(dists_control, dists_aut_control)$p.value)
        )
###Todo, take out sibling pairs in aut_control


plot(density(as.vector(dists_control)), col = "blue", 
     main= "Wunifrac Distances within same phenotype groups \n PCA filtered data", xlab = "Weighted Unifrac Distance")
lines(density(as.vector(dists_aut)), col = "red")
lines(density(dists_aut_control), col = "purple")
