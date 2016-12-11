ps <- readRDS("~/Lab/16S/dada2/phyloseqObj.rds")

ord.nmds.bray <- ordinate(ps, method="RDA", distance="bray")
plot_ordination(ps, ord.nmds.bray, color="treatment", title="Bray NMDS")

ps_aut <- subset_samples(ps, treatment=="Aut")
ps_control <- subset_samples(ps, treatment=="Control")

top20 <- names(sort(taxa_sums(ps_aut), decreasing=TRUE))[1:10]
ps.top20 <- transform_sample_counts(ps_aut, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="SampleID", fill="Genus") 


p_inf_aut <- subset_taxa(ps_aut, Species=="s__parainfluenzae")
p_inf_control <- subset_taxa(ps_control, Species=="s__parainfluenzae")
