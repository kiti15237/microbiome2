library(ggplot2)
library(dada2)
library(optparse)
library(phyloseq)


option_list <- list(
  make_option(c("-p", "--path"), type = "character", default= "~/Lab/data/dada2/raw_fastqs/", help = "path contained fastq files by sample"),
  make_option(c("-m", "--map"), type = "character", default = "~/Lab/data/mapping_combined_ok_corrected3.txt"),
  make_option(c("-o", "--outdir"), type = "character", default="~/Lab/16S/dada2/", help = "out directory"), 
  make_option(c( "--database"), type = "character", default="~/Lab/16S/gg_13_8_train_set_97.fa.gz", help = "out directory")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

path <- opt$path
#path <-"~/Lab/16S/split_allFastq_bySample/"

fns <- list.files(path)
fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
# Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Fully specify the path for the fnFs and fnRs
fnFs <- paste0(path, fnFs)
fnRs <- paste0(path, fnRs)
sample.names


#Plot quality of sequences if you desire. Commented out for speed
#sapply(fnFs, function(file) ggsave(paste(file, "_quality.png", sep=""), dada2::plotQualityProfile(file)))
#sapply(fnRs, function(file) ggsave(paste(file, "_quality.png", sep=""), dada2::plotQualityProfile(file)))


# Make filenames for the filtered fastq files
filtFs <- paste0(path, sample.names, "_F_filt.fastq.gz")
filtRs <- paste0(path, sample.names, "_R_filt.fastq.gz")
# Filter
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    trimLeft=c(10, 10), truncLen=c(140,140),
                    maxN=0, maxEE=2, truncQ=2,
                    compress=TRUE, verbose=TRUE)
}

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

predict <- stats::predict
dadaFs <- dada(derepFs, err=inflateErr(tperr1,3), selfConsist = TRUE)
dadaRs <- dada(derepRs, err=inflateErr(tperr1,3), selfConsist = TRUE)
#sapply(seq(i,length(dadaFs)), function(i) ggsave(paste('/scratch/users/ctataru6/microbiome/dada2_', names(derepFs)[[i], "_errors.png", sep=""), dada2::plotErrors(dadaFs[[i]], nominalQ=TRUE)))



makeFasta <- function(sample){
  counts <- seqtab.nochim[sample,]
  labels <- paste(">", sample, "_", seq(1:sum(counts)), "\n", sep = "")
  counts[is.na(counts)] <- 0
  seqs <- rep(colnames(seqtab.nochim), counts)
  return(paste(labels, seqs, collapse="\n"))
}



minOverlap = 25
maxMismatch = 0
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap = minOverlap, maxMismatch = maxMismatch)
# Inspect the merger data.frame from the first sample
seqtab <- makeSequenceTable(mergers[names(mergers) != "Mock"])
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)

sameLengthSeqs <- seqtab.nochim[nchar(rownames(seqtab.nochim))==233, ]
fasta <- sapply(rownames(sameLengthSeqs), makeFasta)
write(fasta, paste(opt$outdir, "seqs.txt", sep=""))

taxa <- assignTaxonomy(seqtab.nochim, opt$database)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax <- phyloseq::tax_table(taxa)
write.table(tax, paste(opt$outdir, "/taxa_table.txt", sep=""))

otu_table <- phyloseq::otu_table(seqtab.nochim, taxa_are_rows=FALSE)
write.table(t(otu_table), paste(opt$outdir, "/otu_table_raw_", maxMismatch,"mismatch_", minOverlap, "minOverlap.txt", sep=""))

tree <- constructPhyTree(rownames(otu_table))


###################################################
####phyloseq######################################
source("~/Lab/microbiome2/scripts/getTables.R")
mapping <- getMapping() #function from getTables. Some manual path setting

metadata <- sample_data(data.frame( SampleID = mapping$SampleID,
                                    pair=mapping$Pair, 
                                    classifier=mapping$classifier,
                                    age=mapping$age_month_ok, 
                                    batch=factor(mapping$batch),
                                    treatment= factor(mapping$Treatment, levels = c("Aut", "Control"))))
rownames(metadata) <- mapping$SampleID


ps <- phyloseq(otu_table(table, taxa_are_rows=T), 
               sample_data(metadata), 
               tax_table(as.matrix(tax)), 
               phy_tree(tree$tree))
ps
saveRDS(ps, psate(opt$outdir, "/phyloseqObj.rds"))


constructPhyTree <- function(seqs){
  ##################################################################################
  ####Construct Phylogenetic Tree######
  biocLite("msa")
  library(msa)
  mult <- msa(seqs, method="ClustalW", type="dna", order="input")
  
  #The phangorn package is then used to construct a phylogenetic tree. Here we first construct a neighbor-joining tree, and then fit a GTR+G+I maximum likelihood tree using the neighbor-joining tree as a starting point.
  library("phangorn")
  
  phang.align <- msa::msaConvert(mult, type="phangorn::phyDat")
  saveRDS(phang.align, paste(opt$outdic, "/phang_align.rds"))
  phang.align <- as.phyDat(mult, type="DNA")
  attr(phang.align, 'names') <- seqs
  dm <- dist.ml(phang.align)
  
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  treeNJ$tip.label <- seqs
  fit = pml(treeNJ, data=phang.align)
  
  ## negative edges length changed to 0!
  
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
  saveRDS(fitGTR, paste(opt$outdir, "/tree.rds"))
  detach("package:phangorn", unload=TRUE)
  return(fitGTR)
  
}




