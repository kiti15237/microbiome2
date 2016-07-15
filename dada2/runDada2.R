library(optparse)
option_list <- list(
  make_option(c("-p", "--path"), type = "character", default=NULL, help = "path contained fastq files by sample"),
  make_option(c("-o", "--outdir"), type = "character", default=NULL, help = "out directory")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

path <- opt$path
#path <- "~/Lab/data/dada2/raw_fastqs/" # CHANGE ME to the directory containing the fastq files after unzipping.
fns <- list.files(path)
fns
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
library(ggplot2)
library(dada2)
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

dadaFs <- dada(derepFs, err=inflateErr(tperr1,3), selfConsist = TRUE)

dadaRs <- dada(derepRs, err=inflateErr(tperr1,3), selfConsist = TRUE)
#sapply(seq(i,length(dadaFs)), function(i) ggsave(paste('/scratch/users/ctataru6/microbiome/dada2_', names(derepFs)[[i], "_errors.png", sep=""), dada2::plotErrors(dadaFs[[i]], nominalQ=TRUE)))

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap = 25, maxMismatch = 0)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers[names(mergers) != "Mock"])
seqtab_clean <- seqtab[,-which(nchar(colnames(seqtab)) == 207)]
seqtab.nochim <- removeBimeraDenovo(seqtab_clean, verbose=TRUE)

taxa_gg <- assignTaxonomy(seqtab.nochim, "~/Lab/data/dada2/gg_13_8_train_set_97.fa.gz")

mapping<-read.table("C:/Users/ctataru/Documents/Lab/data/mapping_combined_ok_corrected3.txt", header = T, skip=F)
mapping$SampleID<-gsub(".0","",mapping$SampleID, fixed=T)
mapping$SampleID<-gsub(" ","",mapping$SampleID, fixed=T)
mapping <- mapping[order(mapping$SampleID),]
seqtab.nochim_ord <- seqtab.nochim[order(rownames(seqtab.nochim)), ]
mapping$SampleID <- rownames(seqtab.nochim_ord)
write.table(mapping, "~/Lab/data/dada2/mapping_dada2.txt")
samples.out <- rownames(seqtab.nochim)
otu_table <- phyloseq::otu_table(seqtab.nochim, taxa_are_rows=FALSE)

data <- sample_data(data.frame( pair=mapping$Pair, 
                               classifier=mapping$classifier,
                               age=mapping$age_month_ok, 
                               treatment=mapping$Treatment))
rownames(data) <- mapping$SampleID
ps <- phyloseq(otu_table(seqtab.nochim_ord, taxa_are_rows=FALSE), 
               data,
               tax_table(taxa_gg))


taxa_names <- apply(taxa_gg, 1, function(taxa) return(paste0(taxa, collapse=";")))
colnames(otu_table) <- taxa_names
write.table(t(otu_table), paste(opt$dir, "/otu_table_raw_0mismatch.txt", sep=""))

