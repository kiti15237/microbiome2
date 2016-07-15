library(optparse)
option_list <- list(
  make_option(c("-p", "--path"), type = "character", default=NULL, help = "path contained fastq files by sample")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

path <- opt$path
#path <- "/scratch/users/ctataru5/microbiome/dada2/MiSeq_SOP/test/" # CHANGE ME to the directory containing the fastq files after unzipping.
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
                    trimLeft=c(10, 10), truncLen=c(240,160), 
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

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#seqtab <- makeSequenceTable(mergers[names(mergers) != "Mock"])

#dim(seqtab)

# Inspect distribution of sequence lengths
#table(nchar(colnames(seqtab)))

#seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
#dim(seqtab.nochim)

#sum(seqtab.nochim)/sum(seqtab)
