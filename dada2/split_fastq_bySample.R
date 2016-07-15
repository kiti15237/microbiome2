source("http://bioconductor.org/biocLite.R")
biocLite('dada2')
library(optparse)
option_list <- list(
  make_option(c("-i", "--infile"), type = "character", default=NULL, help = "otu table file name")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt$infile)
fastq <- read.delim(opt$infile, header=F)
headers <-fastq[seq(1, to=nrow(fastq), by=4),]
seqs <- as.character(fastq[seq(2, to=nrow(fastq), by=4),])
plus <- as.character(fastq[seq(3, to=nrow(fastq), by=4),])
quals <- as.character(fastq[seq(4, to=nrow(fastq), by=4),])

getSampleName <- function(header){
  ind_dash <- regexpr("_", header)
  sampleName <- substr(header, 2, ind_dash - 1)
  if(sampleName == "180"){
	sampleName <- "180.1"
	print("Sample 180 changed to 181.1")
  }
  if(sampleName == "81"){
	sampleName <- "181.1"
  }
  return(sampleName)
}
writeOut <- function(text, sname){
  if(grepl("R1", opt$infile)){
    write(text, file=paste(sname, "_R1.fastq", sep="", collapse=""))
  }else if(grepl("R2", opt$infile)){
    write(text, file=paste(sname, "_R2.fastq", sep="", collapse=""))
  }
}
sampleNames <- sapply(headers, getSampleName)
print(unique(sampleNames))
newHeaders1 <- sapply(headers, function(header) return(paste("@", sub(".*? ", "", header), sep="")))
newHeaders <- sapply(newHeaders1, function(header) return(sub(" orig(.*)", "", header)))
vectorized <- c(rbind(newHeaders, seqs, plus, quals))
info <- sapply(seq(1,length(vectorized), by=4), function(i) return(paste(vectorized[c(i,i+1,i+2,i+3)], collapse="\n")))
bySample <- split(info, sampleNames)
mapply(writeOut, bySample, names(bySample))

