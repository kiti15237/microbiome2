library(optparse)
option_list <- list(
  make_option(c("-p", "--path"), type = "character", default=NULL, help = "path contained fastq files by sample")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

path <- opt$path
#path <- "~/Lab/data/dada2/raw_fastqs/" # CHANGE ME to the directory containing the fastq files after unzipping.
fns <- list.files(path)
fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
fnFs <- paste0(path, fnFs)
fnRs <- paste0(path, fnRs)
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)

for(i in seq_along(fnFs)) {
  file1 <- fnFs[i]
  file2 <- fnRs[i]
  f1 <- read.delim(file1, header=F)
  f2 <- read.delim(file2, header=F)
  h1 <- f1[seq(1, nrow(f1), by=4),]
  h2 <- f2[seq(1, nrow(f2), by=4),]
  h1 <- sub(" (.*)", "_", h1)
  h2 <- sub(" (.*)", "_", h2)
  
  if(length(setdiff(h1,h2)) != 0 || length(setdiff(h2, h1))!=0){
    seqs1 <- as.character(f1[seq(2, to=nrow(f1), by=4),])
    plus1 <- as.character(f1[seq(3, to=nrow(f1), by=4),])
    quals1 <- as.character(f1[seq(4, to=nrow(f1), by=4),])
    seqs2 <- as.character(f2[seq(2, to=nrow(f2), by=4),])
    plus2 <- as.character(f2[seq(3, to=nrow(f2), by=4),])
    quals2 <- as.character(f2[seq(4, to=nrow(f2), by=4),])
  
    cleanH1 <- h1
    cleanH2 <- h2
    
    if(length(setdiff(h1, h2)) != 0){
      cleanH1 <- h1[-which(h1==setdiff(h1,h2))]
      seqs1 <- seqs1[-which(h1==setdiff(h1,h2))]
      plus2 <- plus1[-which(h1==setdiff(h1,h2))]
      quals2 <- quals1[-which(h1==setdiff(h1,h2))]
      
    }
    if(length(setdiff(h2, h1)) != 0){
      cleanH2 <- h2[-which(h2==setdiff(h2,h1))]
      seqs2 <- seqs2[-which(h2==setdiff(h2,h1))]
      plus2 <- plus2[-which(h2==setdiff(h2,h1))]
      quals2 <- quals2[-which(h2==setdiff(h2,h1))]
    }
    cleanH1 <- sub("_", " 1:N:0:0", cleanH1)
    cleanH2 <- sub("_", " 2:N:0:0", cleanH2)
    
  
    
    vectorized1 <- c(rbind(cleanH1, seqs1, plus1, quals1))
    info1 <- sapply(seq(1,length(vectorized1), by=4), 
                   function(i) return(paste(vectorized1[c(i,i+1,i+2,i+3)], collapse="\n")))
    
    vectorized2 <- c(rbind(cleanH2, seqs2, plus2, quals2))
    info2 <- sapply(seq(1,length(vectorized2), by=4), 
                    function(i) return(paste(vectorized2[c(i,i+1,i+2,i+3)], collapse="\n")))
    
  
    write(info1, file=file1)
    write(info2, file=file2)
    
  }
}
