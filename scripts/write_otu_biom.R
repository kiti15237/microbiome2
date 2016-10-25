library(optparse)
library(tools)
library(RJSONIO)
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default ="~/Lab/data/dada2/output/0mismatch_25overlap/seqtab_otuId.txt" )
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
####Take in a dada2 exported otu table. Should have sequences as row names and samples as column names
path <- dirname(opt$input)
params <- file_path_sans_ext(basename(opt$input))

path = "~/Lab/16S/otuTables/open_oct"
params = "otu_table_pcaFiltered_byDiffs"


table <-  read.csv(paste(path,"/", params,".txt", sep=""), row.names=1, sep="")
#table <- read.csv("~/Lab/data/dada2/differences_noFilter.txt", sep="", row.names=1)
#table <- ddply(table, "OTUID", numcolwise(sum))
#rownames(table) <- table[,1]
#table <- table[,-1]

table <- table[, apply(table, 2, function(col) return(sum(col) > 0))]
colnames(table)<-gsub("X","",colnames(table))
table<- table[ , order(colnames(table))]


#taxa <- read.csv(paste(path, "taxa_", params, ".txt", sep=""), row.names=1, sep="")
#colnames(taxa) <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species' )
#taxa_ord <- taxa[match(rownames(table), rownames(taxa)),]

obj = list()
obj$id = "None"
obj$format = "Biological Observation Matrix 1.0.0"
obj$format_url = "http://biom-format.org"
obj$type = "OTU table"
obj$generated_by = "None"
obj$date = "None"
obj$matrix_type = "dense"
obj$matrix_element_type = "int"
nrows <- nrow(table)
ncols <- length(colnames(table))
obj$shape = list(nrows, ncols)
obj$rows = lapply(rownames(table), function(name) return(list(id=name, metadata=NULL)))
obj$columns = lapply(colnames(table), function(name) return(list(id=name, metadata=NULL)))
obj$data = matrix(as.vector(unlist(table)), nrow=nrows, ncol=ncols)

library(RJSONIO)
library(biom)
write(RJSONIO::toJSON(obj), paste(path,"/", params,  ".json", sep=""))
b <- read_biom(paste(path, "/", params,  ".json", sep=""))
write_biom(biom(b), biom_file = paste(path, "/", params,  ".biom", sep=""))
