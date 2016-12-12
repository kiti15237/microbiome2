library(phyloseq)

css_normalized_otu_path <- "~/Lab/16S/dada2/filtered_otu_table_normCSS.biom"
mapping_path <- "~/Lab/16S/mapping/mapping_file_MMD_11212016_FINAL.txt"

getMapping <- function(){
  mapping <- read.delim(mapping)
  colnames(mapping)[1] = "SampleID"
  mapping$SampleID <- as.character(mapping$SampleID)
  mapping <- mapping[mapping$SampleID != "168.1",] #sequence sucks
  mapping <- mapping[mapping$SampleID != "367",] # sequence sucks
  mapping <- mapping[mapping$SampleID != "173",] #in mapping but no data
  mapping <- mapping[mapping$SampleID != "384",] #in mapping but no data
  mapping <- mapping[mapping$SampleID != "211",] #seqeunce sucks
  mapping <- mapping[mapping$SampleID != "172.saliva",]
  mapping <- mapping[mapping$SampleID != "173.saliva",]
  mapping$classifier[is.na(mapping$classifier)] <-    4
  mapping <- mapping[order(mapping$SampleID),]
  
  return(mapping)
}

getOtuTable <- function(){
  biom <-phyloseq::import_biom(css_normalized_otu_path)
  table <- otu_table(biom, taxa_are_rows = T)
  table <- table[apply(table, 1, function(row) return(sum(row) > 0)), ]
  table <- formatTable(table)
  return(table)
}

getTableQiime <- function(){
  biom <- phyloseq::import_biom("~/Lab/16S/otuTables/open_oct/filtered_otu_table_normCSS.biom")
  table <- otu_table(biom, taxa_are_rows = T)
  table <- table[apply(table, 1, function(row) return(sum(row) > 0)), ]
  table <- formatTable(table)
  return(table)
}
getOtuTableL3 <- function(){
  biom <-import_biom("~/Lab/16S/otuTables/open_oct/phylogeny_levels/filtered_otu_table_normCSS_L3.biom")
  table <- otu_table(biom, taxa_are_rows = T)
  table <- formatTable(table)
  return(table)
}
getOtuTableL4 <- function(){
  biom <-import_biom("~/Lab/16S/otuTables/open_oct/phylogeny_levels/filtered_otu_table_normCSS_L4.biom")
  table <- otu_table(biom, taxa_are_rows = T)
  table <- formatTable(table)
  return(table)
}
getOtuTableL5 <- function(){
  biom <-import_biom("~/Lab/16S/otuTables/open_oct/phylogeny_levels/filtered_otu_table_normCSS_L5.biom")
  table <- otu_table(biom, taxa_are_rows = T)
  table <- formatTable(table)
  return(table)
}
getOtuTableL6 <- function(){
  biom <-import_biom("~/Lab/16S/otuTables/open_oct/phylogeny_levels/filtered_otu_table_normCSS_L6.biom")
  table <- otu_table(biom, taxa_are_rows = T)
  table <- formatTable(table)
  return(table)
}

getPcaFilteredOtuTable <- function(){
  table <- read.csv("~/Lab/16S/otuTables/open_oct/pca_filtered/otu_table_pcaFiltered_byDiffs.txt", row.names=1, sep="")
  table <- formatTable(table)
 # mapping <- mapping[!(mapping$SampleID %in% setdiff(mapping$SampleID, colnames(table))),]
  return(table)
}


formatTable <- function(table){
  colnames(table) <- gsub("X", "", colnames(table))
  table <- table[,order(colnames(table))]
  return(table)
}

getTree <- function(){
  tree = read_tree_greengenes("~/Lab/16S/otuTables/open_oct/rep_set.tre")
  return(tree)
}


splitAutControl <- function(table, mapping){
  table <- table[ , mapping$SampleID != "33"] #thired sibling of pair 2
  mapping <- mapping[ mapping$SampleID != "33", ]
  table <- table[ , mapping$SampleID != "232"] #third sibling of pair 55 w/ autism
  mapping <- mapping[ mapping$SampleID != "232", ]
  
  
  aut <- table[ , mapping$Treatment == "Aut"]
  control <- table[ , mapping$Treatment == "Control"]
  autMap <- mapping[mapping$Treatment == "Aut",]
  cMap <- mapping[mapping$Treatment == "Control",]
  
  aut <- aut[, order(autMap$Pair)]
  control <- control[,order(cMap$Pair)]
  autMap <- autMap[order(autMap$Pair),]
  cMap <- cMap[order(cMap$Pair),]
  
  ####pairing###
  aut <- aut[ ,autMap$Pair %in% cMap$Pair]
  autMap <- autMap[autMap$Pair %in% cMap$Pair,]
  control <- control[, cMap$Pair %in% autMap$Pair]
  cMap <- cMap[cMap$Pair %in% autMap$Pair,]
  
  return(list(aut, control, autMap,cMap))
}
