library(phyloseq)

getMapping <- function(){
  mapping <- read.delim("~/Lab/16S/mapping/mapping_merged_oct.txt")
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
  biom <-import_biom("~/Lab/16S/otuTables/open_oct/otu_table_normCSS.biom")
  table <- otu_table(biom, taxa_are_rows = T)
  table <- formatTable(table)
  return(table)
}

getPcaFilteredOtuTable <- function(){
  biom <-import_biom("~/Lab/16S/otuTables/open_oct/pcaFiltered/otu_table_pcaFiltered_byDiffs.biom")
  table <- otu_table(biom, taxa_are_rows = T)
  table <- formatTable(table)
  return(table)
}


formatTable <- function(table){
  colnames(table) <- gsub("X", "", colnames(table))
  table <- table[, colnames(table)!="168.1"]
  table <- table[, colnames(table)!="367"]
  mapping <- mapping[mapping$SampleID != "211",] #seqeunce sucks
  table <- table[,order(colnames(table))]
  return(table)
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
