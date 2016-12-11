library(dplyr)

library(gbm)

#otu_table <- read.table('otu_table_normCSS.txt')

mapping_merged <- read.csv('~/Lab/16s/mapping/mapping_merged_oct.csv')

otu_table_raw <- read.table('~/Lab/16S/otuTables/open_oct/filtered_otu_table.txt')
otu_table_pca <- read.table('~/Lab/16S/otuTables/open_oct/otu_table_pcaFiltered_byDiffs2.txt')
#colnames(otu_table_raw) <- gsub("X", "", colnames(otu_table_raw))
mapping_merged <- mutate(mapping_merged, X.SampleID = as.character(X.SampleID))

otu_raw <- data.frame(t(otu_table_raw))
otu_pca <- data.frame(t(otu_table_pca))

otu_raw <- otu_raw %>%
  
  mutate(X.SampleID = substring(row.names(otu_raw),2)) %>%
  
  inner_join(select(mapping_merged, X.SampleID, class = Treatment), by = 'X.SampleID') %>%
  
  mutate(class = class == 'Aut')


otu_pca <- otu_pca %>%
  
  mutate(X.SampleID = substring(row.names(otu_pca),2)) %>%
  
  inner_join(select(mapping_merged, X.SampleID, class = Treatment), by = 'X.SampleID') %>%
  
  mutate(class = class == 'Aut')


otu_raw <- select(otu_raw, -which(sapply(otu_raw, function(x) all(x == 0))))
#otu_pca <- select(otu_pca, -which(sapply(otu_pca, function(x) all(x == 0))))

set.seed(10110)

otu_raw <- otu_raw[sample(nrow(otu_raw)), ]
otu_pca <- otu_pca[sample(nrow(otu_pca)), ]

gbm0 <- gbm(class ~ ., data = select(otu_pca, -X.SampleID), train.fraction = 1.0, interaction.depth = 5,
            
            shrinkage = 0.001, n.trees = 2500, bag.fraction = 0.5, cv.folds = 10,
            
            distribution = "bernoulli", verbose = F)

gbm.perf(gbm0, method = 'cv')
