library(e1071)
library(crossval)
library(pROC)




predictSingleSample <- function(table, mapping){
  actual <- c()
  predicted <- c()
  ####leave one out validation
  for (pairId in unique(sort(mapping$Pair))){
    testVectors <- table[, mapping$Pair == pairId ]
    trainTable <- table[,-which(mapping$Pair == pairId)]
    trainMapping <- mapping[-which(mapping$Pair == pairId), ]
    
    otus <- getOtusGreatestMI(trainTable, trainMapping, 30)
    table_filtered <- 1 * (trainTable[rownames(trainTable) %in% otus, ] > 0)
    m <- naiveBayes(t(table_filtered), as.factor(trainMapping$Treatment), laplace = 1)
    
    #Predict cold which class each sample is from
    p <- predict(m, t(testVectors), type="class")
    actual <- c(actual, mapping$Treatment[mapping$Pair == pairId])
    predicted <- c(predicted, p)
    print(pairId)
  }
  
  print(crossval::confusionMatrix(as.factor(actual), as.factor(predicted), negative = as.factor(actual)[2]))
  roc_obj <- roc(actual, predicted)
  
  print(plot.roc(roc_obj))
}


predictPair <- function(aut, control, autMap, cMap){
  actual <- c()
  predicted <- c()
  ####leave one out validation
  for (i in 1:ncol(aut)){
    testVectors <- cbind(aut[,i], control[,i])
    trainTable <- cbind(aut[,-i], control[,-i])
    trainMapping <- rbind(autMap[-i,], cMap[-i,])
    
    #otus <- (rownames(trainTable))
    otus <- getOtusGreatestMI(trainTable, trainMapping, 30)
    table_filtered <- 1 * (trainTable[rownames(trainTable) %in% otus, ] > 0) # changes abundance counts to prescence/absence 
    m <- naiveBayes(t(table_filtered), as.factor(trainMapping$Treatment), laplace = 1)
    p <- predict(m, t(testVectors), type="raw")
    
    #Predict s1 to be autistic or not based on which sample look more autistic, found in column 1
    #Let prediction for s2 follow as the opposite category from s1
    s1_pred <- as.numeric(ifelse(p[1,1] >= p[2,1], 1, 0))
    s2_pred <- as.numeric(ifelse(s1_pred == 1, 0, 1))
    predicted <- c(predicted, c(s1_pred,s2_pred))
    print(i)
  }
  
  actual <- rep(c(1,0), length(predicted) / 2)
  print(crossval::confusionMatrix(as.factor(actual), as.factor(predicted), negative = as.factor(actual)[2]))
  roc_obj <- roc(actual, predicted)
  
  print(plot.roc(roc_obj))
}

source("~/Lab/microbiome2/scripts/getTables.R")
source("~/Lab/16S/scripts/mutualInformation.R")
mapping <- getMapping()
table <- getOtuTable()



temp <- splitAutControl(table, mapping)
aut <- temp[[1]]
control <- temp[[2]]
autMap <- temp[[3]]
cMap <- temp[[4]]


predictSingleSample(table, mapping)
predictPair(aut,control, autMap, cMap)

#We can pretty accurately distinguish between a pair, which one is more autistic. Can't do it cold
#WITH NORMALIZATION THE NUMBERS SHOOT UP
