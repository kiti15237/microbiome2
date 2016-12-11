library(e1071)
library(crossval)
library(pROC)



getTrainTestSet <- function(aut, control, autMap, cMap){
  s <- seq(1, ncol(aut))
  train_seq <- sample(s, 40, replace=F)
  #train_seq <- c(47, 32, 48 ,42 ,19, 26, 44, 33,  5, 15, 41, 18, 45, 14, 43, 25,
  #               21 ,17, 39, 12,  3, 36, 35,  2,  9, 20, 28, 31, 51, 11, 24, 46, 23,  1, 40, 38, 34, 27, 10,  6)
  aut_train <- aut[ , train_seq]
  control_train <- control[ , train_seq]
  autMap_train <- autMap[train_seq, ]
  cMap_train <- cMap[train_seq,]
  
  test_seq <- s[!(s %in% train_seq)]
  #test_seq <- c(4 , 7,  8 ,13, 16, 22, 29, 30, 37, 49, 50)
  aut_test <- aut[ , test_seq]
  control_test <- control[ , test_seq]
  autMap_test <- autMap[test_seq, ]
  cMap_test <- cMap[test_seq, ]
  
  trainTable <- cbind(aut_train, control_train)
  trainMapping <- rbind(autMap_train, cMap_train)
  return(list(trainTable, trainMapping, aut_test, control_test))
}

predictPairs <- function(aut_test, control_test, m){
  predicted <- c()
  for (i in seq(1, ncol(aut_test))){
    p_aut <- predict(m, t(aut_test[,i]), type="raw")
    p_control <- predict(m, t(control_test[,i]), type="raw")
    #If confidence in predicting "Aut" is higher for sample 1 than for sample 2, predict 1 for sample 1, and 0 for 2
    s1_pred <- as.numeric(ifelse(p_aut[1] > p_control[1], 1, 0))
    s2_pred <- as.numeric(ifelse(s1_pred == 1, 0, 1))
    predicted <- c(predicted, c(s1_pred,s2_pred))
  }
  
  actual <- rep(c(1,0), length(predicted) / 2)
  return(list(predicted, actual))
}

predictSingles <- function(aut_test, control_test, m){
  predicted <- c()
  for (i in seq(1, ncol(aut_test))){
    p_aut <- predict(m, t(aut_test[,i]), type="raw")
    p_control <- predict(m, t(control_test[,i]), type="raw")
    predicted <- c(predicted, c(1*(p_aut[1] > p_aut[2]) , 1*(p_control[1] > p_control[2])))
  }
  actual <- rep(c(1,0), length(predicted) / 2)
  print(predicted)
  return(list(predicted, actual))
}



runBayes <- function(table, mapping, aut, control, autMap, cMap){
  roc_objs
  #plot(seq(0,1, by=0.2), seq(0,1, by=0.2), type="n")
  for(i in seq(1,10)){
    temp <- getTrainTestSet(aut, control, autMap, cMap)
    trainTable <- temp[[1]]
    trainMapping <- temp[[2]]
    aut_test <- temp[[3]]
    control_test <- temp[[4]]
    
    #imp_otus <- getOtusGreatestMI(trainTable, trainMapping, 50)
    imp_otus <- getOtusMICutoff(table, mapping, .012)
    table_filtered <- 1 * (trainTable[rownames(trainTable) %in% imp_otus, ]) # changes abundance counts to prescence/absence 
    table_filtered <- cbind(table_filtered, rep(5,nrow(table_filtered)), rep(5, nrow(table_filtered)))
    treatments <- c(as.character(trainMapping$Treatment), "Aut", "Control")
    m <- naiveBayes(t(table_filtered), as.factor(treatments), laplace = 1)
    #predict_training <- predict(m, t(table_filtered))
    #actual_training <- as.factor(treatments)
    
    #Predict on test set
    aut_test_filtered <- aut_test[rownames(aut_test) %in% imp_otus, ]
    control_test_filtered <- control_test[rownames(control_test) %in% imp_otus, ]
    temp <- predictSingles(aut_test_filtered, control_test_filtered, m)
    predicted <- temp[[1]]
    actual <- temp[[2]]
  
    print(crossval::confusionMatrix(as.factor(actual), as.factor(predicted), negative = as.factor(actual)[2]))
   # roc_obj <- roc(actual, predicted)
    #roc_objs <- list(roc_objs, roc_obj)
    #print(lines.roc(roc_obj))
    preds <- prediction(predicted, actual)
    perf <- performance(preds, "acc")
    plot(perf)
  }
  mean(roc_objs)
}


# TO get tax assignments, cross ref with otuTables/open_oct/rep_set_tax_assignments
writeImpOtus <- function(imp_otus){
  tax <- read.delim("~/Lab/16S/otuTables/open_oct/rep_set_tax_assignments.txt", header=FALSE)
  bac <- tax$V2[tax$V1 %in% imp_otus  ]
  l <- sapply(bac, function(string) return(strsplit(as.character(string), ";")))
  #l[[10]] = c(l[[10]], "g__", "s__")
  df <- data.frame(matrix(unlist(l), ncol=7, byrow=T))
  write.csv(df, "~/Lab/16S/figures/imp_otus.csv", row.names = F, quote = F)
  
}



getParameters <- function(table, mapping, aut, control, autMap, cMap){
  #what you want is training and testing error when you are picking parameters
  temp <- getTrainTestSet(aut, control, autMap, cMap)
  trainTable <- temp[[1]]
  trainMapping <- temp[[2]]
  aut_test <- temp[[3]]
  control_test <- temp[[4]]
  
  auc_values_training <- c()
  auc_values_testing <- c()
  laplaceConstants <- seq(1, 10)
  cutoffs <- seq(.001, .02, by= .001)
  combos <- expand.grid(laplaceConstants, cutoffs)
  combos <- paste(combos[,1], combos[,2])
  for(c in cutoffs){
    l = 5
    # = as.numeric(strsplit(combo, " ")[[1]][1])
    #c = as.numeric(strsplit(combo, " ")[[1]][2])
      print(c)
      imp_otus <- getOtusMICutoff(table, mapping, c)
      table_filtered <- 1 * (trainTable[rownames(trainTable) %in% imp_otus, ]) # changes abundance counts to prescence/absence 
      table_filtered <- cbind(table_filtered, rep(l,nrow(table_filtered)), rep(l, nrow(table_filtered)))
      treatments <- c(as.character(trainMapping$Treatment), "Aut", "Control")
      
      m <- naiveBayes(t(table_filtered), as.factor(treatments))
      predict_training <- predict(m, t(table_filtered))
      actual_training <- as.factor(treatments)
      
     
      roc_obj_training <- roc(1*(actual_training=="Aut"), 1*(predict_training=="Aut"))
      auc_values_training <- c(auc_values_training, roc_obj_training$auc)
      
      
      
      #Predict on test set
      aut_test_filtered <- aut_test[rownames(aut_test) %in% imp_otus, ]
      control_test_filtered <- control_test[rownames(control_test) %in% imp_otus, ]
      temp <- predictSingles(aut_test_filtered, control_test_filtered, m)
      predict_testing <- temp[[1]]
      actual_testing <- temp[[2]]
      
     
      roc_obj_testing <- roc(actual_testing, predict_testing)
      auc_values_testing <- c(auc_values_testing, roc_obj_testing$auc)
    
    
  }
  return(list(auc_values_training, auc_values_testing, cutoffs))
  

}

source("~/Lab/microbiome2/scripts/getTables.R")
source("~/Lab/16S/scripts/mutualInformation.R")

mapping <- getMapping()
#mapping <- sample_data(ps)
table <- getOtuTable()
#table <- otu_table(ps)

temp <- splitAutControl(table, mapping)
aut <- temp[[1]]
control <- temp[[2]]
autMap <- temp[[3]]
cMap <- temp[[4]]

runBayes(table, mapping, aut, control, autMap, cMap)
#temp <- getParameters(table, mapping, aut, control, autMap, cMap)
#auc_values_training <- temp[[1]]
#auc_values_testing <- temp[[2]]
#cutoffs <- temp[[3]]
#plot(seq(1:length(combos)) , auc_values, main="Values of laplace constant", ylab="AUC values", xlab="laplace constants")
