#library(optparse)
library(glmnet)
source("~/Lab/microbiome2/scripts/getTables.R")


mapping <- getMapping()
table <- getOtuTable()

#potentially filter out otu's based on low standard dev. across samples.
filtering <- function(table, sampleMin, SDMultiplier, sampleMax){
  #Filter out those KOs that appear in less than 6 samples
  threshold <- sampleMin
  toKeep <- apply(table, 1, function(r) return(sum(r>0) >= threshold))
  table <- table[toKeep, ]
  
  #Filter out those KOS that appear in more than sampleMax samples
  toKeep <- apply(table, 1, function(r) return(sum(r > 0) < sampleMax))
  table <- table[toKeep,]
  
  #Filter out those KOs that have SD less than the mean (or 3 times the mean)
  sd_threshold <- SDMultiplier
  sds <- apply(table, 1, sd)
  table <- table[sds > (SDMultiplier * mean(sds)) , ]
  return(table)
}



plotAccuracy <- function(df){
  plot(df$classifier, col="red", ylim=c(-9, 9))
  par(new=TRUE)
  plot(df$predicted, col="blue", add=T, ylim=c(-9, 9))
}




regression2 <- function(training, testing, glmnet_family, alpha, lambda){
  rownames(testing) <- gsub("X", "", rownames(testing))
  rownames(training) <- gsub("X", "", rownames(training))
  
  if(glmnet_family == "binomial"){
    #1 means affected with autism
    categories_train <- ifelse(mapping$Treatment[match(rownames(training), mapping$SampleID)] =="Aut", 1, 0)
    cvfit <- cv.glmnet(as.matrix(training), as.factor(categories_train), alpha = alpha, family = glmnet_family,nfold=5)
    predictions <- glmnet::predict.cv.glmnet(cvfit, as.matrix(testing), s=lambda, type="class")
    print(predictions)
    
    results <- ifelse(mapping$Treatment[match(rownames(testing), mapping$SampleID)] == "Aut", 1, 0)
    correct <- sum(as.numeric(predictions) == results) / length(results)
  }
  else if(glmnet_family == "poisson"){
    classifier_train <- mapping$classifier[match(rownames(training), mapping$SampleID)]
    classifier_train[classifier_train == 4] <- 0
    classifier_train <- abs(classifier_train)
    cvfit <- cv.glmnet(as.matrix(training), classifier_train, alpha = alpha, family = glmnet_family)
    results <- mapping$classifier[match(rownames(testing), mapping$SampleID)]
    results[results >= 0] <- 0
    results <- abs(results)
    predictions <- glmnet::predict.cv.glmnet(cvfit, as.matrix(testing), s=0.4)
    correct <- sum((predictions - results)^2) / length(results)
  }
  else if (glmnet_family == "gaussian"){
    cvfit <- cv.glmnet(as.matrix(training), mapping$classifier[match(rownames(training), mapping$SampleID)], nfold=4, alpha = alpha, family = glmnet_family)
    results <- mapping$classifier[match(rownames(testing), mapping$SampleID)]
    predictions <- glmnet::predict.cv.glmnet(cvfit, as.matrix(testing), s=lambda)
    print(predictions)
    correct <- sum((predictions - results)^2) / length(results)
  }

  return(correct)
}

iterations <- 100
accuracy <- 0
alpha <- 1
lambda <- 0.1
glmnet_family <- "binomial"
relevantVars <- c()

temp <- splitAutControl(table, mapping)
aut <- temp[[1]]
control <- temp[[2]]
autMap <- temp[[3]]
cMap <- temp[[4]]


aut = t(aut)
control = t(control)

#Rather than doing true k fold validation, I'm picking a random subset of autism and control samples at every iteration.
#training and test sets are constructed to have even numbers of each condition 
for(iter in seq(1:iterations)){
  aut_order = sample(nrow(aut), replace = F)
  aut_shuffle <- aut[aut_order, ] #randomizes order
  autMap_shuffle <- autMap[aut_order, ] #keeps map in that same order
  
  control_order = sample(nrow(control), replace = F)
  control_shuffle <- control[control_order, ] #randomizes order
  cMap_shuffle <- cMap[control_order, ] # keeps map in that same order
  
  aut_train = aut_shuffle[1: (nrow(aut) * 2 / 3), ]
  control_train = control_shuffle[1: (nrow(control) * 2 / 3), ]
  aut_test = aut_shuffle[((nrow(aut) * 2/3) + 1): nrow(aut), ]
  control_test = control_shuffle[((nrow(control) * 2/3) + 1): nrow(control), ]
  
  
  train <- as.data.frame(rbind(aut_train, control_train))
  test <- as.data.frame(rbind(aut_test, control_test))
  results <- regression2(train, test, glmnet_family=glmnet_family, alpha = alpha, lambda = lambda)
  print(results)
  accuracy <- accuracy + results
}

cat("Alpha: ", alpha, "\n")
cat("Method: ", glmnet_family, "\n")
cat("Total Accuracy for : ", accuracy / iterations)

