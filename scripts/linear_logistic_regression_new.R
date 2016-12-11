library(optparse)
library(glmnet)


source("~/Lab/microbiome2/scripts/getTables.R")
mapping <- getMapping()
table<- getOtuTable()


tTable <- t(table)

plotAccuracy <- function(df){
  plot(df$classifier, col="red", ylim=c(-9, 9))
  par(new=TRUE)
  plot(df$predicted, col="blue", add=T, ylim=c(-9, 9))
}

regression2 <- function(training, testing, glmnet_family, alpha){
  rownames(testing) <- gsub("X", "", rownames(testing))
  rownames(training) <- gsub("X", "", rownames(training))
  
  if(glmnet_family == "binomial"){
    categories_train <- ifelse(mapping$classifier[match(rownames(training), mapping$SampleID)] < 0, 1, 0)
    cvfit <- cv.glmnet(as.matrix(training), as.factor(categories_train), alpha = alpha, family = glmnet_family,nfold=4 )
    results <- ifelse(mapping$classifier[match(rownames(testing), mapping$SampleID)]< 0, 1, 0)
    predictions <- glmnet::predict.cv.glmnet(cvfit, as.matrix(testing), s=0.01, type="class")
    correct <- sum(as.numeric(predictions) != results) / length(results)
  }
  else if(glmnet_family == "poisson"){
    classifier_train <- mapping$classifier[match(rownames(training), mapping$SampleID)]
    classifier_train[classifier_train == 4] <- 0
    classifier_train <- abs(classifier_train)
    cvfit <- cv.glmnet(as.matrix(training), classifier_train, alpha = alpha, family = glmnet_family)
    results <- mapping$classifier[match(rownames(testing), mapping$SampleID)]
    results[results == 4] <- 0
    results <- abs(results)
  }
  else if (glmnet_family == "gaussian"){
    cvfit <- cv.glmnet(training, mapping$classifier[match(rownames(training), mapping$SampleID)], nfold=4, alpha = alpha, family = glmnet_family)
    results <- mapping$classifier[match(rownames(testing), mapping$SampleID)]
    predictions <- glmnet::predict.cv.glmnet(cvfit, as.matrix(testing), s=0.5)
    correct <- sum((predictions - results)^2) / length(results)
    df <- data.frame(classifier=results, predicted=as.numeric(predictions))
    #plotAccuracy(df)
  }
  return(correct)
}

iterations <- 20
accuracy <- 0
alpha <- 1
lambda <- 0.6
glmnet_family <- "binomial"
relevantVars <- c()

for(iter in seq(1:iterations)){
  tmp <- tTable[sample(nrow(tTable)),]
  results <- regression2(tmp[1:38, ], tmp[39:48,], glmnet_family=glmnet_family, alpha = alpha)
  print(results)
  accuracy <- accuracy + results
}
#output accuracy
print(table(relevantVars)[order(table(relevantVars))])
cat("Alpha: ", alpha, "\n")
cat("Method: ", glmnet_family, "\n")
cat("Total Accuracy for : ", accuracy / iterations)

