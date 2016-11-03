library(glmnet)
source("~/Lab/microbiome2/scripts/getTables.R")
source("~/Lab/16S/scripts/getTrainTestCohorts.R")

mapping <- getMapping()
table <- getPcaFilteredOtuTable()


letters = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h','i', 'j', 'k', 'l', 'm', 'n')
combos = combn(letters, 3)
labels <- unique(apply(combos, 2, paste, collapse=""))
rownames(table) <- labels[1:nrow(table)]

temp <- splitAutControl(table, mapping)
aut <- temp[[1]]
control <- temp[[2]]
autMap <- temp[[3]]
cMap <- temp[[4]]


aut = t(aut)
control = t(control)

iterations = 5
#Rather than doing true k fold validation, I'm picking a random subset of autism and control samples at every iteration.
#training and test sets are constructed to have even numbers of each condition 
for(iter in seq(1:iterations)){
  l = getTrainTestCohorts(aut, control)
  train <- l[[1]]
  train_y <- l[[2]]
  test <- l[[3]]
  test_y <- l[[4]]
  
  boost <- ada(as.matrix(train), train_y, as.matrix(test), test_y, type="discrete", iter = 40)
  summary(boost)
}