source("~/Lab/microbiome2/scripts/getTables.R")
table <- getPcaFilteredOtuTable()
#table <- data.matrix(table[1:12,])
mapping <- getMapping()

#awkward third sibling

temp <- splitAutControl(table, mapping)
aut <- temp[[1]]
control <- temp[[2]]
autMap <- temp[[3]]
cMap <- temp[[4]]

differences <- control - aut


checkNormal <- function(dist){
  p <- shapiro.test(dist)$p.value
  return(p)
}

ps <- apply(table, 1, checkNormal)
mean(ps)
# mean of 2.415 * 10 ^ -9. Means very confidently NOT normal