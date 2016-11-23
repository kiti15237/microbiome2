



getMI <- function(otuVector, mapping){
  y <- matrix(c(0,0,0,0), nrow=2)
  #Create confustion matrix
  #(1,1) = no Aut, no otu =
  y[1,1] <- sum(mapping$Treatment == "Control" & otuVector == 0)
  #(1,2) = no Aut, otu
  y[1,2] <- sum(mapping$Treatment == "Control" & otuVector > 0)
  #(2,1) = Aut, not otu
  y[2,1] <- sum(mapping$Treatment == "Aut" & otuVector == 0)
  #(2,2) = Aut, otu
  y[2,2] <- sum(mapping$Treatment == "Aut" & otuVector > 0)
  
  return(mi.empirical(y))
}

getOtusGreatestMI <- function(table, mapping, n){
  otu_MI <- apply(table, 1, getMI,mapping)
  otu_MI_sorted <- sort(otu_MI, decreasing = T)
  return(names(otu_MI_sorted[1:n]))
}



source("~/Lab/microbiome2/scripts/getTables.R")
mapping <- getMapping()
table <- getOtuTable()
otus <- getOtusGreatestMI(table, mapping, 30)