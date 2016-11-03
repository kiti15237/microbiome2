

getTrainTestCohorts <- function(aut, control){
  
    aut_order = sample(nrow(aut), replace = F)
    aut_shuffle <- aut[aut_order, ] #randomizes order
    autMap_shuffle <- autMap[aut_order, ] #keeps map in that same order
    
    control_order = sample(nrow(control), replace = F)
    control_shuffle <- control[control_order, ] #randomizes order
    cMap_shuffle <- cMap[control_order, ] # keeps map in that same order
    
    aut_train = aut_shuffle[1: (nrow(aut) * 2 / 3), ]
    aut_train_y = autMap_shuffle$Treatment[1: (nrow(aut) * 2 / 3)]
    
    control_train = control_shuffle[1: (nrow(control) * 2 / 3), ]
    control_train_y = cMap_shuffle$Treatment[1: (nrow(control) * 2 / 3)]
    
    aut_test = aut_shuffle[((nrow(aut) * 2/3) + 1): nrow(aut), ]
    aut_test_y = autMap_shuffle$Treatment[((nrow(aut) * 2/3) + 1): nrow(aut)]
    
    control_test = control_shuffle[((nrow(control) * 2/3) + 1): nrow(control), ]
    control_test_y = cMap_shuffle$Treatment[((nrow(control) * 2/3) + 1): nrow(control)]
    
    train <- as.data.frame(rbind(aut_train, control_train))
    train_y <- c(aut_train_y, control_train_y)
    test <- as.data.frame(rbind(aut_test, control_test))
    test_y <- c(aut_test_y, control_test_y)
    
    return(list(train, train_y, test, test_y))
}
