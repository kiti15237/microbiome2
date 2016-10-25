getTable.R is a fetching script. The path names need to be set appropriately for you local environment

Once path names are correct, any script that will use the otu table or mapping file can simply
1) `source("~/Lab/microbiome2/scripts/getTables.R")`

2) `mapping <- getMapping()`

3) `table <- getOtuTable()`


To get separate tables for aut and control samples, and separate mapping files for each, run: 

`temp <- splitAutControl(table, mapping)`
`aut <- temp[[1]]`
`control <- temp[[2]]`
`autMap <- temp[[3]]`
`cMap <- temp[[4]]`
