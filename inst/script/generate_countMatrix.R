################################################################################
###########################  extdata/countMatrix.tsv  ##########################
################################################################################
##
## Artificially generated count matrix used for examples and tests.
## Size : 20 * 20.
## The number of counts has been chosen so that we have 4 clusters 
## with 5 specific markers each. These numbers have been chosen arbitrarily to 
## have a fast running time for the examples and the tests.
##
## This is the code to generate the extdata/countMatrix.tsv :

## Creation of the matrix : each line in the for loop corresponds to an
## artificial category of cells.
countMatrix <- c()
for (i in seq_len(5)){
    countMatrix <- rbind(countMatrix, c(rep(i,5), rep(5+i,5), rep(10+i,5), 
                                        rep(20+i,5)))
    countMatrix <- rbind(countMatrix, c(rep(20+i,5), rep(i,5), rep(5+i,5), 
                                        rep(10+i,5)))
    countMatrix <- rbind(countMatrix, c(rep(10+i,5), rep(20+i,5), rep(i,5), 
                                        rep(5+i,5)))
    countMatrix <- rbind(countMatrix, c(rep(5+i,5), rep(10+i,5), rep(20+i,5), 
                                        rep(i,5)))
}

## Theses ensembl id have been chosen arbitrarily and regardless of cell or 
## number of counts.
rownames(countMatrix) <- c("Sipa1l2", "Mcts1", "Mov10", "Rtca", "Fmr1", "Sp1", 
                           "Uck1", "Gm2a", "Napsa", "Tsc2", "Rab5b", "Efnb2", 
                           "Apoe", "S100a4", "Eef1e1", "Braf", "Usp32",  
                           "Rapsn", "Ndufa9", "Ap4e1")

## The name of the cells have been chosen arbitrarily.
colnames(countMatrix) <- paste0("c", seq(20))

write.table(countMatrix, file="inst/extdata/countMatrix.tsv",
            sep="\t", quote=FALSE)

