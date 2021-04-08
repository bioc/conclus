################################################################################
###########################   extdata/colData.tsv   ############################
################################################################################
##
## Artificial generated col data used for examples and tests.
## Size : 20 * 3.
##
## For each cell, there are its name, its state (A, B, C or D) and its unique
## cell barcode.
## 
## Cells sharing a similar expression signature received the same state.
## 
## ## This is the code to generate the extdata/colData.tsv :

cell_ID <- paste0("c", seq(20))
state <- rep(c("A", "B", "C", "D"), 5)
cellBarcode <- c("AAAAAAAAAAT", "AAAAAAAAATT", "AAAAAAAATTT", "AAAAAAATTTT",
                 "AAAAAATTTTT", "AAAAATTTTTT", "AAAATTTTTTT", "AAATTTTTTTT",
                 "AATTTTTTTTT", "ATTTTTTTTTT", "GGGGGGGGGGA", "GGGGGGGGGAA",
                 "GGGGGGGGAAA", "GGGGGGGAAAA", "GGGGGGAAAAA", "GGGGGAAAAAA",
                 "GGGGAAAAAAA", "GGGAAAAAAAA", "GGAAAAAAAAA", "GAAAAAAAAAA")

colData <- cbind(cell_ID, state, cellBarcode)
rownames(colData) <- cell_ID

write.table(colData, file="inst/extdata/colData.tsv", quote=FALSE, sep="\t")
