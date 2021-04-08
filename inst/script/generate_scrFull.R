################################################################################
##############################  scrFull.Rdat  ##################################
################################################################################
##
## This object is a typical singlecellRNAseq object that has followed all the 
## steps of CONCLUS and has all slots completed.
##
## This is the code to generate this object:

library(conclus)

## Loading of the count matrix
countMatrixPath <- file.path(system.file("extdata", package = "conclus"),
                                "countMatrix.tsv")
countMatrix <- loadDataOrMatrix(file=countMatrixPath, type="countMatrix",
                                ignoreCellNumber=TRUE)

## Running CONCLUS
scr <- runCONCLUS(outputDirectory=tempdir(),
                    experimentName="LightExperience", countMatrix=countMatrix,
                    species="mouse", tSNENb=1, clusterNumber=4, cores=1,
                    epsilon=c(380, 390, 400), minPoints=c(2,3), 
                    PCs =c(4,5,6,7,8,9,10), perplexities=c(2,3))

## Save the object
save(scr, file="inst/extdata/scrFull.Rdat")
