################################################################################
##############################  scrLight.Rdat  #################################
################################################################################
##
## This object is a typical singlecellRNAseq object created at first with
## mandatory slots and allows to start running CONCLUS.
##
## This is the code to generate this object:

library(conclus)

## Loading of the count matrix
countMatrixPath <- file.path(system.file("extdata", package = "conclus"),
                                "countMatrix.tsv")
countMatrix <- loadDataOrMatrix(file=countMatrixPath, type="countMatrix",
                                ignoreCellNumber=TRUE)

## Creation of the singlecellRNAseq object 
scrLight <- singlecellRNAseq(experimentName = "LightExperience",
                  countMatrix     = countMatrix,
                  species         = "mouse",
                  outputDirectory = tempdir())

## Save the object
save(scrLight, file="inst/extdata/scrLight.Rdat")
