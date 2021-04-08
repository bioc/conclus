################################################################################
################  extdata/expected_normalizedMatrix.Rdat  ######################
################################################################################
##
## This object is a SingleCellExperiment object created by CONCLUS after the
## normalisation.
## 
## This is the code to generate this object:

library(conclus)

## Loading of the count matrix
countMatrixPath <- file.path(system.file("extdata", package = "conclus"),
                                "countMatrix.tsv")
countMatrix <- loadDataOrMatrix(file=countMatrixPath, type="countMatrix",
                                ignoreCellNumber=TRUE)

## Loading of the count matrix
columnsMetaDataPath <- file.path(system.file("extdata", package = "conclus"),
                                "colData.tsv")
columnsMetaData <- loadDataOrMatrix(file=columnsMetaDataPath, type="coldata", 
                                    columnID="cell_ID")

## Creation of the singlecellRNAseq object 
scrLight <- singlecellRNAseq(experimentName = "LightExperience",
                  countMatrix     = countMatrix,
                  species         = "mouse",
                  outputDirectory = "YourOutputDirectory")

## Normalization 
scrNorm <- normaliseCountMatrix(scrLight, coldata=columnsMetaData, info=FALSE)

## Retrieved the SingleCellExperiment object
expectedNormalizedMatrix <- getSceNorm(scrNorm)

## Save the object
file <- "inst/extdata/expected_normalizedMatrix.Rdat"
save(expectedNormalizedMatrix, file=file)