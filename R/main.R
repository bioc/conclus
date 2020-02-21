source("R/AllClasses.R")
source("R/AllGenerics.R")
source("R/methods-accessor-scRNAseq.R")
source("R/methods-replace-scRNAseq.R")
source("R/methods-class-scRNAseq.R")


## Data
outputDirectory <- "./YourOutputDirectory"
experimentName <- "Bergiers"

countMatrix <- read.delim(file.path(system.file("extdata", package="conclus"),		
                                    "Bergiers_counts_matrix_filtered.tsv"), 
                          stringsAsFactors = FALSE)

columnsMetaData <- read.delim(file.path(system.file("extdata",
                                                    package="conclus"), 
                                        "Bergiers_colData_filtered.tsv"))

## Construction
scr <- scRNAseq(experimentName  = experimentName, 
                countMatrix     = countMatrix, 
                colData         = columnsMetaData,
                species         = "mmu",
                outputDirectory = outputDirectory 
)

## Normalization with S4 method
scrS4 <- normaliseCountMatrix(src, colData = columnsMetaData) # return NCM
getNormalizedCountMatrix(scr)    # Slot normalizedCountMatrix empty
getNormalizedCountMatrix(scrS4)  # Slot normalizedCountMatrix full
