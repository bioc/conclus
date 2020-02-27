if (!require("pacman")) install.packages("pacman")

pacman::p_load(ggplot2, pheatmap, zoo, dynamicTreeCut, factoextra, 
               digest, RColorBrewer, doParallel, BiocParallel, scran, 
               monocle, KEGGREST, AnnotationDbi, SingleCellExperiment,
               Cairo, rvest, curl, scater, Matrix, dbscan, fpc, matrixStats,
               dplyr, biomaRt, org.Mm.eg.db, grDevices, S4Vectors, Biobase,
               DataCombine, zoo, rvest, devtools, BiocManager)


source("R/AllClasses.R")
source("R/AllGenerics.R")
source("R/methods-scRNAseq-getters.R")
source("R/methods-scRNAseq-setters.R")
source("R/methods-scRNAseq-normalization.R")
source("R/methods-scRNAseq-tsne.R")


## Data
outputDirectory <- "YourOutputDirectory"
experimentName <- "Bergiers"

countMatrix <- as.matrix(read.delim(file.path("inst/extdata/Bergiers_counts_matrix_filtered.tsv"), 
                          stringsAsFactors = FALSE))

columnsMetaData <- read.delim(
    file.path("inst/extdata/Bergiers_colData_filtered.tsv"))

## Construction
scrS4 <- scRNAseq(experimentName  = experimentName, 
                countMatrix     = countMatrix, 
                colData         = columnsMetaData,
                species         = "mmu",
                outputDirectory = outputDirectory)


## Normalization with S4 method
scrS4Norm <- normaliseCountMatrix(scrS4, colData = columnsMetaData) # return NCM
getNormalizedCountMatrix(scrS4Norm)  # Slot normalizedCountMatrix full
getNormalizedCountMatrix(scrS4)      # Slot normalizedCountMatrix empty

scrS4TSNE <- generateTSNECoordinates(scrS4Norm)

class(getNormalizedCountMatrix(scrS4Norm))
getTSNEResults(scrS4Norm)
getTSNEResults(scrS4TSNE)


