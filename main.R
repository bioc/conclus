if (!require("pacman")) install.packages("pacman")

pacman::p_load(rlist, foreach, ggplot2, pheatmap, zoo, dynamicTreeCut, factoextra, 
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
source("R/methods-scRNAseq-clustering.R")
source("R/methods-scRNAseq-tsne.R")
source("R/methods-scRNAseq-dbscan.R")
source("R/methods-Tsne-accessors.R")

## Data
outputDirectory <- "YourOutputDirectory"
experimentName <- "Bergiers"

countMatrix <- as.matrix(read.delim(file.path("inst/extdata/Bergiers_counts_matrix_filtered.tsv"), stringsAsFactors = FALSE))
columnsMetaData <- read.delim(file.path("inst/extdata/Bergiers_colData_filtered.tsv"))

## Construction
scrS4 <- scRNAseq(experimentName  = experimentName, 
                countMatrix     = countMatrix, 
                colData         = columnsMetaData,
                species         = "mmu",
                outputDirectory = outputDirectory)

## Normalization with S4 method
scrS4Norm <- normaliseCountMatrix(scrS4, colData = columnsMetaData) 

## Test Clustering 
# testClustering(scrS4Norm)

## TSNE
scrS4TSNE <- generateTSNECoordinates(scrS4Norm)

## DBSCAN
scrS4DBSCAN <- runDBSCAN(scrS4TSNE)
# scrS4DBSCAN@dbscanList[[1]]@clustering


