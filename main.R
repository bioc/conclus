if (!require("pacman")) install.packages("pacman")

pacman::p_load(devtools,BiocManager, ggplot2, pheatmap, zoo, dynamicTreeCut, factoextra, 
               digest, RColorBrewer, doParallel, BiocParallel, scran, 
               monocle, KEGGREST, AnnotationDbi, SingleCellExperiment, Cairo, rvest, curl, 
               scater, Matrix, dbscan, fpc, matrixStats, dplyr, biomaRt, org.Mm.eg.db,
               grDevices, S4Vectors, Biobase, DataCombine, zoo, rvest)

library(devtools)
library(BiocManager)

source("R/AllClasses.R")
source("R/AllGenerics.R")
source("R/methods-accessor-scRNAseq.R")
source("R/methods-replace-scRNAseq.R")
source("R/methods-class-scRNAseq.R")



## Data
outputDirectory <- "./YourOutputDirectory"
experimentName <- "Bergiers"

countMatrix <- read.delim(file.path("inst/extdata/Bergiers_counts_matrix_filtered.tsv"), 
                          stringsAsFactors = FALSE)

columnsMetaData <- read.delim(file.path("inst/extdata/Bergiers_colData_filtered.tsv"))

## Construction
scr <- scRNAseq(experimentName  = experimentName, 
                countMatrix     = countMatrix, 
                colData         = columnsMetaData,
                species         = "mmu",
                outputDirectory = outputDirectory)

## Normalization with S4 method
scrS4 <- normaliseCountMatrix(scr, colData = columnsMetaData) # return NCM
getNormalizedCountMatrix(scr)    # Slot normalizedCountMatrix empty
getNormalizedCountMatrix(scrS4)  # Slot normalizedCountMatrix full
