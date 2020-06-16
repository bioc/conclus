if (!require("pacman")) install.packages("pacman")
pacman::p_load(rlist, foreach, ggplot2, pheatmap, zoo, dynamicTreeCut, factoextra,
               digest, RColorBrewer,devtools, BiocManager, BiocParallel, scran, scater,
               monocle, SingleCellExperiment , KEGGREST, AnnotationDbi,
               Cairo, rvest, curl,  Matrix, dbscan, fpc, matrixStats,
               dplyr, biomaRt, org.Mm.eg.db, grDevices, S4Vectors, Biobase,
               DataCombine, zoo, rvest, DataCombine, doParallel, testthat)

source("R/AllGenerics.R")
source("R/AllClasses.R")
source("R/sharedInternals.R")
source("R/getters.R")
source("R/setters.R")
source("R/methods-normalization.R")
source("R/methods-tsne.R")
source("R/methods-dbscan.R")
source("R/methods-clustering.R")
source("R/methods-plot.R")
source("R/methods-export.R")
source("R/methods-markers.R")

## Data
outputDirectory <- "YourOutputDirectory"
experimentName <- "Bergiers"
countMatrix <- as.matrix(read.delim(file.path("inst/extdata/Bergiers_counts_matrix_filtered.tsv"), stringsAsFactors = FALSE))
columnsMetaData <- read.delim(file.path("inst/extdata/Bergiers_colData_filtered.tsv"))

# countMatrix <- as.matrix(read.delim(
#     file.path("tests/testthat/test_data/test_countMatrix.tsv")))


## Construction
scrS4 <- scRNAseq(experimentName  = experimentName, 
                  countMatrix     = countMatrix, 
                  species         = "mouse",
                  outputDirectory = outputDirectory)


## Normalization with S4 method
scrS4Norm <- normaliseCountMatrix(scrS4, colData = columnsMetaData)


## Test Clustering 
# a <- testClustering(scrS4Norm)


## TSNE
scrS4TSNE <- generateTSNECoordinates(scrS4Norm, cores= 15)
# load(file="tests/testthat/test_data/expectedlistTSNEFromConclus.Rdat")
# setTSNEList(scrS4TSNE) <- expectedlistTSNEFromConclus


## DBSCAN
scrS4DBSCAN <- runDBSCAN(scrS4TSNE)
dblist <- getClustering(getDbscanList(scrS4DBSCAN)[[1]])


## clusterCellsInternal
scrS4CCI<- clusterCellsInternal(scrS4DBSCAN, clusterNumber=10, deepSplit=4)


## calculateClustersSimilarity
scrS4Clusters <- calculateClustersSimilarity(scrS4CCI)
getClustersSimilarityMatrix(scrS4Clusters)

## Export
exportResults(scrS4Clusters)

## Marker Genes
scrS4MG <- rankGenes(scrS4Clusters)
getMarkerGenesList(scrS4MG)
markersClusters <- getMarkerGenes(scrS4MG)
infos <- getGenesInfo(gene=markersClusters, species = "mouse", cores = 15)

## Plotting
plotCellSimilarity(scrS4Clusters)
plotClusteredTSNE(scrS4Clusters)
# plotCellHeatmap(scrS4Clusters)  Need geneMarker slot
plotGeneExpression(scrS4Clusters, geneName="Ccl3")
plotClustersSimilarity(scrS4Clusters)


## Save for unit test
# scrFull <- scrS4Clusters
# save(scrFull, file="tests/testthat/test_data/scrFull.Rdat")


