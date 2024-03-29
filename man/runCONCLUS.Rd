% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runCONCLUS.R
\name{runCONCLUS}
\alias{runCONCLUS}
\title{runCONCLUS}
\usage{
runCONCLUS(
        ## General parameters
        outputDirectory, experimentName, countMatrix, species, cores=2,
        clusteringMethod="ward.D2", exportAllResults=TRUE,
        orderClusters=FALSE, clusToAdd=NA, silentPlot=TRUE,

        ## Normalisation parameters
        sizes=c(20,40,60,80,100), rowMetaData=NULL, columnsMetaData = NULL,
        alreadyCellFiltered=FALSE, runQuickCluster=TRUE, info=TRUE,

        ## tSNE parameters
        randomSeed = 42, PCs=c(4, 6, 8, 10, 20, 40, 50),
        perplexities=c(30,40), writeOutputTSne = FALSE,

        ## Dbscan parameters
        epsilon=c(1.3, 1.4, 1.5), minPoints=c(3, 4), writeOutputDbScan=FALSE,

        ## Cell Similarity matrix parameters
        clusterNumber=10, deepSplit=4,

        ## Rank genes parameters
        columnRankGenes="clusters", writeOutputRankGenes=FALSE,

        ## Retrieving top markers parameters
        nTopMarkers=10, removeDuplicates = TRUE, writeTopMarkers=FALSE,

        ## Retrieving genes infos parameters
        groupBy="clusters", orderGenes="initial", getUniprot=TRUE,
        saveInfos=FALSE,

        ## plotCellSimilarity parameters
        colorPalette="default", statePalette="default", writeCSM=FALSE,
        widthCSM=7, heightCSM=6,

        ## plotClusteredTSNE parameters
        savePlotCTSNE=FALSE, widthPlotClustTSNE=6, heightPlotClustTSNE=5,
        tSNENb=NA,

        ## plotCellHeatmap parameters
        meanCentered=TRUE, orderGenesCH=FALSE, savePlotCH=FALSE, widthCH=10,
        heightCH=8.5, clusterCols=FALSE,

        ## plotClustersSimilarity parameters
        savePlotClustSM=FALSE, widthPlotClustSM=7, heightPlotClustSM=5.5)
}
\arguments{
\item{outputDirectory}{Directory to which results should be written. This
needs to be defined even if you choose to not output any results.}

\item{experimentName}{String of the name of the experiment.}

\item{countMatrix}{Matrix containing the raw counts.}

\item{species}{Character string of the species of interest. Shoud be mouse or
human. Other organisms can be added on demand.}

\item{cores}{Maximum number of jobs that CONCLUS can run in parallel. This
parameter is used by ?generateTSNECoordinates, ?runDBSCAN,
?clusterCellsInternal, and ?retrieveGenesInfo. Default=1.}

\item{clusteringMethod}{Clustering method passed to hclust() function. See
?hclust for a list of method. This parameter is used by
?clusterCellsInternal, ?calculateClustersSimilarity, ?plotCellSimilarity,
?plotClusteredTSNE, ?plotCellHeatmap, and ?plotClustersSimilarity.
Default = "ward.D2".}

\item{exportAllResults}{If TRUE, Save all results of CONCLUS. See
?exportResults for details. Default=TRUE.}

\item{orderClusters}{If TRUE, clusters in the cells and clusters similarity
matrix of cells will be ordered by name. Default = FALSE.}

\item{clusToAdd}{If not NA, defines the clustering to be used in theObject.
This is particularly useful when one wants to compare the clustering
performance of different tools. It should be a data frame having two columns
'clusters' and 'cells'. Default=NA.}

\item{silentPlot}{Boolean indicating if the figures should not be output on
the R graphics. Default=TRUE.}

\item{sizes}{Vector of size factors from scran::computeSumFactors() function
used by ?normaliseCountMatrix.}

\item{rowMetaData}{Data frame containing genes informations. Default is NULL.
See ?normaliseCountMatrix.}

\item{columnsMetaData}{Data frame containing cells informations.
Default is NULL. See ?normaliseCountMatrix.}

\item{alreadyCellFiltered}{If TRUE, quality check and filtering will not be
applied during the normalization of the count matrix.
See ?normaliseCountMatrix.}

\item{runQuickCluster}{If TRUE scran::quickCluster() function will
be applied. It usually improves the normalization for medium-size count
matrices. However, it is not recommended for datasets with less than 200
cells and may take too long for datasets with more than 10000 cells.
Default=TRUE. See ?normaliseCountMatrix.}

\item{info}{Logical. If TRUE, additional annotations like ensembl_gene_id,
go_id, name_1006, chromosome_name and gene_biotype are added to the 
row data, for all the genes from the count matrix with ENSEMBL IDs or 
SYMBOL ID. Default: TRUE.}

\item{randomSeed}{Default is 42. Seeds used to generate the tSNE. See
?generateTSNECoordinates.}

\item{PCs}{Vector of first principal components. For example, to take ranges
1:5 and 1:10 write c(5, 10). Default = c(4, 6, 8, 10, 20, 40, 50). See
?generateTSNECoordinates.}

\item{perplexities}{A vector of perplexity (t-SNE parameter). See
?generateTSNECoordinates for details. Default = c(30, 40).}

\item{writeOutputTSne}{If TRUE, write the tsne parameters to the output
directory defined in theObject. Default = FALSE. Ignored if
exportAllResults=TRUE.}

\item{epsilon}{Reachability distance parameter of fpc::dbscan() function.
See Ester et al. (1996) for more details. Default = c(1.3, 1.4, 1.5).}

\item{minPoints}{Reachability minimum no. of points parameter of
fpc::dbscan() function. See Ester et al. (1996) for more details.
Default = c(3, 4).}

\item{writeOutputDbScan}{If TRUE, write the results of the dbScan clustering
to the output directory defined in theObject, in the sub-directory
output_tables. Default = FALSE. Ignored if exportAllResults=TRUE.}

\item{clusterNumber}{Exact number of cluster. Default = NULL
that will determine the number of clusters automatically. 
See ?clusterCellsInternal.}

\item{deepSplit}{Intuitive level of clustering depth. Options are 1, 2, 3, 4.
See ?clusterCellsInternal. Default = 4.}

\item{columnRankGenes}{Name of the column with a clustering result. See
?rankGenes. Default="clusters".}

\item{writeOutputRankGenes}{If TRUE, output one list of marker genes per
cluster in the output directory defined in theObject and in the sub-directory
'marker_genes'. Default=FALSE. Ignored if exportAllResults=TRUE.}

\item{nTopMarkers}{Number of marker genes to retrieve per cluster. See
?retrieveTopClustersMarkers. Default=10.}

\item{removeDuplicates}{If TRUE, duplicated markers are removed from the
lists. See ?retrieveTopClustersMarkers. Default=TRUE.}

\item{writeTopMarkers}{If TRUE, writes one list per cluster in the output
folder defined in theObject, and in the sub-directory
marker_genes/markers_lists. Default=FALSE. Ignored if exportAllResults=TRUE.}

\item{groupBy}{A column in the input table used for grouping the genes in
the output tables. This option is useful if a table contains genes from
different clusters. See ?retrieveGenesInfo. Default = "clusters".}

\item{orderGenes}{If "initial" then the order of genes will not be changed.
The other option is "alphabetical". See ?retrieveGenesInfo.
Default="initial".}

\item{getUniprot}{Boolean, whether to get information from UniProt or not.
See ?retrieveGenesInfo. Default = TRUE.}

\item{saveInfos}{If TRUE, save the genes infos table in the directory
defined in theObject (?getOutputDirectory) and in the sub-directory
'marker_genes/saveGenesInfo'. Default=FALSE. Ignored if
exportAllResults=TRUE.}

\item{colorPalette}{A vector of colors for clusters. This parameter is used
by all plotting methods. Default = "default". See ?plotClustersSimilarity
for details.}

\item{statePalette}{A vector of colors for states or conditions. This
parameter is used by all plotting functions except ?plotClusteredTSNE.
See ?plotClustersSimilarity for details.}

\item{writeCSM}{If TRUE, the cells similarity heatmap is saved in the
directory defined in theObject (?getOutputDirectory) and in the sub-directory
"pictures". Default=FALSE. Ignored if exportAllResults=TRUE.}

\item{widthCSM}{Width of the plot in the pdf file. See ?pdf for more details.
Default = 7.}

\item{heightCSM}{Height of the plot in the pdf file. See ?pdf for more
details. Default = 6.}

\item{savePlotCTSNE}{If TRUE, the heatmap of the clustered tSNE is saved in
the directory defined in theObject (?getOutputDirectory) and in the
sub-directory "pictures/tSNE_pictures". Default=FALSE. Ignored if
exportAllResults=TRUE.}

\item{widthPlotClustTSNE}{Width of the clustered tSNE plot in the pdf file.
See ?pdf for more details. Default = 6.}

\item{heightPlotClustTSNE}{Height of the clustered tSNE plot in the pdf file.
See ?pdf for more details. Default = 5.}

\item{tSNENb}{Give the number of the tSNE to plot. If NA, all tSNE solutions
are plotted (14 tSNE by default). Default=NA.}

\item{meanCentered}{Boolean indicating if mean centering should be applied
to the expression matrix. See ?plotCellHeatmap. Default = TRUE.}

\item{orderGenesCH}{Boolean, should the heatmap be structured by gene. See
?plotCellHeatmap. Default=FALSE.}

\item{savePlotCH}{If TRUE save the cell heatmap in pdf format.
The heatmap is saved in the output directory defined in theObject
(?getOutputDirectory) and in the sub-directory 'pictures'. Default=FALSE.
Ignored if exportAllResults=TRUE.}

\item{widthCH}{Width of the cell heatmap saved in ?pdf. Default = 10.}

\item{heightCH}{Height of the cell heatmap saved in ?pdf. Default = 8.5.}

\item{clusterCols}{If TRUE, the columns representing the clusters are also
taken into account in the hierarchical clustering of the cell heatmap.
Default=FALSE.}

\item{savePlotClustSM}{If TRUE, save the cluster similarity heatmap in pdf
format. The heatmap is saved in the output directory defined in theObject
(?getOutputDirectory) and in the sub-directory 'pictures'. Default=FALSE.
Ignored if exportAllResults=TRUE.}

\item{widthPlotClustSM}{Width of the clusters similarity heatmap in the pdf
file. See ?pdf for more details. Default = 7.}

\item{heightPlotClustSM}{Height of the clusters similarity heatmap in the pdf
file. See ?pdf for more details. Default = 5.5.}
}
\value{
A \code{scRNAseq} object containing the similarity matrices and the
marker genes.
}
\description{
This function is a wrapper to run the whole CONCLUS workflow.
See details.
}
\details{
CONCLUS is a tool for robust clustering and positive marker features
selection of single-cell RNA-seq (sc-RNA-seq) datasets. Of note, CONCLUS does
not cover the preprocessing steps of sequencing files obtained following
next-generation sequencing.

CONCLUS is organized into the following steps:

1) Generation of multiple t-SNE plots with a range of parameters including
different selection of genes extracted from PCA. \cr
2) Use the Density-based spatial clustering of applications with noise
(DBSCAN) algorithm for idenfication of clusters in each generated t-SNE plot.
3) All DBSCAN results are combined into a cell similarity matrix. \cr
4) The cell similarity matrix is used to define "CONSENSUS" clusters
conserved accross the previously defined clustering solutions. \cr
5) Identify marker genes for each concensus cluster. cr

This wrapper function performs the following steps:

1) Building the single-cell RNA-Seq object. See ?scRNAseq-class. \cr
2) Performing the normalization. See ?normaliseCountMatrix. \cr
3) Calculating all tSNEs. See ?generateTSNECoordinates. \cr
4) Clustering with DbScan. See ?runDBSCAN. \cr
5) Computing the cells similarity matrix. See ?clusterCellsInternal. \cr
6) Computing the clusters similarity matrix. If clusToAdd is not NA, add
the provided clustering. See ?calculateClustersSimilarity and
?addClustering.  \cr
7) Ranking genes. See ?rankGenes.  \cr
8) Getting marker genes. See ?retrieveTopClustersMarkers.  \cr
9) Getting genes info. See ?retrieveGenesInfo.  \cr
10) Plot the cell similarity matrix. See ?plotCellSimilarity.  \cr
11) Plot clustered tSNE. See ?plotClusteredTSNE.  \cr
12) Plot the cell heatmap. See ?plotCellHeatmap.  \cr
13) Plot the clusters similarity heatmap. See ?plotClustersSimilarity.  \cr
14) Exporting all results to outputDirectory if exportAllResults=TRUE.
See ?exportAllResults.  \cr
15) Return an object containing all the results provided by CONCLUS.  \cr

If exportAllResults=TRUE, in your "outputDirectory", the sub-folder
pictures contains all tSNE with dbscan coloration (sub-folder
tSNE_pictures), the cell similarity matrix
(Test_cells_correlation_X_clusters.pdf), the cell heatmap
(Test_clustersX_meanCenteredTRUE_orderClustersFALSE_orderGenesFALSE
markrsPerCluster.pdf`), and the cluster similarity matrix
(`Test_clusters_similarity_10_clusters.pdf`). You will also find in the
sub-folder `Results`:

+ `1_MatrixInfo`: The normalized count matrix and its meta-data for both
rows and columns. \cr
+ `2_TSNECoordinates`: The tSNE coordinates for each parameter of principal
components (PCs) and perplexities.  \cr
+ `3_dbScan`: The different clusters given by DBscan according to different
parameters. Each file gives a cluster number for each cell.  \cr
+ `4_CellSimilarityMatrix`: The matrix underlying the cells similarity
heatmap. \cr
+ `5_ClusterSimilarityMatrix`: The matrix underlying the clusters similarity
heatmap.  \cr
+ `6_ConclusResult`: A table containing the result of the consensus
clustering. This table contains two columns: clusters-cells. \cr
+ `7_fullMarkers`: Files containing markers for each cluster, defined by the
consensus clustering. \cr
+ `8_TopMarkers`: Files containing the top 10 markers for each cluster. \cr
+ `9_genesInfos`: Files containing gene information for the top markers
defined in the previous folder. \cr
}
\examples{
experimentName <- "Bergiers"
outputDirectory <- "YourOutputDirectory"
species <- "mouse"

## Load the count matrix
countmatrixPath <- system.file("extdata/countMatrix.tsv", package="conclus")
countMatrix <- loadDataOrMatrix(file=countmatrixPath, type="countMatrix",
                                ignoreCellNumber=TRUE)

## Load the coldata
coldataPath <- system.file("extdata/colData.tsv", package="conclus")
columnsMetaData <- loadDataOrMatrix(file=coldataPath, type="coldata",
columnID="cell_ID")

## Use runCONCLUS
## These parameters are tweaked to fit our example data and reduce
## computing time, please consider using the default parameters or
## adjusted to your dataset.
scr <- runCONCLUS(outputDirectory, experimentName, countMatrix, species,
        columnsMetaData=columnsMetaData, perplexities=c(2,3), tSNENb=1,
        PCs=c(4,5,6,7,8,9,10), epsilon=c(380, 390, 400), minPoints=c(2,3),
        clusterNumber=2)

## Remove the results
unlink(outputDirectory, recursive=TRUE)

}
\author{
Nicolas Descostes
}
