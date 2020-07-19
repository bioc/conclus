#' .mkOutlierScoreDf
#'
#' @param mat Matrix row : dbscan solutions, column : cells
#'
#' @return Data frame of cells with associated outliers score
.mkOutlierScoreDf <- function(mat){
    
    outlierScoreDf <- as.data.frame(colnames(mat))
    colnames(outlierScoreDf) <- "cellName"
    outlierScoreDf <- dplyr::mutate(outlierScoreDf, outlierScore=NA)
    
    for(i in 1:ncol(mat)){
        vec <- mat[, i]
        outlierScoreDf$outlierScore[i] <- length(vec[vec == 0])
    }
    
    outlierScoreDf$outlierScorePer <- outlierScoreDf$outlierScore / nrow(mat)
    return(outlierScoreDf)
}



#' .excludeOutliers
#' 
#' @description Exclude outliers cells from Dbscan clustering by creating
#' a dataframe with outliers score and applying the threshold to remove them
#' 
#'
#' @param theObject An Object of class scRNASeq for which 
#' ?runDBSCAN was used.
#' @param threshold Threshold to remove outliers cells.
#' @param minPoints Reachability minimum no. of points parameter of 
#' fpc::dbscan() function. See Ester et al. (1996) for more details. 
#' Default = c(3, 4)
#' @param epsilon Reachability distance parameter of fpc::dbscan() function.
#' See Ester et al. (1996) for more details. Default = c(1.3, 1.4, 1.5)
#'
#' @return
#' @export
#'
#' @examples
.excludeOutliers <- function(theObject, threshold=0.3, minPoints, epsilon){
    # exclude outliers based on DBSCAN clustering
    # outliers are the cells which cannot be assigned
    # to any final cluster
    
    dbscanList <- getDbscanList(theObject)
    sceObject <- getSceNorm(theObject)
    
    
    ## Transform dbscan list to matrix
    l <- lapply(dbscanList, function(element){
        clustering <- getClustering(element)
        return(clustering)
    })
    dbscanMatrix <- do.call(rbind, l)
    outlierInfo <- .mkOutlierScoreDf(dbscanMatrix)
    colData <- 
        SummarizedExperiment::colData(
            sceObject)[SummarizedExperiment::colData(sceObject)$cellName
                       %in% colnames(dbscanMatrix), ]
    
    if(is.vector(colData)){
        colData <- S4Vectors::DataFrame(cellName=colData, row.names=colData)
    }
    
    colData$outlierScorePer <- NULL
    colData$outlierScore <- NULL
    
    colData <- merge(colData, outlierInfo,
                     by.x="cellName", by.y="cellName",
                     all.x=TRUE, all.y=TRUE, sort=FALSE)
    rownames(colData) <- colData$cellName
    
    numberOfCellsBefore <- dim(colData)[1]
    print(threshold)
    sceObject <- sceObject[, colData$outlierScorePer < threshold]
    colData <- colData[colData$outlierScorePer < threshold, ]
    numberOfCellsAfter <- dim(colData)[1]
    
    dbscanMatrix <- dbscanMatrix[, outlierInfo$outlierScorePer < threshold]
    SummarizedExperiment::colData(sceObject) <- colData
    
    message(
        paste(numberOfCellsBefore - numberOfCellsAfter,
              "outliers were excluded from the SingleCellExperiment object.\n"))
    
    setSceNorm(theObject) <- sceObject
    
    return(theObject)
}







# .checkRunConclus <- function(deleteOutliers, manualClusteringObject){
#     
#     if(!is.logical(deleteOutliers))
#         stop("'deleteOutliers' should be boolean.")
#         
#     if(!is.na(manualClusteringObject) && 
#        class(manualClusteringObject) != "scRNAseq")
#         stop("'manualClusteringObject' should be an scRNAseq object of CONCLUS.")
#     
#     clusters <- colData(getSceNorm(manualClusteringObject))$clusters
#     
#     if (is.null(clusters))
#         stop(paste("'manualClusteringObject' should got with addClustering ",
#                    "method."))  
# }

#' runCONCLUS
#'
#' @description This function performs the core CONCLUS workflow. 
#' 
#' @details This function performs the following steps:
#'  1) It generates PCA and t-SNE coordinates
#'  2) runs DBSCAN 
#'  3) Calculates similarity matrices of cells and clusters
#'  4) Assigns cells to clusters
#'  5) Searches for positive markers for each cluster
#'  6) Saves plots and tables into outputDirectory.
#'  
#'  columnsMetaData -- Dataframe containing three columns:
#'  cellName, state, and cellBarcode.
#'  Not used if manualClusteringObject is defined.
#'
#' colorPalette/statePalette -- A vector of colors for clusters/states or 
#''default' value. If 'default' is selected, the number of clusters is limited
#' to 16. 
#' If an error message is thrown, re-run the function with your own color vector.  
#'
#' manualClusteringObject -- After running once runCONCLUS, one could which to
#' modify the obtained clusters manually. This is achieved with the function
#' addClusteringManually'. The result of 'addClusteringManually' should be passed
#' to the manualClusteringObject parameter to re-run CONCLUS 
#' on the new defined clusters.
#'
#' 
#' @param outputDirectory CONCLUS will create this directory if it doesn't exist
#' and store there all output files.
#' @param experimentName Prefix used for output files.
#' @param countMatrix The count matrix with raw values or UMI
#' @param columnsMetaData {A data frame with information about cells. Not used
#' if manualClusteringObject is defined. See details. Default = NA.
#' @param species Currently limited to human and mouse. Possible values are 
#' either 'mmu' or 'human'. Not used if manualClusteringObject is defined. 
#' Default = NA.
#' @param colorPalette A vector of colors for clusters. Default = "default",
#'  see details.
#' @param statePalette A vector of colors for states or conditions. 
#' Default = "default", See details.
#' @param clusteringMethod Clustering method passed to hclust() function.
#'  See ?hclust for a list of method. Default = "ward.D2"
#' @param epsilon Reachability distance parameter of fpc::dbscan() function.
#'  See Ester et al. (1996) for more details. Default = c(1.3, 1.4, 1.5)
#' @param minPoints Reachability minimum no. of points parameter of 
#' fpc::dbscan() function. See Ester et al. (1996) for more details. 
#' Default = c(3, 4)
#' @param PCs {a vector of first principal components. For example, to take
#'  ranges 1:5 and 1:10 write c(5, 10). Default = c(4, 6, 8, 10, 20, 40, 50)
#' @param perplexities Numeric scalar defining the perplexity parameter,
#'  see ‘?Rtsne’ for more details. Default = c(30, 40)
#' @param randomSeed Random seed for reproducibility. Default = 42.
#' @param clusterNumber Exact number of cluster. Default = 0 that will determine
#' the number of clusters automatically.
#' @param deepSplit Intuitive level of clustering depth.
#' Options are 1, 2, 3, 4. Default = 4
#' @param preClustered Boolean precising if DBSCAN is run to calculate similarity
#'  matrices. Should be TRUE if manualClusteringObject is defined. 
#'  Default="FALSE"
#' @param orderClusters If True, clusters in the similarity matrix of cells will
#'  be ordered by name. Default = FALSE
#' @param cores maximum number of jobs that CONCLUS can run in parallel.
#'  Default = 1
#' @param plotPDFcellSim if TRUE, the similarity matrix of cells will be saved
#'  in pdf format; png if FALSE. FALSE is recommended for count matrices
#'   with more than 2500 cells due to large pdf file size. Default = TRUE
#' @param deleteOutliers Boolean indicating if whether cells which were often 
#' defined as outliers by dbscan must be deleted. It will require recalculating
#'  of the similarity matrix of cells. Default = FALSE.
#' @param removeDuplicates If TRUE, duplicated markers are removed from the
#'  lists. Default=TRUE.
#' @param tSNEalreadyGenerated TRUE if you already ran CONCLUS ones and have
#'  t-SNE coordinated saved. Default = FALSE
#' @param tSNEresExp experimentName of t-SNE coordinates which you want to use.
#'  This argument allows copying and pasting t-SNE coordinates between different
#'   CONCLUS runs without renaming the files. Default = ""
#' @param manualClusteringObject Result of the function addClusteringManually.
#'  Default = NA. See details
#'
#' @return \code{scRNAseq} object containing the similarity matrices and the 
#' marker genes. Write also marker genes in file for each clusters, tSNE and 
#' heatmaps.
#' 
#' @examples
#' 
#' experimentName <- "Bergiers"
#' countMatrix <- as.matrix(read.delim(file.path(
#' "tests/testthat/test_data/test_countMatrix.tsv")))
#' outputDirectory <- "./"
#' columnsMetaData <- read.delim(
#' file.path("extdata/Bergiers_colData_filtered.tsv"))
#' species <- "mouse"
#' 
#' runCONCLUS(outputDirectory=outputDirectory, countMatrix=countMatrix,
#'            columnsMetaData=columnsMetaData, species=species)
#' 
#' 
#' 
#' @seealso \code{normaliseCountMatrix}, \code{generateTSNECoordinates}
#' \code{runDBSCAN}, \code{clusterCellsInternal}, 
#' \code{calculateClustersSimilarity},  \code{rankGenes}
#' \code{retrieveTopClustersMarkers}, \code{retrieveGenesInfo}, 
#' \code{saveGenesInfo}, \code{plotClusteredTSNE}, \code{plotCellSimilarity},
#' \code{plotClustersSimilarity}, \code{exportMatrix},
#' @export
#' @author Ilyess RACHEDI
runCONCLUS <- function(
		## General parameters
		utputDirectory, experimentName, countMatrix, species, cores=1, 
		clusteringMethod="ward.D2", exportAllResults=TRUE, orderClusters=FALSE,
		
		## Normalisation parameters
		sizes=c(20,40,60,80,100), rowMetaData=NULL, columnsMetaData = NULL,
		alreadyCellFiltered=FALSE, runQuickCluster=TRUE,
		
		## tSNE parameters
		randomSeed = 42, PCs=c(4, 6, 8, 10, 20, 40, 50), perplexities=c(30,40),
		writeOutputTSne = FALSE,
		
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
		
		## plotCellHeatmap parameters
		meanCentered=TRUE, orderGenesCH=FALSE, savePlotCH=FALSE, widthCH=10,
		heightCH=8.5, clusterCols=FALSE, heightPlotCustSM=5.5,
		
		## plotClustersSimilarity parameters
		savePlotClustSM=FALSE, widthPlotCustSM=7, 
		
                       preClustered = FALSE,  
                       plotPDFcellSim = TRUE, deleteOutliers = TRUE,
                       tSNEalreadyGenerated = FALSE, tSNEresExp = "",
                       manualClusteringObject = NA){
    
	
	## Verify parameters
	!!
    # .checkRunConclus(deleteOutliers,manualClusteringObject)
	
	if(exportAllResults){
		
		writeOutputTSne <- FALSE
		writeOutputDbScan <- FALSE
		writeOutputRankGenes <- FALSE
		writeTopMarkers <- FALSE
		saveInfos <- FALSE
		writeCSM <- FALSE
		savePlotCTSNE <- FALSE
		savePlotCH <- FALSE
		savePlotClustSM <- FALSE
	}		
			
	message("## Building the single-cell RNA-Seq object ##")
    scr <- scRNAseq(experimentName = experimentName,
                    countMatrix     = countMatrix,
                    species         = species,
                    outputDirectory = outputDirectory)
    
	## Processing
	
    message("## Performing the normalization ##")
    scrNorm <- normaliseCountMatrix(scr, sizes=sizes, rowdata=rowMetaData, 
			coldata=columnsMetaData, alreadyCellFiltered=alreadyCellFiltered,
			runQuickCluster=runQuickCluster)
	
	message("## Calculating all tSNEs ##")
    scrTsne <- generateTSNECoordinates(scrNorm, randomSeed=randomSeed, 
			cores=cores, PCs=PCs, perplexities=perplexities, 
			writeOutput=writeOutputTSne)
	
	message("## Clustering with DbScan ##")
    scrDbscan <- runDBSCAN(scrTsne, cores=cores, epsilon=epsilon, 
			minPoints=minPoints, writeOutput=writeOutputDbScan)

    message("## Computing the cells similarity matrix ##")
    scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=clusterNumber,
                                   deepSplit=deepSplit, cores=cores,
                                   clusteringMethod=clusteringMethod)

	message("## Computing the clusters similarity matrix ##")
    scrCSM <- calculateClustersSimilarity(scrCCI, 
			clusteringMethod=clusteringMethod)

	## Markers
	
	message("## Ranking genes ##")
    scrS4MG <- rankGenes(scrCSM, column=columnRankGenes, 
			writeMarkerGenes=writeOutputRankGenes)
	
    message("## Getting marker genes ##")
    scrFinal <- retrieveTopClustersMarkers(scrS4MG, nTop=nTopMarkers, 
			removeDuplicates=removeDuplicates, writeMarkerGenes=writeTopMarkers)
    
	message("## Getting genes info ##")
    scrInfos <- retrieveGenesInfo(scrFinal, cores=cores, groupBy=groupBy,
			orderGenes=orderGenes, getUniprot=getUniprot, saveInfos=saveInfos)
    
	## Plotting
	
	message("## Plot the cell similarity matrix ##")
	plotCellSimilarity(scrInfos, colorPalette=colorPalette, 
			statePalette=statePalette, clusteringMethod=clusteringMethod,
			orderClusters=orderClusters, savePlot=writeCSM, widthCSM=7, 
			heightCSM=6)
	
	message("## Plot clustered tSNE ##")
	plotClusteredTSNE(scrInfos, colorPalette=colorPalette, PCs=PCs, 
			perplexities=perplexities, columnName=columnRankGenes, 
			savePlot=savePlotCTSNE, width=widthPlotClustTSNE, 
			height=heightPlotClustTSNE)

    message("## Plot the cell heatmap ##")
	plotCellHeatmap(scrInfos, meanCentered=meanCentered, 
			colorPalette=colorPalette, statePalette=statePalette, 
			clusteringMethod=clusteringMethod, orderClusters=orderClusters, 
			orderGenes=orderGenesCH, savePlot=savePlotCH, width=widthCH, 
			height=heightCH, clusterCols=clusterCols)
	
	message("## Plot the clusters similarity heatmap ##")
	plotClustersSimilarity(scrInfos, colorPalette=colorPalette, 
			statePalette=statePalette, clusteringMethod=clusteringMethod, 
			savePlot=savePlotClustSM, width=widthPlotCustSM, 
			height=heightPlotCustSM)
			  
	!! see the original code for the rest
	
	
	if(exportAllResults)
		exportResults(scrInfos, saveAll=TRUE)
	
    return(scr)
}
