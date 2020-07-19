#' runCONCLUS
#'
#' @description This function is a wrapper to run the whole CONCLUS workflow. 
#' See details.
#'
#' @usage 
#' runCONCLUS(
#' 		## General parameters
#' 		utputDirectory, experimentName, countMatrix, species, cores=1, 
#' 		clusteringMethod="ward.D2", exportAllResults=TRUE, orderClusters=FALSE,
#' 		clusToAdd=NA,
#' 		
#' 		## Normalisation parameters
#' 		sizes=c(20,40,60,80,100), rowMetaData=NULL, columnsMetaData = NULL,
#' 		alreadyCellFiltered=FALSE, runQuickCluster=TRUE,
#' 		
#' 		## tSNE parameters
#' 		randomSeed = 42, PCs=c(4, 6, 8, 10, 20, 40, 50), perplexities=c(30,40),
#' 		writeOutputTSne = FALSE,
#' 		
#' 		## Dbscan parameters
#' 		epsilon=c(1.3, 1.4, 1.5), minPoints=c(3, 4), writeOutputDbScan=FALSE,
#' 		
#' 		## Cell Similarity matrix parameters
#' 		clusterNumber=10, deepSplit=4,
#' 		
#' 		## Rank genes parameters
#' 		columnRankGenes="clusters", writeOutputRankGenes=FALSE,
#' 		
#' 		## Retrieving top markers parameters
#' 		nTopMarkers=10, removeDuplicates = TRUE, writeTopMarkers=FALSE,
#' 		
#' 		## Retrieving genes infos parameters
#' 		groupBy="clusters", orderGenes="initial", getUniprot=TRUE, 
#' 		saveInfos=FALSE, 
#' 		
#' 		## plotCellSimilarity parameters
#' 		colorPalette="default", statePalette="default", writeCSM=FALSE, 
#' 		widthCSM=7, heightCSM=6,
#' 		
#' 		## plotClusteredTSNE parameters
#' 		savePlotCTSNE=FALSE, widthPlotClustTSNE=6, heightPlotClustTSNE=5,
#' 		
#' 		## plotCellHeatmap parameters
#' 		meanCentered=TRUE, orderGenesCH=FALSE, savePlotCH=FALSE, widthCH=10,
#' 		heightCH=8.5, clusterCols=FALSE, heightPlotCustSM=5.5,
#' 		
#' 		## plotClustersSimilarity parameters
#' 		savePlotClustSM=FALSE, widthPlotCustSM=7)
#' 
#' 
#' @param 
#'
#'
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
		clusToAdd=NA,
		
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
		savePlotClustSM=FALSE, widthPlotCustSM=7){
    
	
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
			
	message("## Building the single-cell RNA-Seq object (step 1/13) ##")
    scr <- scRNAseq(experimentName = experimentName,
                    countMatrix     = countMatrix,
                    species         = species,
                    outputDirectory = outputDirectory)
    
	## Processing
	
    message("## Performing the normalization (step 2/13) ##")
	message("\t Note: The connection to biomaRt can take a while sometimes.")
    scrNorm <- normaliseCountMatrix(scr, sizes=sizes, rowdata=rowMetaData, 
			coldata=columnsMetaData, alreadyCellFiltered=alreadyCellFiltered,
			runQuickCluster=runQuickCluster)
	
	message("## Calculating all tSNEs (step 3/13) ##")
    scrTsne <- generateTSNECoordinates(scrNorm, randomSeed=randomSeed, 
			cores=cores, PCs=PCs, perplexities=perplexities, 
			writeOutput=writeOutputTSne)
	
	message("## Clustering with DbScan (step 4/13) ##")
    scrDbscan <- runDBSCAN(scrTsne, cores=cores, epsilon=epsilon, 
			minPoints=minPoints, writeOutput=writeOutputDbScan)

    message("## Computing the cells similarity matrix (step 5/13) ##")
    scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=clusterNumber,
                                   deepSplit=deepSplit, cores=cores,
                                   clusteringMethod=clusteringMethod)

	message("## Computing the clusters similarity matrix (step 6/13) ##")
    scrCSM <- calculateClustersSimilarity(scrCCI, 
			clusteringMethod=clusteringMethod)

	if(!is.na(clusToAdd)){
		message("Adding the provided clustering manually.")
		scrCSM <- addClustering(scrCSM, clusToAdd=clusToAdd)
	}
		
	## Markers
	
	message("## Ranking genes (step 7/13) ##")
    scrS4MG <- rankGenes(scrCSM, column=columnRankGenes, 
			writeMarkerGenes=writeOutputRankGenes)
	
    message("## Getting marker genes (step 8/13) ##")
    scrFinal <- retrieveTopClustersMarkers(scrS4MG, nTop=nTopMarkers, 
			removeDuplicates=removeDuplicates, writeMarkerGenes=writeTopMarkers)
    
	message("## Getting genes info (step 9/13) ##")
	message("\t Note: The connection to biomaRt can take a while sometimes.")
    scrInfos <- retrieveGenesInfo(scrFinal, cores=cores, groupBy=groupBy,
			orderGenes=orderGenes, getUniprot=getUniprot, saveInfos=saveInfos)
    
	## Plotting
	
	message("## Plot the cell similarity matrix (step 10/13) ##")
	plotCellSimilarity(scrInfos, colorPalette=colorPalette, 
			statePalette=statePalette, clusteringMethod=clusteringMethod,
			orderClusters=orderClusters, savePlot=writeCSM, widthCSM=7, 
			heightCSM=6)
	
	message("## Plot clustered tSNE (step 11/13) ##")
	plotClusteredTSNE(scrInfos, colorPalette=colorPalette, PCs=PCs, 
			perplexities=perplexities, columnName=columnRankGenes, 
			savePlot=savePlotCTSNE, width=widthPlotClustTSNE, 
			height=heightPlotClustTSNE)

    message("## Plot the cell heatmap (step 12/13) ##")
	plotCellHeatmap(scrInfos, meanCentered=meanCentered, 
			colorPalette=colorPalette, statePalette=statePalette, 
			clusteringMethod=clusteringMethod, orderClusters=orderClusters, 
			orderGenes=orderGenesCH, savePlot=savePlotCH, width=widthCH, 
			height=heightCH, clusterCols=clusterCols)
	
	message("## Plot the clusters similarity heatmap (step 13/13) ##")
	plotClustersSimilarity(scrInfos, colorPalette=colorPalette, 
			statePalette=statePalette, clusteringMethod=clusteringMethod, 
			savePlot=savePlotClustSM, width=widthPlotCustSM, 
			height=heightPlotCustSM)
	
	
	if(exportAllResults){
		message("Exporting all results to ", outputDirectory)
		exportResults(scrInfos, saveAll=TRUE)
	}
		
	message("Done.")
	
    return(scrInfos)
}
