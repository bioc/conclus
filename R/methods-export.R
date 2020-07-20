#' .exportGenesInfos
#'
#' @description 
#' Export the genes information tables per cluster.
#'
#' @param theObject scRNASeq object on which the methods ?retrieveGenesInfo 
#' was run.
#' @param outputInfos Output directory defined as dataDirectory/9_genesInfos. 
#' 
#' @keywords internal
#' @noRd
.exportGenesInfos <- function(theObject, outputInfos){

	infos <- getGenesInfos(theObject)
	infosList <- split(infos, infos$clusters)
	
	invisible(lapply(infosList, function(clusterDF){
						
						clustNB <- unique(clusterDF$clusters)
						
						if(!isTRUE(all.equal(length(clustNB), 1)))
							stop("Error in saveGenesInfo. Contact the",
									" developper.")
						
						outputFile <- file.path(outputInfos, 
								paste0("genes_info_clust", clustNB, 
										".csv"))
						write.csv(clusterDF, file=outputFile, 
								row.names=FALSE)
						
					}))
	
	message("Genes infos saved.")
}


#' .exportMarkers
#'
#' @description 
#' Export the full markers lists (?rankGenes) or the top markers lists 
#' (?retrieveTopClustersMarkers).
#'
#' @param theObject scRNASeq object on which the methods ?rankGenes or 
#' ?retrieveTopClustersMarkers were run.
#' @param outputDir Output directory defined as dataDirectory/7_fullMarkers or 
#' dataDirectory/8_TopMarkers. dataDirectory is directly retrieved from the 
#' scRNASeq object.
#' @param listType String displayed in the final message to indicate the type of
#'  information retrieved.
#' @param top Default=FALSE. If TRUE, exports the top markers otherwise the 
#' complete lists. 
#' 
#' @keywords internal
#' @noRd
.exportMarkers <- function(theObject, outputDir, listType="Full marker lists", 
		top=FALSE){
	
	if(top){
		
		markers <- getClustersMarkers(theObject)
		firstElement <- markers[1,"geneName"]
		slotName <- "clustersMarkers"
		funName <- "retrieveTopClustersMarkers"
		markers <- split(markers, markers$clusters)
	}else{
		
		markers <- getMarkerGenesList(theObject)
		firstElement <- markers[[1]][1,"Gene"]
		slotName <- "markerGenesList"
		funName <- "rankGenes"
	}
	
	
	if(isTRUE(all.equal(firstElement, "gene1")))
		stop("The 'scRNAseq' object that you're using with 'exportResults' ",
				"function doesn't have its '", slotName, "' slot updated. ",
				"Please use '", funName,"' on the object before.")
		
	invisible(mapply(function(currentMarkers, clustNb){
						
						fileName <- paste0("markers_clust", clustNb, ".csv")
						write.csv(currentMarkers, file=file.path(outputDir, 
										fileName))
					}, markers, seq_len(length(markers))))
	
	
	message(listType, " saved.")
	
}

#' .conclusResult
#'
#' @description 
#' Export the result of the concensus clustering that is defined in the row 
#' and col metadata of the sceNorm slot of the scRNASeq object.
#'
#' @param theObject scRNASeq object on which the method 
#' ?calculateClustersSimilarity was run.
#' @param experimentName Name of the experiment retrieved in the corresponding 
#' slot of the scRNASeq object.
#' @param sceObject SingleCellExperiment object retrieved with ?getSceNorm. See
#'  ?SingleCellExperiment::SingleCellExperiment for more details. 
#' @param outputDir Output directory defined as dataDirectory/6_ConclusResult. 
#' dataDirectory is directly retrieved from the scRNASeq object.
#' 
#' @keywords internal
#' @noRd
.conclusResult <- function(theObject, sceObject, outputDir){
	
	## Check that the cluster similarity matrix was computed
	## This step is necessary to  update the 'clusters' column
	## of the col metadata.
	matrix <- getClustersSimilarityMatrix(theObject)
	
	if(all(dim(matrix) == c(1,1)))
		stop("The 'scRNAseq' object that you're using with 'exportResults' ",
				"function doesn't have its columns metadata ",
				"updated. Please use 'calculateClustersSimilarity' on the ",
				"object before.")
	
	colMetaDf <- SummarizedExperiment::colData(sceObject)
	results <- data.frame(clusters=colMetaDf$clusters, 
			cells=colMetaDf$cellName)
	
	fileName <- paste0(experimentName,"_clusters_table.tsv")
	write.table(results, file=file.path(outputDir, fileName), quote=FALSE, 
			sep="\t", row.names=FALSE, col.names=TRUE)
	message("Clusters table saved.")
}


#' .exportCSM
#'
#' @description 
#' Export the cells and clusters similarity matrices that were obtained with 
#' ?clusterCellsInternal and ?calculateClustersSimilarity respectively.
#'
#' @param theObject scRNASeq object on which the methods mentioned above were 
#' run.
#' @param experimentName Name of the experiment retrieved in the corresponding 
#' slot of the scRNASeq object.
#' @param outputDir Output directory defined as 
#' dataDirectory/4_CellSimilarityMatrix or 
#' dataDirectory/5_ClusterSimilarityMatrix. dataDirectory is directly retrieved 
#' from the scRNASeq object.
#' @param cell If TRUE, exports the cells similarity matrix otherwise 
#' exports the clusters similarity matrix. Default=TRUE.
#' 
#' @keywords internal
#' @noRd
.exportCSM <- function(theObject, experimentName, outputDir, cell=TRUE){
	
	if(cell){
		
		matrix <- getCellsSimilarityMatrix(theObject)
		matType <- "cellsSimilarityMatrix"
		funName <- "clusterCellsInternal"
		fileName <- "_cells_similarity_matrix.csv"
		
	}else{
		matrix <- getClustersSimilarityMatrix(theObject)
		matType <- "clustersSimilarityMatrix"
		funName <- "calculateClustersSimilarity"
		fileName <- "_clusters_similarity_matrix.csv" 
	}
		
	if(all(dim(matrix) == c(0,0)) || all(dim(matrix) == c(1,1)))
		stop("The 'scRNAseq' object that you're using with 'exportResults' ",
				"function doesn't have its '", matType,"' slot ",
				"updated. Please use '", funName,"' on the object ",
				"before.")
		
	fileName <- paste0(experimentName, fileName)
	write.csv(matrix, file=file.path(outputDir, fileName))
	message(matType, " saved.")
}


#' .exportDBScan
#'
#' @description 
#' Export the dbscan clustering results for each combination of epsilon and 
#' minPts. These results were obtained with ?runDBSCAN.
#'
#' @param theObject scRNASeq object on which ?runDBSCAN was run.
#' @param outputDir Output directory defined as dataDirectory/3_dbScan. 
#' dataDirectory is directly retrieved from the scRNASeq object.
#' 
#' @keywords internal
#' @noRd
.exportDBScan <- function(theObject, outputDir){
	
	dbscanList <- getDbscanList(theObject)
	
	if(isTRUE(all.equal(length(getName(dbscanList[[1]])), 0)))
		stop("The 'scRNAseq' object that you're using with 'exportResults' ",
				"method doesn't have its 'dbscanList' slot updated. Please use",
				" 'runDBSCAN' on the object before.")
	
	invisible(lapply(dbscanList, function(currentDBScan){
						
						name <- getName(currentDBScan)
						eps <- getEpsilon(currentDBScan)
						mPt <- getMinPoints(currentDBScan)
						fileName <- paste(name, eps, mPt, sep="_")
						fileName <- paste0(fileName, ".tsv")
						write.table(getClustering(currentDBScan), 
								file=file.path(outputDir, fileName),
								quote=FALSE, row.names=TRUE, col.names=TRUE)
					}))
	
	message("dbScan clustering saved.")
}


#' .exportTsne
#'
#' @description 
#' Export the tsne coordinates for each combination of PCs and perplexities.
#' These results were obtained with ?generateTSNECoordinates.
#'
#' @param theObject scRNASeq object on which ?generateTSNECoordinates was run.
#' @param outputDir Output directory defined as dataDirectory/2_TSNECoordinates. 
#' dataDirectory is directly retrieved from the scRNASeq object.
#' 
#' @keywords internal
#' @noRd
.exportTsne <- function(theObject, outputDir){
	
	tsneList <- getTSNEList(theObject)
	
	if(isTRUE(all.equal(length(getName(tsneList[[1]])), 0)))
		stop("The 'scRNAseq' object that you're using with 'exportResults' ",
				"method doesn't have its 'tSNEList' slot updated. Please use ",
				"'generateTSNECoordinates' on the object before.")
	
	invisible(lapply(tsneList, function(currentTSNE, outputDir){
						
						fileName <- paste0(getName(currentTSNE), ".tsv")
						write.table(getCoordinates(currentTSNE), 
								file=file.path(outputDir, fileName), 
								quote=FALSE, sep="\t", row.names=TRUE, 
								col.names=TRUE)
					},outputDir))
	
	message("Tsne coordinates saved.")
	
}


#' .exportNormInfo
#'
#' @description 
#' Export the normalized matrix, the columns, and the rows metadata to the 
#' sub-directory '1_MatrixInfo'. These results were obtained with 
#' ?normaliseCountMatrix.
#'
#' @param theObject scRNASeq object on which ?normaliseCountMatrix was run.
#' @param outputDir Output directory defined as dataDirectory/1_MatrixInfo. 
#' dataDirectory is directly retrieved from the scRNASeq object.
#' @param experimentName Name of the experiment retrieved in the corresponding 
#' slot of the scRNASeq object.
#' @param sceObject SingleCellExperiment object retrieved with ?getSceNorm. See
#'  ?SingleCellExperiment::SingleCellExperiment for more details. 
#' @param f Function defining the type of information retrieved from the 
#' sceObject. Can be the normalized count matrix with f=Biobase::exprs, the row
#'  metadata with f=SummarizedExperiment::rowData, or the col metadata with 
#' f=SummarizedExperiment::colData.
#' @param fileSuffix Suffix to add to the output files.
#' @param objectName String used in the message to indicate what object was 
#' saved.
#' 
#' @keywords internal
#' @noRd
.exportNormInfo <- function(theObject, outputDir, experimentName, sceObject, f,
		fileSuffix, objectName){
	
	if(all(dim(sceObject) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with ",
				"'exportResults' method doesn't have its ",
				"'sceNorm' slot updated. Please use 'normaliseCountMatrix'",
				" on the object before.")
	
	fileName <- paste0(experimentName, fileSuffix)
	write.table(f(sceObject), file=file.path(outputDir, fileName),
			quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
	message(objectName, " saved.")
}
	

.createFolder <- function(outputDir){
	
	if(!file.exists(outputDir))
		dir.create(outputDir,  showWarnings=FALSE)
}


#' exportResults
#'
#' @description 
#' Export all the results of Conclus to a Results sub-directory.
#'
#' @usage 
#' exportResults(theObject, saveClusteringResults=TRUE, 
#' saveNormalizedMatrix=FALSE, saveColData=FALSE, saveRowData=FALSE, 
#' saveTsne=FALSE, saveDBScan=FALSE, saveCellsSimilarityMatrix=FALSE,
#' saveClustersSimilarityMatrix=FALSE, saveFullMarkers=FALSE, 
#' saveTopMarkers=FALSE, saveGenesInfos=FALSE, saveAll=FALSE)
#' 
#' @param theObject An Object of class scRNASeq for which different steps of 
#' CONCLUS was applied to. The number of steps to run depends on what is 
#' wanted to be saved.
#' @param saveClusteringResults Default=TRUE. Save the final clustering results
#' giving the corresponding cluster number to each cell. The method 
#' ?calculateClustersSimilarity should have been run on the object. It is saved
#'  in the sub-directory 6_ConclusResult.
#' @param saveAll Default=FALSE. Save all results of CONCLUS (see details). The
#' last step run on the scRNASeq object should be ?retrieveGenesInfo.
#' @param saveNormalizedMatrix Default=FALSE. Save the normalized count matrix
#' as a csv file. It is obtained with ?normaliseCountMatrix. The matrix is saved
#'  to the sub-directory '1_MatrixInfo'.
#' @param saveColData Default=FALSE. Save the columns metadata of the normalized
#'  count matrix as a tsv file. These data are obtained with 
#' ?normaliseCountMatrix or were given as input of the method. These data are 
#' saved in the sub-directory '1_MatrixInfo'.
#' @param saveRowData Default=FALSE. Save the raw metadata of the normalized 
#' count matrix as a tsv file. These data are obtained with 
#' ?normaliseCountMatrix. They are saved in the sub-directory '1_MatrixInfo'.   
#' @param saveTsne Default=FALSE. Save the tsne coordinates for each combination
#'  of PCs and perplexities as tsv files. These coordinates were obtained with
#'  ?generateTSNECoordinates. They are saved in the sub-directory 
#' '2_TSNECoordinates'. 
#' @param saveDBScan Default=FALSE. Save the clustering results of dbscan as tsv
#'  files. The number of clustering solutions is 
#' PCs*perplexity*epsilon*minPoints (see ?runDBSCAN, 84 solutions by default). 
#' These are saved in the sub-directory '3_dbScan'.
#' @param saveCellsSimilarityMatrix Default=FALSE. Save the cells similarity 
#' matrix that was obtained with ?clusterCellsInternal. This matrix is saved in
#'  the sub-directory '4_CellSimilarityMatrix'.
#' @param saveClustersSimilarityMatrix Default=FALSE. Save the cluster 
#' similarity matrix that was obtained with ?calculateClustersSimilarity. It is
#'  saved in the sub-directory '5_ClusterSimilarityMatrix'.
#' @param saveFullMarkers Default=FALSE. Save the lists of markers that were 
#' obtained with ?rankGenes to the sub-directory 7_fullMarkers.
#' @param saveTopMarkers Default=FALSE. Save the top markers per clusters as csv
#'  files in the sub-directory '8_TopMarkers'. See ?retrieveTopClustersMarkers 
#' for more details. 
#' @param saveGenesInfos Default=FALSE. Save the genes information for each 
#' cluster as csv files in the sub-directory '9_genesInfos'. See 
#' ?retrieveGenesInfo for more details.
#' 
#' @aliases calculateClustersSimilarity
#'  
#' @author 
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
#' 
#' @rdname 
#' exportResults-scRNAseq
#' 
#' 
#' @examples
#' experimentName <- "Bergiers"
#' countMatrix <- as.matrix(read.delim(file.path(
#' "tests/testthat/test_data/test_countMatrix.tsv")))
#' outputDirectory <- "./"
#' columnsMetaData <- read.delim(
#' file.path("extdata/Bergiers_colData_filtered.tsv"))
#' 
#' ## Create the initial object
#' scr <- scRNAseq(experimentName = experimentName, 
#'                 countMatrix     = countMatrix, 
#'                 species         = "mouse",
#'                 outputDirectory = outputDirectory)
#' 
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#' 
#' ## Compute the tSNE coordinates
#' scrTsne <- generateTSNECoordinates(scrNorm, cores=5)
#' 
#' ## Perform the clustering with dbScan
#' scrDbscan <- runDBSCAN(scrTsne, cores=5)
#' 
#' ## Compute the cell similarity matrix
#' scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=10, cores=4)
#' 
#' ## Calculate clusters similarity
#' scrCSM <- calculateClustersSimilarity(scrCCI)
#' 
#' ## Ranking genes
#' scrS4MG <- rankGenes(scrCSM)
#' 
#' ## Getting marker genes
#' scrFinal <- retrieveTopClustersMarkers(scrS4MG, removeDuplicates = F)
#' 
#' ## Getting genes info
#' scrInfos <- retrieveGenesInfo(scrFinal, cores=5)
#' 
#' ## Saving all results
#' exportResults(scrInfos, saveAll=TRUE)
#' @exportMethod
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData
#' @importFrom Biobase exprs


setMethod(
		
		f = "exportResults",
		
		signature = "scRNAseq",
		
		definition = function(theObject, saveClusteringResults=TRUE, 
				saveAll=FALSE, saveNormalizedMatrix=FALSE, saveColData=FALSE, 
				saveRowData=FALSE, saveTsne=FALSE, saveDBScan=FALSE, 
				saveCellsSimilarityMatrix=FALSE, 
				saveClustersSimilarityMatrix=FALSE, saveFullMarkers=FALSE,
				saveTopMarkers=FALSE, saveGenesInfos=FALSE){
			
			
			if(saveAll){
				
				message("Saving all results.")
				saveNormalizedMatrix <- TRUE
				saveColData <- TRUE 
				saveRowData <- TRUE
				saveTsne <- TRUE
				saveDBScan <- TRUE 
				saveCellsSimilarityMatrix <- TRUE 
				saveClustersSimilarityMatrix <- TRUE
				saveFullMarkers <- TRUE
				saveTopMarkers <- TRUE
				saveGenesInfos <- TRUE
			}
			
			dataDirectory  <- getOutputDirectory(theObject)
			experimentName <- getExperimentName(theObject)
			sceObject <- getSceNorm(theObject)
			outputDir <- file.path(dataDirectory, "Results")
			.createFolder(outputDir)
			
			##########
			## Saving results after normalization
			#########
			
			if(saveNormalizedMatrix || saveRowData || saveColData){
				
				outputNorm <- file.path(outputDir, "1_MatrixInfo")
				.createFolder(outputNorm)
				
				## Export Normalized expression matrix
				if(saveNormalizedMatrix)
					.exportNormInfo(theObject, outputNorm, experimentName, 
							sceObject, Biobase::exprs, "_expression_matrix.tsv", 
							"Normalized expression matrix")
				
				## Export RowData
				if(saveRowData)
					.exportNormInfo(theObject, outputNorm, experimentName, 
							sceObject, SummarizedExperiment::rowData, 
							"_rowData.tsv", "RowData")
				
				## Export ColData
				if(saveColData)
					.exportNormInfo(theObject, outputNorm, experimentName, 
							sceObject, SummarizedExperiment::colData, 
							"_colData.tsv", "ColData")
			}
			
			
			##########
			## Saving results after tSNE calculation
			#########
			
			if(saveTsne){
				
				outputTSNE <- file.path(outputDir, "2_TSNECoordinates")
				.createFolder(outputTSNE)
				
				## Export all tSNE coordinates
				.exportTsne(theObject, outputTSNE)
			}
			
			
			##########
			## Saving results after running dbScan
			#########
			
			if(saveDBScan){
				
				outputDBSCAN <- file.path(outputDir, "3_dbScan")
				.createFolder(outputDBSCAN)
				
				## Export all clustering results given by dbscan
				.exportDBScan(theObject, outputDBSCAN)
			}
			
			##########
			## Saving results cells and clusters similarity matrices
			#########
			
			## Export CellsSimilarityMatrix
			if(saveCellsSimilarityMatrix){
				
				outputCellSM <- file.path(outputDir, "4_CellSimilarityMatrix")
				.createFolder(outputCellSM)
				.exportCSM(theObject, experimentName, outputCellSM)
			}
				
				
			## Export ClustersSimilarityMatrix
			if(saveClustersSimilarityMatrix){
				
				outClustSM <- file.path(outputDir, "5_ClusterSimilarityMatrix")
				.createFolder(outClustSM)
				.exportCSM(theObject, experimentName, outClustSM, cell=FALSE)
			}
				
					
			##########
			## Saving a table cluster-cell being the result of the consensus
			## clustering
			#########
			
			if (saveClusteringResults){
				
				outputClust <- file.path(outputDir, "6_ConclusResult")
				.createFolder(outputClust)
				
				## Export Clustering results
				.conclusResult(theObject, sceObject, outputClust)
			}
			
			
			##########
			## Saving the full markers lists and the top markers list
			##########
			
			if(saveFullMarkers){
				
				outputFull <- file.path(outputDir, "7_fullMarkers")
				.createFolder(outputFull)
				
				## Export the full marker lists
				.exportMarkers(theObject, outputFull)
			}
			
			if(saveTopMarkers){
				
				outputTop <- file.path(outputDir, "8_TopMarkers")
				.createFolder(outputTop)
				
				## Export top markers
				.exportMarkers(theObject, outputTop, "Top markers", top=TRUE)
			}
			
			###########
			## Saving genes infos
			###########
			
			if(saveGenesInfos){
				
				outputInfos <- file.path(outputDir, "9_genesInfos")
				.createFolder(outputInfos)
				
				## Export genes infos
				.exportGenesInfos(theObject, outputInfos)
			}
		})
