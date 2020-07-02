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

setMethod(
		
		f = "exportResults",
		
		signature = "scRNAseq",
		
		definition = function(theObject, saveClusteringResults=TRUE, 
				saveNormalizedMatrix=FALSE, saveColData=FALSE, 
				saveRowData=FALSE, saveTsne=FALSE, saveDBScan=FALSE, 
				saveCellsSimilarityMatrix=FALSE, 
				saveClustersSimilarityMatrix=FALSE, saveFullMarkers=FALSE,
				saveTopMarkers=FALSE, saveGenesInfos=FALSE, saveAll=FALSE){
			
			if(saveAll){
				
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
				saveGenesInfo(theObject, outputInfos)
				message("Genes infos saved.")
			}
		})
