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
		
	if(all(dim(matrix) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with 'exportResults' ",
				"function doesn't have its '", matType,"' slot ",
				"updated. Please use '", funName,"' on the object ",
				"before.")
		
	fileName <- paste0(experimentName, fileName)
	write.csv(matrix, file=file.path(outputDir, fileName))
	message(matType, " saved.")
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
		
		definition = function(theObject, saveNormalizedMatrix=TRUE,
				saveColData=TRUE, saveRowData=TRUE, saveTsne=TRUE,
				saveCellsSimilarityMatrix=TRUE, 
				saveClustersSimilarityMatrix=TRUE, saveWorkspace=TRUE,
				saveClusteringResults=TRUE){
			
			dataDirectory  <- getOutputDirectory(theObject)
			experimentName <- getExperimentName(theObject)
			sceObject <- getSceNorm(theObject)
			outputDir <- file.path(dataDirectory, "output_tables")
			.createFolder(outputDir)
			
			##########
			## Saving results after normalization
			#########
			
			outputNorm <- file.path(outputDir, "MatrixInfo")
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
			
			##########
			## Saving results after tSNE calculation
			#########
			
			outputTSNE <- file.path(outputDir, "TSNECoordinates")
			.createFolder(outputTSNE)
			
			## Export all tSNE coordinates
			if(saveTsne)
				.exportTsne(theObject, outputTSNE)
			
			
			##########
			## Saving results after running dbScan
			#########
			
			outputDBSCAN <- file.path(outputDir, "dbScan")
			.createFolder(outputDBSCAN)
			
			## Export all clustering results given by dbscan
			
			
			
			
			## Export CellsSimilarityMatrix
			if(saveCellsSimilarityMatrix)
				.exportCSM(theObject, experimentName, outputDir)
				
			## Export ClustersSimilarityMatrix
			if(saveClustersSimilarityMatrix)
				.exportCSM(theObject, experimentName, outputDir, cell=FALSE)
					
			## Save workspace
			if(saveWorkspace){
				fileName <- paste0(experimentName, "_full_workspace.RData")
				save.image(file=file.path(outputDir, fileName))
				message("Workspace saved.")
			}
			
			## Export Clustering results
			if (saveClusteringResults){
				tableData <- S4Vectors::DataFrame(
						clusters=SummarizedExperiment::colData(sceObject)$clusters,
						row.names=SummarizedExperiment::colData(sceObject)$cellName)
				
				write.table(tableData,
						file=file.path(dataDirectory, "output_tables",
								paste0(experimentName,"_", 
										"clusters_table.tsv")),
						sep="\t", quote=FALSE)
				message("Clusters table saved.")
			}
			
		})
