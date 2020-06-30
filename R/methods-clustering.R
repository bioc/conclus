##################
## testClustering
##################

# similar function as plotDistanceGraph,
# but with already known epsilon value
.plotDistanceGraphWithEpsilon <- function(tSNEData, minNeighbours=5,
		epsilon=1.2){
	
	dbscan::kNNdistplot(tSNEData, k=minNeighbours)
	abline(h=epsilon, lty=2)
}


# plots test DBSCAN on one of the pictures
# for being ensured that clustering will
# probably work successful
.plotTestClustering <- function(tSNEData, minNeighbours=5,epsilon=1.2){
	
	dbscanList <- fpc::dbscan(tSNEData, eps=epsilon, MinPts=minNeighbours)
	p <- factoextra::fviz_cluster(dbscanList, tSNEData, ellipse=TRUE, 
			geom="point", legend="bottom")
	return(p)
}



.checkParamsTests <- function(sceObject, dbscanEpsilon, minPts, perplexities, 
		PCs){
	
	
	if(all(dim(sceObject) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with ",
				"'testClustering' method doesn't have its ",
				"'sceNorm' slot updated. Please use 'normaliseCountMatrix'",
				" on the object before.")
	
	## Check parameters
	if(!is.numeric(dbscanEpsilon))
		stop("'dbscanEpsilon' parameter should be an integer.")
	
	## Check minPts argument
	if(!is.numeric(minPts))
		stop("'minPts' parameter should be an integer")
	
	## Check perplexities argument
	if(!is.numeric(perplexities))
		stop("'perplexities' parameter should be a vector of numeric.")
	
	## Check PCs argument
	if(!is.numeric(PCs))
		stop("'PCs' parameter should be a vector of numeric.")
	
	## Check randomSeed argument
	if(!is.numeric(randomSeed))
		stop("'randomSeed' parameter should be an integer.")
}


.printTSNE <- function(writeOutput, dataDirectory, width, height, onefile, tSNE, 
		...){
	
	if(writeOutput){
		message("Saving results tSNE.")
		pdf(file.path(dataDirectory, "test_clustering", "test_tSNE.pdf"),
				width=width, height=height, onefile=onefile, ...)
		print(tSNE)
		dev.off()
	}else
		print(tSNE)
	
}


.printDist <- function(writeOutput, dataDirectory, width, height, onefile, tSNE,
		dbscanEpsilon, minPts, ...){
	
	if(writeOutput){
		
		message("Saving results distance graph.")
		fileDist <- "distance_graph.pdf"
		pdf(file.path(dataDirectory, "test_clustering", fileDist), 
				width=width, height=height, onefile=onefile, ...)
		.plotDistanceGraphWithEpsilon(tSNE$data, epsilon=dbscanEpsilon,
				minNeighbours=minPts)
		dev.off()
	}else{
		
		dev.new()
		.plotDistanceGraphWithEpsilon(tSNE$data, epsilon=dbscanEpsilon,
				minNeighbours=minPts)
	}
	
}


.printDBScan <- function(writeOutput, dataDirectory, width, height, onefile, 
		tSNE, epsilon, minPts, ...){
	
	p <- .plotTestClustering(tSNE$data, epsilon=dbscanEpsilon,
			minNeighbours=minPts)
	
	if(writeOutput){
		
		message("Saving dbscan results.")
		fileClust <- "test_clustering.pdf"
		pdf(file.path(dataDirectory, "test_clustering", fileClust),
				width=width, height=height, onefile=onefile, ...)
		print(p)
		dev.off()
	}else{
		dev.new()
		print(p)
	}
}	



setMethod(
		
		f="testClustering",
		
		signature="scRNAseq",
		
		definition=function(theObject, dbscanEpsilon=1.4, minPts=5, 
				perplexities=c(30), PCs=c(4), randomSeed=42, width=7, height=7,
				onefile=FALSE, cores=1, writeOutput=FALSE, ...){
			
			validObject(theObject)
			
			## Get the values from slots
			sceObject <- getSceNorm(theObject)
			dataDirectory <- getOutputDirectory(theObject)
			experimentName <- getExperimentName(theObject)
			
			initialisePath(dataDirectory)
			dir.create(file.path(dataDirectory, "test_clustering"), 
					showWarnings=FALSE)
			
			## Checking parameters
			.checkParamsTests(sceObject, dbscanEpsilon, minPts, perplexities, 
					PCs)
			
			message("Generating TSNE.")
			
			## 1. Generating 2D tSNE plots
			tSNE <- .getTSNEresults(theObject, Biobase::exprs(sceObject),
					cores=cores, PCs=PCs, perplexities=perplexities, 
					randomSeed=randomSeed)
			.printTSNE(writeOutput, dataDirectory, width, height, onefile, tSNE, 
					...)
			
			## 2. Clustering with dbscan
			.printDist(writeOutput, dataDirectory, width, height, onefile, tSNE,
					dbscanEpsilon, minPts, ...)
			.printDBScan(writeOutput, dataDirectory, width, height, onefile, 
					tSNE, epsilon, minPts, ...)
		})



##################
## clusterCellsInternal
##################


## This function calculates how many time a pair of cells were assigned to
## a cluster by dbscan.

.mkSimMat <- function(dbscanList, cores=14){
	
	## This foreach gives, for each dbscan result, a matrix containing 1 and
	## 0. 1 means that a pair of cells was allocated to one cluster, 
	## 0 means that the pair of cells was allocated to no cluster.
	
	myCluster <- parallel::makeCluster(cores, type="PSOCK")
	doParallel::registerDoParallel(myCluster)
	simMatsList <- foreach::foreach(i=seq_len(length(dbscanList)), .export=
							c("getClustering")) %dopar% {
				
				## Get the clustering
				clustering <- getClustering(dbscanList[[i]])
				
				## Get the number of cells in the current clustering 
				nrow <- ncol(clustering)
				ncol <- ncol(clustering)
				
				## Create a square matrix 
				simMat <- matrix(0, ncol=ncol, nrow=nrow)
				colnames(simMat) <- colnames(clustering)
				rownames(simMat) <- colnames(clustering)
				
				## Get cluster assignments
				clusters <- unique(as.vector(clustering))
				
				## For each cluster number (clusters), we compute a similarity
				## matrix with 1. In the resulting list (l), there is one matrix
				## for each cluster number. 
				l <- lapply(clusters[clusters!=0], 
						function(cluster, simMat, clustering){
							
							## Get the cells in the same cluster
							selCol <- colnames(clustering)[
									clustering == cluster]
							## Add +1 for each cell seen in the same cluster
							simMat[rownames(simMat) %in% selCol,
									colnames(simMat) %in% selCol] <- 1
							return(simMat)
						}, simMat, clustering) 
				
				## Sum all the matrix to have a single one showing cells 
				## that were allocated to one cluster (1) and cells that were
				## not allocated to one cluster (0).
				simMat <- Reduce("+", l)}
	parallel::stopCluster(myCluster)
	
	## Sum the similiratity matrices to obtain a single one giving howw many 
	## times a pair of cells was clustered together accross all the dbscan 
	## parameters 
	simMat <- Reduce('+', simMatsList)
	simMat <- simMat / length(dbscanList)
	stopifnot(isSymmetric(simMat))
	
	return(simMat)
}

.checkClusteringMethod <- function(clusteringMethod){
	
	clusteringMethods <- c("ward.D", "ward.D2", "single", "complete", "average",
			"mcquitty", "median", "centroid")
	availableMethods <- paste(clusteringMethods, collapse="; ")
	if(!clusteringMethod %in% clusteringMethods)
		stop("'clusteringMethod' should be one of: ", availableMethods)
	
}


.checkParamsInternal <- function(sceObject, dbscanList, clusterNumber, 
		deepSplit, cores, clusteringMethod){
	
	## Check if the normalized matrix is correct
	if(all(dim(sceObject) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with ",
				"'clusterCellsInternal' function doesn't have its 'sceNorm' ",
				"slot updated. Please use 'normaliseCountMatrix' on the object",
				" before.")
	
	## Check if the dbscan list is correct
	if (length(dbscanList) <= 1)
		stop("The 'scRNAseq' object that you're using with ",
				"'clusterCellsInternal' function doesn't have its 'dbscanList'",
				" slot updated. Please use 'runDBSCAN' on the object before.")
	
	## Check the cluster number
	if(!is.numeric(clusterNumber))
		stop("'clusterNumber' parameter should be a numeric.")
	
	## Check the deepSplit
	if(!is.numeric(deepSplit))
		stop("'deepSplit' parameter should be a numeric.")
	
	## Check the cores
	if(!is.numeric(cores))
		stop("'cores' parameter should be a numeric.")
	
	## Check the clustering method
	.checkClusteringMethod(clusteringMethod) 
	
}


setMethod(
		
		f = "clusterCellsInternal",
		
		signature = "scRNAseq",
		
		definition = function(theObject, clusterNumber=0, deepSplit=4, cores=1,
				clusteringMethod="ward.D2"){
			
			## Check if the Object is valid
			validObject(theObject)
			
			sceObject  <- getSceNorm(theObject)
			dbscanList <- getDbscanList(theObject)
			
			.checkParamsInternal(sceObject, dbscanList, clusterNumber, deepSplit, 
					cores, clusteringMethod)
			
			message("Calculating cells similarity matrix.")
			cellsSimilarityMatrix <- .mkSimMat(dbscanList, cores=cores)
			
			
			distanceMatrix <- as.dist(sqrt((1-cellsSimilarityMatrix)/2))
			clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
			
			if(clusterNumber == 0){
				message(paste0("Assigning cells to clusters. DeepSplit = ",
								deepSplit))
				clusters <- unname(cutreeDynamic(clusteringTree,
								distM=as.matrix(distanceMatrix),
								verbose=0,
								deepSplit=deepSplit))
			} else {
				message(paste0("Assigning cells to ", clusterNumber, " clusters."))
				clusters <- cutree(clusteringTree, k=clusterNumber)
			}
			
			SummarizedExperiment::colData(sceObject)$clusters <- factor(clusters)
			
			setCellsSimilarityMatrix(theObject) <- cellsSimilarityMatrix
			setSceNorm(theObject) <- sceObject
			return(theObject)
		})




##################
## calculateClustersSimilarity
##################


### This function returns a matrix with "protocells" representing clusters.
### Values show how much two "protocells" are similar.
### 1 if clusters are perfectly similar, 0 if totally different.

.computeClusMat <- function(clustersNames, mat, clusters){
	
	resultMedList <- lapply(clustersNames, function(currentClustName, 
					fullmat, clusts){
				
				return(matrixStats::rowMedians(fullmat[, clusts == 
												currentClustName]))
			}, mat, clusters)	
	medMat <- do.call(cbind, resultMedList)		
	return(medMat)
}


.mkSimMed <- function(simMat, clusters, clustersNames){
	
	clusMed <- matrix(ncol=length(unique(clusters)), nrow=nrow(simMat))
	clusMed <- .computeClusMat(clustersNames, simMat, clusters)
	clusMed <- t(clusMed)
	
	simMed <- matrix(ncol=length(unique(clusters)),
			nrow=length(unique(clusters)))
	simMed <- .computeClusMat(clustersNames, clusMed, clusters)
	colnames(simMed) <- clustersNames
	rownames(simMed) <- clustersNames
	
	return(simMed)
}

.checkParamsClusterSim <- function(clusteringMethod, cellsSimilarityMatrix, 
		sceObject){
	
	## Check the clustering method
	.checkClusteringMethod(clusteringMethod)
	
	## Check if the normalized matrix is correct
	if(all(dim(sceObject) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with ",
				"'calculateClustersSimilarity' function doesn't have its ",
				"'sceNorm' slot updated. Please use 'normaliseCountMatrix' on",
				" the object before.")
	
	## Check if the normalized matrix is filtered
	if(!("clusters" %in% names(colData(sceObject))))
		stop("The 'scRNAseq' object that you're using with ",
				"'calculateClustersSimilarity' function doesn't have a correct",
				" 'sceNorm' slot. This slot should be a 'SingleCellExperiment'",
				" object containing 'clusters' column in its colData. Please ",
				"check if you correctly used 'clusterCellsInternal' on the ",
				"object.")
	
	## Check the cell similarity matrix
	if(all(dim(cellsSimilarityMatrix) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with ",
				"'calculateClustersSimilarity' function doesn't have its ",
				"'cellsSimilarityMatrix' slot updated by clusterCellsInternal.",
				" Please use 'clusterCellsInternal' on the object before.")
}


setMethod(
		
		f = "calculateClustersSimilarity",
		
		signature = "scRNAseq",
		
		definition = function(theObject, clusteringMethod = "ward.D2"){
			
			## Check if the Object is valid
			validObject(theObject)
			
			cellsSimilarityMatrix <- getCellsSimilarityMatrix(theObject)
			sceObject  <- getSceNorm(theObject)
			
			.checkParamsClusterSim(clusteringMethod, cellsSimilarityMatrix, 
					sceObject)
			
			clusters <- SummarizedExperiment::colData(sceObject)$clusters
			clustersNumber <- length(unique(clusters))
			clustersNames <- levels(clusters)
			
			## Calculating cluster similarity.
			clustersSimilarityMatrix <- .mkSimMed(as.matrix(cellsSimilarityMatrix), 
					clusters, clustersNames)
			
			## Plotting matrix
			distanceMatrix <- as.dist(sqrt((1 - clustersSimilarityMatrix)/2))
			clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
			
			clustersSimOrdered <- data.frame(clusterNames=clustersNames,
					clusterIndexes=1:clustersNumber)
			rownames(clustersSimOrdered) <- clustersSimOrdered$clusterIndexes
			clustersSimOrdered <- clustersSimOrdered[clusteringTree$order, ]
			
			setClustersSimilarityMatrix(theObject)  <- clustersSimilarityMatrix
			setClustersSimiliratyOrdered(theObject) <-
					as.factor(clustersSimOrdered$clusterNames)
			return(theObject)
		})


##################
## addClusteringManually
##################

setMethod(
		
		f = "addClusteringManually",
		
		signature = "scRNAseq",
		
		definition = function(theObject, fileName, columnName = "clusters"){
			
			## Check if the Object is valid
			validObject(theObject)
			
			## Check if the normalized matrix is correct
			sceObject  <- getSceNorm(theObject)
			
			if (all(dim(sceObject) == c(0,0)))
				stop("The 'scRNAseq' object that you're using with ",
						"'addClusteringManually' function doesn't have its ",
						"'sceNorm' slot updated. Please use ",
						"'normaliseCountMatrix' on the object before.")
			
			!! not good use the slot !!
			tableData <- read.table(file.path(dataDirectory, "output_tables",
							paste0(experimentName, "_", fileName)), sep="\t")
			
			## Already check in validObject
			dataDirectory  <- getOutputDirectory(theObject)
			experimentName <- getExperimentName(theObject)
			
			colDf <- SummarizedExperiment::colData(sceObject)
			
			if(all(rownames(colDf) %in% rownames(tableData))){
				
				if(isTRUE(all.equal(ncol(tableData), 1)))
					colDf$clusters <- factor(tableData[rownames(colDf), ])
				else
					colDf$clusters <- factor(
							tableData[rownames(colDf), ][ columnName])
				
				setSceNorm(theObject) <- sceObject
				return(theObject)
				
			}else{
				
				msg <- paste("Rownames in colDf are not equal to rownames in",
						"table. Returning SCE object with cells intersecting ",
						"with clusters_table.\n",
						sep="")
				message(msg)
				sceObject <- sceObject[, colnames(sceObject) %in%
								intersect(colnames(sceObject),
										rownames(tableData))]
				tableData$randomColumn <- NA
				tableData <- tableData[rownames(tableData) %in%
								intersect(colnames(sceObject),
										rownames(tableData)), ]
				colDf$clusters <- factor(
						tableData[rownames(colDf),][,columnName])
				setSceNorm(theObject) <- sceObject
				
				return(theObject)
			}
		})
