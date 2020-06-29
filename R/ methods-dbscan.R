## This function calculates dbscan for all t-SNE from tSNEList with all
## combinations of parameters from epsilon and minPoints.
## It does not set random seed. It allows to vary this parameter automatically.
## It returns a matrix where columns are iterations.
## The number of iterations is equal to 
## (tSNEList)*length(epsilon)*length(minPoints)


.mkDbscan <- function(tSNEList, cores=14, epsilon=c(1.2, 1.5, 1.8), 
		minPoints=c(15, 20)){
	
	myCluster <- parallel::makeCluster(cores, type="PSOCK")
	doParallel::registerDoParallel(myCluster)
	dbscanResults <- foreach::foreach(i=rep(rep(1:length(tSNEList),
									each=length(minPoints)),
							length(epsilon)),
					eps=rep(epsilon,
							each=length(tSNEList)*
									length(minPoints)),
					MinPts=rep(minPoints,
							length(tSNEList)*
									length(epsilon)),
					.combine='cbind',
					.export = 'getCoordinates'
			) %dopar% {
				fpc::dbscan(
						getCoordinates(tSNEList[[i]]),
						eps=eps,
						MinPts=MinPts)$cluster
			}
	parallel::stopCluster(myCluster)
	return(dbscanResults)
}


.checkParamDbScan <- function(sceObject, tSNEList, cores, epsilon, minPoints, 
		writeOutput){
	
	## Check if the normalized count matrix is correct
	if(all(dim(sceObject) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with 'runDBSCAN' ",
				"function doesn't have its 'sceNorm' slot updated. Please ",
				"use 'normaliseCountMatrix' on the object before.")
	
	## Check if the tsne list is correct
	if(length(tSNEList) <= 1)
		stop("The 'scRNAseq' object that you're using with 'runDBSCAN' ",
				"function doesn't have its 'tSNEList' slot updated. Please",
				" use 'generateTSNECoordinates' on the object before.")
	
	## Check cores argument
	if(!is.numeric(cores))
		stop("'cores' parameter should be an integer")
	
	## Check epsilon argument
	if(!is.numeric(epsilon))
		stop("'epsilon' parameter should be a vector of numeric.")
	
	## Check minPoints argument
	if(!is.numeric(minPoints))
		stop("'minPoints' parameter should be a vector of numeric.")
	
	## Check writeOutput argument
	if(!is.logical(writeOutput))
		stop("'writeOutput' parameter should be a boolean.")
	
}



.writeDBScanResults <- function(theObject, dbscanResults){
	
	## Get the values of slots
	dataDirectory  <- getOutputDirectory(theObject)
	experimentName <- getExperimentName(theObject)
	outputDataDirectory <- "output_tables"
	filename <- paste0(experimentName, "_dbscan_results.tsv")
	outputfile <- file.path(dataDirectory, outputDataDirectory, filename)
	write.table(dbscanResults, outputfile, quote=FALSE, sep="\t")
}



#' runDBSCAN
#'
#' @description 
#' Run clustering iterations with selected parameters using DBSCAN.
#'
#' @usage 
#' runDBSCAN(theObject, cores=1, epsilon=c(1.3, 1.4, 1.5), minPoints=c(3, 4), 
#' writeOutput=FALSE)
#' 
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized and the tSNE coordinates were calculated. 
#' See ?normaliseCountMatrix and ?generateTSNECoordinates.
#' @param cores Maximum number of jobs that CONCLUS can run in parallel. 
#' Default is 1.
#' @param epsilon Reachability distance parameter of fpc::dbscan() function. 
#' See Ester et al. (1996) for more details. Default = c(1.3, 1.4, 1.5)
#' @param minPoints Reachability minimum no. of points parameter of 
#' fpc::dbscan() function. See Ester et al. (1996) for more details. 
#' Default = c(3, 4) 
#' @param writeOutput If TRUE, write the results of the dbScan clustering to 
#' the output directory defined in theObject, in the sub-directory 
#' output_tables. Default = FALSE.
#' 
#' @aliases runDBSCAN
#' 
#' @details
#' Following the calculation of t-SNE coordinates, DBSCAN is run with a range 
#' of epsilon and MinPoints values which will yield a total of 84 clustering 
#' solutions (PCs x perplexities x MinPoints x epsilon). *minPoints* is the 
#' minimum cluster size which you assume to be meaningful for your experiment 
#' and *epsilon* is the radius around the cell where the algorithm will try to 
#' find *minPoints* dots. Optimal *epsilon* must lay one the knee of the k-NN 
#' function as shown in the "test_clustering/distance_graph.pdf" (**!! Might need to be changed**).
#' 
#' @author 
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
#' 
#' @rdname 
#' runDBSCAN-scRNAseq
#' 
#' @return 
#' An object of class scRNASeq with its dbscanList slot updated. Also writes 
#' the clustering results in "dataDirectory/output_tables" subfolder if the 
#' parameter writeOutput is TRUE.
#' 
#' @examples
#' experimentName <- "Bergiers"
#' countMatrix <- matrix(sample(seq_len(40), size=4000, replace = TRUE), 
#' nrow=20, ncol=200)
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
#' @seealso
#' normaliseCountMatrix generateTSNECoordinates
#' 
#' @exportMethod
setMethod(
		
		f = "runDBSCAN",
		
		signature = "scRNAseq",
		
		definition = function(theObject, cores=1, epsilon=c(1.3, 1.4, 1.5), 
				minPoints=c(3, 4), writeOutput=FALSE){
			
			
			## Check the Object
			validObject(theObject)
			
			## Check function parameters
			sceObject <- getSceNorm(theObject)
			tSNEList <- getTSNEList(theObject)
			.checkParamDbScan(sceObject, tSNEList, cores, epsilon, minPoints, 
					writeOutput)
			
			## Taking only cells from the sceObject
			tSNEList <- lapply(tSNEList, function(currentTsne, sceObj){
						
						tmp <- getCoordinates(currentTsne)
						setCoordinates(currentTsne) <- tmp[colnames(sceObj), ]
						return(currentTsne)
					}, sceObject)
			
			## Calculating dbscan combinaisons
			dbscanResults <- .mkDbscan(tSNEList=tSNEList, cores=cores, 
					epsilon=epsilon, minPoints=minPoints)
			dbscanResults <- t(dbscanResults)
			colnames(dbscanResults) <-
					SummarizedExperiment::colData(sceObject)$cellName
			
			## Output file
			if(writeOutput)
				.writeDBScanResults(theObject, dbscanResults)
			
			## Creation of a list of Dbscan objects
			totalLength <- length(tSNEList) * length(minPoints)
			epsilonCombinaison <- rep(epsilon, each=totalLength)
			minPtsCombinaison  <- rep(minPoints, 
					length(tSNEList)*length(epsilon))
			rowDbscanList <- split(dbscanResults, seq_len(nrow(dbscanResults)))
			rowNamesVec <- paste0("clust.", seq_len(nrow(dbscanResults)))
			dbscanObjNameVec <- paste0("Clustering_", 
					seq_len(nrow(dbscanResults)))
			
			dbscanList <- mapply(function(rowName, dbscanObjName, epsComb, 
							minPts, rowDbscan){
						
						clustering <- t(rowDbscan)
						colnames(clustering) <- colnames(dbscanResults)
						rownames(clustering) <- rowName
						dbscanObj<- Dbscan(name= dbscanObjName, epsilon=epsComb,
								minPoints=minPts, clustering=clustering)
						return(dbscanObj)
						
					}, rowNamesVec, dbscanObjNameVec, epsilonCombinaison,
					minPtsCombinaison, rowDbscanList, SIMPLIFY = FALSE, 
					USE.NAMES = FALSE)
			
			setDbscanList(theObject) <- dbscanList
			
			return(theObject)
		})
