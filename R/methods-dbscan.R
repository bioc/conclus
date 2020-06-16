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
	minPtsCombinaison  <- rep(minPoints, length(tSNEList) * length(epsilon))
	rowDbscanList <- split(dbscanResults, seq_len(nrow(dbscanResults)))
	rowNamesVec <- paste0("clust.", seq_len(nrow(dbscanResults)))
	dbscanObjNameVec <- paste0("Clustering_", seq_len(nrow(dbscanResults)))
	
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
