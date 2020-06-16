.checkParamTsne <- function(theObject, randomSeed, cores, PCs, perplexities, 
		writeOutput){
	
	## Test the normalized count matrix
	sceObject <- getSceNorm(theObject)
	
	if(all(dim(sceObject) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with ",
				"'generateTSNECoordinates' function doesn't have its ",
				"'sceNorm' slot updated. Please use 'normaliseCountMatrix'",
				" on the object before.")
	
	## Check parameters
	if(!is.numeric(randomSeed))
		stop("'randomSeed' parameter should be an integer.")
	
	## Check cores argument
	if(!is.numeric(cores))
		stop("'cores' parameter should be an integer")
	
	## Check PCs argument
	if(!is.numeric(PCs))
		stop("'PCs' parameter should be a vector of numeric.")
	
	## Check perplexities argument
	if(!is.numeric(perplexities))
		stop("'perplexities' parameter should be a vector of numeric.")
	
	## Check writeOutput argument
	if(!is.logical(writeOutput))
		stop("'writeOutput' parameter should be a boolean.")
		
	return(sceObject)
}



.buildingTsneObjects <- function(PCs, perplexities, experimentName, TSNEres){
	
	return(sapply(seq_len(length(PCs)*length(perplexities)), function(i, PCA, 
							perp){
						
						name <- paste0(experimentName, '_tsne_coordinates_', i, 
								"_" , PCA[i], "PCs_", perp[i], "perp")
						
						tSNEObj <- Tsne(name = name, 
								coordinates = as.matrix(TSNEres[1, i][[1]]), 
								perplexity = perp[i],  
								pc  = PCA[i])
						
						return(tSNEObj)
						
					}, rep(PCs, length(perplexities)), 
					rep(perplexities, each=length(PCs))))
	
}



.writeOutputTsne <- function(theObject, PCs, perplexities, experimentName, 
		TSNEres){
	
	dataDirectory <- getOutputDirectory(theObject)
	tSNEDirectory <- "tsnes"
	
	invisible(sapply(seq_len(length(PCs)*length(perplexities)), function(i, PCA, 
							perp){
						
						name <- paste0(experimentName, '_tsne_coordinates_', i, 
								"_" , PCA[i], "PCs_", perp[i], "perp")
						
						write.table(TSNEres[1, i][[1]], 
								file = file.path(dataDirectory, tSNEDirectory,
										paste0(name, ".tsv")), quote=FALSE, 
								sep='\t')
						
						
					}, rep(PCs, length(perplexities)), 
					rep(perplexities, each=length(PCs))))
	
	initialisePath(dataDirectory)
	outputDir <- file.path(dataDirectory, tSNEDirectory)
	filesList <- list.files(outputDir, pattern="_tsne_coordinates_")
	deletedFiles <- sapply(filesList, function(fileName) 
				file.remove(file.path(outputDir, fileName)))
	saveRDS(TSNEres, file=file.path(dataDirectory, "output_tables",
					paste0(experimentName, "_tSNEResults.rds")))
	
}



setMethod(
		
    f = "generateTSNECoordinates",
	
    signature = "scRNAseq",
	
    definition = function(theObject, randomSeed=42, cores=1, 
			PCs=c(4, 6, 8, 10, 20, 40, 50), perplexities=c(30,40), 
			writeOutput = FALSE){
		
        ## Check the Object
        validObject(theObject)
        
		## Check method parameters
		sceObject <- .checkParamTsne(theObject, randomSeed, cores, PCs, 
				perplexities, writeOutput)
        
		message("Running TSNEs using ", cores, " cores.")
		TSNEres <- .getTSNEresults(theObject, Biobase::exprs(sceObject),
				cores=cores, PCs=PCs, perplexities=perplexities, 
				randomSeed=randomSeed)
		
		message("Building TSNEs objects.")
		experimentName <- getExperimentName(theObject)
		setTSNEList(theObject) <- .buildingTsneObjects(PCs, perplexities, 
				experimentName, TSNEres)
		
		  
        if(writeOutput)
			.writeOutputTsne(theObject, PCs, perplexities, experimentName, 
					TSNEres)
		
		return(theObject)
    })
