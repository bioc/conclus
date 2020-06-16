.checkObject <- function(obj, val, getElementFunction, className){
	
	element <- getElementFunction(obj)
	
	if(isTRUE(all.equal(className, "Dbscan")))
		nbElement <- ncol(element)
	else
		nbElement <- nrow(element)
	
	if(isTRUE(all.equal(nbElement, val)) && is(obj) == className)
		return(TRUE)
	else
		return(FALSE)
}

checkList <- function(theobjectList, getFunction, className){
	
	if(!isTRUE(all.equal(className, "Dbscan")) &&
			!isTRUE(all.equal(className, "Tsne")))
		stop("In checkList, className should be Dbscan or Tsne.")
	
	if(isTRUE(all.equal(className, "Dbscan")))
		val <- ncol(getFunction(theobjectList[[1]]))
	else
		val <- nrow(getFunction(theobjectList[[1]]))
	
	vec <- sapply(theobjectList, .checkObject, val, getFunction, 
			className)
	
	if(!all(vec)){
		theobjectList[!vec]
		names <- sapply(theobjectList , getName)
		names <- paste(names , collapse="; ")
		stop("The elements in ", className, "List slot don't have the same ",
				"number of cells or the same class")
	}
}


checkMarkerGenesList <- function(markerGeneobjectlist, 
		clustersSimiliratyOrdered=NULL){
	
	expecteColumn <- c("Gene", "mean_log10_fdr", "n_05", "score")
	
	vec <- sapply(markerGeneobjectlist, function(df, expcol){
				
				nameColumn <- colnames(df)
				
				if(all(expcol %in% nameColumn))
					return(TRUE)
				else
					return(FALSE)
			}, expecteColumn)	
	
	if(!all(vec))
		stop("'markerGenesList' slot should contain a list of ",
				"dataframes with at least following columns : 'Gene', ",
				"'mean_log10_fdr', ", "'n_05', 'score'.\n")
	
	
	if(!is.null(clustersSimiliratyOrdered) &&
			!is.na(markerGeneobjectlist[[1]]$score))
		if(!isTRUE(all.equal(length(markerGeneobjectlist),
						length(clustersSimiliratyOrdered))))
			stop("'markerGenesList' should contain as many dataframes as ",
					"clusters found. Number of dataframes :",
					length(markerGeneobjectlist), " and the number of cluters ",
					"found is :", length(clustersSimiliratyOrdered), ".")
}


checkClusterMarkers <- function(clusterMarkers, clustersSimiliratyOrdered=NULL){
	
	
	expectedColumn <- c("geneName", "clusters")
	
	if(!all(expectedColumn %in% colnames(clusterMarkers)))
		stop("The clusterMarkers data frame should have the columns 'geneName'",
				" and 'clusters'")
	
	
	nbClustMark <- length(unique(clusterMarkers$clusters))
	nbClust <- length(clustersSimiliratyOrdered)
	
	if(nrow(clusterMarkers) > 1 && !isTRUE(all.equal(nbClust, nbClustMark)))
		stop("clusterMarkers should have the same number of clusters than the",
				" number of clusters found. Nb clusters for markers: ", 
				nbClustMark, ". Nb of clusters: ", nbClust)
}