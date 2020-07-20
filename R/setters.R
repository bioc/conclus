################################################################################
########################  Setters for scRNAseq class  ##########################
################################################################################

#' @exportMethod setExperimentName
setReplaceMethod(
    f = "setExperimentName",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@experimentName <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setCountMatrix
setReplaceMethod(
    f = "setCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@countMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setSceNorm
setReplaceMethod(
    f = "setSceNorm",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@sceNorm <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setSpecies
setReplaceMethod(
    f = "setSpecies",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@species <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setOutputDirectory
setReplaceMethod(
    f = "setOutputDirectory",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@outputDirectory <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setTSNEList
setReplaceMethod(
    f = "setTSNEList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@tSNEList <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setDbscanList
setReplaceMethod(
    f = "setDbscanList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@dbscanList <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setCellsSimilarityMatrix
setReplaceMethod(
    f = "setCellsSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@cellsSimilarityMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setClustersSimilarityMatrix
setReplaceMethod(
    f = "setClustersSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@clustersSimilarityMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setClustersSimiliratyOrdered
setReplaceMethod(
    f = "setClustersSimiliratyOrdered",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@clustersSimiliratyOrdered <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setMarkerGenesList
setReplaceMethod(
    f = "setMarkerGenesList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@markerGenesList <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setClustersMarkers
setReplaceMethod(
		f="setClustersMarkers",
		signature="scRNAseq",
		definition = function(theObject, value){
			theObject@clustersMarkers<- value
			validObject(theObject)
			return(theObject)
		})


#' @exportMethod setGenesInfos
setReplaceMethod(
		f="setGenesInfos",
		signature="scRNAseq",
		definition = function(theObject, value){
			theObject@genesInfos<- value
			validObject(theObject)
			return(theObject)
		})



################################################################################
############################  Setters for Tsne class  ##########################
################################################################################

#' @exportMethod setName
setReplaceMethod(
    f = "setName",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@name <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setPC
setReplaceMethod(
    f = "setPC",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@pc <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setPerplexity
setReplaceMethod(
    f = "setPerplexity",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@perplexity <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setCoordinates
setReplaceMethod(
    f = "setCoordinates",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@coordinates <- value
        validObject(theObject)
        return(theObject)
    })


################################################################################
###########################  Setters for Dbscan class  #########################
################################################################################

#' @exportMethod setName
setReplaceMethod(
    f = "setName",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@name <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setEpsilon
setReplaceMethod(
    f = "setEpsilon",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@epsilon <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setMinPoints
setReplaceMethod(
    f = "setMinPoints",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@minPoints <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setClustering
setReplaceMethod(
    f = "setClustering",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@clustering <- value
        validObject(theObject)
        return(theObject)
    })


