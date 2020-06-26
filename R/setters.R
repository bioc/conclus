################################################################################
########################  Setters for scRNAseq class  ##########################
################################################################################

#' @rdname setExperimentName-scRNAseq
#' @aliases setExperimentName<-,setExperimentName-method
setReplaceMethod(
    f = "setExperimentName",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@experimentName <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setCountMatrix-scRNAseq
#' @aliases setCountMatrix<-,setCountMatrix-method
setReplaceMethod(
    f = "setCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@countMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setSceNorm-scRNAseq
#' @aliases setSceNorm<-,setSceNorm-method
setReplaceMethod(
    f = "setSceNorm",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@sceNorm <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setSpecies-scRNAseq
#' @aliases setSpecies<-,setSpecies-method
setReplaceMethod(
    f = "setSpecies",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@species <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setOutputDirectory-scRNAseq
#' @aliases setOutputDirectory<-,setOutputDirectory-method
setReplaceMethod(
    f = "setOutputDirectory",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@outputDirectory <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setTSNEList-scRNAseq
#' @aliases setTSNEList<-,setTSNEList-method
setReplaceMethod(
    f = "setTSNEList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@tSNEList <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setDbscanList-scRNAseq
#' @aliases setDbscanList<-,setDbscanList-method
setReplaceMethod(
    f = "setDbscanList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@dbscanList <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setCellsSimilarityMatrix-scRNAseq
#' @aliases setCellsSimilarityMatrix<-,setCellsSimilarityMatrix-method
setReplaceMethod(
    f = "setCellsSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@cellsSimilarityMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setClustersSimilarityMatrix-scRNAseq
#' @aliases setClustersSimilarityMatrix<-,setClustersSimilarityMatrix-method
setReplaceMethod(
    f = "setClustersSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@clustersSimilarityMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setCountMatrix-scRNAseq
#' @aliases setClustersSimiliratyOrdered<-,setClustersSimiliratyOrdered-method
setReplaceMethod(
    f = "setClustersSimiliratyOrdered",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@clustersSimiliratyOrdered <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setMarkerGenesList-scRNAseq
#' @aliases setMarkerGenesList<-,setMarkerGenesList-method
setReplaceMethod(
    f = "setMarkerGenesList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@markerGenesList <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setClustersMarkers-scRNAseq
#' @aliases setClustersMarkers<-,setClustersMarkers-method
setReplaceMethod(
		f="setClustersMarkers",
		signature="scRNAseq",
		definition = function(theObject, value){
			theObject@clustersMarkers<- value
			validObject(theObject)
			return(theObject)
		})


#' @rdname setGenesInfos-scRNAseq
#' @aliases setGenesInfos<-,setGenesInfos-method
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


#' @rdname setTsneName-Tsne
#' @aliases setTsneName-method,setName-method
setReplaceMethod(
    f = "setTsneName",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@name <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setPC-Tsne
#' @aliases setPC<-,setPC-method
setReplaceMethod(
    f = "setPC",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@pc <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setPerplexity-Tsne
#' @aliases setPerplexity<-,setPerplexity-method
setReplaceMethod(
    f = "setPerplexity",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@perplexity <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setCoordinates-Tsne
#' @aliases setCoordinates<-,setCoordinates-method
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


#' @rdname setDbscanName-Dbscan
#' @aliases setDbscanName-method,setName-method
setReplaceMethod(
    f = "setDbscanName",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@name <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setEpsilon-Dbscan
#' @aliases setEpsilon<-,setEpsilon-method
setReplaceMethod(
    f = "setEpsilon",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@epsilon <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setMinPoints-Dbscan
#' @aliases setMinPoints<-,setMinPoints-method
setReplaceMethod(
    f = "setMinPoints",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@minPoints <- value
        validObject(theObject)
        return(theObject)
    })


#' @rdname setClustering-Dbscan
#' @aliases setClustering<-,setClustering-method
setReplaceMethod(
    f = "setClustering",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@clustering <- value
        validObject(theObject)
        return(theObject)
    })


