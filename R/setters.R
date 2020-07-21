################################################################################
########################  Setters for scRNAseq class  ##########################
################################################################################

#' @exportMethod setExperimentName<-
#' @noRd
## Described in scRNAseq-class
setReplaceMethod(
    f = "setExperimentName",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@experimentName <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setCountMatrix<-
#' @noRd
## Described in scRNAseq-class
setReplaceMethod(
    f = "setCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@countMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setSceNorm<-
#' @noRd
## Described in scRNAseq-class
setReplaceMethod(
    f = "setSceNorm",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@sceNorm <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setSpecies<-
#' @noRd
## Described in scRNAseq-class
setReplaceMethod(
    f = "setSpecies",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@species <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setOutputDirectory<-
#' @noRd
## Described in scRNAseq-class
setReplaceMethod(
    f = "setOutputDirectory",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@outputDirectory <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setTSNEList<-
#' @noRd
## Described in scRNAseq-class
setReplaceMethod(
    f = "setTSNEList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@tSNEList <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setDbscanList<-
#' @noRd
## Described in scRNAseq-class
setReplaceMethod(
    f = "setDbscanList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@dbscanList <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setCellsSimilarityMatrix<-
#' @noRd
## Described in scRNAseq-class
setReplaceMethod(
    f = "setCellsSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@cellsSimilarityMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setClustersSimilarityMatrix<-
#' @noRd
## Described in scRNAseq-class
setReplaceMethod(
    f = "setClustersSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@clustersSimilarityMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setClustersSimiliratyOrdered<-
#' @noRd
## Described in scRNAseq-class
setReplaceMethod(
    f = "setClustersSimiliratyOrdered",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@clustersSimiliratyOrdered <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setMarkerGenesList<-
#' @noRd
## Described in scRNAseq-class
setReplaceMethod(
    f = "setMarkerGenesList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@markerGenesList <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setClustersMarkers<-
#' @noRd
## Described in scRNAseq-class
setReplaceMethod(
		f="setClustersMarkers",
		signature="scRNAseq",
		definition = function(theObject, value){
			theObject@clustersMarkers<- value
			validObject(theObject)
			return(theObject)
		})


#' @exportMethod setGenesInfos<-
#' @noRd
## Described in scRNAseq-class
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

#' @exportMethod setName<-
#' @noRd
## Described in Tsne-class
setReplaceMethod(
    f = "setName",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@name <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setPC<-
#' @noRd
## Described in Tsne-class
setReplaceMethod(
    f = "setPC",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@pc <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setPerplexity<-
#' @noRd
## Described in Tsne-class
setReplaceMethod(
    f = "setPerplexity",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@perplexity <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setCoordinates<-
#' @noRd
## Described in Tsne-class
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

#' @exportMethod setName<-
#' @noRd
## Described in Dbscan-class
setReplaceMethod(
    f = "setName",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@name <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setEpsilon<-
#' @noRd
## Described in Dbscan-class
setReplaceMethod(
    f = "setEpsilon",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@epsilon <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setMinPoints<-
#' @noRd
## Described in Dbscan-class
setReplaceMethod(
    f = "setMinPoints",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@minPoints <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setClustering<-
#' @noRd
## Described in Dbscan-class
setReplaceMethod(
    f = "setClustering",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@clustering <- value
        validObject(theObject)
        return(theObject)
    })


