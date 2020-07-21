################################################################################
########################  Setters for scRNAseq class  ##########################
################################################################################

#' @usage
#' setExperimentName(theObject) <- value
#' 
#' @description
#' Update the experiment name slot with value.
#' 
#' @param theObject A scRNA-seq object to update.
#' @param value The value to update the slot with.
#' 
#' @rdname setters
#' @aliases setExperimentName<-
#' @title setters
#' @name setters
#' 
#' @exportMethod setExperimentName<-
setReplaceMethod(
    f = "setExperimentName",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@experimentName <- value
        validObject(theObject)
        return(theObject)
    })

#' @usage
#' setCountMatrix(theObject) <- value
#' 
#' @description
#' Update the countMatrix slot with value.
#' 
#' @rdname setters
#' @aliases setCountMatrix<-
#' @exportMethod setCountMatrix<-
setReplaceMethod(
    f = "setCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@countMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setSceNorm<-
#' @describeIn scRNAseq-class Set the SingleCellExperiment object used.
setReplaceMethod(
    f = "setSceNorm",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@sceNorm <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setSpecies<-
#' @describeIn scRNAseq-class Set the species. 
setReplaceMethod(
    f = "setSpecies",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@species <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setOutputDirectory<-
#' @describeIn scRNAseq-class Set the path of the output directory.
setReplaceMethod(
    f = "setOutputDirectory",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@outputDirectory <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setTSNEList<-
#' @describeIn scRNAseq-class Set the list of Tsne objects.
setReplaceMethod(
    f = "setTSNEList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@tSNEList <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setDbscanList<-
#' @describeIn scRNAseq-class Set the list of Dbscan objects. 
setReplaceMethod(
    f = "setDbscanList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@dbscanList <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setCellsSimilarityMatrix<-
#' @param theObject An object of class scRNAseq.
#' @param value A numeric matrix. 
#' @describeIn scRNAseq-class Set the cell similarity matrix.
setReplaceMethod(
    f = "setCellsSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@cellsSimilarityMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setClustersSimilarityMatrix<-
#' @describeIn scRNAseq-class Set the cluster similarity matrix.
setReplaceMethod(
    f = "setClustersSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@clustersSimilarityMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setClustersSimiliratyOrdered<-
#' @describeIn scRNAseq-class Set the clusters ordered by similarity. 
setReplaceMethod(
    f = "setClustersSimiliratyOrdered",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@clustersSimiliratyOrdered <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setMarkerGenesList<-
#' @describeIn scRNAseq-class Set the list of marker genes by clusters. 
setReplaceMethod(
    f = "setMarkerGenesList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@markerGenesList <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setClustersMarkers<-
#' @describeIn scRNAseq-class Set the most significant markers by clusters. 
setReplaceMethod(
		f="setClustersMarkers",
		signature="scRNAseq",
		definition = function(theObject, value){
			theObject@clustersMarkers<- value
			validObject(theObject)
			return(theObject)
		})


#' @exportMethod setGenesInfos<-
#' @describeIn scRNAseq-class Set a data.frame containing informations about 
#' the marker genes.
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
#' @describeIn Tsne-class Set the name of the tSNE. 
setReplaceMethod(
    f = "setName",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@name <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setPC<-
#' @describeIn Tsne-class Set the PC parameter.
setReplaceMethod(
    f = "setPC",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@pc <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setPerplexity<-
#' @describeIn Tsne-class Set the perplexity parameter.
setReplaceMethod(
    f = "setPerplexity",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@perplexity <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setCoordinates<-
#' @describeIn Tsne-class Set the matrix of tSNE coordinates.
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
#' @describeIn Dbscan-class Set the name of the Dbscan. 
setReplaceMethod(
    f = "setName",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@name <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setEpsilon<-
#' @describeIn Dbscan-class Set the epsilon used.
setReplaceMethod(
    f = "setEpsilon",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@epsilon <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setMinPoints<-
#' @describeIn Dbscan-class Set the minPoints used.
setReplaceMethod(
    f = "setMinPoints",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@minPoints <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setClustering<-
#' @describeIn Dbscan-class Set the matrix of Dbscan clustering. 
setReplaceMethod(
    f = "setClustering",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@clustering <- value
        validObject(theObject)
        return(theObject)
    })

