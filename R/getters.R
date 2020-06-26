################################################################################
########################  Getters of scRNAseq class ############################
################################################################################


#' @rdname getExperimentName-scRNAseq
#' @aliases getExperimentName,getExperimentName-method
setMethod(
    f = "getExperimentName",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@experimentName)
    })


#' @rdname getCountMatrix-scRNAseq
#' @aliases getCountMatrix,getCountMatrix-method
setMethod(
    f = "getCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@countMatrix)
    })


#' @rdname getSceNorm-scRNAseq
#' @aliases getSceNorm,getSceNorm-method
setMethod(
    f = "getSceNorm",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@sceNorm)
    })


#' @aliases getSpecies,getSpecies-method
#' @rdname getSpecies-scRNAseq
setMethod(
    f = "getSpecies",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@species)
    })



#' @rdname getOutputDirectory-scRNAseq
#' @aliases getOutputDirectory,getOutputDirectory-method
setMethod(
    f = "getOutputDirectory",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@outputDirectory)
    })



#' @rdname getTSNEList-scRNAseq
#' @aliases getTSNEList,getTSNEList-method
setMethod(
    f = "getTSNEList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@tSNEList)
    })


#' @rdname getDbscanList-scRNAseq
#' @aliases getDbscanList,getDbscanList-method
setMethod(
    f = "getDbscanList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@dbscanList)
    })

#' @rdname getCellsSimilarityMatrix-scRNAseq
#' @aliases getCellsSimilarityMatrix,getCellsSimilarityMatrix-method
setMethod(
    f = "getCellsSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@cellsSimilarityMatrix)
    })

#' @rdname getClustersSimilarityMatrix-scRNAseq
#' @aliases getClustersSimilarityMatrix,getClustersSimilarityMatrix-method
setMethod(
    f = "getClustersSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@clustersSimilarityMatrix)
    })

#' @rdname getClustersSimilarityMatrix-scRNAseq
#' @aliases getClustersSimilarityMatrix,getClustersSimilarityMatrix-method
setMethod(
    f = "getClustersSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@clustersSimiliratyOrdered)
    })

#' @rdname getMarkerGenesList-scRNAseq
#' @aliases getMarkerGenesList,getExperimentName-method
setMethod(
    f = "getMarkerGenesList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@markerGenesList)
    })


#' @rdname getClustersMarkers-scRNAseq
#' @aliases getClustersMarkers,getClustersMarkers-method
setMethod(
    f="getClustersMarkers",
    signature="scRNAseq",
    definition = function(theObject){
        return(theObject@clustersMarkers)
    })


#' @rdname getGenesInfos-scRNAseq
#' @aliases getGenesInfos,getGenesInfos-method
setMethod(
    f="getGenesInfos",
    signature="scRNAseq",
    definition = function(theObject){
        return(theObject@genesInfos)
    })


################################################################################
############################  Getter of Tsne Class  ############################
################################################################################

#' @rdname getTsneName-Tsne
#' @aliases getTsneName,getTsneName-method
setMethod(
    f = "getTsneName",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@name)
    })


#' @rdname getPerplexity-Tsne
#' @aliases getPerplexity,getPerplexity-method
setMethod(
    f = "getPerplexity",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@perplexity)
    })


#' @rdname getPC-Tsne
#' @aliases getPC,getPC-method
setMethod(
    f = "getPC",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@pc)
    })


#' @rdname getCoordinates-Tsne
#' @aliases getCoordinates,getCoordinates-method
setMethod(
    f = "getCoordinates",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@coordinates)
    })



################################################################################
###########################  Getter of Dbscan Class  ###########################
################################################################################

#' @rdname getDbscanName-Dbscan
#' @aliases getDbscanName,getDbscanName-method
setMethod(
    f = "getDbscanName",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@name)
    })



#' @rdname getEpsilon-Dbscan
#' @aliases getEpsilon,getEpsilon-method
setMethod(
    f = "getEpsilon",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@epsilon)
    })


#' @rdname getMinPoints-Dbscan
#' @aliases getMinPoints,getExperimentName-method
setMethod(
    f = "getMinPoints",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@minPoints)
    })


#' @rdname getClustering-Dbscan
#' @aliases getClustering,getClustering-method
setMethod(
    f = "getClustering",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@clustering)
    })
