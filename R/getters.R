################################################################################
########################  Getters of scRNAseq class ############################
################################################################################

#' @exportMethod getExperimentName
#' @noRd
## Described in scRNAseq-class
setMethod(
    f = "getExperimentName",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@experimentName)
    })


#' @exportMethod getCountMatrix
#' @noRd
## Described in scRNAseq-class
setMethod(
    f = "getCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@countMatrix)
    })


#' @exportMethod getSceNorm
#' @noRd
## Described in scRNAseq-class
setMethod(
    f = "getSceNorm",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@sceNorm)
    })


#' @exportMethod getSpecies
#' @noRd
## Described in scRNAseq-class
setMethod(
    f = "getSpecies",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@species)
    })


#' @exportMethod getOutputDirectory
#' @noRd
## Described in scRNAseq-class
setMethod(
    f = "getOutputDirectory",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@outputDirectory)
    })


#' @exportMethod getTSNEList
#' @noRd
## Described in scRNAseq-class
setMethod(
    f = "getTSNEList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@tSNEList)
    })


#' @exportMethod getDbscanList
#' @noRd
## Described in scRNAseq-class
setMethod(
    f = "getDbscanList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@dbscanList)
    })

#' @exportMethod getCellsSimilarityMatrix
#' @noRd
## Described in scRNAseq-class
setMethod(
    f = "getCellsSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@cellsSimilarityMatrix)
    })


#' @exportMethod getClustersSimilarityMatrix
#' @noRd
## Described in scRNAseq-class
setMethod(
    f = "getClustersSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@clustersSimilarityMatrix)
    })


#' @exportMethod getClustersSimiliratyOrdered
#' @noRd
## Described in scRNAseq-class
setMethod(
		f = "getClustersSimiliratyOrdered",
		signature = "scRNAseq",
		definition = function(theObject){
			return(theObject@clustersSimiliratyOrdered)
		})


#' @exportMethod getMarkerGenesList
#' @noRd
## Described in scRNAseq-class
setMethod(
    f = "getMarkerGenesList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@markerGenesList)
    })


#' @exportMethod getClustersMarkers
#' @noRd
## Described in scRNAseq-class
setMethod(
    f="getClustersMarkers",
    signature="scRNAseq",
    definition = function(theObject){
        return(theObject@clustersMarkers)
    })


#' @exportMethod getGenesInfos
#' @noRd
## Described in scRNAseq-class
setMethod(
    f="getGenesInfos",
    signature="scRNAseq",
    definition = function(theObject){
        return(theObject@genesInfos)
    })


################################################################################
############################  Getter of Tsne Class  ############################
################################################################################

#' @exportMethod getName
#' @noRd
## Described in Tsne-class
setMethod(
    f = "getName",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@name)
    })


#' @exportMethod getPerplexity
#' @noRd
## Described in Tsne-class
setMethod(
    f = "getPerplexity",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@perplexity)
    })


#' @exportMethod getPC
#' @noRd
## Described in Tsne-class
setMethod(
    f = "getPC",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@pc)
    })


#' @exportMethod getCoordinates
#' @noRd
## Described in Tsne-class
setMethod(
    f = "getCoordinates",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@coordinates)
    })



################################################################################
###########################  Getter of Dbscan Class  ###########################
################################################################################

#' @exportMethod getName
#' @noRd
## Described in Dbscan-class
setMethod(
    f = "getName",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@name)
    })


#' @exportMethod getEpsilon
#' @noRd
## Described in Dbscan-class
setMethod(
    f = "getEpsilon",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@epsilon)
    })


#' @exportMethod getMinPoints
#' @noRd
## Described in Dbscan-class
setMethod(
    f = "getMinPoints",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@minPoints)
    })


#' @exportMethod getClustering
#' @noRd
## Described in Dbscan-class
setMethod(
    f = "getClustering",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@clustering)
    })
