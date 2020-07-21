################################################################################
########################  Getters of scRNAseq class ############################
################################################################################

#' @exportMethod getExperimentName
#' @describeIn scRNAseq-class Get the name of the experiment.
setMethod(
    f = "getExperimentName",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@experimentName)
    })


#' @exportMethod getCountMatrix
#' @describeIn scRNAseq-class Get the count matrix.
setMethod(
    f = "getCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@countMatrix)
    })


#' @exportMethod getSceNorm
#' @describeIn scRNAseq-class Get the SingleCellExperiment object used.
setMethod(
    f = "getSceNorm",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@sceNorm)
    })


#' @exportMethod getSpecies
#' @describeIn scRNAseq-class Get the species.
setMethod(
    f = "getSpecies",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@species)
    })


#' @exportMethod getOutputDirectory
#' @describeIn scRNAseq-class Get the path of the output directory.
setMethod(
    f = "getOutputDirectory",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@outputDirectory)
    })


#' @exportMethod getTSNEList
#' @describeIn scRNAseq-class Get the list of Tsne objects.
setMethod(
    f = "getTSNEList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@tSNEList)
    })


#' @exportMethod getDbscanList
#' @describeIn scRNAseq-class Get the list of Dbscan objects.
setMethod(
    f = "getDbscanList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@dbscanList)
    })

#' @exportMethod getCellsSimilarityMatrix
#' @describeIn scRNAseq-class Get the cell similarity matrix.
setMethod(
    f = "getCellsSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@cellsSimilarityMatrix)
    })


#' @exportMethod getClustersSimilarityMatrix
#' @describeIn scRNAseq-class Get the cluster similarity matrix.
setMethod(
    f = "getClustersSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@clustersSimilarityMatrix)
    })


#' @exportMethod getClustersSimiliratyOrdered
#' @describeIn scRNAseq-class Get the clusters ordered by similarity.
setMethod(
		f = "getClustersSimiliratyOrdered",
		signature = "scRNAseq",
		definition = function(theObject){
			return(theObject@clustersSimiliratyOrdered)
		})


#' @exportMethod getMarkerGenesList
#' @describeIn scRNAseq-class Get the list of marker genes by clusters
setMethod(
    f = "getMarkerGenesList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@markerGenesList)
    })


#' @exportMethod getClustersMarkers
#' @describeIn scRNAseq-class Get the most significant markers by clusters.
setMethod(
    f="getClustersMarkers",
    signature="scRNAseq",
    definition = function(theObject){
        return(theObject@clustersMarkers)
    })


#' @exportMethod getGenesInfos
#' @describeIn scRNAseq-class Get informations about marker genes.
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
#' @describeIn Tsne-class Get the name of the tSNE. 
setMethod(
    f = "getName",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@name)
    })


#' @exportMethod getPerplexity
#' @describeIn Tsne-class Get the perplexity used.
setMethod(
    f = "getPerplexity",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@perplexity)
    })


#' @exportMethod getPC
#' @describeIn Tsne-class Get the PC used.
setMethod(
    f = "getPC",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@pc)
    })


#' @exportMethod getCoordinates
#' @describeIn Tsne-class Get the matrix of tSNE coordinates. 
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
#' @describeIn Dbscan-class Get the name of the Dbscan.
setMethod(
    f = "getName",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@name)
    })


#' @exportMethod getEpsilon
#' @describeIn Dbscan-class Get the epsilon used.
setMethod(
    f = "getEpsilon",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@epsilon)
    })


#' @exportMethod getMinPoints
#' @describeIn Dbscan-class Get the MinPoint used.
setMethod(
    f = "getMinPoints",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@minPoints)
    })


#' @exportMethod getClustering
#' @describeIn Dbscan-class Get the matrix of DBSCAN clustering.
setMethod(
    f = "getClustering",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@clustering)
    })
