################################################################################
########################  Getters of scRNAseq class ############################
################################################################################

#' @usage
#' getExperimentName(theObject)
#' 
#' @description
#' Retrieve data of a slot of a scRNA-seq object.
#' 
#' @param theObject A scRNA-seq object. See description or ?scRNAseq.
#' 
#' @description
#' getExperimentName: Get the name of the experiment.
#' 
#' @rdname getters-scRNAseq
#' @aliases getExperimentName
#' @title getters-scRNAseq
#' @name getters-scRNAseq
#' 
#' @exportMethod getExperimentName
setMethod(
    f = "getExperimentName",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@experimentName)
    })


#' @usage
#' getCountMatrix(theObject)
#' 
#' @description
#' getCountMatrix: Get the count matrix.
#' 
#' @rdname getters-scRNAseq
#' @aliases getCountMatrix
#' 
#' @exportMethod getCountMatrix
setMethod(
    f = "getCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@countMatrix)
    })

#' @usage
#' getSceNorm(theObject)
#' 
#' @description
#' getSceNorm: Get the SingleCellExperiment object used.
#' 
#' @rdname getters-scRNAseq
#' @aliases getSceNorm
#' 
#' @exportMethod getSceNorm
setMethod(
    f = "getSceNorm",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@sceNorm)
    })


#' @usage
#' getSpecies(theObject)
#' 
#' @description
#' getSpecies: Get the species.
#' 
#' @rdname getters-scRNAseq
#' @aliases getSpecies
#' 
#' @exportMethod getSpecies
setMethod(
    f = "getSpecies",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@species)
    })


#' @usage
#' getOutputDirectory(theObject)
#' 
#' @description
#' getOutputDirectory: Get the path of the output directory.
#' 
#' @rdname getters-scRNAseq
#' @aliases getOutputDirectory
#' 
#' @exportMethod getOutputDirectory
setMethod(
    f = "getOutputDirectory",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@outputDirectory)
    })


#' @usage
#' getTSNEList(theObject)
#' 
#' @description
#' getTSNEList: Get the list of Tsne objects.
#' 
#' @rdname getters-scRNAseq
#' @aliases getTSNEList
#' 
#' @exportMethod getTSNEList
setMethod(
    f = "getTSNEList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@tSNEList)
    })


#' @usage
#' getDbscanList(theObject)
#' 
#' @description
#' getDbscanList: Get the list of Dbscan objects.
#' 
#' @rdname getters-scRNAseq
#' @aliases getDbscanList
#' 
#' @exportMethod getDbscanList
setMethod(
    f = "getDbscanList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@dbscanList)
    })

#' @usage
#' getCellsSimilarityMatrix(theObject)
#' 
#' @description
#' getCellsSimilarityMatrix: Get the cell similarity matrix.
#' 
#' @rdname getters-scRNAseq
#' @aliases getCellsSimilarityMatrix
#' 
#' @exportMethod getCellsSimilarityMatrix
setMethod(
    f = "getCellsSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@cellsSimilarityMatrix)
    })


#' @usage
#' getClustersSimilarityMatrix(theObject)
#' 
#' @description
#' getClustersSimilarityMatrix: Get the cluster similarity matrix.
#' 
#' @rdname getters-scRNAseq
#' @aliases getClustersSimilarityMatrix
#' 
#' @exportMethod getClustersSimilarityMatrix
setMethod(
    f = "getClustersSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@clustersSimilarityMatrix)
    })


#' @usage
#' getClustersSimiliratyOrdered(theObject)
#' 
#' @description
#' getClustersSimiliratyOrdered: Get the clusters ordered by similarity.
#' 
#' @rdname getters-scRNAseq
#' @aliases getClustersSimiliratyOrdered
#' 
#' @exportMethod getClustersSimiliratyOrdered
setMethod(
		f = "getClustersSimiliratyOrdered",
		signature = "scRNAseq",
		definition = function(theObject){
			return(theObject@clustersSimiliratyOrdered)
		})


#' @usage
#' getMarkerGenesList(theObject)
#' 
#' @description
#' getMarkerGenesList: Get the list of marker genes by clusters.
#' 
#' @rdname getters-scRNAseq
#' @aliases getMarkerGenesList
#' 
#' @exportMethod getMarkerGenesList
setMethod(
    f = "getMarkerGenesList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@markerGenesList)
    })


#' @usage
#' getClustersMarkers(theObject)
#' 
#' @description
#' getClustersMarkers: Get the most significant markers by clusters.
#' 
#' @rdname getters-scRNAseq
#' @aliases getClustersMarkers
#' 
#' @exportMethod getClustersMarkers
setMethod(
    f="getClustersMarkers",
    signature="scRNAseq",
    definition = function(theObject){
        return(theObject@clustersMarkers)
    })


#' @usage
#' getGenesInfos(theObject)
#' 
#' @description
#' getGenesInfos: Get informations about marker genes.
#' 
#' @rdname getters-scRNAseq
#' @aliases getGenesInfos
#' 
#' @exportMethod getGenesInfos
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
