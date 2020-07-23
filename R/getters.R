################################################################################
########################  Getters of scRNAseq class ############################
################################################################################


#' @description
#' Retrieve data of a slot of a scRNA-seq, Tsne or Dbscan object.
#' 
#' @param theObject A scRNA-seq, Tsne or Dbscan object. See description or 
#' ?scRNAseq, ?Tsne, ?Dbscan.
#'  
#' @rdname getters
#' @name getters
#' @title getters
NULL


#' @usage
#' getExperimentName(theObject)
#' 
#' @description
#' getExperimentName: Get the name of the experiment.
#' 
#' @param theObject A scRNA-seq object. See description or ?scRNAseq.
#'  
#' @rdname getters
#' @aliases getExperimentName
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
#' @rdname getters
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
#' @rdname getters
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
#' @rdname getters
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
#' @rdname getters
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
#' @rdname getters
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
#' @rdname getters
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
#' @rdname getters
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
#' @rdname getters
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
#' @rdname getters
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
#' @rdname getters
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
#' @rdname getters
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
#' @rdname getters
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

#' @rdname getters
#' @exportMethod getName,Tsne-method
setMethod(
		f = "getName",
		signature = c("Tsne"),
		definition = function(theObject){
			return(theObject@name)
		})

#' @usage
#' getPerplexity(theObject)
#' 
#' @description
#' getPerplexity: Get the perplexity used.
#' 
#' @rdname getters
#' @aliases getPerplexity
#' 
#' @exportMethod getPerplexity
setMethod(
    f = "getPerplexity",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@perplexity)
    })


#' @usage
#' getPC(theObject)
#' 
#' @description
#' getPC: Get the PC used.
#' 
#' @rdname getters
#' @aliases getPC
#' 
#' @exportMethod getPC
setMethod(
    f = "getPC",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@pc)
    })


#' @usage
#' getCoordinates(theObject)
#' 
#' @description
#' getCoordinates: Get the matrix of tSNE coordinates.
#' 
#' @rdname getters
#' @aliases getCoordinates
#' 
#' @exportMethod getCoordinates 
setMethod(
    f = "getCoordinates",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@coordinates)
    })



################################################################################
###########################  Getter of Dbscan Class  ###########################
################################################################################


#' @usage
#' getName(theObject)
#'
#' @description
#' getName: Get the name of the tSNE or Dbscan object.
#'  
#' @rdname getters
#' @aliases getName
#' 
#' @exportMethod getName,Dbscan-method 
setMethod(
		f = "getName",
		signature = c("Dbscan"),
		definition = function(theObject){
			return(theObject@name)
		})


#' @usage
#' getEpsilon(theObject)
#' 
#' @description
#' getEpsilon: Get the epsilon used.
#' 
#' @rdname getters
#' @aliases getEpsilon
#' 
#' @exportMethod getEpsilon
setMethod(
    f = "getEpsilon",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@epsilon)
    })

#' @usage
#' getMinPoints(theObject)
#' 
#' @description
#' getMinPoints: Get the MinPoint used.
#' 
#' @rdname getters
#' @aliases getMinPoints
#' 
#' @exportMethod getMinPoints
setMethod(
    f = "getMinPoints",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@minPoints)
    })


#' @usage
#' getClustering(theObject)
#' 
#' @description
#' getClustering: Get the matrix of DBSCAN clustering.
#' 
#' @rdname getters
#' @aliases getClustering
#' 
#' @exportMethod getClustering
setMethod(
    f = "getClustering",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@clustering)
    })



