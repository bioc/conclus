################################################################################
########################  Getters of scRNAseq class ############################
################################################################################


#' @usage
#'  
#' 
#' @description
#' Retrieve data of a slot of a scRNA-seq object.
#' 
#' @param theObject A scRNA-seq object. See description or ?scRNAseq.
#'  
#' @rdname gettersscRNAseq
#' @aliases gettersscRNAseq
#' @name gettersscRNAseq
#' @title gettersscRNAseq
#'  
#' @exportMethod gettersscRNAseq
setMethod(
		f = "getterscRNAseq",
		signature = "scRNAseq",
		definition = function(theObject){
			return(theObject)
		})


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
#' @rdname gettersscRNAseq
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
#' @rdname gettersscRNAseq
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
#' @rdname gettersscRNAseq
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
#' @rdname gettersscRNAseq
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
#' @rdname gettersscRNAseq
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
#' @rdname gettersscRNAseq
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
#' @rdname gettersscRNAseq
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
#' @rdname gettersscRNAseq
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
#' @rdname gettersscRNAseq
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
#' @rdname gettersscRNAseq
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
#' @rdname gettersscRNAseq
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
#' @rdname gettersscRNAseq
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
#' @rdname gettersscRNAseq
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

#' @usage
#' getName(theObject)
#' 
#' @description
#' Retrieve the data of a slot of a Tsne-class object.
#' 
#' @param theObject A Tsne object. See description or ?Tsne.
#' 
#' @description
#' getName: Get the name of the tSNE.
#' 
#' @rdname getters-Tsne
#' @aliases getName
#' @title getters-Tsne
#' @name getters-Tsne
#' 
#' @exportMethod getName 
setMethod(
    f = "getName",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@name)
    })


#' @usage
#' getPerplexity(theObject)
#' 
#' @description
#' getPerplexity: Get the perplexity used.
#' 
#' @rdname getters-Tsne
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
#' @rdname getters-Tsne
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
#' @rdname getters-Tsne
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
#' Retrieve data of a slot of a Dbscan-class object.
#' 
#' @param theObject A Dbscan object. See description or ?Dbscan.
#' 
#' @description
#' getName: Get the name of the Dbscan.
#' 
#' @rdname getters-Dbscan
#' @aliases getNameDbscan
#' @title getters-Dbscan
#' @name getters-Dbscan
#' 
#' @exportMethod getName
setMethod(
    f = "getName",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@name)
    })


#' @usage
#' getEpsilon(theObject)
#' 
#' @description
#' getEpsilon: Get the epsilon used.
#' 
#' @rdname getters-Dbscan
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
#' @rdname getters-Dbscan
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
#' @rdname getters-Dbscan
#' @aliases getClustering
#' 
#' @exportMethod getClustering
setMethod(
    f = "getClustering",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@clustering)
    })
