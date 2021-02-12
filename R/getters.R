################################################################################
########################  Getters of scRNAseq class ############################
################################################################################


#' @description
#' Retrieve the data of the slots of a scRNA-seq, Tsne or Dbscan object.
#'
#' @param theObject A scRNA-seq, Tsne or Dbscan object. See description or
#' ?scRNAseq, ?Tsne, ?Dbscan.
#'
#' @rdname getters
#'
#' @examples
#' ## Load the count matrix
#' countmatrixPath <- system.file("extdata/countMatrix.tsv", package="conclus")
#' countMatrix <- loadDataOrMatrix(file=countmatrixPath, type="countMatrix")
#'
#' ## Load the coldata
#' coldataPath <- system.file("extdata/colData.tsv", package="conclus")
#' columnsMetaData <- loadDataOrMatrix(file=coldataPath, type="coldata",
#'                                     columnID="cell_ID")
#'
#' ## Create the initial object
#' scr <- singlecellRNAseq(experimentName = "Bergiers",
#'                 countMatrix     = countMatrix,
#'                 species         = "mouse",
#'                 outputDirectory = "YourOutputDirectory")
#'
#' experimentName <- getExperimentName(scr)
#' countMatrix <- getCountMatrix(scr)
#' species <- getSpecies(scr)
#' outputDirectory <- getOutputDirectory(scr)
#'
#' @name getters
#' @title getters
#' @author Ilyess RACHEDI
NULL


#' @usage
#' getExperimentName(theObject)
#'
#' @return
#' getExperimentName: Get the name of the experiment (scRNA-seq).
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
#' @return
#' getCountMatrix: Get the count matrix (scRNA-seq).
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
#' @return
#' getSceNorm: Get the SingleCellExperiment object used (scRNA-seq).
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
#' @return
#' getSpecies: Get the species (scRNA-seq).
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
#' @return
#' getOutputDirectory: Get the path of the output directory (scRNA-seq).
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
#' @return
#' getTSNEList: Get the list of Tsne objects (scRNA-seq).
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
#' @return
#' getDbscanList: Get the list of Dbscan objects (scRNA-seq).
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
#' @return
#' getCellsSimilarityMatrix: Get the cell similarity matrix (scRNA-seq).
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
#' @return
#' getClustersSimilarityMatrix: Get the cluster similarity matrix (scRNA-seq).
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
#' getClustersSimilarityOrdered(theObject)
#'
#' @return
#' getClustersSimilarityOrdered: Get the clusters ordered by similarity
#' (scRNA-seq).
#'
#' @rdname getters
#' @aliases getClustersSimilarityOrdered
#'
#' @exportMethod getClustersSimilarityOrdered
setMethod(
        f = "getClustersSimilarityOrdered",
        signature = "scRNAseq",
        definition = function(theObject){
            return(theObject@clustersSimiliratyOrdered)
        })


#' @usage
#' getMarkerGenesList(theObject)
#'
#' @return
#' getMarkerGenesList: Get the list of marker genes by clusters  (scRNA-seq).
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
#' @return
#' getClustersMarkers: Get the most significant markers by clusters (scRNA-seq).
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
#' @return
#' getGenesInfos: Get informations about marker genes  (scRNA-seq).
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
setMethod(
        f = "getName",
        signature = c("Tsne"),
        definition = function(theObject){
            return(theObject@name)
        })

#' @usage
#' getPerplexity(theObject)
#'
#' @return
#' getPerplexity: Get the perplexity used (Tsne).
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
#' @return
#' getPC: Get the PC used (Tsne).
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
#' @return
#' getCoordinates: Get the matrix of tSNE coordinates  (Tsne).
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
#' @return
#' getName: Get the name of the tSNE or Dbscan object (Dbscan).
#'
#' @rdname getters
#' @aliases getName
#'
#' @exportMethod getName
setMethod(
        f = "getName",
        signature = c("Dbscan"),
        definition = function(theObject){
            return(theObject@name)
        })


#' @usage
#' getEpsilon(theObject)
#'
#' @return
#' getEpsilon: Get the epsilon used (Dbscan).
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
#' @return
#' getMinPoints: Get the MinPoint used (Dbscan).
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
#' @return
#' getClustering: Get the matrix of DBSCAN clustering (Dbscan).
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
