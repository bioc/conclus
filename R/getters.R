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
#' 
#' @examples
#' ## Load the count matrix
#' countmatrixPath <- system.file("extdata/test_countMatrix.tsv", 
#'                             package="conclus")
#' countMatrix <- loadCountMatrix(file=countmatrixPath, header=TRUE, dec=".",
#'                                 sep='\t')
#' 
#' ## Load the coldata
#' coldataPath <- system.file("extdata/test_colData_filtered.tsv", 
#'                             package="conclus")
#' columnsMetaData <- loadColdata(file=coldataPath, columnCell="cell_ID",
#'                                 header=TRUE, dec=".", sep='\t')
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
NULL


#' @usage
#' getExperimentName(theObject)
#'
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
#' getName: Get the name of the tSNE or Dbscan object.
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
#' @return
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
#' @return
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



