################################################################################
########################  Setters for scRNAseq class  ##########################
################################################################################


#' @description
#' Update a slot of a scRNA-seq, Tsne or Dbscan object.
#'
#' @param theObject A scRNA-seq, Tsne or Dbscan object to update. See
#' description or ?scRNAseq, ?Tsne or ?Dbscan.
#' @param value The value to update the slot with. See ?scRNAseq, ?Tsne or
#' ?Dbscan.
#'
#' @rdname setters
#'
#' @examples
#' ## Load the count matrix
#' countmatrixPath <- system.file("extdata/example_countMatrix.tsv",
#'                             package="conclus")
#' countMatrix <- loadDataOrMatrix(file=countmatrixPath, type="countMatrix")
#'
#' ## Load the coldata
#' coldataPath <- system.file("extdata/example_colData.tsv",
#'                             package="conclus")
#' columnsMetaData <- loadDataOrMatrix(file=coldataPath, type="coldata",
#' columnID="cell_ID")
#'
#' ## Create the initial object
#' scr <- singlecellRNAseq(experimentName = "Bergiers",
#'                 countMatrix     = countMatrix,
#'                 species         = "mouse",
#'                 outputDirectory = "YourOutputDirectory")
#'
#' setExperimentName(scr) <- "newName"
#' setCountMatrix(scr) <- countMatrix[seq_len(15), seq_len(100)]
#' setSpecies(scr) <- "human"
#' setOutputDirectory(scr) <- "newPath"
#'
#' @title setters
#' @name setters
NULL


#' @usage
#' setExperimentName(theObject) <- value
#'
#' @return
#' setExperimentName: Update the experiment name slot with a character
#' string (scRNA-seq).
#'
#' @rdname setters
#' @aliases setExperimentName<- setExperimentName
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
#' @return
#' setCountMatrix: Update the countMatrix slot with a matrix of numeric
#' (scRNA-seq).
#'
#' @rdname setters
#' @aliases setCountMatrix<- setCountMatrix
#' @exportMethod setCountMatrix<-
setReplaceMethod(
    f = "setCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@countMatrix <- value
        validObject(theObject)
        return(theObject)
    })

#' @usage
#' setSceNorm(theObject) <- value
#'
#' @return
#' setSceNorm: Update the normalized countMatrix slot with SingleCellExperiment
#' object (scRNA-seq).
#'
#' @rdname setters
#' @aliases setSceNorm<- setSceNorm
#' @exportMethod setSceNorm<-
setReplaceMethod(
    f = "setSceNorm",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@sceNorm <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setSpecies(theObject) <- value
#'
#' @return
#' setSpecies: Update the species slot with a character string. Value should be
#' mouse or human. Other organisms can be added on demand  (scRNA-seq).
#'
#' @rdname setters
#' @aliases setSpecies<- setSpecies
#' @exportMethod setSpecies<-
setReplaceMethod(
    f = "setSpecies",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@species <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setOutputDirectory(theObject) <- value
#'
#' @return
#' setOutputDirectory: Update the outputDirectory slot with a character string.
#' Value should be a valid path (scRNA-seq).
#'
#' @rdname setters
#' @aliases setOutputDirectory<- setOutputDirectory
#' @exportMethod setOutputDirectory<-
setReplaceMethod(
    f = "setOutputDirectory",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@outputDirectory <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setTSNEList(theObject) <- value
#'
#' @return
#' setTSNEList: Update the tSNEList slot with a list of tSNE objects. See
#' ?Tsne-class (scRNA-seq).
#'
#' @rdname setters
#' @aliases setTSNEList<- setTSNEList
#' @exportMethod setTSNEList<-
setReplaceMethod(
    f = "setTSNEList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@tSNEList <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setDbscanList(theObject) <- value
#'
#' @return
#' setDbscanList: Update the dbscanList slot with a list of dbscan objects. See
#' ?Dbscan-class (scRNA-seq).
#'
#' @rdname setters
#' @aliases setDbscanList<- setDbscanList
#' @exportMethod setDbscanList<-
setReplaceMethod(
    f = "setDbscanList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@dbscanList <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setCellsSimilarityMatrix(theObject) <- value
#'
#' @return
#' setCellsSimilarityMatrix: Update the cellsSimilarityMatrix slot with a
#' numeric matrix (scRNA-seq).
#'
#' @rdname setters
#' @aliases setCellsSimilarityMatrix<- setCellsSimilarityMatrix
#' @exportMethod setCellsSimilarityMatrix<-
setReplaceMethod(
    f = "setCellsSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@cellsSimilarityMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setClustersSimilarityMatrix(theObject) <- value
#'
#' @return
#' setClustersSimilarityMatrix: Update the clustersSimilarityMatrix slot with a
#' numeric matrix (scRNA-seq).
#'
#' @rdname setters
#' @aliases setClustersSimilarityMatrix<- setClustersSimilarityMatrix
#' @exportMethod setClustersSimilarityMatrix<-
setReplaceMethod(
    f = "setClustersSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@clustersSimilarityMatrix <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setClustersSimiliratyOrdered(theObject) <- value
#'
#' @return
#' setClustersSimiliratyOrdered: Update the clustersSimilarityOrdered slot with
#' a numeric factor (scRNA-seq).
#'
#' @rdname setters
#' @aliases setClustersSimiliratyOrdered<- setClustersSimiliratyOrdered
#' @exportMethod setClustersSimiliratyOrdered<-
setReplaceMethod(
    f = "setClustersSimiliratyOrdered",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@clustersSimiliratyOrdered <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setMarkerGenesList(theObject) <- value
#'
#' @return
#' setMarkerGenesList: Update the markerGenesList slot with a list of data
#' frame. The data frame structure should be:
#' data.frame(Gene = c("gene1"), mean_log10_fdr = c(NA), n_05 = c(NA),
#' score = c(NA))  (scRNA-seq).
#'
#' @rdname setters
#' @aliases setMarkerGenesList<- setMarkerGenesList
#' @exportMethod setMarkerGenesList<-
setReplaceMethod(
    f = "setMarkerGenesList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@markerGenesList <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setClustersMarkers(theObject) <- value
#'
#' @return
#' setClustersMarkers: Update the clustersMarkers slot with a data frame.
#' The data frame structure should be:
#' data.frame(geneName="gene1", clusters=NA). (scRNA-seq)
#'
#' @rdname setters
#' @aliases setClustersMarkers<- setClustersMarkers
#' @exportMethod setClustersMarkers<-
setReplaceMethod(
        f="setClustersMarkers",
        signature="scRNAseq",
        definition = function(theObject, value){
            theObject@clustersMarkers<- value
            validObject(theObject)
            return(theObject)
        })


#' @usage
#' setGenesInfos(theObject) <- value
#'
#' @return
#' setGenesInfos: Update the genesInfos slot with a data frame.
#' The data frame structure should be:
#' data.frame(uniprot_gn_symbol=c("symbol"), clusters="1",
#'                 external_gene_name="gene", go_id="GO1,GO2",
#'                 mgi_description="description",
#'                 entrezgene_description="descr",
#'                 gene_biotype="gene", chromosome_name="1", Symbol="symbol",
#'                 ensembl_gene_id="ENS", mgi_id="MGI", entrezgene_id="1",
#'                 uniprot_gn_id="ID"). (scRNA-seq)
#'
#' @rdname setters
#' @aliases setGenesInfos<- setGenesInfos
#' @exportMethod setGenesInfos<-
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


#' @rdname setters
setReplaceMethod(
        f = "setName",
        signature = c("Tsne"),
        definition = function(theObject, value){
            theObject@name <- value
            validObject(theObject)
            return(theObject)
        })

#' @usage
#' setPC(theObject) <- value
#'
#' @return
#' setPC: Update the pc slot with a vector of numeric (Tsne).
#'
#' @rdname setters
#' @aliases setPC<- setPC
#' @exportMethod setPC<-
setReplaceMethod(
    f = "setPC",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@pc <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setPerplexity(theObject) <- value
#'
#' @return
#' setPerplexity: Update the perplexity slot with a vector of numeric (Tsne).
#'
#' @rdname setters
#' @aliases setPerplexity<- setPerplexity
#' @exportMethod setPerplexity<-
setReplaceMethod(
    f = "setPerplexity",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@perplexity <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setCoordinates(theObject) <- value
#'
#' @return
#' setCoordinates: Update the coordinates slot with a matrix of numeric (Tsne).
#'
#' @rdname setters
#' @aliases setCoordinates<- setCoordinates
#' @exportMethod setCoordinates<-
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

#' @usage
#' setName(theObject) <- value
#'
#' @return
#' setName: Update the Tsne or Dbscan name slot with a character string
#' (Dbscan).
#'
#' @rdname setters
#' @aliases setName<-
#'
#' @exportMethod setName<-
setReplaceMethod(
        f = "setName",
        signature = c("Dbscan"),
        definition = function(theObject, value){
            theObject@name <- value
            validObject(theObject)
            return(theObject)
        })

#' @usage
#' setEpsilon(theObject) <- value
#'
#' @return
#' setEpsilon: Update the epsilon slot with a vector of numeric (Dbscan).
#'
#' @rdname setters
#' @aliases setEpsilon<- setEpsilon
#'
#' @exportMethod setEpsilon<-
setReplaceMethod(
    f = "setEpsilon",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@epsilon <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setMinPoints(theObject) <- value
#'
#' @return
#' setMinPoints: Update the minPoints slot with a vector of numeric (Dbscan).
#'
#' @rdname setters
#' @aliases setMinPoints<- setMinPoints
#'
#' @exportMethod setMinPoints<-
setReplaceMethod(
    f = "setMinPoints",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@minPoints <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setClustering(theObject) <- value
#'
#' @return
#' setClustering: Update the clustering slot with a matrix of numeric (Dbscan).
#'
#' @rdname setters
#' @aliases setClustering<- setClustering
#'
#' @exportMethod setClustering<-
setReplaceMethod(
    f = "setClustering",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@clustering <- value
        validObject(theObject)
        return(theObject)
    })
