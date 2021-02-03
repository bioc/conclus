################################################################################
########################  Setters for scRNAseq class  ##########################
################################################################################

#' @title setters
#'
#' @description
#' Update a slot of a scRNA-seq, Tsne or Dbscan object.
#'
#' @param theObject A scRNA-seq, Tsne or Dbscan object to update. See
#' description or ?scRNAseq, ?Tsne or ?Dbscan.
#' @param value The value to update the slot with. See ?scRNAseq, ?Tsne or
#' ?Dbscan.
#'
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
#' @name setters
#' @encoding UTF-8
#' @docType methods
#' @aliases
#' setExperimentName<-
#' setCountMatrix<-
#' setSpecies<-
#' setOutputDirectory<-
#' setSceNorm<-
#' setTSNEList<-
#' setDbscanList<-
#' setCellsSimilarityMatrix<-
#' setClustersSimilarityMatrix<-
#' setClustersSimiliratyOrdered<-
#' setMarkerGenesList<-
#' setClustersMarkers<-
#' setGenesInfos<-
#' setPC<-
#' setPerplexity<-
#' setCoordinates<-
#' setName<-
#' setEpsilon<-
#' setMinPoints<-
#' setClustering<-
#'
#' @usage
#' setExperimentName(theObject) <- value
#' setCountMatrix(theObject)(theObject) <- value
#' setSpecies(theObject) <- value
#' setOutputDirectory(theObject) <- value
#' setSceNorm(theObject) <- value
#' setTSNEList(theObject) <- value
#' setDbscanList(theObject) <- value
#' setCellsSimilarityMatrix(theObject) <- value
#' setClustersSimilarityMatrix(theObject) <- value
#' setClustersSimiliratyOrdered(theObject) <- value
#' setMarkerGenesList(theObject) <- value
#' setClustersMarkers(theObject) <- value
#' setGenesInfos(theObject) <- value
#' setPC(theObject) <- value
#' setPerplexity(theObject) <- value
#' setCoordinates(theObject) <- value
#' setName(theObject) <- value
#' setEpsilon(theObject) <- value
#' setMinPoints(theObject) <- value
#' setClustering(theObject) <- value
#'
#' @return
#' setExperimentName: Update the experiment name slot with a character
#' string (scRNA-seq).
#' setCountMatrix: Update the countMatrix slot with a matrix of numeric
#' (scRNA-seq).
#' setSceNorm: Update the normalized countMatrix slot with SingleCellExperiment
#' object (scRNA-seq).
#' setSpecies: Update the species slot with a character string. Value should be
#' mouse or human. Other organisms can be added on demand  (scRNA-seq).
#' setOutputDirectory: Update the outputDirectory slot with a character string.
#' Value should be a valid path (scRNA-seq).
#' setTSNEList: Update the tSNEList slot with a list of tSNE objects. See
#' ?Tsne-class (scRNA-seq).
#' setDbscanList: Update the dbscanList slot with a list of dbscan objects. See
#' ?Dbscan-class (scRNA-seq).
#' setCellsSimilarityMatrix: Update the cellsSimilarityMatrix slot with a
#' numeric matrix (scRNA-seq).
#' setClustersSimilarityMatrix: Update the clustersSimilarityMatrix slot with a
#' numeric matrix (scRNA-seq).
#' setClustersSimiliratyOrdered: Update the clustersSimilarityOrdered slot with
#' a numeric factor (scRNA-seq).
#' setMarkerGenesList: Update the markerGenesList slot with a list of data
#' frame. The data frame structure should be:
#' data.frame(Gene = c("gene1"), mean_log10_fdr = c(NA), n_05 = c(NA),
#' score = c(NA))  (scRNA-seq).
#' setClustersMarkers: Update the clustersMarkers slot with a data frame.
#' The data frame structure should be:
#' data.frame(geneName="gene1", clusters=NA). (scRNA-seq)
#' setGenesInfos: Update the genesInfos slot with a data frame.
#' The data frame structure should be:
#' data.frame(uniprot_gn_symbol=c("symbol"), clusters="1",
#'                 external_gene_name="gene", go_id="GO1,GO2",
#'                 mgi_description="description",
#'                 entrezgene_description="descr",
#'                 gene_biotype="gene", chromosome_name="1", Symbol="symbol",
#'                 ensembl_gene_id="ENS", mgi_id="MGI", entrezgene_id="1",
#'                 uniprot_gn_id="ID"). (scRNA-seq)
#' setPC: Update the pc slot with a vector of numeric (Tsne).
#' setPerplexity: Update the perplexity slot with a vector of numeric (Tsne).
#' setCoordinates: Update the coordinates slot with a matrix of numeric (Tsne).
#' setName: Update the Tsne or Dbscan name slot with a character string
#' (Dbscan).
#' setEpsilon: Update the epsilon slot with a vector of numeric (Dbscan).
#' setMinPoints: Update the minPoints slot with a vector of numeric (Dbscan).
#' setClustering: Update the clustering slot with a matrix of numeric (Dbscan).
NULL



#' @export
setReplaceMethod(
    f = "setExperimentName",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@experimentName <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
    f = "setCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@countMatrix <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
    f = "setSceNorm",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@sceNorm <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
    f = "setSpecies",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@species <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
    f = "setOutputDirectory",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@outputDirectory <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
    f = "setTSNEList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@tSNEList <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
    f = "setDbscanList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@dbscanList <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
    f = "setCellsSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@cellsSimilarityMatrix <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
    f = "setClustersSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@clustersSimilarityMatrix <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
    f = "setClustersSimiliratyOrdered",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@clustersSimiliratyOrdered <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
    f = "setMarkerGenesList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        theObject@markerGenesList <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
        f="setClustersMarkers",
        signature="scRNAseq",
        definition = function(theObject, value){
            theObject@clustersMarkers<- value
            validObject(theObject)
            return(theObject)
        })



#' @export
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


#' @export
setReplaceMethod(
        f = "setName",
        signature = c("Tsne"),
        definition = function(theObject, value){
            theObject@name <- value
            validObject(theObject)
            return(theObject)
        })


#' @export
setReplaceMethod(
    f = "setPC",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@pc <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
    f = "setPerplexity",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@perplexity <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
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


#' @export
setReplaceMethod(
        f = "setName",
        signature = c("Dbscan"),
        definition = function(theObject, value){
            theObject@name <- value
            validObject(theObject)
            return(theObject)
        })


#' @export
setReplaceMethod(
    f = "setEpsilon",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@epsilon <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
    f = "setMinPoints",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@minPoints <- value
        validObject(theObject)
        return(theObject)
    })



#' @export
setReplaceMethod(
    f = "setClustering",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@clustering <- value
        validObject(theObject)
        return(theObject)
    })
