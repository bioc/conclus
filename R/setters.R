################################################################################
########################  Setters for scRNAseq class  ##########################
################################################################################


#' @description
#' Update a slot of a scRNA-seq object.
#' 
#' @param theObject A scRNA-seq object to update. See description or ?scRNAseq.
#' @param value The value to update the slot with. See description or ?scRNAseq.
#' 
#' @rdname settersscRNAseq
#' @title settersscRNAseq
#' @name settersscRNAseq
NULL


#' @usage
#' setExperimentName(theObject) <- value
#' 
#' @description
#' setExperimentName: Update the experiment name slot with a character string.
#'  
#' @rdname settersscRNAseq
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
#' @description
#' setCountMatrix: Update the countMatrix slot with a matrix of numeric.
#' 
#' @rdname settersscRNAseq
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
#' @description
#' setSceNorm: Update the normalized countMatrix slot with SingleCellExperiment 
#' object.
#' 
#' @rdname settersscRNAseq
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
#' @description
#' setSpecies: Update the species slot with a character string. Value should be
#'  mouse or human.
#' 
#' @rdname settersscRNAseq
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
#' @description
#' setOutputDirectory: Update the outputDirectory slot with a character string.
#' Value should be a valid path.
#' 
#' @rdname settersscRNAseq
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
#' @description
#' setTSNEList: Update the tSNEList slot with a list of tSNE objects. See 
#' ?Tsne-class.
#' 
#' @rdname settersscRNAseq
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
#' @description
#' setDbscanList: Update the dbscanList slot with a list of dbscan objects. See 
#' ?Dbscan-class.
#' 
#' @rdname settersscRNAseq
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
#' @description
#' setCellsSimilarityMatrix: Update the cellsSimilarityMatrix slot with a 
#' numeric matrix.
#' 
#' @rdname settersscRNAseq
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
#' @description
#' setClustersSimilarityMatrix: Update the clustersSimilarityMatrix slot with a 
#' numeric matrix.
#' 
#' @rdname settersscRNAseq
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
#' @description
#' setClustersSimiliratyOrdered: Update the clustersSimilarityOrdered slot with 
#' a numeric factor.
#' 
#' @rdname settersscRNAseq
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
#' @description
#' setMarkerGenesList: Update the markerGenesList slot with a list of data 
#' frame. The data frame structure should be: 
#' data.frame(Gene = c("gene1"), mean_log10_fdr = c(NA), n_05 = c(NA), 
#' score = c(NA))
#' 
#' @rdname settersscRNAseq
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
#' @description
#' setClustersMarkers: Update the clustersMarkers slot with a data frame. 
#' The data frame structure should be: 
#' data.frame(geneName="gene1", clusters=NA)
#' 
#' @rdname settersscRNAseq
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
#' @description
#' setGenesInfos: Update the genesInfos slot with a data frame. 
#' The data frame structure should be: 
#' data.frame(uniprot_gn_symbol=c("symbol"), clusters="1",
#' 				external_gene_name="gene", go_id="GO1,GO2", 
#' 				mgi_description="description", entrezgene_description="descr",
#' 				gene_biotype="gene", chromosome_name="1", Symbol="symbol",
#' 				ensembl_gene_id="ENS", mgi_id="MGI", entrezgene_id="1",
#' 				uniprot_gn_id="ID")
#' 
#' @rdname settersscRNAseq
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


#' @description
#' Update a slot of a Tsne-class object.
#' 
#' @param theObject A Tsne object to update. See description or ?Tsne.
#' @param value The value to update the slot with. See description or ?Tsne.
#' 
#' @rdname settersTsne
#' @title settersTsne
#' @name settersTsne
NULL


#' @usage
#' setName(theObject) <- value
#'  
#' @description
#' setName: Update the Tsne name slot with a character string.
#' 
#' @rdname settersTsne
#' @aliases setName<- setName
#' 
#' @exportMethod setName<- 
setReplaceMethod(
    f = "setName",
    signature = "Tsne",
    definition = function(theObject, value){
        theObject@name <- value
        validObject(theObject)
        return(theObject)
    })

#' @usage
#' setPC(theObject) <- value
#' 
#' @description
#' setPC: Update the pc slot with a vector of numeric.
#' 
#' @rdname settersTsne
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
#' @description
#' setPerplexity: Update the perplexity slot with a vector of numeric.
#' 
#' @rdname settersTsne
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
#' @description
#' setCoordinates: Update the coordinates slot with a matrix of numeric.
#' 
#' @rdname settersTsne
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


#' @description
#' Update a slot of a Dbscan-class object.
#' 
#' @param theObject A Dbscan object to update. See description or ?Dbscan.
#' @param value The value to update the slot with. See description or ?Dbscan.
#' 
#' @rdname settersDbscan
#' @title settersDbscan
#' @name settersDbscan
NULL

#' @usage
#' setName(theObject) <- value
#' 
#' @description
#' setName: Update the Dbscan name slot with a character string.
#' 
#' @rdname settersDbscan
#' @aliases setName<- setName
#' 
#' @exportMethod setName<- 
setReplaceMethod(
    f = "setName",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@name <- value
        validObject(theObject)
        return(theObject)
    })


#' @usage
#' setEpsilon(theObject) <- value
#' 
#' @description
#' setEpsilon: Update the epsilon slot with a vector of numeric.
#' 
#' @rdname settersDbscan
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
#' @description
#' setMinPoints: Update the minPoints slot with a vector of numeric.
#' 
#' @rdname settersDbscan
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
#' @description
#' setClustering: Update the clustering slot with a matrix of numeric.
#' 
#' @rdname settersDbscan
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

