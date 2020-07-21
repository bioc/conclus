################################################################################
########################  Setters for scRNAseq class  ##########################
################################################################################

#' @usage
#' setExperimentName(theObject) <- value
#' 
#' @description
#' Update a slot of a scRNA-seq object.
#' 
#' @param theObject A scRNA-seq object to update. See description or ?scRNAseq.
#' @param value The value to update the slot with. See description or ?scRNAseq.
#' 
#' @description
#' setExperimentName: Update the experiment name slot with a character string.

#' @rdname setters-scRNAseq
#' @aliases setExperimentName<- setExperimentName
#' @title setters-scRNAseq
#' @name setters-scRNAseq
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
#' @rdname setters-scRNAseq
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
#' @rdname setters-scRNAseq
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
#' @rdname setters-scRNAseq
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
#' @rdname setters-scRNAseq
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
#' @rdname setters-scRNAseq
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
#' @rdname setters-scRNAseq
#' @aliases setDbscanList<- setDbscanList
#' @exportMethod setDbscanList<-
#' @describeIn scRNAseq-class Set the list of Dbscan objects. 
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
#' @rdname setters-scRNAseq
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
#' @rdname setters-scRNAseq
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
#' @rdname setters-scRNAseq
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
#' @rdname setters-scRNAseq
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
#' @rdname setters-scRNAseq
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
#' @rdname setters-scRNAseq
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

#' @usage
#' setName(theObject) <- value
#' 
#' @description
#' Update a slot of a Tsne-class object.
#' 
#' @param theObject A Tsne object to update. See description or ?Tsne.
#' @param value The value to update the slot with. See description or ?Tsne.
#' 
#' @description
#' setName: Update the Tsne name slot with a character string.

#' @rdname setters-Tsne
#' @aliases setName<- setName
#' @title setters-Tsne
#' @name setters-Tsne
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
#' @rdname setters-Tsne
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
#' @rdname setters-Tsne
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
#' @rdname setters-Tsne
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

#' @exportMethod setName<-
#' @describeIn Dbscan-class Set the name of the Dbscan. 
setReplaceMethod(
    f = "setName",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@name <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setEpsilon<-
#' @describeIn Dbscan-class Set the epsilon used.
setReplaceMethod(
    f = "setEpsilon",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@epsilon <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setMinPoints<-
#' @describeIn Dbscan-class Set the minPoints used.
setReplaceMethod(
    f = "setMinPoints",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@minPoints <- value
        validObject(theObject)
        return(theObject)
    })


#' @exportMethod setClustering<-
#' @describeIn Dbscan-class Set the matrix of Dbscan clustering. 
setReplaceMethod(
    f = "setClustering",
    signature = "Dbscan",
    definition = function(theObject, value){
        theObject@clustering <- value
        validObject(theObject)
        return(theObject)
    })

