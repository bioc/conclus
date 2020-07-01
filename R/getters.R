################################################################################
########################  Getters of scRNAseq class ############################
################################################################################


setMethod(
    f = "getExperimentName",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@experimentName)
    })


setMethod(
    f = "getCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@countMatrix)
    })


setMethod(
    f = "getSceNorm",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@sceNorm)
    })


setMethod(
    f = "getSpecies",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@species)
    })



setMethod(
    f = "getOutputDirectory",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@outputDirectory)
    })



setMethod(
    f = "getTSNEList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@tSNEList)
    })


setMethod(
    f = "getDbscanList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@dbscanList)
    })


setMethod(
    f = "getCellsSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@cellsSimilarityMatrix)
    })


setMethod(
    f = "getClustersSimilarityMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@clustersSimilarityMatrix)
    })


setMethod(
		f = "getClustersSimiliratyOrdered",
		signature = "scRNAseq",
		definition = function(theObject){
			return(theObject@clustersSimiliratyOrdered)
		})


setMethod(
    f = "getMarkerGenesList",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@markerGenesList)
    })


setMethod(
    f="getClustersMarkers",
    signature="scRNAseq",
    definition = function(theObject){
        return(theObject@clustersMarkers)
    })


setMethod(
    f="getGenesInfos",
    signature="scRNAseq",
    definition = function(theObject){
        return(theObject@genesInfos)
    })


################################################################################
############################  Getter of Tsne Class  ############################
################################################################################

setMethod(
    f = "getName",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@name)
    })



setMethod(
    f = "getPerplexity",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@perplexity)
    })


setMethod(
    f = "getPC",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@pc)
    })


setMethod(
    f = "getCoordinates",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@coordinates)
    })



################################################################################
###########################  Getter of Dbscan Class  ###########################
################################################################################

setMethod(
    f = "getName",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@name)
    })




setMethod(
    f = "getEpsilon",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@epsilon)
    })


setMethod(
    f = "getMinPoints",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@minPoints)
    })


setMethod(
    f = "getClustering",
    signature = "Dbscan",
    definition = function(theObject){
        return(theObject@clustering)
    })
