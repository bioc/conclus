## Getters of scRNAseq class

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
    f = "getNormalizedCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@normalizedCountMatrix)
    })


setMethod(
    f = "getColData",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@colData)
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
    f = "getClusteringResults",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@clusteringResults)
    })

