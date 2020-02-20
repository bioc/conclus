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
    f = "getSceObject",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@sceObject)
    })


setMethod(
    f = "getOutputDirectory",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@outputDirectory)
    })


setMethod(
    f = "getPCs",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@PCs)
    })


setMethod(
    f = "getPerplexities",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@perplexities)
    })


setMethod(
    f = "getRandomSeed",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@randomSeed)
    })


setMethod(
    f = "getTSNEResults",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@tSNEResults)
    })




setMethod(
    f = "getClusteringResults",
    signature = "scRNAseq",
    definition = function(theObject){
        return(theObject@clusteringResults)
    })

