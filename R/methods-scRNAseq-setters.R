## Setters for scRNAseq class

setReplaceMethod(
    f = "setExperimentName",
    signature = "scRNAseq",
    definition = function(theObject, value){
        if (!is.character(value)) stop("Experiment Name should be character")
        validObject(theObject)
        theObject@experimentName <- value
        return(theObject)
    })


setReplaceMethod(
    f = "setCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        if (!is.data.frame(value)) stop("Count Matrix should be a data frame")
        validObject(theObject)
        theObject@countMatrix <- value
        return(theObject)
    })


setReplaceMethod(
    f = "setNormalizedCountMatrix",
    signature = "scRNAseq",
    definition = function(theObject, value){
        if (!class(value) == "SingleCellExperiment") 
            stop("Normalized Count Matrix should be SingleCellExperiment")
        validObject(theObject)
        theObject@normalizedCountMatrix <- value
        return(theObject)
    })


setReplaceMethod(
    f = "setColData",
    signature = "scRNAseq",
    definition = function(theObject, value){
        if (!class(value) == "standardGeneric") 
            stop("colData should be standardGeneric")
        validObject(theObject)
        theObject@colData <- value
        return(theObject)
    })


setReplaceMethod(
    f = "setSpecies",
    signature = "scRNAseq",
    definition = function(theObject, value){
        if (!is.character(value)) stop("Species should be character")
        validObject(theObject)
        theObject@species <- value
        return(theObject)
    })


setReplaceMethod(
    f = "setOutputDirectory",
    signature = "scRNAseq",
    definition = function(theObject, value){
        if (!is.character(value)) 
            stop("outputDirectory should be SingleCellExperiment")
        validObject(theObject)
        theObject@outputDirectory <- value
        return(theObject)
    })


setReplaceMethod(
    f = "setTSNEList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        if (!is.list(value)) stop("tSNEList should be list")
        validObject(theObject)
        theObject@tSNEList <- value
        return(theObject)
    })


setReplaceMethod(
    f = "setDbscanList",
    signature = "scRNAseq",
    definition = function(theObject, value){
        if (!is.list(value)) stop("dbscanList should be list")
        validObject(theObject)
        theObject@dbscanList <- value
        return(theObject)
    })


setReplaceMethod(
    f = "setClusteringResults",
    signature = "scRNAseq",
    definition = function(theObject, value){
        if (!is.matrix(value)) stop("clusteringResults should be matrix")
        validObject(theObject)
        theObject@clusteringResults <- value
        return(theObject)
    })

