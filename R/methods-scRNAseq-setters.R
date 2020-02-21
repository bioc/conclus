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
    f = "setSceObject",
    signature = "scRNAseq",
    definition = function(theObject, value){
        if (!class(value) == "SingleCellExperiment") 
            stop("sceObject should be SingleCellExperiment")
        validObject(theObject)
        theObject@sceObject <- value
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
    f = "setPCs",
    signature = "scRNAseq",
    definition = function(theObject, value){
        if (!is.double(value)) stop("PCs should be double")
        validObject(theObject)
        theObject@PCs <- value
        return(theObject)
    })


setReplaceMethod(
    f = "setPerplexities",
    signature = "scRNAseq",
    definition = function(theObject, value){
        if (!is.double(value)) stop("Perplexities should be double")
            validObject(theObject)
            theObject@perplexities <- value
            return(theObject)
    })


setReplaceMethod(
    f = "setTSNEResults",
    signature = "scRNAseq",
    definition = function(theObject, value){
        if (!is.matrix(value)) stop("tSNEResults should be matrix")
        validObject(theObject)
        theObject@tSNEResults <- value
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


        

        
        
        





