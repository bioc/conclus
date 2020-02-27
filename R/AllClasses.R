library(SingleCellExperiment)

scRNAseq <- setClass(
    # Set the name for the class
    "scRNAseq",
    
    # Define the slots
    slots = c(
        experimentName = "character",
        countMatrix = "matrix", # matrix
        normalizedCountMatrix = "SingleCellExperiment",
        colData = "data.frame",
        species = "character",
        outputDirectory = "character",
        PCs = "numeric",
        perplexities = "numeric",
        tSNEResults = "matrix", # liste sceObject
        clusteringResults = "matrix"
    ),
    
        # Set the default values for the slots. (optional)
    prototype = list(
        PCs=c(4, 6, 8, 10, 20, 40, 50),
        perplexities=c(30, 40)
        ),
        
        # Make a function that can test to see if the data is consistent.
        # This is not called if you have an initialize function defined!
    validity = function(object)
    {
        if(ncol(countMatrix) < 100) {
            return("Not enough cells in the count matrix")
        }
        
        # Vérifiez les types des slots
        return(TRUE)
    }
)

TSNE <- setClass(
    # Set the name for the class
    "TSNE",

    # Define the slots
    slots = c(
        normalizedCountMatrix = "SingleCellExperiment",
        PC = "numeric",
        perplexity = "numeric"
    ),

    # Set the default values for the slots. (optional)
    # prototype = list(
    #     PCs=c(4, 6, 8, 10, 20, 40, 50),
    #     perplexities=c(30, 40),
    #     randomSeed = 42
    # ),
    #
    # Make a function that can test to see if the data is consistent.
    # This is not called if you have an initialize function defined!
    validity = function(object)
    {
        if(ncol(countMatrix) < 100) {
            return("Not enough cells in the count matrix")
        }

        # Vérifiez les types des slots
        return(TRUE)
    }
)

