library(SingleCellExperiment)

scRNAseq <- setClass(
    # Set the name for the class
    "scRNAseq",
    
    # Define the slots
    slots = c(
        experimentName = "character",
        countMatrix = "matrix", 
        normalizedCountMatrix = "SingleCellExperiment",
        colData = "data.frame",
        species = "character",
        outputDirectory = "character",
        tSNEList = "list",
        dbscanList = "list",
        cellsSimilarityMatrix = "matrix",
        clustersSimilarityMatrix = "matrix",
        clusters = "SingleCellExperiment"
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

Tsne <- setClass(
    # Set the name for the class
    "Tsne",

    # Define the slots
    slots = c(
        name = "character",
        pc = "numeric",
        perplexity = "numeric",
        coordinates = "list"
        )
)

Dbscan <- setClass(
    # Set the name for the class
    "Dbscan",
    
    # Define the slots
    slots = c(
        name = "character",
        clustering = "matrix",
        epsilon    = "numeric",
        minPoints  = "numeric"
    )
)


