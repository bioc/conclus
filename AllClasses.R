setwd("/g/lancrin/People/Ilyess/conclus_S4/R")

scRNAseq <- setClass(
    # Set the name for the class
    "scRNAseq",
    
    # Define the slots
    slots = c(
        experimentName = "character",
        countMatrix = "data.frame", # matrix
        normalizedCountMatrix = "SingleCellExperiment",
        colData = "data.frame",
        species = "character",
        sceObject = "SingleCellExperiment",
        outputDirectory = "character",
        PCs = "numeric",
        perplexities = "numeric",
        randomSeed = "numeric",
        tSNEResults = "matrix", # liste sceObject
        clusteringResults = "matrix"
    ),
    
        # Set the default values for the slots. (optional)
    prototype = list(
        PCs=c(4, 6, 8, 10, 20, 40, 50),
        perplexities=c(30, 40),
        randomSeed = 42
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
    
 



        # Set the inheritance for this class
    # contains = 
# )


# scRNAseq <- setClass(
#     # Set the name for the class
#     "Tsne",
#     
#     # Define the slots
#     slots = c(
#         PC = "numeric",
#         perplexity = "numeric",
#         randomSeed = "numeric",
#         tSNEResults = "matrix",
#         cores = "numeric",
#         clusteringResults = " ", 
#         
#         
#         
#     ),
#     # Set the default values for the slots. (optional)
#     prototype=
#         
#         # Make a function that can test to see if the data is consistent.
#         # This is not called if you have an initialize function defined!
#         validity=
#         
#         # Set the inheritance for this class
#         contains = 
# )














# 
# 
# ## normaliseCountMatrix
# setGeneric(name="normaliseCountMatrix",
#            def=function(getCountMatrix, getSpecies, 
#                         getColData)
#                
#            {
#                standardGeneric("normaliseCountMatrix")
#            }
# )
# 
# setMethod(f="normaliseCountMatrix",
#           signature="scRNAseq",
#           definition=function(getCountMatrix, getSpecies, 
#                               getColData)
#           {
# 
#           }
# )
# 
# ## testClustering
# setGeneric(name="testClustering",
#            def=function(getCountMatrix, getSpecies, 
#                         getColData)
#                
#            {
#                standardGeneric("normaliseCountMatrix")
#            }
# )
# 
# setMethod(f="testClustering",
#           signature="scRNAseq",
#           definition=function(getSceObject, getOutputDirectory, 
#                               getExperimentName)
#           {
#               
#           }
# )
# 


