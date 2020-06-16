setMethod(
    f="exportResults",
    signature="scRNAseq",
    definition=function(theObject,
                        saveCellsSimilarityMatrix=TRUE,
                        saveClustersSimilarityMatrix=TRUE,
                        saveNormalizedMatrix=TRUE,
                        saveColData=TRUE,
                        saveRowData=TRUE,
                        saveWorkspace=TRUE,
                        saveClusteringResults=TRUE
                        ){
  
        dataDirectory  <- getOutputDirectory(theObject)
        experimentName <- getExperimentName(theObject)
        outputDataDirectory <- "output_tables"

        ## Export CellsSimilarityMatrix
        if (saveCellsSimilarityMatrix){
            matrix <- getCellsSimilarityMatrix(theObject)
            if (all(dim(matrix) == c(0,0))){
                msg <- paste("The 'scRNAseq' object that you're using with ",
                             "'exportResults' function doesn't have ",
                             "its 'cellsSimilarityMatrix' slot updated.\n",
                             "Please use 'clusterCellsInternal' on the object ",
                             "before.\n",
                             sep="")
                stop(msg)
            }
            
            fileName <- paste0(experimentName, "_",
                               "cells_similarity_matrix", ".csv")
            write.table(matrix,
                        file=file.path(dataDirectory,
                                       "output_tables", fileName),
                        sep=",")
            message("Cells similarity matrix saved.")
        }
        
        ## Export ClustersSimilarityMatrix
        if (saveClustersSimilarityMatrix){
            matrix <- getClustersSimilarityMatrix(theObject)
            if (all(dim(matrix) == c(0,0))){
                msg <- paste("The 'scRNAseq' object that you're using with ",
                             "'exportResults' function doesn't have ",
                             "its 'clustersSimilarityMatrix' slot updated.",
                             "\nPlease use 'clusterCellsInternal' on ",
                             "the object before.\n",
                             sep="")
                stop(msg)
            }
            
            fileName <- paste0(experimentName, "_", 
                               "clusters_similarity_matrix", ".csv")
            write.table(matrix,
                        file=file.path(dataDirectory,
                                       "output_tables", fileName),
                        sep=",")
            message("Clusters similarity matrix saved.")
        }
        
        ## Check the SingleCellExperiment object
        sceObject <- getSceNorm(theObject)
        if  (!("clusters" %in% names(colData(sceObject)))){
            msg <- paste("The 'scRNAseq' object that you're using with ",
                         "'calculateClustersSimilarity' function doesn't  ",
                         "have a correct 'sceNorm' slot.\n",
                         "This slot should be a 'SingleCellExperiment' ",
                         "object containing 'clusters' column in its ",
                         "colData.\nPlease check if you correctly used ",
                         "'clusterCellsInternal' on the object.\n",
                         sep="")
            stop(msg)
            
        } else {
            ## Export Normalized expression matrix
            if (saveNormalizedMatrix){
                write.table(Biobase::exprs(sceObject),
                            file=file.path(dataDirectory,
                                           outputDataDirectory,
                                           paste0(experimentName, "_",
                                                  "expression_matrix.tsv")),
                            sep="\t",
                            row.names=TRUE, quote=FALSE, col.names=TRUE)
              message("Normalized expression matrix saved.")
            }
        
            ## Export RowData
            if (saveRowData){
                write.table(rowData(sceObject),
                            file=file.path(dataDirectory,
                                           outputDataDirectory,
                                           paste0(experimentName, "_",
                                                  "rowData.tsv")), sep="\t",
                            row.names=TRUE, quote=FALSE, col.names=TRUE)
              message("RowData saved.")
            }
            
            ## Export ColData
            if (saveColData){
                write.table(SummarizedExperiment::colData(sceObject),
                            file=file.path(dataDirectory,
                                           outputDataDirectory,
                                           paste0(experimentName, "_",
                                                  "colData.tsv")), sep="\t",
                            row.names=TRUE, quote=FALSE, col.names=TRUE)
              message("ColData saved.")
            }
        }
        
        ## Save workspace
        if (saveWorkspace){
            save.image(file=file.path(dataDirectory,
                                      outputDataDirectory,
                                      paste0(experimentName, "_",
                                             "full_workspace.RData")))
          message("Workspace saved.")
        }
        
        ## Export Clustering results
        if (saveClusteringResults){
            tableData <- S4Vectors::DataFrame(
                clusters=SummarizedExperiment::colData(sceObject)$clusters,
                row.names=SummarizedExperiment::colData(sceObject)$cellName)
            
            write.table(tableData,
                        file=file.path(dataDirectory, "output_tables",
                                         paste0(experimentName,"_", 
                                                "clusters_table.tsv")),
                        sep="\t", quote=FALSE)
            message("Clusters table saved.")
        }

})
