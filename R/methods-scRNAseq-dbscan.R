### This function calculates dbscan for all t-SNE from tSNEList with all
### combinations of paramenters from epsilon and minPoints
### it does not set random seed. It allows to vary this parameter automatically
### it returns a matrix where columns are iterations
### number of iterations is equal to
### ncol is (tSNEList)*length(epsilon)*length(epsilon)

# tSNEList = getTSNEList(scrS4TSNE)
.mkDbscan <- function(tSNEList,
                      cores = 14,
                      epsilon = c(1.2, 1.5, 1.8),
                      minPoints = c(15, 20)){ 
    myCluster <- parallel::makeCluster(cores, # number of cores to use
                                       type = "PSOCK") # type of cluster
    doParallel::registerDoParallel(myCluster)
    dbscanResults <- foreach::foreach(i=rep(rep(1:length(tSNEList),
                                                each=length(minPoints)),
                                            length(epsilon)),
                                      eps=rep(epsilon,
                                              each=length(tSNEList)*
                                                  length(minPoints)),
                                      MinPts=rep(minPoints,
                                                 length(tSNEList)*
                                                     length(epsilon)),
                                      .combine='cbind') %dopar% {
                                          source("R/AllGenerics.R")
                                          source("R/methods-Tsne-accessors.R")
                                          fpc::dbscan(
                                              getCoordinates(tSNEList[[i]]),
                                              eps=eps,
                                              MinPts=MinPts)$cluster
                                      }
    parallel::stopCluster(myCluster)
    return(dbscanResults)
}



runDBSCAN <- function(theObject,
                      dataDirectory=getOutputDirectory(theObject),
                      cores=1,
                      epsilon=c(1.3, 1.4, 1.5),
                      minPoints=c(3, 4)){
    
    outputDataDirectory <- "output_tables"
    normalizedMatrix <- getNormalizedCountMatrix(theObject)
    tSNEList <- getTSNEList(theObject)
    experimentName <- getExperimentName(theObject)
    epsilonCombinaison <-rep(epsilon, each=length(tSNEList) * length(minPoints))
    minPtsCombinaison <- rep(minPoints, length(tSNEList) * length(epsilon))
    
    
    ## Taking only cells from the normalizedMatrix
    for(i in 1:length(tSNEList)){
      tmp <- getCoordinates(tSNEList[[i]])
      setCoordinates(tSNEList[[i]]) <- tmp[colnames(normalizedMatrix), ]
    }

    

    ## Calculating dbscan combinaisons 
    dbscanResults <- .mkDbscan(tSNEList = tSNEList,
                               cores = cores,
                               epsilon = epsilon,
                               minPoints = minPoints)
    dbscanResults <- t(dbscanResults)
    colnames(dbscanResults) <- 
      SummarizedExperiment::colData(normalizedMatrix)$cellName
    
    
    ## Output file
    write.table(dbscanResults, 
                file.path(dataDirectory,
                          outputDataDirectory,
                          paste0(experimentName,
                                 "_dbscan_results.tsv")),
                quote = FALSE,
                sep = "\t")
    
    
    ## Creation of a list of Dbscan objects
    dbscanList <- list()
    for (i in 1:(nrow(dbscanResults))){
      clustering <- t(dbscanResults[i, ])
      rownames(clustering) <- paste("clust.", i, sep = "")
      dbscanObj <- Dbscan(name       = paste("Clustering n°", i, sep =""),
                          epsilon    = epsilonCombinaison[i],
                          minPoints  = minPtsCombinaison[i],
                          clustering = clustering)
      dbscanList <- c(dbscanList, dbscanObj) 
      
    }
    
    setDbscanList(theObject) <- dbscanList
    rm(tmp)
    return(theObject)
}
