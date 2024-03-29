#' .mkDbscan
#'
#' @description
#' Use doParallel to calculate clusters with dbScan algorithm.
#'
#' @param tSNEList List of tSNE coordinates obtained with
#' ?generateTSNECoordinates.
#' @param cores The number of cores to use for parallelisation. Default = 1.
#' @param epsilon Reachability distance parameter of fpc::dbscan() function.
#' See Ester et al. (1996) for more details. Default = c(1.3, 1.4, 1.5)
#' @param minPoints Reachability minimum no. of points parameter of
#' fpc::dbscan() function. See Ester et al. (1996) for more details.
#' Default = c(3, 4)
#'
#' @details
#' This function calculates dbscan for all t-SNE from tSNEList with all
#' combinations of parameters from epsilon and minPoints. It does not set
#' random seed. It allows to vary this parameter automatically. It returns a
#' matrix where columns are iterations. The number of iterations is equal to
#' (tSNEList)*length(epsilon)*length(minPoints).
#'
#' @keywords internal
#'
#' @importFrom fpc dbscan
#' @return Returns a matrix of the combinations of dbscan results
#' @noRd
.mkDbscan <- function(tSNEList, cores=cores, epsilon=c(1.3, 1.4, 1.5),
        minPoints=c(3, 4)){
    myCluster <- parallel::makeCluster(cores, type="PSOCK")
    doParallel::registerDoParallel(myCluster)
    dbscanResults <- foreach::foreach(iMkDbscan=rep(rep(
                                    seq_len(length(tSNEList)),
                                    each=length(minPoints)), length(epsilon)),
                    epsMkDbscan=rep(epsilon,
                            each=length(tSNEList)*length(minPoints)),
                    MinPtsMkDbscan=rep(minPoints,
                            length(tSNEList)*length(epsilon)),
                    .combine='cbind', .export = 'getCoordinates') %dopar% {
                fpc::dbscan(getCoordinates(tSNEList[[iMkDbscan]]),
                        eps=epsMkDbscan,
                        MinPts=MinPtsMkDbscan)$cluster}
    parallel::stopCluster(myCluster)
    return(dbscanResults)
}


#' .checkParamDbScan
#'
#' @description
#' Check the parameters of the method runDBSCAN.
#'
#' @param sceObject Normalized count matrix retrieved with ?getSceNorm.
#' @param tSNEList List of tSNE retrieved with getTSNEList.
#' @param cores Maximum number of jobs that CONCLUS can run in parallel.
#' Default is 1.
#' @param epsilon Reachability distance parameter of fpc::dbscan() function.
#' See Ester et al. (1996) for more details. Default = c(1.3, 1.4, 1.5)
#' @param minPoints Reachability minimum no. of points parameter of
#' fpc::dbscan() function. See Ester et al. (1996) for more details.
#' Default = c(3, 4)
#' @param writeOutput If TRUE, write the results of the dbScan clustering to
#' the output directory defined in theObject, in the sub-directory
#' output_tables. Default = FALSE.
#'
#' @keywords internal
#' @noRd
.checkParamDbScan <- function(sceObject, tSNEList, cores, epsilon, minPoints,
        writeOutput){

    ## Check if the normalized count matrix is correct
    if(all(dim(sceObject) == c(0,0)))
        stop("The 'scRNAseq' object that you're using with 'runDBSCAN' ",
                "function doesn't have its 'sceNorm' slot updated. Please ",
                "use 'normaliseCountMatrix' on the object before.")

    ## Check if the tsne list is correct
    if(length(tSNEList) <= 1)
        stop("The 'scRNAseq' object that you're using with 'runDBSCAN' ",
                "function doesn't have its 'tSNEList' slot updated. Please",
                " use 'generateTSNECoordinates' on the object before.")

    ## Check cores argument
    if(!is.numeric(cores))
        stop("'cores' parameter should be an integer")

    ## Check epsilon argument
    if(!is.numeric(epsilon))
        stop("'epsilon' parameter should be a vector of numeric.")

    ## Check minPoints argument
    if(!is.numeric(minPoints))
        stop("'minPoints' parameter should be a vector of numeric.")

    ## Check writeOutput argument
    if(!is.logical(writeOutput))
        stop("'writeOutput' parameter should be a boolean.")

}


#' .writeDBScanResults
#'
#' @description
#' Export the dbScan clustering results to an output folder if the parameter
#' \code{writeOutput} of the method \code{runDBSCAN} is TRUE.
#'
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized and the tSNE coordinates were calculated.
#' See ?normaliseCountMatrix and ?generateTSNECoordinates.
#' @param dbscanResults Result of the function fpc::dbscan that is called in the
#' internal function .mkDbscan.
#'
#' @keywords internal
#'
#' @return Nothing. Write dbScan clustering results to the output directory in
#' the sub-directory output_tables.
#' @noRd
.writeDBScanResults <- function(theObject, dbscanResults){

    ## Get the values of slots
    dataDirectory  <- getOutputDirectory(theObject)
    experimentName <- getExperimentName(theObject)
    outputDataDirectory <- "output_tables"
    filename <- paste0(experimentName, "_dbscan_results.tsv")
    outputfile <- file.path(dataDirectory, outputDataDirectory, filename)
    write.table(dbscanResults, outputfile, quote=FALSE, sep="\t")
}


#' .dbscanComb
#'
#' @description
#' Compute the clustering combinations using dbscan.
#'
#' @param tSNEList List of tSNE retrieved with getTSNEList.
#' @param cores Maximum number of jobs that CONCLUS can run in parallel.
#' Default is 1.
#' @param epsilon Reachability distance parameter of fpc::dbscan() function.
#' See Ester et al. (1996) for more details. Default = c(1.3, 1.4, 1.5)
#' @param minPoints Reachability minimum no. of points parameter of
#' fpc::dbscan() function. See Ester et al. (1996) for more details.
#' Default = c(3, 4)
#' @param sceObject Normalized count matrix retrieved with ?getSceNorm.
#'
#' @keywords internal
#' @importFrom SummarizedExperiment colData
#' @return The results of the dbscan clustering.
#' @noRd
.dbscanComb <- function(tSNEList, cores, epsilon, minPoints, sceObject){

    dbscanResults <- .mkDbscan(tSNEList=tSNEList, cores=cores,
            epsilon=epsilon, minPoints=minPoints)
    dbscanResults <- t(dbscanResults)
    colnames(dbscanResults) <-
            SummarizedExperiment::colData(sceObject)$cellName

    return(dbscanResults)

}


#' .createDbscanList
#'
#' @description
#' Create a list of dbscan object to pass to the main object.
#'
#' @param tSNEList List of tSNE retrieved with getTSNEList.
#' @param minPoints Reachability minimum no. of points parameter of
#' fpc::dbscan() function. See Ester et al. (1996) for more details.
#' Default = c(3, 4)
#' @param epsilon Reachability distance parameter of fpc::dbscan() function.
#' See Ester et al. (1996) for more details. Default = c(1.3, 1.4, 1.5)
#' @param dbscanResults The result of the function .dbscanComb.
#'
#' @keywords internal
#' @return The results of the dbscan clustering.
#' @noRd
.createDbscanList <- function(tSNEList, minPoints, epsilon, dbscanResults){

    totalLength <- length(tSNEList) * length(minPoints)
    epsilonCombinaison <- rep(epsilon, each=totalLength)
    minPtsCombinaison  <- rep(minPoints,
            length(tSNEList)*length(epsilon))
    rowDbscanList <- split(dbscanResults, seq_len(nrow(dbscanResults)))
    rowNamesVec <- paste0("clust.", seq_len(nrow(dbscanResults)))
    dbscanObjNameVec <- paste0("Clustering_", seq_len(nrow(dbscanResults)))

    dbscanList <- mapply(function(rowName, dbscanObjName, epsComb,
                    minPts, rowDbscan){

                clustering <- t(rowDbscan)
                colnames(clustering) <- colnames(dbscanResults)
                rownames(clustering) <- rowName
                dbscanObj<- DbscanCluster(name= dbscanObjName, epsilon=epsComb,
                        minPoints=minPts, clustering=clustering)
                return(dbscanObj)

            }, rowNamesVec, dbscanObjNameVec, epsilonCombinaison,
            minPtsCombinaison, rowDbscanList, SIMPLIFY = FALSE,
            USE.NAMES = FALSE)

    return(dbscanList)
}



#' .retrieveClustersNumberK
#'
#' @description
#' Use the list of dbscan object to find which clusters number would be the 
#' most suitable
#'
#' @param dbscanList List of dbscan solutions retrieved with getDbscanList.
#'
#' @keywords internal
#' @return Suggested clusters number to use in clusterCellsInternal(), 
#' corresponding to the median value of the found clusters number among the
#' dbscan solutions.
#' @importFrom stats median
#' @noRd
.retrieveClustersNumberK <- function(dbscanList){

    ## Get for each dbscan the cluster number, represented by the highest number
    ## in the clustering
    l <- vapply(dbscanList, function(obj){
                max(obj@clustering)}, numeric(1))

    ## Put the result in a decreasing table
    tab <- sort(decreasing = TRUE, table(l))
    mat <- rbind(as.numeric(names(tab)), tab)
    rownames(mat) <- c("Number of clusters k :", "Count :")
    colnames(mat) <- seq(tab)

    msg <- paste("The following matrix shows how many times a number of",
                "clusters 'k' has been found among the dbscan solutions :\n")
    cat(msg)
    print(mat)
    
    cat("\nStatistics about number of clusters 'k' among dbscan solutions:\n")
    print(summary(l))
    
    median_k <- round(median(l))
    
    msg <- paste0("\nSuggested clusters number to use in " , 
        "clusterCellsInternal() : clusterNumber=", median_k, ".")
    cat(msg)
    
    return(median_k)
}


#' runDBSCAN
#'
#' @description
#' Run clustering iterations with selected parameters using DBSCAN.
#'
#' @usage
#' runDBSCAN(theObject, cores=2, epsilon=c(1.3, 1.4, 1.5), minPoints=c(3, 4),
#' writeOutput=FALSE)
#'
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized and the tSNE coordinates were calculated.
#' See ?normaliseCountMatrix and ?generateTSNECoordinates.
#' @param cores Maximum number of jobs that CONCLUS can run in parallel.
#' Default is 1.
#' @param epsilon Reachability distance parameter of fpc::dbscan() function.
#' See Ester et al. (1996) for more details. Default = c(1.3, 1.4, 1.5)
#' @param minPoints Reachability minimum no. of points parameter of
#' fpc::dbscan() function. See Ester et al. (1996) for more details.
#' Default = c(3, 4)
#' @param writeOutput If TRUE, write the results of the dbScan clustering to
#' the output directory defined in theObject, in the sub-directory
#' output_tables. Default = FALSE.
#'
#' @aliases runDBSCAN
#'
#' @details
#' Following the calculation of t-SNE coordinates, DBSCAN is run with a range
#' of epsilon and MinPoints values which will yield a total of 84 clustering
#' solutions (PCs x perplexities x MinPoints x epsilon). minPoints is the
#' minimum cluster size which you assume to be meaningful for your experiment
#' and epsilon is the radius around the cell where the algorithm will try to
#' find minPoints dots. Optimal epsilon must lay one the knee of the k-NN
#' function as shown in the "test_clustering/distance_graph.pdf".
#'
#' @author
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
#'
#' @rdname
#' runDBSCAN-scRNAseq
#'
#' @return
#' An object of class scRNASeq with its dbscanList slot updated. Also writes
#' the clustering results in "dataDirectory/output_tables" subfolder if the
#' parameter writeOutput is TRUE.
#'
#' @examples
#' ## Object scr containing the results of previous steps
#' load(system.file("extdata/scrFull.Rdat", package="conclus"))
#'
#' ## Perform the clustering with dbScan
#' ## These parameters are tweaked to fit our example data and reduce
#' ## computing time, please consider using the default parameters or
#' ## adjusted to your dataset.
#' scr <- runDBSCAN(scr, epsilon=c(380, 390, 400), minPoints=c(2,3), cores=2)
#'
#' @exportMethod runDBSCAN
#'
#' @importFrom methods validObject
#'
#' @seealso
#' normaliseCountMatrix generateTSNECoordinates
#'
setMethod(

        f = "runDBSCAN",

        signature = "scRNAseq",

        definition = function(theObject, cores=2, epsilon=c(1.3, 1.4, 1.5),
                minPoints=c(3, 4), writeOutput=FALSE){


            ## Check the Object
            validObject(theObject)

            ## Check function parameters
            sceObject <- getSceNorm(theObject)
            tSNEList <- getTSNEList(theObject)
            .checkParamDbScan(sceObject, tSNEList, cores, epsilon, minPoints,
                    writeOutput)

            ## Taking only cells from the sceObject
            tSNEList <- lapply(tSNEList, function(currentTsne, sceObj){

                        tmp <- getCoordinates(currentTsne)
                        setCoordinates(currentTsne) <- tmp[colnames(sceObj), ]
                        return(currentTsne)
                    }, sceObject)

            ## Calculating dbscan combinaisons
            dbscanResults <- .dbscanComb(tSNEList, cores, epsilon, minPoints,
                    sceObject)

            ## Output file
            if(writeOutput)
                .writeDBScanResults(theObject, dbscanResults)

            ## Creation of a list of Dbscan objects
            dbscanList <- .createDbscanList(tSNEList, minPoints, epsilon,
                    dbscanResults)

            setDbscanList(theObject) <- dbscanList
            
            ## Get a suggested clusters number to use
            k <- .retrieveClustersNumberK(dbscanList)
            setSuggestedClustersNumber(theObject) <- k
            
            return(theObject)
        })
