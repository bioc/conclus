##################
## testClustering
##################

#' .plotDistanceGraphWithEpsilon
#'
#' @description
#' Generate a diagnostic plot to determine the epsilon of dbscan.
#'
#' @param tSNEData Data column of the object tSNE obtained with .getTSNEresults.
#' @param minNeighbours MinPoints parameter of ?testClustering.
#' @param epsilon dbscanEpsillon parameter of ?testClustering.
#'
#' @keywords internal
#'
#' @importFrom dbscan kNNdistplot
#' @importFrom graphics abline
#' @noRd
.plotDistanceGraphWithEpsilon <- function(tSNEData, minNeighbours=5,
        epsilon=1.2){

    dbscan::kNNdistplot(tSNEData, k=minNeighbours)
    abline(h=epsilon, lty=2)
}


#' .plotTestClustering
#'
#' @description
#' Plot the result of the dbscan classification.
#'
#' @param tSNEData Data column of the object tSNE obtained with .getTSNEresults.
#' @param minNeighbours MinPoints parameter of ?testClustering.
#' @param epsilon dbscanEpsillon parameter of ?testClustering.
#'
#' @keywords internal
#'
#' @importFrom fpc dbscan
#' @importFrom factoextra fviz_cluster
#' @return A ggplot object
#' @noRd
.plotTestClustering <- function(tSNEData, minNeighbours=5,epsilon=1.2){

    dbscanList <- fpc::dbscan(tSNEData, eps=epsilon, MinPts=minNeighbours)
    p <- factoextra::fviz_cluster(dbscanList, tSNEData, ellipse=TRUE,
            geom="point", legend="bottom")
    return(p)
}



.checkParamsTests <- function(sceObject, dbscanEpsilon, minPts, perplexities,
        PCs, randomSeed){


    if(all(dim(sceObject) == c(0,0)))
        stop("The 'scRNAseq' object that you're using with ",
                "'testClustering' method doesn't have its ",
                "'sceNorm' slot updated. Please use 'normaliseCountMatrix'",
                " on the object before.")

    ## Check parameters
    if(!is.numeric(dbscanEpsilon))
        stop("'dbscanEpsilon' parameter should be an integer.")

    ## Check minPts argument
    if(!is.numeric(minPts))
        stop("'minPts' parameter should be an integer")

    ## Check perplexities argument
    if(!is.numeric(perplexities))
        stop("'perplexities' parameter should be a vector of numeric.")

    ## Check PCs argument
    if(!is.numeric(PCs))
        stop("'PCs' parameter should be a vector of numeric.")

    ## Check randomSeed argument
    if(!is.numeric(randomSeed))
        stop("'randomSeed' parameter should be an integer.")

    ## Check number of elements
    if(!isTRUE(all.equal(length(dbscanEpsilon),1)) ||
            !isTRUE(all.equal(length(minPts),1)) ||
            !isTRUE(all.equal(length(perplexities), 1)) ||
            !isTRUE(all.equal(length(PCs), 1)) ||
            !isTRUE(all.equal(length(randomSeed), 1)))
        stop("dbscanEpsilon, minPts, perplexities, PCs, and randomSeed ",
                "should be a single value.")
}

#' .printTSNE
#'
#' @description
#' Generate the tSNE plot.
#'
#' @param writeOutput If TRUE, write the results of the test to the output
#' directory defined in theObject in the sub-directory 'test_clustering'.
#' Default = FALSE.
#' @param dataDirectory Path to the output folder in which the pdf files will be
#' written.
#' @param width Width of the pdf file. Default=7. See ?pdf for details.
#' @param height Height of the pdf file. Default=7. See ?pdf for details.
#' @param tSNE Result of the function .getTSNEresults.
#' @param fileTSNE Name of the pdf file. Default="test_tSNE.pdf"
#' @param silent If TRUE, do not plot graphics. Default=FALSE.
#' @param ... Options for generating the pdf files. See ?pdf for a list.
#'
#' @keywords internal
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @noRd

.printTSNE <- function(writeOutput, dataDirectory, width, height, tSNE,
        fileTSNE, silent, ...){

    if(writeOutput){
        message("Saving results tSNE.")
        pdf(file.path(dataDirectory, "test_clustering", fileTSNE),
                width=width, height=height, ...)
        print(tSNE)
        dev.off()
    }
    
    if(!silent)
        print(tSNE)
    
    return(tSNE)
}


#' .printDist
#'
#' @description
#' Generate a diagnostic plot to determine the optimal epsilon parameter.
#'
#' @param writeOutput If TRUE, write the results of the test to the output
#' directory defined in theObject in the sub-directory 'test_clustering'.
#' Default = FALSE.
#' @param dataDirectory Path to the output folder in which the pdf files will be
#' written.
#' @param width Width of the pdf file. Default=7. See ?pdf for details.
#' @param height Height of the pdf file. Default=7. See ?pdf for details.
#' @param tSNE Result of the function .getTSNEresults.
#' @param dbscanEpsilon Single value for the distance parameter of dbscan.
#' Default = 1.4. See ?runDBSCAN for more details.
#' @param minPoints Single value for the minimum no. of points parameter of
#' dbscan. Default = 5. See ?runDBSCAN for more details.
#' @param fileDist Name of the pdf file. Default="distance_graph.pdf"
#' @param plotKNN If TRUE, output the kNN plot on graphics. Default=TRUE.
#' @param ... Options for generating the pdf files. See ?pdf for a list.
#'
#' @details
#' See ?testClustering for more details.
#'
#' @keywords internal
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom grDevices dev.new
#' @noRd

.printDist <- function(writeOutput, dataDirectory, width, height, tSNE,
        dbscanEpsilon, minPts, fileDist, plotKNN, ...){

    if(writeOutput){

        message("Saving results distance graph.")
        pdf(file.path(dataDirectory, "test_clustering", fileDist),
                width=width, height=height, ...)
        .plotDistanceGraphWithEpsilon(tSNE$data, epsilon=dbscanEpsilon,
                minNeighbours=minPts)
        dev.off()
    }
    
    if(plotKNN){

        dev.new()
        .plotDistanceGraphWithEpsilon(tSNE$data, epsilon=dbscanEpsilon,
                minNeighbours=minPts)
    }

}

#' .printDBScan
#'
#' @description
#' Generate the dbscan clustering on the tSNE plot.
#'
#' @param writeOutput If TRUE, write the results of the test to the output
#' directory defined in theObject in the sub-directory 'test_clustering'.
#' Default = FALSE.
#' @param dataDirectory Path to the output folder in which the pdf files will be
#' written.
#' @param width Width of the pdf file. Default=7. See ?pdf for details.
#' @param height Height of the pdf file. Default=7. See ?pdf for details.
#' @param tSNE Result of the function .getTSNEresults.
#' @param dbscanEpsilon Single value for the distance parameter of dbscan.
#' Default = 1.4. See ?runDBSCAN for more details.
#' @param minPoints Single value for the minimum no. of points parameter of
#' dbscan. Default = 5. See ?runDBSCAN for more details.
#' @param fileClust Name of the pdf file. Default="test_clustering.pdf"
#' @param silent If TRUE, do not plot graphics. Default=FALSE.
#' @param ... Options for generating the pdf files. See ?pdf for a list.
#'
#' @keywords internal
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom grDevices dev.new
#'
#' @noRd

.printDBScan <- function(writeOutput, dataDirectory, width, height, tSNE,
        dbscanEpsilon, minPts, fileClust, silent, ...){

    p <- .plotTestClustering(tSNE$data, epsilon=dbscanEpsilon,
            minNeighbours=minPts)

    if(writeOutput){

        message("Saving dbscan results.")
        pdf(file.path(dataDirectory, "test_clustering", fileClust),
                width=width, height=height, ...)
        print(p)
        dev.off()
    }
    
    if(!silent){
        dev.new()
        print(p)
    }
    
    return(p)
}


#' testClustering
#'
#' @description
#' This function generates a single clustering iteration of CONCLUS to check
#' whether the chosen parameters of tSNE and dbscan are suitable for your data.
#'
#' @usage
#' testClustering(theObject, dbscanEpsilon=1.4, minPts=5,
#'                 perplexities=30, PCs=4, randomSeed=42, width=7, height=7,
#'                 cores=2, writeOutput=FALSE, fileTSNE="test_tSNE.pdf",
#'                 fileDist="distance_graph.pdf",
#'                 fileClust="test_clustering.pdf", silent=FALSE, plotKNN=TRUE, 
#'                  ...)
#'
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized. See ?normaliseCountMatrix.
#' @param dbscanEpsilon Single value for the distance parameter of dbscan.
#' Default = 1.4. See ?runDBSCAN for more details.
#' @param minPts Single value for the minimum no. of points parameter of
#' dbscan. Default = 5. See ?runDBSCAN for more details.
#' @param perplexities A single value of perplexity (t-SNE parameter).
#' Default = 30. See ?generateTSNECoordinates for details.
#' @param PCs Single value of first principal components. Default=4. See
#' ?generateTSNECoordinates for details.
#' @param randomSeed  Default is 42. Seeds used to generate the tSNE.
#' @param width Width of the pdf file. Default=7. See ?pdf for details.
#' @param height Height of the pdf file. Default=7. See ?pdf for details.
#' @param cores Maximum number of jobs that CONCLUS can run in parallel.
#' Default is 1.
#' @param writeOutput If TRUE, write the results of the test to the output
#' directory defined in theObject in the sub-directory 'test_clustering'.
#' Default = FALSE.
#' @param fileTSNE Name of the pdf file for tSNE. Default="test_tSNE.pdf".
#' @param fileDist Name of the pdf file for NN distance.
#' Default="distance_graph.pdf"
#' @param fileClust Name of the pdf file for dbscan.
#' Default="test_clustering.pdf"
#' @param silent If TRUE, do not plot graphics. Default=FALSE.
#' @param plotKNN If TRUE, output the kNN plot on graphics. Default=TRUE.
#' @param ... Options for generating the pdf files. See ?pdf for a list.
#' 
#' @return A ggplot object of the tSNE and the dbscan clustering.
#'
#' @aliases testClustering
#'
#' @details
#' The TestClustering function runs one clustering round out of the 84 (default)
#' rounds that CONCLUS normally performs. This step can be useful to determine
#' if the default DBSCAN parameters are suitable for your dataset. By default,
#' they are dbscanEpsilon = c(1.3, 1.4, 1.5) and minPts = c(3,4). If the dashed
#' horizontal line in the k-NN distance plot lays on the "knee" of the curve,
#' it means that optimal epsilon is equal to the intersection of the line to the
#' y-axis. In our example, optimal epsilon is 1.4 for 5-NN distance where 5
#' corresponds to MinPts.
#'
#' In the "test_clustering" folder under outputDirectory, the three plots will
#' be saved where one corresponds to the "distance_graph.pdf", another one to
#' "test_tSNE.pdf", and the last one will be saved as "test_clustering.pdf".
#'
#' @author
#' Ilyess RACHEDI, based on code by Konstantin CHUKREV and Nicolas DESCOSTES.
#'
#' @rdname
#' testClustering-scRNAseq
#'
#' @examples
#' ## Load the count matrix
#' countmatrixPath <- system.file("extdata/test_countMatrix.tsv", 
#'                             package="conclus")
#' countMatrix <- loadCountMatrix(file=countmatrixPath, header=TRUE, dec=".",
#'                                 sep='\t')
#' 
#' ## Load the coldata
#' coldataPath <- system.file("extdata/test_colData_filtered.tsv", 
#'                             package="conclus")
#' columnsMetaData <- loadColdata(file=coldataPath, columnCell="cell_ID",
#'                                 header=TRUE, dec=".", sep='\t')
#' 
#' ## Create the initial object
#' scr <- singlecellRNAseq(experimentName = "Bergiers",
#'                 countMatrix     = countMatrix,
#'                 species         = "mouse",
#'                 outputDirectory = "YourOutputDirectory")
#'
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#'
#' ## Test the clustering with default parameters (dbscanEpsilon=1.4, minPts=5)
#' testClustering(scrNorm)
#'
#' ## Test the clustering writing pdfs to test_clustering folder
#' testClustering(scrNorm, writeOutput=TRUE)
#'
#' ## Removing the written results
#' unlink("YourOutputDirectory/", recursive = TRUE)
#'
#' @seealso
#' normaliseCountMatrix runDBSCAN pdf
#'
#' @exportMethod testClustering
#' @importFrom Biobase exprs
#' @importFrom methods validObject

setMethod(

        f="testClustering",

        signature="scRNAseq",

        definition=function(theObject, dbscanEpsilon=1.4, minPts=5,
                perplexities=30, PCs=4, randomSeed=42, width=7, height=7,
                cores=2, writeOutput=FALSE, fileTSNE="test_tSNE.pdf",
                fileDist="distance_graph.pdf", fileClust="test_clustering.pdf",
                silent=FALSE, plotKNN=TRUE, ...){
            
            if(!writeOutput && silent)
                warning("writeOutput=FALSE and silent=TRUE, nothing will ",
                        "happen.", immediate. =TRUE)
            
            validObject(theObject)

            ## Get the values from slots
            sceObject <- getSceNorm(theObject)
            dataDirectory <- getOutputDirectory(theObject)
            experimentName <- getExperimentName(theObject)

            ## Checking parameters
            .checkParamsTests(sceObject, dbscanEpsilon, minPts, perplexities,
                    PCs, randomSeed)

            message("Generating TSNE.")

            if(writeOutput){
                
                initialisePath(dataDirectory)
                dir.create(file.path(dataDirectory, "test_clustering"),
                        showWarnings=FALSE)
            }
            
            p <- list()
            
            ## 1. Generating 2D tSNE plots
            tSNE <- .getTSNEresults(expressionMatrix=Biobase::exprs(sceObject),
                    cores=cores, PCs=PCs, perplexities=perplexities,
                    randomSeed=randomSeed)
            p[[1]] <- .printTSNE(writeOutput, dataDirectory, width, height, 
                    tSNE, fileTSNE, silent, ...)

            ## 2. Clustering with dbscan
            .printDist(writeOutput, dataDirectory, width, height, 
                    tSNE, dbscanEpsilon, minPts, fileDist, plotKNN, ...)
            p[[2]] <- .printDBScan(writeOutput, dataDirectory, width, height, 
                    tSNE, dbscanEpsilon, minPts, fileClust, silent, ...)
            
            return(p)
        })



##################
## clusterCellsInternal
##################


#' .computeSimMat
#'
#' @description
#' For each cluster number (clusters), computes a similarity matrix with 1. In
#' the resulting list (l), there is one matrix for each cluster number.
#'
#' @param clusters Unique list of clusters names.
#' @param simMat Empty matrix containing only 0 that will be filled by the
#' function.
#' @param clustering dbscan clustering results.
#'
#' @keywords internal
#' @noRd
#' @return A list of similarity matrix.
.computeSimMat <- function(clusters, simMat, clustering){
    
    l <- lapply(clusters[clusters!=0],
            function(cluster, simMat, clustering){
                
                ## Get the cells in the same cluster
                selCol <- colnames(clustering)[
                        clustering == cluster]
                ## Add +1 for each cell seen in the same cluster
                simMat[rownames(simMat) %in% selCol,
                        colnames(simMat) %in% selCol] <- 1
                return(simMat)
            }, simMat, clustering)
    return(l)
}


#' .mkSimMat
#'
#' @description
#' This function calculates how many time a pair of cells were assigned to
#' a cluster by dbscan.
#'
#' @param dbscanList List of dbscan results given by ?getDbscanList.
#' @param cores Maximum number of jobs that CONCLUS can run in parallel.
#' Default is 1.
#'
#' @keywords internal
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom utils globalVariables
#' @return A cell similarity matrix
#' @noRd
.mkSimMat <- function(dbscanList, cores){

    ## This foreach gives, for each dbscan result, a matrix containing 1 and
    ## 0. 1 means that a pair of cells was allocated to one cluster,
    ## 0 means that the pair of cells was allocated to no cluster.

    myCluster <- parallel::makeCluster(cores, type="PSOCK")
    doParallel::registerDoParallel(myCluster)
    simMatsList <- foreach::foreach(iMkSimMat=seq_len(length(dbscanList)),
                    .export=c("getClustering", ".computeSimMat")) %dopar% {

                ## Get the clustering
                clustering <- getClustering(dbscanList[[iMkSimMat]])

                ## Get the number of cells in the current clustering
                nrow <- ncol(clustering)
                ncol <- ncol(clustering)

                ## Create a square matrix
                simMat <- matrix(0, ncol=ncol, nrow=nrow)
                colnames(simMat) <- colnames(clustering)
                rownames(simMat) <- colnames(clustering)

                ## Get cluster assignments
                clusters <- unique(as.vector(clustering))
                l <- .computeSimMat(clusters, simMat, clustering)

                ## Sum all the matrix to have a single one showing cells
                ## that were allocated to one cluster (1) and cells that were
                ## not allocated to one cluster (0).
                simMat <- Reduce("+", l)}
    parallel::stopCluster(myCluster)

    ## Sum the similiratity matrices to obtain a single one giving howw many
    ## times a pair of cells was clustered together accross all the dbscan
    ## parameters
    simMat <- Reduce('+', simMatsList)
    simMat <- simMat / length(dbscanList)
    stopifnot(isSymmetric(simMat))

    return(simMat)
}

.checkClusteringMethod <- function(clusteringMethod){

    clusteringMethods <- c("ward.D", "ward.D2", "single", "complete", "average",
            "mcquitty", "median", "centroid")
    availableMethods <- paste(clusteringMethods, collapse="; ")
    if(!clusteringMethod %in% clusteringMethods)
        stop("'clusteringMethod' should be one of: ", availableMethods)

}


.checkParamsInternal <- function(sceObject, dbscanList, clusterNumber,
        deepSplit, cores, clusteringMethod){

    ## Check if the normalized matrix is correct
    if(all(dim(sceObject) == c(0,0)))
        stop("The 'scRNAseq' object that you're using with ",
                "'clusterCellsInternal' function doesn't have its 'sceNorm' ",
                "slot updated. Please use 'normaliseCountMatrix' on the object",
                " before.")

    ## Check if the dbscan list is correct
    if (length(dbscanList) <= 1)
        stop("The 'scRNAseq' object that you're using with ",
                "'clusterCellsInternal' function doesn't have its 'dbscanList'",
                " slot updated. Please use 'runDBSCAN' on the object before.")

    ## Check the cluster number
    if(!is.numeric(clusterNumber))
        stop("'clusterNumber' parameter should be a numeric.")

    ## Check the deepSplit
    if(!is.numeric(deepSplit))
        stop("'deepSplit' parameter should be a numeric.")

    ## Check the cores
    if(!is.numeric(cores))
        stop("'cores' parameter should be a numeric.")

    ## Check the clustering method
    .checkClusteringMethod(clusteringMethod)

}


.reinitializeObject <- function(theObject){
    
    setClustersSimilarityMatrix(theObject) <-  matrix(nrow = 1, ncol = 1, 
            dimnames = list("1", "1"), data = 1)
    setMarkerGenesList(theObject) <- list(data.frame(Gene = c("gene1"), 
                    mean_log10_fdr = c(NA), n_05 = c(NA), score = c(NA)))
    setClustersMarkers(theObject) <- data.frame(geneName="gene1", clusters=NA)
    setGenesInfos(theObject) <- data.frame(uniprot_gn_symbol=c("symbol"), 
            clusters="1", external_gene_name="gene", go_id="GO1,GO2", 
            mgi_description="description", entrezgene_description="descr",
            gene_biotype="gene", chromosome_name="1", Symbol="symbol",
            ensembl_gene_id="ENS", mgi_id="MGI", entrezgene_id="1",
            uniprot_gn_id="ID")
    setClustersSimiliratyOrdered(theObject) <- factor(1)
    
    return(theObject)
}


#' clusterCellsInternal
#'
#' @description
#' Returns consensus clusters by using hierarchical clustering on the similarity
#' matrix of cells.
#'
#' @usage
#' clusterCellsInternal(theObject, clusterNumber=0, deepSplit=4, cores=2,
#'                 clusteringMethod="ward.D2")
#'
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized (see ?normaliseCountMatrix), tSNE were calculated (see
#' ?generateTSNECoordinates), and dbScan was run (see ?runDBSCAN),
#' @param clusterNumber Exact number of cluster. Default = 0 that will determine
#' the number of clusters automatically.
#' @param deepSplit Intuitive level of clustering depth. Options are 1, 2, 3, 4.
#' Default = 4
#' @param cores Maximum number of jobs that CONCLUS can run in parallel.
#' Default is 1.
#' @param clusteringMethod Clustering method passed to hclust() function. See
#' ?hclust for a list of method. Default = "ward.D2".
#'
#' @aliases clusterCellsInternal
#'
#' @author
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
#'
#' @rdname
#' clusterCellsInternal-scRNAseq
#'
#' @return
#' An object of class scRNASeq with its cellsSimilarityMatrix and sceNorm slots
#' updated.
#'
#' @examples
#' ## Load the count matrix
#' countmatrixPath <- system.file("extdata/test_countMatrix.tsv", 
#'                             package="conclus")
#' countMatrix <- loadCountMatrix(file=countmatrixPath, header=TRUE, dec=".",
#'                                 sep='\t')
#' 
#' ## Load the coldata
#' coldataPath <- system.file("extdata/test_colData_filtered.tsv", 
#'                             package="conclus")
#' columnsMetaData <- loadColdata(file=coldataPath, columnCell="cell_ID",
#'                                 header=TRUE, dec=".", sep='\t')
#' 
#' ## Create the initial object
#' scr <- singlecellRNAseq(experimentName = "Bergiers",
#'                 countMatrix     = countMatrix,
#'                 species         = "mouse",
#'                 outputDirectory = "YourOutputDirectory")
#'
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#'
#' ## Compute the tSNE coordinates
#' scrTsne <- generateTSNECoordinates(scrNorm, cores=2)
#'
#' ## Perform the clustering with dbScan
#' scrDbscan <- runDBSCAN(scrTsne, cores=2)
#'
#' ## Compute the cell similarity matrix
#' scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=10, cores=2)
#'
#' @seealso plotCellSimilarity
#'
#' @exportMethod clusterCellsInternal
#' @importFrom SummarizedExperiment colData
#' @importFrom stats hclust
#' @importFrom stats as.dist
#' @importFrom stats cutree
#' @importFrom utils capture.output
#' @importFrom methods validObject

setMethod(

        f = "clusterCellsInternal",

        signature = "scRNAseq",

        definition = function(theObject, clusterNumber=0, deepSplit=4, cores=2,
                clusteringMethod="ward.D2"){

            ## Check if the Object is valid
            validObject(theObject)
            .reinitializeObject(theObject)

            sceObject  <- getSceNorm(theObject)
            dbscanList <- getDbscanList(theObject)

            .checkParamsInternal(sceObject, dbscanList, clusterNumber,
                    deepSplit, cores, clusteringMethod)

            message("Calculating cells similarity matrix.")
            cellsSimilarityMatrix <- .mkSimMat(dbscanList, cores=cores)


            distMatrix <- as.dist(sqrt((1-cellsSimilarityMatrix)/2))
            clusteringTree <- hclust(distMatrix, method=clusteringMethod)

            if(clusterNumber == 0){
                message(paste0("Assigning cells to clusters. DeepSplit = ",
                                deepSplit))
                clusters <- unname(dynamicTreeCut::cutreeDynamic(clusteringTree,
                                distM=as.matrix(distMatrix),
                                verbose=0,
                                deepSplit=deepSplit))
            } else {
                message(paste0("Assigning cells to ", clusterNumber,
                                " clusters."))
                clusters <- cutree(clusteringTree, k=clusterNumber)
            }

            SummarizedExperiment::colData(sceObject)$clusters <-
                    factor(clusters)

            setCellsSimilarityMatrix(theObject) <- cellsSimilarityMatrix
            setSceNorm(theObject) <- sceObject
            message(paste0(capture.output(
                table(SummarizedExperiment::colData(sceObject)$clusters,
                dnn=list("Cells distribution by clusters: "))),collapse = "\n"))

            return(theObject)
        })




##################
## calculateClustersSimilarity
##################


#' .computeClusMat
#'
#' @description
#' This function returns a matrix with "protocells" representing clusters.
#' Values show how much two "protocells" are similar. 1 if clusters are
#' perfectly similar, 0 if totally different. It is based on the median.
#'
#' @param clustersNames Names of the different clusters typically of the form
#' clusterX with X being the cluster number.
#' @param mat Matrix containing the cells similarity or the clusters median.
#' @param clusters The different clusters for each cell given by
#' SummarizedExperiment::colData(sceObject)$clusters
#'
#' @keywords internal
#' @return A median matrix.
#' @noRd

.computeClusMat <- function(clustersNames, mat, clusters){

    resultMedList <- lapply(clustersNames, function(currentClustName,
                    fullmat, clusts){
                return(matrixStats::rowMedians(as.matrix(fullmat[, clusts ==
                                                currentClustName])))
            }, mat, clusters)
    medMat <- do.call(cbind, resultMedList)
    return(medMat)
}

#' .mkSimMed
#'
#' @description
#' Generate the clusters similarity matrix.
#'
#' @param simMat The cellsSimilarityMatrix.
#' @param clusters The different clusters for each cell given by
#' SummarizedExperiment::colData(sceObject)$clusters
#' @param clustersNames Names of the different clusters typically of the form
#' clusterX with X being the cluster number.
#'
#' @keywords internal
#' @return A similarity matrix
#' @noRd

.mkSimMed <- function(simMat, clusters, clustersNames){

    clusMed <- matrix(ncol=length(unique(clusters)), nrow=nrow(simMat))
    clusMed <- .computeClusMat(clustersNames, simMat, clusters)
    clusMed <- t(clusMed)

    simMed <- matrix(ncol=length(unique(clusters)),
            nrow=length(unique(clusters)))
    simMed <- .computeClusMat(clustersNames, clusMed, clusters)
    colnames(simMed) <- clustersNames
    rownames(simMed) <- clustersNames

    return(simMed)
}

.checkParamsClusterSim <- function(clusteringMethod, cellsSimilarityMatrix,
        sceObject){

    ## Check the clustering method
    .checkClusteringMethod(clusteringMethod)

    ## Check if the normalized matrix is correct
    if(all(dim(sceObject) == c(0,0)))
        stop("The 'scRNAseq' object that you're using with ",
                "'calculateClustersSimilarity' function doesn't have its ",
                "'sceNorm' slot updated. Please use 'normaliseCountMatrix' on",
                " the object before.")

    ## Check if the normalized matrix is filtered
    if(!("clusters" %in% names(colData(sceObject))))
        stop("The 'scRNAseq' object that you're using with ",
                "'calculateClustersSimilarity' function doesn't have a correct",
                " 'sceNorm' slot. This slot should be a 'SingleCellExperiment'",
                " object containing 'clusters' column in its colData. Please ",
                "check if you correctly used 'clusterCellsInternal' on the ",
                "object.")

    ## Check the cell similarity matrix
    if(all(dim(cellsSimilarityMatrix) == c(1,1)))
        stop("The 'scRNAseq' object that you're using with ",
                "'calculateClustersSimilarity' function doesn't have its ",
                "'cellsSimilarityMatrix' slot updated by clusterCellsInternal.",
                " Please use 'clusterCellsInternal' on the object before.")
}


#' calculateClustersSimilarity
#'
#' @description
#' Having computed cells similarity, pools information into clusters.
#'
#' @usage
#' calculateClustersSimilarity(theObject, clusteringMethod = "ward.D2")
#'
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized (see ?normaliseCountMatrix), tSNE were calculated (see
#' ?generateTSNECoordinates), dbScan was run (see ?runDBSCAN), and cells were
#' clustered (see ?clusterCellsInternal).
#' @param clusteringMethod Clustering method passed to hclust() function. See
#' ?hclust for a list of method. Default = "ward.D2".
#'
#' @aliases calculateClustersSimilarity
#'
#' @author
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
#'
#' @rdname
#' calculateClustersSimilarity-scRNAseq
#'
#' @return
#' An object of class scRNASeq with its clustersSimilarityMatrix and
#' clustersSimiliratyOrdered slots updated.
#'
#' @examples
#' ## Load the count matrix
#' countmatrixPath <- system.file("extdata/test_countMatrix.tsv", 
#'                             package="conclus")
#' countMatrix <- loadCountMatrix(file=countmatrixPath, header=TRUE, dec=".",
#'                                 sep='\t')
#' 
#' ## Load the coldata
#' coldataPath <- system.file("extdata/test_colData_filtered.tsv", 
#'                             package="conclus")
#' columnsMetaData <- loadColdata(file=coldataPath, columnCell="cell_ID",
#'                                 header=TRUE, dec=".", sep='\t')
#' 
#' ## Create the initial object
#' scr <- singlecellRNAseq(experimentName = "Bergiers",
#'                 countMatrix     = countMatrix,
#'                 species         = "mouse",
#'                 outputDirectory = "YourOutputDirectory")
#'
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#'
#' ## Compute the tSNE coordinates
#' scrTsne <- generateTSNECoordinates(scrNorm, cores=2)
#'
#' ## Perform the clustering with dbScan
#' scrDbscan <- runDBSCAN(scrTsne, cores=2)
#'
#' ## Compute the cell similarity matrix
#' scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=10, cores=2)
#'
#' ## Calculate clusters similarity
#' scrCSM <- calculateClustersSimilarity(scrCCI)
#'
#' @seealso
#' plotClustersSimilarity
#'
#' @exportMethod calculateClustersSimilarity
#' @importFrom SummarizedExperiment colData
#' @importFrom stats as.dist
#' @importFrom stats hclust
#' @importFrom methods validObject

setMethod(

        f = "calculateClustersSimilarity",

        signature = "scRNAseq",

        definition = function(theObject, clusteringMethod = "ward.D2"){

            ## Check if the Object is valid
            validObject(theObject)

            cellsSimilarityMatrix <- getCellsSimilarityMatrix(theObject)
            sceObject  <- getSceNorm(theObject)

            .checkParamsClusterSim(clusteringMethod, cellsSimilarityMatrix,
                    sceObject)

            clusters <- SummarizedExperiment::colData(sceObject)$clusters
            clustersNumber <- length(unique(clusters))
            clustersNames <- levels(clusters)

            ## Calculating cluster similarity.
            clustersSimilarityMatrix <- .mkSimMed(
                    as.matrix(cellsSimilarityMatrix), clusters, clustersNames)

            ## Plotting matrix
            distMatrix <- as.dist(sqrt((1 - clustersSimilarityMatrix)/2))
            clusteringTree <- hclust(distMatrix, method=clusteringMethod)

            clustersSimOrdered <- data.frame(clusterNames=clustersNames,
                    clusterIndexes=seq_len(clustersNumber))
            rownames(clustersSimOrdered) <- clustersSimOrdered$clusterIndexes
            clustersSimOrdered <- clustersSimOrdered[clusteringTree$order, ]

            setClustersSimilarityMatrix(theObject)  <- clustersSimilarityMatrix
            setClustersSimiliratyOrdered(theObject) <-
                    as.factor(clustersSimOrdered$clusterNames)
            return(theObject)
        })



##################
## retrieveTableClustersCells
##################
        


#' retrieveTableClustersCells
#'
#' @description
#' Having computed clusterCellsInternal, retrieve to what cluster each cell 
#' belongs. The output data.frame can be passed to the method ?addClustering.
#'
#' @usage
#' retrieveTableClustersCells(theObject)
#'
#' @param theObject An Object of class scRNASeq for which the cells were 
#' clustered internally. See ?clusterCellsInternal.
#'
#' @aliases retrieveTableClustersCells
#'
#' @author
#' Nicolas DESCOSTES.
#'
#' @rdname
#' retrieveTableClustersCells-scRNAseq
#'
#' @return
#' A data frame containing two columns 'clusters' and 'cells' indicating the 
#' result of the consensus clustering at the cellular level.
#'
#' @examples
#' ## Load the count matrix
#' countmatrixPath <- system.file("extdata/test_countMatrix.tsv", 
#'                             package="conclus")
#' countMatrix <- loadCountMatrix(file=countmatrixPath, header=TRUE, dec=".",
#'                                 sep='\t')
#' 
#' ## Load the coldata
#' coldataPath <- system.file("extdata/test_colData_filtered.tsv", 
#'                             package="conclus")
#' columnsMetaData <- loadColdata(file=coldataPath, columnCell="cell_ID",
#'                                 header=TRUE, dec=".", sep='\t')
#' 
#' ## Create the initial object
#' scr <- singlecellRNAseq(experimentName = "Bergiers",
#'                 countMatrix     = countMatrix,
#'                 species         = "mouse",
#'                 outputDirectory = "YourOutputDirectory")
#'
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#'
#' ## Compute the tSNE coordinates
#' scrTsne <- generateTSNECoordinates(scrNorm, cores=2)
#'
#' ## Perform the clustering with dbScan
#' scrDbscan <- runDBSCAN(scrTsne, cores=2)
#'
#' ## Compute the cell similarity matrix
#' scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=10, cores=2)
#'
#' ## Retrieving the table clusters-cells.
#' cellClustDf <- retrieveTableClustersCells(scrCCI)
#'
#' @seealso
#' addClustering
#'
#' @exportMethod retrieveTableClustersCells
#' @importFrom SummarizedExperiment colData

setMethod(
        
        f = "retrieveTableClustersCells",
        
        signature = "scRNAseq",
        
        definition = function(theObject){
            
            ## Check if the Object is valid
            validObject(theObject)
            
            ## Retrieve the clustering result in theObject
            sceObject  <- getSceNorm(theObject)
            colDf <- SummarizedExperiment::colData(sceObject)
            
            if(!any(names(colDf) == "clusters"))
                stop("clusterCellsInternal should be performed before ",
                        "retrieving this information.")
            
            colDf <- data.frame(clusters=colDf$clusters, cells=colDf$cellName)
            
            return(colDf)
        })

        
        
##################
## addClustering
##################

.checkParamsAddClustering <- function(theObject, clusToAdd, colDf){

    ## Check that the cluster similarity matrix was computed
    ## This step is necessary to  update the 'clusters' column
    ## of the col metadata.
    matrix <- getClustersSimilarityMatrix(theObject)

    if(all(dim(matrix) == c(1,1)))
        stop("The 'scRNAseq' object that you're using with 'exportResults' ",
                "function doesn't have its columns metadata ",
                "updated. Please use 'calculateClustersSimilarity' on the ",
                "object before.")

    if(!isTRUE(all.equal(ncol(clusToAdd), 2)))
        stop("The file given to filePathAdd  should contain two columns ",
                "'clusters' and 'cells'. Instead it contains: ",
                paste(colnames(clusToAdd), collapse="-"))

    if(!all(colnames(clusToAdd) %in% c("clusters", "cells")))
        stop("The file given to filePathAdd  should contain two columns ",
                "'clusters' and 'cells' Instead it contains: ",
                paste(colnames(clusToAdd), collapse="-"))

    if(!all(colDf$cells %in% clusToAdd$cells))
        stop("The cells column in theObject clustering results contains cells ",
                "names that are not the same then the ones of the cluster to ",
                "add. Make sure that the cells names of the cluster to add ",
                " are the same.")
}


#' addClustering
#'
#' @description
#' This method enables to add a clustering to the existing object in order to
#' change the coloration of the t-sne. It is particularly useful to compare the
#' performance of different tools.
#'
#' @usage
#' addClustering(theObject, filePathAdd=NA, clusToAdd=NA)
#'
#' @param theObject An Object of class scRNASeq for which
#' ?calculateClustersSimilarity was used.
#' @param filePathAdd Default=NA. Path to the file containing the clustering to
#' replace in the object. It should be made of two columns 'clusters' and
#' 'cells'. This should be left to NA if clusToAdd is used. Default=NA.
#' @param clusToAdd Data frame having two columns 'clusters' and 'cells'
#' containing the clustering to be used in theObject. Should be left to NA if
#' filePathAdd is used. Default=NA.
#'
#' @aliases addClustering
#'
#' @author
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
#'
#' @rdname
#' addClustering-scRNAseq
#'
#' @return
#' An object of class scRNASeq with its column name metadata updated.
#'
#' @examples
#' ## Load the count matrix
#' countmatrixPath <- system.file("extdata/test_countMatrix.tsv", 
#'                             package="conclus")
#' countMatrix <- loadCountMatrix(file=countmatrixPath, header=TRUE, dec=".",
#'                                 sep='\t')
#' 
#' ## Load the coldata
#' coldataPath <- system.file("extdata/test_colData_filtered.tsv", 
#'                             package="conclus")
#' columnsMetaData <- loadColdata(file=coldataPath, columnCell="cell_ID",
#'                                 header=TRUE, dec=".", sep='\t')
#' 
#' ## Load the clusters table with the new cluster to add
#' clustAddTab <- read.delim(
#' system.file("extdata/Bergiers_clusters_table.tsv", package="conclus"))
#'
#' ## Create the initial object
#' scr <- singlecellRNAseq(experimentName = "Bergiers",
#'                 countMatrix     = countMatrix,
#'                 species         = "mouse",
#'                 outputDirectory = "YourOutputDirectory")
#'
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#'
#' ## Compute the tSNE coordinates
#' scrTsne <- generateTSNECoordinates(scrNorm, cores=2)
#'
#' ## Perform the clustering with dbScan
#' scrDbscan <- runDBSCAN(scrTsne, cores=2)
#'
#' ## Compute the cell similarity matrix
#' scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=10, cores=2)
#'
#' ## Calculate clusters similarity
#' scrCSM <- calculateClustersSimilarity(scrCCI)
#'
#' ## Ranking genes
#' scrS4MG <- rankGenes(scrCSM)
#' 
#' ## Getting marker genes
#' scrFinal <- retrieveTopClustersMarkers(scrS4MG, removeDuplicates = FALSE)
#'
#' ## Getting genes info
#' scrInfos <- retrieveGenesInfo(scrFinal, cores=2)
#' 
#' ## Retrieving the table indicating to which cluster each cell belongs
#' clustCellsDf <- retrieveTableClustersCells(scrInfos)
#' 
#' ## Replace “10” by “9” to merge 9/10
#' clustCellsDf$clusters[which(clustCellsDf$clusters == 10)] <- 9
#' 
#' ## Modifying the object to take into account the new classification
#' scrUpdated <- addClustering(scrInfos, clusToAdd=clustCellsDf)
#' 
#'
#' @exportMethod addClustering
#' @importFrom SummarizedExperiment colData
#' @importFrom utils read.table
#' @importFrom methods validObject

setMethod(

        f = "addClustering",

        signature = "scRNAseq",

        definition = function(theObject, filePathAdd=NA, clusToAdd=NA){

            ## Check if the Object is valid
            validObject(theObject)

            ## Retrieve the clustering to add
            if(isTRUE(all.equal(length(clusToAdd), 1)) && is.na(clusToAdd))
                if(!is.na(filePathAdd))
                    clusToAdd <- read.table(filePathAdd, header=TRUE, sep="\t")
                else
                    stop("Either filePathAdd or clustToAdd should be given.")

            ## Retrieve the clustering result in theObject
            sceObject  <- getSceNorm(theObject)
            colDf <- SummarizedExperiment::colData(sceObject)
            colDf <- data.frame(clusters=colDf$clusters, cells=colDf$cellName)


            ## Check Parameters
            .checkParamsAddClustering(theObject, clusToAdd, colDf)

            ## Modify sceNorm
            colDf$clusters <- clusToAdd[colDf$cells, "clusters"]
            colDf$clusters <- as.factor(as.vector(colDf$clusters))
            clNb <- length(unique(colDf$clusters))
            SummarizedExperiment::colData(sceObject)$clusters <- colDf$clusters
            setSceNorm(theObject) <- sceObject
            
            ## Re-initialize values before re-computing markers
            theObject <- .reinitializeObject(theObject)
            
            ## Recompute markers
            message("Computing new markers..")
            theObject <- calculateClustersSimilarity(theObject)
            theObject <- rankGenes(theObject)
            theObject <- retrieveTopClustersMarkers(theObject, 
                    removeDuplicates=FALSE)
            theObject <- retrieveGenesInfo(theObject)
            
            return(theObject)
        })
