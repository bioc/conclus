.checkParamTsne <- function(theObject, randomSeed, cores, PCs, perplexities,
        writeOutput){

    ## Test the normalized count matrix
    sceObject <- getSceNorm(theObject)

    if(all(dim(sceObject) == c(0,0)))
        stop("The 'scRNAseq' object that you're using with ",
                "'generateTSNECoordinates' function doesn't have its ",
                "'sceNorm' slot updated. Please use 'normaliseCountMatrix'",
                " on the object before.")

    ## Check parameters
    if(!is.numeric(randomSeed))
        stop("'randomSeed' parameter should be an integer.")

    ## Check cores argument
    if(!is.numeric(cores))
        stop("'cores' parameter should be an integer")

    ## Check PCs argument
    if(!is.numeric(PCs))
        stop("'PCs' parameter should be a vector of numeric.")

    ## Check perplexities argument
    if(!is.numeric(perplexities))
        stop("'perplexities' parameter should be a vector of numeric.")

    ## Check writeOutput argument
    if(!is.logical(writeOutput))
        stop("'writeOutput' parameter should be a boolean.")

    return(sceObject)
}


#' .buildingTsneObjects
#'
#' @description
#' This internal function builds a list of tSNE object. The function is used in
#' generateTSNECoordinates to update the tsneList slot of theObject with
#' \code{setTSNEList(theObject) <-}.
#'
#' @param PCs Vector of first principal components. For example, to take ranges
#' 1:5 and 1:10 write c(5, 10). Default = c(4, 6, 8, 10, 20, 40, 50)
#' @param perplexities A vector of perplexity (t-SNE parameter).
#' Default = c(30, 40).
#' @param experimentName Name of the analysis accessed from a scRNASeq object
#' with \code{getExperimentName}.
#' @param TSNEres Result of the function scater::runTSNE that is called in the
#' internal function .getTSNEresults.
#' @keywords internal
#'
#' @return Returns a list of Tsne objects.
#' @noRd
.buildingTsneObjects <- function(PCs, perplexities, experimentName, TSNEres){
    

    vec <- unlist(lapply(seq_len(length(PCs)*length(perplexities)),
                        function(i, PCA, perp){

        name <- paste0(experimentName, '_tsne_coordinates_', i,
                        "_" , PCA[i], "PCs_", perp[i], "perp")

        tSNEObj <- TsneCluster(name = name,
                        coordinates = as.matrix(TSNEres[1, i][[1]]),
                        perplexity = perp[i],
                        pc  = PCA[i])

        return(tSNEObj)

    }, rep(PCs, length(perplexities)),
    rep(perplexities, each=length(PCs))))



    return(vec)

}


#' .writeOutputTsne
#'
#' @description
#' Export the tSNE coordinates to an output folder if the parameter
#' \code{writeOutput} of the method \code{generateTSNECoordinates} is TRUE.
#'
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized. See ?normaliseCountMatrix.
#' @param PCs Vector of first principal components. For example, to take ranges
#' 1:5 and 1:10 write c(5, 10). Default = c(4, 6, 8, 10, 20, 40, 50)
#' @param perplexities A vector of perplexity (t-SNE parameter).
#' Default = c(30, 40).
#' @param experimentName Name of the analysis accessed from a scRNASeq object
#' with \code{getExperimentName}.
#' @param TSNEres Result of the function scater::runTSNE that is called in the
#' internal function .getTSNEresults.
#'
#' @keywords internal
#'
#' @return Nothing. Write tSNE coordinates to the output directory in the sub-
#' directory tsnes.
#' @noRd
.writeOutputTsne <- function(theObject, PCs, perplexities, experimentName,
        TSNEres){
    
    dataDirectory <- getOutputDirectory(theObject)
    initialisePath(dataDirectory)
    tSNEDirectory <- "tsnes"
    outputDir <- file.path(dataDirectory, tSNEDirectory)
    
    if(!file.exists(outputDir))
        dir.create(outputDir, showWarnings=FALSE, recursive = TRUE)
    
    invisible(lapply(seq_len(length(PCs)*length(perplexities)), function(i, PCA,
                            perp){

                        name <- paste0(experimentName, '_tsne_coordinates_', i,
                                "_" , PCA[i], "PCs_", perp[i], "perp")

                        write.table(TSNEres[1, i][[1]],
                                file = file.path(dataDirectory, tSNEDirectory,
                                        paste0(name, ".tsv")), quote=FALSE,
                                sep='\t')


                    }, rep(PCs, length(perplexities)),
                    rep(perplexities, each=length(PCs))))

    filesList <- list.files(outputDir, pattern="_tsne_coordinates_")
    deletedFiles <- lapply(filesList, function(fileName)
                file.remove(file.path(outputDir, fileName)))
    saveRDS(TSNEres, file=file.path(dataDirectory, "output_tables",
                    paste0(experimentName, "_tSNEResults.rds")))

}


#' generateTSNECoordinates
#'
#' @description
#' The function generates several t-SNE coordinates based on given perplexity
#' and ranges of PCs. The final number of t-SNE plots is
#' length(PCs)*length(perplexities).
#'
#' @usage
#' generateTSNECoordinates(theObject, randomSeed=42, cores=2,
#'                 PCs=c(4, 6, 8, 10, 20, 40, 50), perplexities=c(30,40),
#'                 writeOutput = FALSE)
#'
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized. See ?normaliseCountMatrix.
#' @param randomSeed  Default is 42. Seeds used to generate the tSNE.
#' @param cores Maximum number of jobs that CONCLUS can run in parallel.
#' Default is 1.
#' @param PCs Vector of first principal components. For example, to take ranges
#' 1:5 and 1:10 write c(5, 10). Default = c(4, 6, 8, 10, 20, 40, 50)
#' @param perplexities A vector of perplexity (t-SNE parameter). See details.
#' Default = c(30, 40)
#' @param writeOutput If TRUE, write the tsne parameters to the output directory
#' defined in theObject. Default = FALSE.
#'
#' @details
#' Generates an object of fourteen (by default) tables with tSNE coordinates.
#' Fourteen because it will vary seven values of principal components
#' *PCs=c(4, 6, 8, 10, 20, 40, 50)* and two values of perplexity
#' *perplexities=c(30, 40)* in all possible combinations. The chosen values of
#' PCs and perplexities can be changed if necessary. We found that this
#' combination works well for sc-RNA-seq datasets with 400-2000 cells. If you
#' have 4000-9000 cells and expect more than 15 clusters, we recommend to take
#' more first PCs and higher perplexity, for example,
#' *PCs=c(8, 10, 20, 40, 50, 80, 100)* and *perplexities=c(200, 240)*. For
#' details about perplexities parameter see ‘?Rtsne’.
#'
#' @rdname generateTSNECoordinates-scRNAseq
#' @aliases generateTSNECoordinates
#'
#' @return
#' An object of class scRNASeq with its tSNEList slot updated. Also writes
#' coordinates in "dataDirectory/tsnes" subfolder if the parameter writeOutput
#' is TRUE.
#'
#' @examples
#' ## Load the count matrix
#' countmatrixPath <- system.file("extdata/test_countMatrix.tsv", 
#'                             package="conclus")
#' countMatrix <- loadDataOrMatrix(file=countmatrixPath, type="countMatrix")
#' 
#' ## Load the coldata
#' coldataPath <- system.file("extdata/test_colData_filtered.tsv", 
#'                             package="conclus")
#' columnsMetaData <- loadDataOrMatrix(file=coldataPath, type="coldata",
#' 										columnID="cell_ID")
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
#' @seealso
#' normaliseCountMatrix
#'
#' @exportMethod generateTSNECoordinates
#' @importFrom methods validObject
#'
#' @author
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.

setMethod(

        f = "generateTSNECoordinates",

        signature = "scRNAseq",

        definition = function(theObject, randomSeed=42, cores=2,
                PCs=c(4, 6, 8, 10, 20, 40, 50), perplexities=c(30,40),
                writeOutput = FALSE){

            ## Check the Object
            validObject(theObject)

            ## Check method parameters
            sceObject <- .checkParamTsne(theObject, randomSeed, cores, PCs,
                    perplexities, writeOutput)

            message("Running TSNEs using ", cores, " cores.")
            TSNEres <- .getTSNEresults(
                    expressionMatrix=Biobase::exprs(sceObject),
                    cores=cores, PCs=PCs, perplexities=perplexities,
                    randomSeed=randomSeed)

            message("Building TSNEs objects.")
            experimentName <- getExperimentName(theObject)
            setTSNEList(theObject) <- .buildingTsneObjects(PCs, perplexities,
                    experimentName, TSNEres)


            if(writeOutput)
                .writeOutputTsne(theObject, PCs, perplexities, experimentName,
                        TSNEres)

            return(theObject)
        })
