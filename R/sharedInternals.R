#' initialisePath
#'
#' Create the output directory
#'
#' @param dataDirectory Directory path
#'
#' @keywords internal
#' @noRd
initialisePath <- function(dataDirectory){

    graphsDirectory <- file.path(dataDirectory, "pictures")
    markerGenesDirectory <- file.path(dataDirectory, "marker_genes")
    tSNEDirectory <- file.path(dataDirectory, "tsnes")
    outputDataDirectory <- file.path(dataDirectory, "output_tables")
    tSNEPicturesDirectory <- file.path(dataDirectory, "pictures",
                                        "tSNE_pictures")


    if(!file.exists(dataDirectory))
        dir.create(dataDirectory, showWarnings=FALSE)

    if(!file.exists(graphsDirectory))
        dir.create(graphsDirectory, showWarnings=FALSE)

    if(!file.exists(tSNEPicturesDirectory))
        dir.create(tSNEPicturesDirectory, showWarnings=FALSE)

    if(!file.exists(markerGenesDirectory))
        dir.create(markerGenesDirectory, showWarnings=FALSE)

    if(!file.exists(tSNEDirectory))
        dir.create(tSNEDirectory, showWarnings=FALSE)

    if(!file.exists(outputDataDirectory))
        dir.create(outputDataDirectory, showWarnings=FALSE)
}


#' createDirectory
#'
#' Create output directory and one sub output directory inside.
#'
#' @param dataDirectory Path of the output directory
#' @param directory Name of the output sub-directory
#'
#' @keywords internal
#' @noRd
createDirectory <- function(dataDirectory, directory){

    newDir <- file.path(dataDirectory, directory)

    if(!file.exists(dataDirectory))
        dir.create(dataDirectory, showWarnings=FALSE)

    if(!file.exists(newDir))
        dir.create(newDir, showWarnings=FALSE)
}

#' .getTSNEresults
#'
#' Use doParallel to calculate several combinaison of Tsne
#'
#' @param theObject The raw count matrix
#' @param expressionMatrix The normalized count matrix
#' @param cores The number of cores to use for parallelisation
#' @param PCs Vector of principal components to create combinations of tSNE
#' @param perplexities Vector of perplexities to create combinations of tSNE
#' @param randomSeed Seeds used to create the tSNE
#'
#' @keywords internal
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom stats prcomp
#' @importFrom doParallel registerDoParallel
#' @import Rtsne
#' @return Returns the combinations of tSNES
#' @noRd
.getTSNEresults <- function(expressionMatrix, cores, PCs, perplexities,
        randomSeed){

    PCAData <- prcomp(t(expressionMatrix))$x
    myCluster <- parallel::makeCluster(cores, type = "PSOCK")
    doParallel::registerDoParallel(myCluster)
    
    tSNECoordinates <- foreach::foreach(
                    PCAGetTSNEresults=rep(PCs, length(perplexities)),
                    perpGetTSNEresults=rep(perplexities, each=length(PCs)),
                    .combine='cbind',
                    .packages="SingleCellExperiment") %dopar% {

                listsce <- list(logcounts=t(PCAData[,
                                        seq_len(PCAGetTSNEresults)]))
                sce <- SingleCellExperiment::SingleCellExperiment(
                        assays=listsce)

                set.seed(randomSeed)

                tsneCoord <- scater::runTSNE(sce, scale=FALSE,
                perplexity=perpGetTSNEresults)
        scater::plotTSNE(tsneCoord)
        }

    parallel::stopCluster(myCluster)
    message("Calculated ", length(PCs)*length(perplexities), " 2D-tSNE plots.")
    return(tSNECoordinates)
}


#' .pickDefaultPalette
#'
#' Create a default color palette. Used in .choosePalette
#'
#' @param clustersNumber The number of clusters
#' @param colorPalette26 A palette of 26 colors
#'
#' @keywords internal
#'
#' @return Returns a palette with a number of colors corresponding to the
#' number of clusters
#' @noRd
.pickDefaultPalette <- function(clustersNumber, colorPalette26){

    if(clustersNumber < 13)
        return(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                        "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6",
                        "#6A3D9A", "#FFFF99", "#B15928")[
                        seq_len(clustersNumber)])
    else
        return(colorPalette26[seq_len(clustersNumber)])
}


#' .choosePalette
#'
#' Create a default color palette. Used in .choosePalette
#'
#' @param clustersNumber The number of clusters
#' @param colorPalette Either a palette of colors or the "default" value
#'
#' @keywords internal
#'
#' @return If colorPalette is not equal to 'default', returns the exact same
#' colorPalette, otherwise returns a default color palette with a number of
#' colors equal to the number of clusters.
#' @noRd
.choosePalette <- function(colorPalette, clustersNumber){

    colorPalette26 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                        "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6",
                        "#6A3D9A", "#FFFF99", "#B15928", "darkgoldenrod1",
                        "coral1", "deeppink",
                        "indianred", "coral4", "darkmagenta",
                        "darkcyan", "mediumorchid", "plum2", "gray73",
                        "cadetblue3", "khaki",
                        "violetred3", "brown3")

    if(isTRUE(all.equal(colorPalette, "default")) && clustersNumber > 26)
        stop("The default option is limited to 26 colors, please provide your ",
                "own color vector.")

    if(!isTRUE(all.equal(colorPalette,"default"))  &&
            clustersNumber > length(colorPalette))
        stop("The number of clusters is greater than the number of given ",
                "colors.")

    if(isTRUE(all.equal(colorPalette, "default")))
        return(.pickDefaultPalette(clustersNumber, colorPalette26))

    return(colorPalette)
}



#' .tryUseMart
#'
#' This function retrieves an instance of ensembl biomaRt.
#'
#' @param biomart Name of the database.
#' @param dataset Name of the dataset. Currently mmusculus_gene_ensembl or
#' hsapiens_gene_ensembl.
#'
#' @keywords internal
#'
#' @return A Mart object.
#' @importFrom biomaRt useEnsembl
#' @noRd

.tryUseMart <- function(biomart="ensembl", dataset){

    c <- 1

    repeat{
    message("# Attempt ", c, "/5 # ",
            "Connection to Ensembl ... ")
    ensembl <- try(useEnsembl(biomart, dataset=dataset), silent=FALSE)

    if(isTRUE(is(ensembl, "try-error"))){
        c <- c + 1
        error_type <- attr(ensembl, "condition")
        message(error_type$message)

        if(c > 5){
            stop(error_type$message)
            # stop("There is a problem of connexion to Ensembl for ",
            #         "now. Please retry later.")
        }

    }else{
        message("Connected with success.")
        return(ensembl)
        }
    }

}


#' .tryGetBM
#'
#' This function retrieves the user specified attributes from the BioMart
#' database one is connected to, with five tries to succeed in case of
#' connection problem.
#'
#' @param attributes A vector of attributes you want to retrieve.
#' @param ensembl Object of class Mart, created with  useEnsembl or useEnsembl
#' function
#' @param values Values of the filter/
#' @param filters Filters used in the query.
#'
#' @keywords internal
#'
#' @return A data.frame with attributes.
#' @noRd
.tryGetBM <- function(attributes, ensembl, values=NULL, filters=NULL){

    c <- 1

    repeat{

    message("# Attempt ", c, "/5 # ",
            "Retrieving information about genes from biomaRt ...")


    if (is.null(values) && is.null(filters))
        res <- try(getBM(attributes=attributes, mart=ensembl), silent=FALSE)
    else
        res <- try(getBM(attributes=attributes, mart=ensembl, values=values,
                            filters=filters), silent=FALSE)

    if(isTRUE(is(res, "try-error"))){
        c <- c + 1
        error_type <- attr(res, "condition")
        message(error_type$message)

        if(c > 5){
            stop(error_type$message)
            # stop("There is a problem of connexion to Ensembl for ",
            #         "now. Please retry later.")
        }

    }else{
        message("Information retrieved with success.")
        return(res)
        }
    }

}
