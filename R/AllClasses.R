## All classes defined for the scRNAseq package


################################################################################
############################### Tsne class  ####################################
################################################################################

#' The Tsne class
#'
#' @description
#' S4 class containing the features to plot tSNEs. This constructor is internal
#' and is used by the method generateTSNECoordinates.
#'
#' @rdname Tsne-class
#' @aliases Tsne-class Tsne
#'
#' @slot name A 'character' string representing the name of the tSNE
#'             coordinates.
#' @slot pc A 'numeric' value representing the number of principal
#'         components used by CONCLUS to perfom a PCA before
#'         calculating the tSNE.
#' @slot perplexity  A 'numeric' vector. Default: c(30, 40)
#' @slot coordinates A 'numeric' matrix that contains the coordinates
#' of one tSNE solution.
#'
#' @details
#'     Tsne is a vector of principal components (PC) and perplexity that are
#' the parameters necessary to reduce the dimensionality of the data in the
#' the form of a t-distributed stochastic neighbor embedding (t-SNE).
#' For details about perplexities parameter see ‘?Rtsne’.
#'
#' @section Constructor:
#'     Tsne(name = "character", pc = "numeric", perplexity = "numeric",
#'         coordinates = "matrix")
#'
#'     name:        Empty character string or the name of the tSNE. \cr
#'     pc:          Empty 'numeric' number of PCs. \cr
#'     perplexity:  Empty 'numeric' perplexity values. \cr
#'     coordinates: Empty 'numeric' "matrix" or matrix of coordinates. \cr
#'
#' @section Accessors:
#'
#' In the following snippets, x is a Tsne object.
#'
#' getName(x):    Get the name of the tSNE. \cr
#' getPC(x):          Get the PC used. \cr
#' getPerplexity(x):  Get the perplexity used. \cr
#' getCoordinates(x): Get the matrix of tSNE coordinates. \cr
#'
#' @section Subsetting:
#'
#'     In the following snippets, x is a Tsne object.
#'
#'     setName(x) <- value:    Set the name of the tSNE. \cr
#'     setPC(x) <- value:          Set the PC parameter. \cr
#'     setPerplexity(x) <- value:  Set the perplexity parameter. \cr
#'     setCoordinates(x) <- value: Set the matrix of tSNE coordinates. \cr
#'
#' @author Ilyess Rachedi and Nicolas Descostes
#' @seealso generateTSNECoordinates

Tsne <- setClass(
    "Tsne",
    slots = c(
        name = "character",
        pc = "numeric",
        perplexity = "numeric",
        coordinates = "matrix"
    ),

    validity = function(object){

        coordinates <- getCoordinates(object)

        if(isFALSE(all.equal(ncol(coordinates), 2)) ||
            isFALSE(identical(colnames(coordinates), c("X", "Y"))))
            stop("Coordinates should be a matrix with two columns X and Y.")

        if(isFALSE(is.numeric(coordinates)))
            stop("Coordinates should be a matrix of numeric values.")
    })



################################################################################
############################### Dbscan class ###################################
################################################################################

#' The Dbscan class
#'
#' @description
#' S4 class containing the features to plot DBSCAN. This constructor is internal
#' and is used by the method runDBSCAN.
#'
#' @rdname Dbscan-class
#' @aliases Dbscan Dbscan-class
#'
#' @slot name A character string representing the name of the Dbscan
#' clustering.
#' @slot epsilon A numeric vector. The epsilon is the distance to \cr
#' consider two points belonging to the same cluster.
#' Default = c(1.3, 1.4, 1.5).
#' @slot minPoints A numeric value. The minPoints is the minimum number
#' of points to construct a cluster.
#' @slot clustering A matrix that contains the result of one DBSCAN
#' clustering solution.
#'
#' @section Constructor:
#'
#'     Dbscan(name = "character", epsilon = "numeric", minPoints = "numeric",
#'         clustering = "matrix")
#'
#'     name:       Empty character string or the name of the tSNE. \cr
#'     epsilon:    Empty 'numeric' representing the epsilon. \cr
#'     minPoints:  Empty 'numeric' representing the minPoints value. \cr
#'     clustering: Empty 'numeric' "matrix" or matrix of clustering. \cr
#'
#'
#' @section Accessors:
#'
#'     In the following snippets, x is a Dbscan object.
#'
#'     getName(x): Get the name of the Dbscan. \cr
#'     getEpsilon(x):    Get the epsilon used. \cr
#'     getMinPoints(x):  Get the MinPoint used. \cr
#'     getClustering(x): Get the matrix of DBSCAN clustering. \cr
#'
#'
#' @section Subsetting:
#'
#'     In the following snippets, x is a Dbscan object.
#'
#'     setName(x) <- value: Set the name of the Dbscan. \cr
#'     setEpsilon(x) <- value:    Set the epsilon used. \cr
#'     setMinPoints(x) <- value:  Set the minPoints used. \cr
#'     setClustering(x) <- value: Set the matrix of Dbscan clustering. \cr
#'
#' @author Ilyess Rachedi and Nicolas Descostes
#' @seealso runDBSCAN

Dbscan <- setClass(
    "Dbscan",
    slots = c(
        name = "character",
        epsilon    = "numeric",
        minPoints  = "numeric",
        clustering = "matrix"),

    validity = function(object){

        clustering<- getClustering(object)

        if(isFALSE(is.numeric(clustering)))
            stop("'Clustering' slot should be a matrix of integer values.")

    })



################################################################################
############################## scRNAseq class ##################################
################################################################################


.testExperimentNameSlot <- function(object){

    experimentName <- getExperimentName(object)

    if(isTRUE(all.equal(length(experimentName),0)) ||
            isTRUE(all.equal(experimentName,"")))
        stop("'experimentName' slot is empty. Please fill it.")

    if(!is.character(experimentName) | grepl(" ", experimentName))
        stop("Experiment name should contain a single string ",
                "describing the experiment, '", experimentName,
                "' is not correct.")
}


.testCountMatrixSlot <- function(object){

    countMatrix <- getCountMatrix(object)

    if(isFALSE(is.numeric(countMatrix)))
        stop("The count matrix is empty or does not contain whole numbers. ",
                "Please check your count matrix.\n")

    if(isFALSE(is.character(rownames(countMatrix))))
        stop("The name of the lines should be character class. ",
                "Please check your count matrix.\n")

    if(isFALSE(is.character(colnames(countMatrix))))
        stop("The name of the columns should be character class. ",
                "Please check your count matrix.\n")
}

.testsceNormSlot <- function(object){

    sceNorm <- getSceNorm(object)

    if(!is(sceNorm, "SingleCellExperiment"))
        stop("Normalized count matrix should be a SingleCellExperiment ",
                "object and not a '", is(sceNorm), "'.")
}

.testSpeciesSlot <- function(object){

    species <- getSpecies(object)

    if(isTRUE(all.equal(length(species), 0)))
        stop("The species is empty. Please fill it.")

    if (!is.element(species, c("mouse","human")))
        stop("species should be 'mouse' or 'human'. '", species,
                "' is currently not supported.\n")

}

.testOutputDirectorySlot <- function(object){

    outputDirectory <- getOutputDirectory(object)

    if (isTRUE(all.equal(length(outputDirectory), 0)) ||
            isTRUE(all.equal(outputDirectory, "")))
        stop("'outputDirectory' slot is empty. Please fill it.")

    if(!is.character(outputDirectory) | grepl(" ", outputDirectory))
        stop("'outputDirectory' should be a conform folder path:",
                "'", outputDirectory, "' is not.")
}

.testtSNEListSlot <- function(object){

    tSNEList <- getTSNEList(object)

    if(isTRUE(all.equal(length(tSNEList), 0)))
        stop("tSNEList is empty. This should be a list of tSNE objects.\n")

    if(isFALSE(all(vapply(tSNEList, is, class2 = "Tsne",
                            FUN.VALUE = logical(1)))))
        stop("tSNEList should be a list of Tsne objects.")

    invisible(checkList(tSNEList, getCoordinates, "Tsne"))
}

.testDbscanSlot <- function(object){

    dbscanList <- getDbscanList(object)

    if(isTRUE(all.equal(length(dbscanList), 0)))
        stop("dbscanList is empty. This should be a list of dbScan ",
                "objects.\n")

    if(isFALSE(all(
                    vapply(dbscanList, is, class2 = "Dbscan",
                            FUN.VALUE = logical(1)))))
        stop("dbscanList should be a list of Dbscan objects.")


    invisible(checkList(dbscanList, getClustering, "Dbscan"))
}

.testCellsSimilarityMatrixSlot <- function(object){

    cellsSimilarityMatrix <- getCellsSimilarityMatrix(object)

    if(is.null(rownames(cellsSimilarityMatrix)) ||
            is.null(colnames(cellsSimilarityMatrix)))
        stop("'cellsSimilarityMatrix' should have column and row names ",
                "corresponding to cell names.")

    if(isFALSE(identical(colnames(cellsSimilarityMatrix),
                    rownames(cellsSimilarityMatrix))))
        stop("'cellsSimilarityMatrix' should be a square matrix with ",
                "identical names in rows and columns.")

    if(isFALSE(all(
                    vapply(cellsSimilarityMatrix, is, class2 = "numeric",
                            FUN.VALUE = logical(1)))))
        stop("'cellsSimilarityMatrix' should contain only numeric values.")
}

.testClustersSimilarityMatrixSlot <- function(object){

    clustersSimilarityMatrix <- getClustersSimilarityMatrix(object)

    if(is.null(rownames(clustersSimilarityMatrix)) ||
            is.null(colnames(clustersSimilarityMatrix)))
        stop("'clustersSimilarityMatrix' should have column and row names ",
                "corresponding to cluster names.")

    if(isFALSE(identical(colnames(clustersSimilarityMatrix),
                    rownames(clustersSimilarityMatrix))))
        stop("'clustersSimilarityMatrix' should be a square matrix with ",
                "identical names in rows and colums.")

    if(isFALSE(all(
                    vapply(clustersSimilarityMatrix, is, class2 = "numeric",
                            FUN.VALUE = logical(1)))))
        stop("'clustersSimilarityMatrix' should contain only numeric ",
                "values.")
}

.testClustersSimiliratyOrderedSlot <- function(object){

    clustersSimiliratyOrdered <- getClustersSimilarityOrdered(object)
    clustersSimilarityMatrix <- getClustersSimilarityMatrix(object)

    if(isTRUE(all.equal(nrow(clustersSimiliratyOrdered), 0)) &&
            isTRUE(all.equal(ncol(clustersSimiliratyOrdered), 0)))
        stop("'clustersSimiliratyOrdered' is empty. It should be a matrix.")

    if (all(!clustersSimiliratyOrdered %in%
                    rownames(clustersSimilarityMatrix)))
        stop("'clustersSimiliratyOrdered' slot should contain the same ",
                "clusters as 'clustersSimilarityMatrix'.")
}

.testGetMarkerGenesListSlot <- function(object){

    markerGenesList <- getMarkerGenesList(object)
    clustersSimiliratyOrdered <- getClustersSimilarityOrdered(object)

    if(isTRUE(all.equal(length(markerGenesList), 0)))
        stop("markerGenesList is empty. This should be a list of dataframe")

    invisible(checkMarkerGenesList(markerGenesList,
                    clustersSimiliratyOrdered))
}


.testTopMarkersSlot <- function(object){

    topMarkers <- getTopMarkers(object)
    clustersSimiliratyOrdered <- getClustersSimilarityOrdered(object)

    if(isTRUE(all.equal(length(topMarkers), 0)))
        stop("topMarkers is empty. This should be a dataframe")

    invisible(checkTopMarkers(topMarkers,
                    clustersSimiliratyOrdered))
}


.testGenesInfosSlot <- function(object){
    
    genesInfos <- getGenesInfos(object)
    clustersSimiliratyOrdered <- getClustersSimilarityOrdered(object)
    species <- getSpecies(object)
    
    if(isTRUE(all.equal(length(genesInfos), 0)))
        stop("genesInfos is empty. This should be a dataframe")
    
    invisible(checkGenesInfos(genesInfos, species,
                                clustersSimiliratyOrdered))
    
}


#' The scRNAseq class
#'
#' @description
#' S4 class and the main class used by CONCLUS containing the results of
#' the different steps to analyse rare cell populations.
#'
#' @slot experimentName 'character' string representing the name of the
#' experiment.
#' @slot countMatrix An 'integer matrix' representing the raw count matrix
#' with reads or unique molecular identifiers (UMIs).
#' @slot sceNorm Object of class SingleCellExperiment that contains the
#' colData giving informations about cells and the rowData giving
#' informations about genes. It also contains the normalized count matrix.
#' @slot species 'character' string representing the species of interest.
#' Currently limited to "mouse" and "human". Other organisms can be
#' added on demand.
#' @slot outputDirectory A 'character' string of the path to the root
#' output folder.
#' @slot tSNEList List of 'Tsne' objects representing the different tSNE
#' coordinates generated by CONCLUS.
#' @slot dbscanList List of 'Dbscan' objects representing the different
#' Dbscan clustering generated by CONCLUS.
#' @slot suggestedClustersNumber A number got from the dbscan list representing
#' a suggested clusters number to use in clusterCellsInternal().
#' @slot cellsSimilarityMatrix A numeric Matrix defining how many times
#' two cells have been associated to the same cluster across the 84
#' solutions (by default) of clustering.
#' @slot clustersSimilarityMatrix A numeric matrix comparing the
#' robustness of the consensus clusters.
#' @slot clustersSimiliratyOrdered A factor representing the clusters
#' ordered by similarity.
#' @slot markerGenesList List of data.frames. Each data frame contains
#' the ranked genes of one cluster.
#' @slot topMarkers A data frame containing the top 10 (by default)
#' marker genes of each clusters.
#' @slot genesInfos A data frame containing informations of the markers
#' genes for each clusters.
#' 
#' @rdname scRNAseq-class
#' @aliases scRNAseq-class
#'
#' @section Constructor:
#'
#' singlecellRNAseq(experimentName = "character", countMatrix = "matrix",
#' species = "character", outputDirectory = "character")
#'
#' experimentName:  String of the name of the experiment.
#'
#' countMatrix:     Matrix containing the raw counts.
#'
#' species:         'character' string representing the species of interest.
#'                 Shoud be mouse or human. Other organisms can be added on
#'                 demand.
#'
#' outputDirectory: 'character' string representing the path to the output
#' directory.
#'
#'
#' @section Accessors:
#'
#'     In the following snippets, x is a scRNAseq object.
#'
#'     getExperimentName(x):            Get the name of the experiment. \cr
#'     getCountMatrix(x):               Get the count matrix. \cr
#'     getSceNorm(x):                   Get the SingleCellExperiment object
#' used  \cr
#'     getSpecies(x):                   Get the species. \cr
#'     getOutputDirectory(x):           Get the path of the output directory.
#' \cr
#'     getTSNEList(x):                  Get the list of Tsne objects. \cr
#'     getDbscanList(x):                Get the list of Dbscan objects. \cr
#'     getSuggestedClustersNumber(x):   Get the suggested clusters number. \cr
#'     getCellsSimilarityMatrix(x):     Get the cell similarity matrix. \cr
#'     getClustersSimilarityMatrix(x):  Get the cluster similarity matrix. \cr
#'     getClustersSimilarityOrdered(x): Get the clusters ordered by
#' similarity. \cr
#'     getMarkerGenesList(x):           Get the list of marker genes by
#' clusters. \cr
#'     getTopMarkers(x):           Get the most significant markers by
#'                                     clusters into a data.frame. \cr
#'     getGenesInfos(x):                Get a data frame containing informations
#'                                     about marker genes. \cr     
#'     
#' @section Subsetting:
#'
#'     In the following snippets, x is a scRNAseq object.
#'
#'     setExperimentName(x):            Set the name of the experiment. \cr
#'     setCountMatrix(x):               Set the count matrix. \cr
#'     setSceNorm(x):                   Set the SingleCellExperiment object
#' used. \cr
#'     setSpecies(x):                   Set the species. \cr
#'     setOutputDirectory(x):           Set the path of the output directory.
#' \cr
#'     setTSNEList(x):                  Set the list of Tsne objects. \cr
#'     setDbscanList(x):                Set the list of Dbscan objects. \cr
#'     setCellsSimilarityMatrix(x):     Set the cell similarity matrix. \cr
#'     setClustersSimilarityMatrix(x):  Set the cluster similarity matrix. \cr
#'     setClustersSimiliratyOrdered(x): Set the clusters ordered by
#' similarity. \cr
#'     setMarkerGenesList(x):           Set the list of marker genes by
#' clusters \cr
#'     setTopMarkers(x):           Set the most significant markers by
#'                                     clusters. \cr
#'     setGenesInfos(x):                Set a data.frame containing informations
#'                                         about the marker genes. \cr
#'
#' @exportClass scRNAseq
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods new
#' @seealso singlecellRNAseq
#' @author Ilyess Rachedi and Nicolas Descostes

scRNAseq <- setClass(
    "scRNAseq",
    slots = c(
        experimentName = "character",
        countMatrix = "matrix",
        sceNorm = "SingleCellExperiment",
        species = "character",
        outputDirectory = "character",
        tSNEList = "list",
        dbscanList = "list",
        suggestedClustersNumber = "numeric",
        cellsSimilarityMatrix = "matrix",
        clustersSimilarityMatrix = "matrix",
        clustersSimiliratyOrdered = "factor",
        markerGenesList = "list",
        topMarkers = "data.frame",
        genesInfos = "data.frame"
    ),
    prototype = list(
        sceNorm = SingleCellExperiment::SingleCellExperiment(),
        tSNEList = list(new("Tsne")),
        dbscanList = list(new("Dbscan")),
        suggestedClustersNumber = 1,
        cellsSimilarityMatrix =  matrix(nrow = 1, ncol = 1,
                                dimnames = list("c1", "c1"), data = 1),
        clustersSimilarityMatrix = matrix(nrow = 1, ncol = 1,
                                    dimnames = list("1", "1"), data = 1),
        clustersSimiliratyOrdered = factor(1),
        markerGenesList = list(data.frame(Gene = c("gene1"),
                                        mean_log10_fdr = c(NA),
                                        n_05 = c(NA), score = c(NA))),
        topMarkers = data.frame(geneName="gene1", clusters=NA),

        genesInfos = data.frame(uniprot_gn_symbol=c("symbol"), clusters="1",
                external_gene_name="gene", go_id="GO1,GO2",
                mgi_description="description", entrezgene_description="descr",
                gene_biotype="gene", chromosome_name="1", Symbol="symbol",
                ensembl_gene_id="ENS", mgi_id="MGI", entrezgene_id="1",
                uniprot_gn_id="ID")
    ),

    validity = function(object) {

        .testExperimentNameSlot(object)
        .testCountMatrixSlot(object)
        .testsceNormSlot(object)
        .testSpeciesSlot(object)
        .testOutputDirectorySlot(object)
        .testtSNEListSlot(object)
        .testDbscanSlot(object)
        .testCellsSimilarityMatrixSlot(object)
        .testClustersSimilarityMatrixSlot(object)
        .testClustersSimiliratyOrderedSlot(object)
        .testGetMarkerGenesListSlot(object)
        .testTopMarkersSlot(object)
        .testGenesInfosSlot(object)
    }
)