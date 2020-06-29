#' initialisePath
#'
#' Create the output directory
#' \code{\link{AnnotationDbi}} to retrieve annotations.
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
        dir.create(dataDirectory, showWarnings=F)
    
    if(!file.exists(graphsDirectory))
        dir.create(graphsDirectory, showWarnings=F)
    
    if(!file.exists(tSNEPicturesDirectory))
        dir.create(tSNEPicturesDirectory, showWarnings=F)
    
    if(!file.exists(markerGenesDirectory))
        dir.create(markerGenesDirectory, showWarnings=F)
    
    if(!file.exists(tSNEDirectory))
        dir.create(tSNEDirectory, showWarnings=F)
    
    if(!file.exists(outputDataDirectory))
        dir.create(outputDataDirectory, showWarnings=F)
}


#' createDirectory
#'
#' Create output directory and one sub output directory inside.
#'
#' @param dataDirectory Path of the output directory
#' @param directory Name of the output directory inside dataDirectory
#'
#' @keywords internal
#' @noRd
createDirectory <- function(dataDirectory, directory){
    
    newDir <- file.path(dataDirectory, directory)
    
    if(!file.exists(dataDirectory))
        dir.create(dataDirectory, showWarnings=F)
    
    if(!file.exists(newDir))
        dir.create(newDir, showWarnings=F)
}

#' .getTSNEresults
#'
#' Use doParallel to calculate several combinaison of Tsne
#'
#' @param theObject The raw count matrix
#' @param expressionMatrix The studied species
#' @param cores rowData of class data.frame, it contains gene names of the
#' @param PCs rowData of class data.frame, it contains gene names of the
#' @param perplexities The studied species
#' @param randomSeed rowData of class data.frame, it contains gene names of the

#' @keywords internal
#' 
#' @importFrom stats prcomp
#' @importFrom parallel makeCluster stopCluster
#' @import doParallel
#' @return Returns the combinations of tSNES
#' @noRd
.getTSNEresults <- function(theObject, expressionMatrix, cores, PCs, 
		perplexities, randomSeed){
	
    PCAData <- prcomp(t(expressionMatrix))$x
    myCluster <- parallel::makeCluster(cores, type = "PSOCK")
    doParallel::registerDoParallel(myCluster)
	
    tSNECoordinates <- foreach::foreach(PCA=rep(PCs, length(perplexities)),
					perp=rep(perplexities, each=length(PCs)), .combine='cbind',
					.packages="SingleCellExperiment") %dopar% {
				
				listsce <- list(logcounts=t(PCAData[, 1:PCA]))
				sce <- SingleCellExperiment(assays=listsce)
				
				tsneCoord <- scater::runTSNE(sce, scale_features=FALSE,
                perplexity=perp, rand_seed=randomSeed, theme_size=13,
                return_SCESet=FALSE)
            scater::plotTSNE(tsneCoord)
        }

    parallel::stopCluster(myCluster)
    message("Calculated ", length(PCs)*length(perplexities), " 2D-tSNE plots.")
    return(tSNECoordinates)
}


.pickDefaultPalette <- function(clustersNumber, colorPalette26){
	
	if(clustersNumber < 13)
		return(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
						"#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", 
						"#6A3D9A", "#FFFF99", "#B15928")[1:clustersNumber])
	else
		return(colorPalette26[seq_len(length(clustersNumber))]) 
}


.choosePalette <- function(colorPalette, clustersNumber){
    
    colorPalette26 <- c("yellow", "darkgoldenrod1", "coral1", "deeppink",
                         "indianred", "coral4", "darkblue", "darkmagenta",
                         "darkcyan", "mediumorchid", "plum2", "gray73", 
						 "cadetblue3", "khaki", "darkolivegreen1", "lightgreen",
						 "limegreen", "darkolivegreen4", "green", "#CC79A7", 
						 "violetred3", "brown3", "darkslategray1", "gray51", 
						 "slateblue2", "blue")
	
    if(isTRUE(all.equal(colorPalette, "default")) && clustersNumber > 26)
        stop("The default option is limited to 26 colors, please provide your ",
				"own color vector.")
    
    if(!isTRUE(all.equal(colorPalette,"default"))  && 
			clustersNumber > length(colorPalette))
        stop("The number of clusters is greater than the number of given ",
				"colors.")
    
    if(isTRUE(all.equal(colorPalette, "default")))
        return(.pickDefaultPalette(clustersNumber))
    
    return(colorPalette)
}


