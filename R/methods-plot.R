################################################################################
## plotCellSimilarity
################################################################################

#' .generateAnnotationColors
#'
#' @description 
#' Generate annotation_colors, using by pheatmap(). It's list for specifying 
#' annotation_row and annotation_col track colors manually. It is possible to 
#' define the colors for only some of the features.
#' 
#' @param colData A data frame with information about cells. 
#' @param colorPaletteParameter A vector of colors for clusters. Default = "default",
#'  see details.
#' @param statePalette A vector of colors for states or conditions. See details.
#' 
#' @keywords internal
#' @return list for specifying annotation_row and annotation_col 
#' track colors manually.
#' @include sharedInernals.R
#' @noRd
.generateAnnotationColors <- function(colData, colorPaletteParameter,
                                      statePalette){
    
    clusters <- levels(colData$clusters)
    states <- unique(colData$state)
    clusterNumber <- length(unique(colData$clusters))
    
    colorAnnotationClusters <- .choosePalette(colorPaletteParameter,
                                              clusterNumber)
    colorAnnotationState <- .choosePalette(statePalette, length(states))
    names(colorAnnotationState) <- states
    names(colorAnnotationClusters) <- clusters
    
    return(list(state=colorAnnotationState, clusters=colorAnnotationClusters))
}


#' .orderCellsInCluster
#'
#' @description Order cells according to clustering results.
#' Uses for ordering matrix to further plot it with pheatmap()
#'
#' @param cluster Label of cluster
#' @param colData A data frame with information about cells. 
#' @param mtx count matrix expression
#' @param clusteringMethod Clustering method passed to hclust() function. 
#' See ?hclust for a list of method. Default = "ward.D2"
#' 
#' @keywords internal
#' @return Cells ordered by clustering results.
#' @importFrom stats hclust
#' @noRd
.orderCellsInCluster <- function(cluster, colData, mtx, 
                                 clusteringMethod="ward.D2"){
    
    
    
    cells <- colData[colData$clusters == cluster, ]$cellName
    if(length(cells) > 2){
        tree <- hclust(dist(t(mtx[, cells])), method=clusteringMethod)
        return(cells[tree$order])
    }else
        return(cells)
}


#' .checkParamCellSimilarity
#' 
#' @description Check parameters of plotCellSimilarity
#'
#' @param theObject A scRNAseq object with the cluster similarity matrix got 
#' with calculateClustersSimilarity method.
#' @param orderClusters If True, clusters in the similarity matrix of cells will
#'  be ordered by name. Default = FALSE.
#' @param plotPDF If TRUE export to pdf, if FALSE export to png.
#' @param returnPlot Boolean indicating if the pHeatmap object should  be
#' returned by the function. Default = FALSE.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to TRUE, and forced to true if file is a pipe.
#'
#' @keywords internal
#' @noRd
.checkParamCellSimilarity <-function(theObject, orderClusters, plotPDF,
                                     returnPlot, width, height, onefile,
                                     showRowNames, showColnames, fontsize,
                                     fontsizeRow, widthPNG, heightPNG){
    
    
    ## Verify the object contains clustersSimilarityMatrix
    clustersSimilarityMatrix <- getClustersSimilarityMatrix(theObject)
    clusters <- colData(getSceNorm(theObject))$clusters
    nbrClusters <- length(unique(clusters))
    
    if (!isTRUE((ncol(clustersSimilarityMatrix) == nbrClusters)))
        stop(paste("You have to calculate the cluster similarity matrix",
             "before plotting."))  
    
    ## Verify orderClusters
    if (!is.logical(orderClusters))
        stop("orderClusters should be a boolean.")
    
    ## Verify plotPDF
    if (!is.logical(plotPDF))
        stop("plotPDF should be a boolean.")  
    
    ## Verify returnPlot
    if (!is.logical(returnPlot))
        stop("returnPlot should be a boolean.")  
    
    ## Verify width
    if (!is.numeric(width))
        stop("width should be a numeric.")
    
    ## Verify height
    if (!is.numeric(height))
        stop("height should be a numeric.")  
    
    ## Verify onefile
    if (!is.logical(onefile))
        stop("onefile should be a boolean.")  
    
    ## Verify showRowNames
    if (!is.logical(showRowNames))
        stop("showRowNames should be a boolean.")  
    
    ## Verify showColnames
    if (!is.logical(showColnames))
        stop("showColnames should be a boolean.")  
    
    ## Verify fontsize
    if (!is.numeric(fontsize))
        stop("fontsize should be a numeric.")
    
    ## Verify fontsizeRow
    if (!is.numeric(fontsizeRow))
        stop("fontsizeRow should be a numeric.")
    
    ## Verify widthPNG
    if (!is.numeric(widthPNG))
        stop("widthPNG should be a numeric.")
    
    ## Verify heightPNG
    if (!is.numeric(heightPNG))
        stop("heightPNG should be a numeric.") 
}


#' plotCellSimilarity
#' 
#' @description This function plots similarity matrix as a heatmap, so one can 
#' see similarity between parts of different clusters.
#'
#' @param theObject A scRNAseq object with the cluster similarity matrix got 
#' with calculateClustersSimilarity method.
#' @param colorPalette A vector of colors for clusters. Default = "default",
#'  see details.
#' @param statePalette A vector of colors for states or conditions. See details.
#' @param clusteringMethod Clustering method passed to hclust() function. 
#' See ?hclust for a list of method. Default = "ward.D2"
#' @param orderClusters If True, clusters in the similarity matrix of cells will
#'  be ordered by name. Default = FALSE.
#' @param plotPDF If TRUE export to pdf, if FALSE export to png.
#' @param returnPlot Boolean indicating if the pHeatmap object should  be
#' returned by the function. Default = FALSE.
#' @param width Width of the plot. Default = 7.
#' @param height Height of the plot. Default = 6.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to TRUE, and forced to true if file is a pipe.
#' @param showRowNames pheatmap parameter. Boolean specifying if column names 
#' are be shown.
#' @param showColnames pheatmap parameter. Boolean specifying if column names
#' are be shown.
#' @param fontsize pheatmap parameter. Base fontsize for the plot
#' @param fontsizeRow pheatmap parameter. Fontsize for rownames 
#' @param widthPNG Width of the png
#' @param heightPNG Height of the png
#' @param ... Additional parameters to pass to pdf() or png(), 
#' and pheatmap() functions.
#' 
#' @examples
#' experimentName <- "Bergiers"
#' countMatrix <- as.matrix(read.delim(file.path(
#' "tests/testthat/test_data/test_countMatrix.tsv")))
#' outputDirectory <- "./"
#' columnsMetaData <- read.delim(
#' file.path("extdata/Bergiers_colData_filtered.tsv"))
#' 
#' ## Create the initial object
#' scr <- scRNAseq(experimentName = experimentName, 
#'                 countMatrix     = countMatrix, 
#'                 species         = "mouse",
#'                 outputDirectory = outputDirectory)
#' 
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#' 
#' ## Compute the tSNE coordinates
#' scrTsne <- generateTSNECoordinates(scrNorm, cores=5)
#' 
#' ## Perform the clustering with dbScan
#' scrDbscan <- runDBSCAN(scrTsne, cores=5)
#' 
#' ## Compute the cell similarity matrix
#' scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=10, cores=4)
#' 
#' ## Calculate clusters similarity
#' scrCSM <- calculateClustersSimilarity(scrCCI)
#' 
#' ## Plot the heatmap of the similarity matrix
#' plotCellSimilarity(scrCSM)
#' 
#' 
#' @seealso calculateClustersSimilarity 
#' @exportMethod 
#' @importFrom SummarizedExperiment colData
#' @importFrom stats hclust
#' @importFrom pheatmap pheatmap
#' @author 
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
setMethod(
    
    f = "plotCellSimilarity",
    
    signature = "scRNAseq",
    
    definition = function(theObject, colorPalette="default", 
                          statePalette="default", clusteringMethod="ward.D2", 
                          orderClusters=FALSE, plotPDF=TRUE, returnPlot=FALSE, width=7, 
                          height=6, onefile=FALSE, showRowNames=FALSE, showColnames=FALSE,
                          fontsize=7.5, fontsizeRow=0.03, widthPNG=800, heightPNG=750,
                          ...){
        
        ## Verify parameters
        validObject(theObject)
        .checkParamCellSimilarity(theObject, orderClusters, plotPDF,
                                  returnPlot, width, height, onefile,
                                  showRowNames, showColnames, fontsize,
                                  fontsizeRow, widthPNG, heightPNG)
        sceObject <- getSceNorm(theObject)
        cellsSimilarityMatrix <- getCellsSimilarityMatrix(theObject)
        dataDirectory  <- getOutputDirectory(theObject)
        experimentName <- getExperimentName(theObject)
        
        graphsDirectory <- "pictures"
        colData <- SummarizedExperiment::colData(sceObject)
        clustersNumber <- length(unique(colData$clusters))
        
        if(orderClusters){
            # Ordering expressionMatrixrix
            newOrder <- sapply(levels(colData$clusters), function(cluster){ 
                .orderCellsInCluster(cluster, colData, expressionMatrix,
                                     clusteringMethod)})
            newOrder <- unname(unlist(newOrder))
            cellsSimilarityMatrix <- cellsSimilarityMatrix[newOrder, newOrder]
            clusterCols <- FALSE
            clusterRows <- FALSE
            
        } else {
            distanceMatrix <- as.dist(sqrt((1-cellsSimilarityMatrix)/2))
            clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
            clusterCols <- clusteringTree
            clusterRows <- clusteringTree
        }
        
        annotationColors <- .generateAnnotationColors(colData, colorPalette, 
                                                      statePalette)
        
        columnsToPlot <- switch(is.null(colData$state) + 1, c("clusters", 
                                                              "state"), c("clusters"))
        
        annotationCol <- as.data.frame(colData[columnsToPlot])
        nbCol <- ncol(cellsSimilarityMatrix)
        nbRow <- nrow(cellsSimilarityMatrix)
        mainTitle <- paste0("Cells similarity matrix ", nbCol, " columns, ", 
                            nbRow, " rows.")
        pheatmapObject <- pheatmap::pheatmap(cellsSimilarityMatrix,
                                             show_colnames=showColnames,
                                             show_rownames=showRowNames,
                                             annotation_col=annotationCol,
                                             annotation_colors=annotationColors,
                                             fontsize_row=fontsizeRow,
                                             cluster_cols=clusterCols,
                                             cluster_rows=clusterRows,
                                             fontsize=fontsize,
                                             main=mainTitle)
        
        message("\nSaving a heatmap with the cell similarity matrix.")
        
        subdir <- file.path(dataDirectory, graphsDirectory)
        if(!file.exists(subdir))
            dir.create(subdir, showWarnings=F, recursive = TRUE)
        
        if(plotPDF)
            pdf(file=file.path(dataDirectory, graphsDirectory, 
                               paste(experimentName, "cells_correlation", 
                                     clustersNumber, "clusters.pdf", sep="_")),
                width=width, height=height, onefile=onefile)
        else{
            message("Plot type is not pdf. Saving in png.")
            filePath <- file.path(dataDirectory, graphsDirectory, 
                                  paste(experimentName, "cells_correlation", clustersNumber,
                                        "clusters.png", sep="_"))
            png(filename= filePath, width=widthPNG, height=heightPNG, 
                type="cairo")
        }
        
        grid::grid.newpage()
        grid::grid.draw(pheatmapObject$gtable)
        dev.off()
        
        if(returnPlot)
            return(pheatmapObject)
    })


################################################################################
## plotClusteredTSNE
################################################################################

#' .checkParamPlotTSNE
#' 
#' @description Check parameters of plotClusteredTSNE
#'
#' @param theObject A scRNAseq object with the cluster similarity matrix got 
#' with calculateClustersSimilarity method.
#' @param PCs A vector of first principal components. For example, to take 
#' ranges 1:5 and 1:10 write c(5, 10). Default = c(4, 6, 8, 10, 20, 40, 50)
#' @param perplexities Numeric scalar defining the perplexity parameter, 
#' see ‘?Rtsne’ for more details. Default = c(30, 40)
#' @param columnName Name of the column to plot on t-SNE dimensions.
#'  Possible values are clusters, noColor, or state.
#' @param returnPlot Boolean indicating if the pHeatmap object should  be
#' returned by the function. Default = FALSE.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to TRUE, and forced to true if file is a pipe.
#' 
#' @keywords internal
#' @noRd
.checkParamPlotTSNE <- function(theObject, PCs, perplexities, columnName,
                                returnPlot, width, height, onefile){
    
    ## Verify the object contains clustersSimilarityMatrix
    clustersSimilarityMatrix <- getClustersSimilarityMatrix(theObject)
    clusters <- colData(getSceNorm(theObject))$clusters
    nbrClusters <- length(unique(clusters))
    
    if (!isTRUE((ncol(clustersSimilarityMatrix) == nbrClusters)))
        stop(paste("You have to calculate the cluster similarity matrix",
             "before plotting."))  
    
    ## Verify PCs 
    if(!is.numeric(PCs))
        stop("'PCs' parameter should be a vector of numeric.")
    
    ## Verify perplexities
    if(!is.numeric(perplexities))
        stop("'perplexities' parameter should be a vector of numeric.")
    
    ## Verify returnPlot
    if (!is.logical(returnPlot))
        stop("returnPlot should be a boolean.")  
    
    ## Verify columnName
    if(!isTRUE(all.equal(columnName, "clusters")) && 
       !isTRUE(all.equal(columnName, "noColor")) &&
       !isTRUE(all.equal(columnName,"state")))
        stop("columnName should be: clusters, noColor, or state.")
    
    ## Verify width
    if (!is.numeric(width))
        stop("width should be a numeric.")
    
    ## Verify height
    if (!is.numeric(height))
        stop("height should be a numeric.")  
    
    ## Verify onefile
    if (!is.logical(onefile))
        stop("onefile should be a boolean.")  
}			



#' .createTSNEDir
#'
#' @description Create subfolder for TSNE plot and the colorPalette
#' 
#' @param theObject A scRNAseq object with the cluster similarity matrix got 
#' with calculateClustersSimilarity method.
#' @param columnName Name of the column to plot on t-SNE dimensions.
#'  Possible values are clusters, noColor, or state.
#' @param sceObject A SingleCellExperiment object with your experiment.
#' @param colorPalette A vector of colors for clusters. Default = "default",
#'  see details.
#'
#' @return List containint path of the subdirectory for tSNE plots and the color
#' palette
#' @importFrom SummarizedExperiment colData
#' @include sharedInernals.R
#' @keywords internal
#' @noRd
.createTSNEDir <- function(theObject, columnName, sceObject, colorPalette){
    
    dataDirectory  <- getOutputDirectory(theObject)
    graphsDirectory     <- "pictures"
    graphsTSNEDirectory <- "tSNE_pictures"
    
    if(isTRUE(all.equal(columnName, "noColor")))
        numberElements <- NULL
    else{
        
        nb <- unique(SummarizedExperiment::colData(sceObject)[, columnName])
        numberElements <- length(nb)
        colorPalette <- .choosePalette(colorPalette, numberElements)
    }
    
    
    outputDir <- file.path(dataDirectory, graphsDirectory, 
                           graphsTSNEDirectory, paste("tSNE", numberElements, columnName,
                                                      sep="_"))
    
    if(!file.exists(outputDir))
        dir.create(outputDir, showWarnings=F, recursive = TRUE)
    
    return(list(outputDir, colorPalette))
}			


#' plotClusteredTSNE
#' 
#' @description Plot t-SNE generated with different PCs and perplexities.
#' It can also use a coloring scheme by clusters or states.
#'
#' @usage plotClusteredTSNE(sceObject, dataDirectory, experimentName,
#' tSNEresExp = "", colorPalette = "default", PCs = c(4, 6, 8, 10, 20,40, 50), 
#' perplexities = c(30, 40), columnName = "clusters", returnPlot = FALSE,
#'  width = 6, height = 5, ...)
#'  
#' @param theObject A scRNAseq object with the cluster similarity matrix got 
#' with calculateClustersSimilarity method.
#' @param colorPalette A vector of colors for clusters. Default = "default",
#'  see details.
#' @param PCs A vector of first principal components. For example, to take 
#' ranges 1:5 and 1:10 write c(5, 10). Default = c(4, 6, 8, 10, 20, 40, 50)
#' @param perplexities Numeric scalar defining the perplexity parameter, 
#' see ‘?Rtsne’ for more details. Default = c(30, 40)
#' @param columnName Name of the column to plot on t-SNE dimensions.
#'  Possible values are clusters, noColor, or state.
#' @param returnPlot Boolean indicating if the pHeatmap object should  be
#' returned by the function. Default = FALSE.
#' @param width Width of the plot. Default = 6.
#' @param height Height of the plot. Default = 5.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' @param ... Additional parameters to pass to pdf() functions.
#' 
#'
#' @examples
#' experimentName <- "Bergiers"
#' countMatrix <- as.matrix(read.delim(file.path(
#' "tests/testthat/test_data/test_countMatrix.tsv")))
#' outputDirectory <- "./"
#' columnsMetaData <- read.delim(
#' file.path("extdata/Bergiers_colData_filtered.tsv"))
#' 
#' ## Create the initial object
#' scr <- scRNAseq(experimentName = experimentName, 
#'                 countMatrix     = countMatrix, 
#'                 species         = "mouse",
#'                 outputDirectory = outputDirectory)
#' 
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#' 
#' ## Compute the tSNE coordinates
#' scrTsne <- generateTSNECoordinates(scrNorm, cores=5)
#' 
#' ## Perform the clustering with dbScan
#' scrDbscan <- runDBSCAN(scrTsne, cores=5)
#' 
#' ## Compute the cell similarity matrix
#' scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=10, cores=4)
#' 
#' ## Calculate clusters similarity
#' scrCSM <- calculateClustersSimilarity(scrCCI)
#' 
#' ## Plot the heatmap of the similarity matrix
#' plotClusteredTSNE(scrCSM)
#' 
#' @seealso calculateClustersSimilarity 
#' @importFrom ggplot2 ggplot
#' @importFrom SummarizedExperiment colData
#' @keywords internal
#' @noRd
setMethod(
    
    f = "plotClusteredTSNE",
    
    signature = "scRNAseq",
    
    definition = function(theObject, colorPalette="default",
                          PCs=c(4, 6, 8, 10, 20, 40, 50), perplexities=c(30, 40),
                          columnName="clusters", returnPlot=FALSE, width=6,
                          height=5, onefile=FALSE, ...){
        
        ## Verify parameters
        validObject(theObject)
        .checkParamPlotTSNE(theObject, PCs, perplexities, columnName,
                            returnPlot, width, height, onefile)
        
        ## Creating output folder			
        sceObject <- getSceNorm(theObject)
        resultDir <- .createTSNEDir(theObject, columnName, sceObject, 
                                    colorPalette)
        outputdir <- resultDir[[1]]
        colorPalette <- resultDir[[2]]
        
        # plots picture based on t-SNE coordinates from
        # generateTSNECoordinates() and clustering results
        # from clusterCellsInternal() or runClustering()
        
        ### Plot all precalculated pSNEs to show your clusters
        PCA <- rep(PCs, length(perplexities))
        perp <- rep(perplexities, each=length(PCs))
        tSNEList <- getTSNEList(theObject)
        totalLength <- length(PCs)*length(perplexities)
        
        if(!isTRUE(all.equal(length(tSNEList), totalLength)))
            stop("The number of elements of TSNEList is not equal to ",
                 "PCs x perplexities. Contact the developper.")
        
        tSNEplots <- lapply(tSNEList, function(currentTsne, sceObj, outputDir,
                                               widthpdf, heightpdf, onefilepdf){
            
            tSNEres <- as.data.frame(getCoordinates(currentTsne))
            coordinatesName <- getName(currentTsne)
            colDatDf <- SummarizedExperiment::colData(sceObj)
            
            ## Get coordinates of filtered cells
            cellVec <- colDatDf$cellName
            tSNEres <-tSNEres[rownames(tSNEres) %in% cellVec, ]
            
            ## Add new columns
            if(!isTRUE(all.equal(columnName, "noColor")))
                tSNEres[columnName] <- factor(colDatDf[, columnName])
            
            ## Create the pdf
            pdf(paste0(outputDir, "/", coordinatesName, ".pdf"),
                width=widthpdf, height=heightpdf, 
                onefile=onefilepdf, ...)
            
            if(isTRUE(all.equal(columnName, "noColor")))
                tmp <- ggplot2::ggplot(tSNEres, 
                                       aes_string(x=names(tSNEres)[1],
                                                  y=names(tSNEres)[2])) +
                geom_point(size=I(1)) + theme_bw()
            
            else {
                tmp <- ggplot2::ggplot(tSNEres, aes_string(
                    x=names(tSNEres)[1],
                    y=names(tSNEres)[2],
                    color=columnName)) +
                    geom_point(size=I(1)) +
                    scale_color_manual(values=colorPalette) + 
                    theme_bw()
            }
            
            print(tmp)
            dev.off()
            
            return(tmp)
            
        }, sceObject, outputdir, width, height, onefile)
        
        if(returnPlot){
            return(tSNEplots)
        }
        rm(PCA, perp)
    })


#################
## plotCellHeatmap
#################

#' .orderGenesInCluster
#'
#' @description Order cells according to clustering results.
#' Uses for ordering matrix to further plot it with pheatmap()
#' 
#' @param cluster Label of cluster
#' @param markersClusters Data frame containing two columns geneName and 
#' clusters.
#' @param mtx count matrix expression
#' @param clusteringMethod Clustering method passed to hclust() function. 
#' See ?hclust for a list of method. Default = "ward.D2"
#' 
#' @keywords internal
#' @return Genes ordered
#' @importFrom stats hclust
#' @noRd
.orderGenesInCluster <- function(cluster, markersClusters, mtx, 
                                 clusteringMethod="ward.D2"){
    
    genes <- markersClusters[markersClusters$clusters == cluster, ]$geneName
    if(length(genes) > 2){
        tree <- hclust(dist(mtx[genes, ]), method=clusteringMethod)
        return(genes[tree$order])
    } else {
        return(genes)
    }
}

#' .plotCellHeatmap
#' 
#'
#' @param markersClusters Data frame containing two columns geneName and 
#' clusters.
#' @param sceObject A SingleCellExperiment object with your experiment.
#' @param dataDirectory Output directory of a given CONCLUS run
#' @param experimentName Prefix used for output files
#' @param fileName Name of the output file
#' @param meanCentered Boolean, should mean centering be applied to the 
#' expression data or not. Default = TRUE.
#' @param colorPalette A vector of colors for clusters. Default = "default",
#'  see details.
#' @param statePalette A vector of colors for states or conditions. See details.
#' @param clusteringMethod Clustering method passed to hclust() function. 
#' See ?hclust for a list of method. Default = "ward.D2"
#' @param orderClusters If True, clusters in the similarity matrix of cells will
#'  be ordered by name. Default = FALSE.
#' @param similarity Boolean, should the heatmap be structured by similarities.
#'  Default = FALSE.
#' @param orderGenes Boolean, should the heatmap be structured by gene.
#'  Default = FALSE.
#' @param returnPlot If TRUE export to pdf, if FALSE export to png.
#' @param saveHeatmapTable Boolean, whether to save the expression matrix 
#' used for heatmap into a .csv file or not. 
#' The file will be saved into 'dataDirectory/output_tables' with the same name
#'  as the .pdf plot.
#' @param width Width of the plot. Default = 7.
#' @param height Height of the plot. Default = 6.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to TRUE, and forced to true if file is a pipe.
#' @param clusterCols Boolean.
#' @param showColnames pheatmap parameter. Boolean specifying if column names
#' are be shown.
#' @param fontsize pheatmap parameter. Base fontsize for the plot
#' @param fontsizeRow pheatmap parameter. Fontsize for rownames 
#' 
#' @keywords internal
#' @importFrom SummarizedExperiment colData
#' @importFrom Biobase exprs
#' @importFrom stats hclust
#' @importFrom pheatmap pheatmap
#' @noRd
.plotCellHeatmap <- function(markersClusters, sceObject, dataDirectory,
                             experimentName, fileName, meanCentered=TRUE, colorPalette="default",
                             statePalette="default", clusteringMethod="ward.D2", orderClusters=FALSE,
                             similarity=FALSE, orderGenes=FALSE, returnPlot=FALSE, 
                             saveHeatmapTable=FALSE, width=10, height=8.5, onefile=FALSE, 
                             clusterCols=FALSE, showColnames=FALSE, fontsize=7.5, fontsizeRow=8){
    
    # plots correlation between cells between clusters
    colData <- SummarizedExperiment::colData(sceObject)
    exprsTmp <- Biobase::exprs(sceObject)
    rowTmp <- rownames(exprsTmp)
    expressionMatrix <- exprsTmp[rowTmp %in% markersClusters$geneName, ]
    
    if(meanCentered){
        meanRows <- rowSums(expressionMatrix) / ncol(expressionMatrix)
        expressionMatrix <- expressionMatrix - meanRows
    }
    
    if(!orderClusters & orderGenes){
        message("Genes cannot be ordered without clusters.
						Returning heatmap with similarity=TRUE")
        similarity <- TRUE
    }
    
    if(orderClusters){
        # Ordering expressionMatrixrix
        newOrder <- unname(unlist(sapply(levels(colData$clusters), 
                                         function(cluster)
                                             .orderCellsInCluster(cluster, colData,
                                                                  expressionMatrix, clusteringMethod))))
        expressionMatrix <- expressionMatrix[, newOrder]
        clusterCols <- FALSE
        
        if(orderGenes){
            newOrder <- unname(unlist(sapply(levels(colData$clusters), 
                                             function(cluster)
                                                 .orderGenesInCluster(cluster, 
                                                                      markersClusters, 
                                                                      expressionMatrix, 
                                                                      clusteringMethod))))
            expressionMatrix <- expressionMatrix[newOrder, ]
            clusterRows <- FALSE
        }
        
    } else if(similarity){
        
        cellsSimilarityMatrix <- read.delim(file.path(dataDirectory, 
                                                      "output_tables", paste0(experimentName, 
                                                                              "_cellsSimilarityMatrix.csv")), 
                                            stringsAsFactors=FALSE, header=TRUE, sep=",")
        
        clustersSimOrdered <- calculateClustersSimilarity(cellsSimilarityMatrix,
                                                          sceObject,
                                                          clusteringMethod)[[2]]
        newOrder <- unname(unlist(sapply(clustersSimOrdered, function(cluster)
            .orderCellsInCluster(cluster, colData,
                                 expressionMatrix, 
                                 clusteringMethod))))
        expressionMatrix <- expressionMatrix[, newOrder]
        clusterCols <- FALSE
        
        if(orderGenes){
            newOrder <- unname(unlist(sapply(clustersSimOrdered, 
                                             function(cluster)
                                                 .orderGenesInCluster(cluster, 
                                                                      markersClusters, 
                                                                      expressionMatrix,
                                                                      clusteringMethod))))
            expressionMatrix <- expressionMatrix[newOrder, ]
            clusterRows <- FALSE
        }
    } else if(!orderClusters){
        distanceMatrix <- dist(t(expressionMatrix))
        clusterCols <- hclust(distanceMatrix, method="ward.D2")
    }
    
    if(!orderGenes)
        clusterRows <- hclust(dist(expressionMatrix), method="ward.D2")
    
    annotationColors <- .generateAnnotationColors(colData, colorPalette,
                                                  statePalette)
    
    columnsToPlot <- switch(is.null(colData$state) + 1, c("clusters", "state"),
                            c("clusters"))
    
    if(is.null(colData$clusters)){
        
        annCol <- switch(is.null(colData$state) + 1, 
                         as.data.frame(colData["state"]), NA)
        annColors <- switch(is.null(colData$state) + 1, 
                            annotationColors[names(annotationColors) == 
                                                 "state"], NA)
    } else {
        annCol <- as.data.frame(colData[columnsToPlot])
        annColors <- annotationColors
    }
    
    color <- colorRampPalette(c("#023b84", "#4b97fc", "#c9d9ef", "#FEE395",
                                "#F4794E", "#D73027", "#a31008", "#7a0f09"))(100)		
    
    pheatmapObject <- pheatmap::pheatmap(expressionMatrix, 
                                         show_colnames=showColnames, annotation_col=annCol, 
                                         annotation_colors=annColors, fontsize_row=fontsizeRow, 
                                         cluster_cols=clusterCols, cluster_rows=clusterRows,
                                         color=color, fontsize=fontsize)
    
    
    subdir <- file.path(dataDirectory, "pictures")
    if(!file.exists(subdir))
        dir.create(subdir, showWarnings=F, recursive = TRUE)
    
    pdf(file.path(dataDirectory, "pictures", 
                  paste0(experimentName, "_", fileName, ".pdf")), width=width,
        height=height, onefile=onefile)
    grid::grid.newpage()
    grid::grid.draw(pheatmapObject$gtable)
    dev.off()
    
    if(saveHeatmapTable)
        exportMatrix(expressionMatrix, dataDirectory, experimentName, fileName)
    
    if(returnPlot)
        return(pheatmapObject)
}



#' .checkParamCellHeatmap 
#'
#' @description checks parameters of plotCellHeatmap
#' 
#' @param theObject A scRNAseq object with the cluster similarity matrix got 
#' with calculateClustersSimilarity method.
#' @param fileName
#' @param meanCentered Boolean, should mean centering be applied to the 
#' expression data or not. Default = TRUE.
#' @param orderClusters If True, clusters in the similarity matrix of cells will
#'  be ordered by name. Default = FALSE.
#' @param similarity Boolean, should the heatmap be structured by similarities.
#'  Default = FALSE.
#' @param orderGenes Boolean, should the heatmap be structured by gene.
#'  Default = FALSE.
#' @param returnPlot If TRUE export to pdf, if FALSE export to png.
#' @param saveHeatmapTable Boolean, whether to save the expression matrix 
#' used for heatmap into a .csv file or not. 
#' The file will be saved into 'dataDirectory/output_tables' with the same name
#' as the .pdf plot.
#' @param width Width of the plot. Default = 10.
#' @param height Height of the plot. Default = 8.5.
#' 
#' @keywords internal
#' @noRd
.checkParamCellHeatmap <- function(theObject, fileName, 
                                   meanCentered, orderClusters, similarity,
                                   orderGenes, returnPlot, saveHeatmapTable,
                                   width, height){
    
    ## Verify the object 
    markersClusters <- getClustersMarkers(theObject)
    
    if (!isTRUE((nrow(markersClusters) > 1)))
        stop(paste("You have to calculate the cluster markers before plotting.",
                   "Please see retrieveTopClustersMarkers method."))
    
    ## Verify fileName
    if(!is.character(fileName) || grepl("/", fileName , fixed = TRUE))
        stop("fileName should be a string, no path.")
    
    ## Verify meanCentered
    if (!is.logical(meanCentered))
        stop("meanCentered should be a boolean.")
    
    ## Verify orderClusters
    if (!is.logical(orderClusters))
        stop("orderClusters should be a boolean.")
    
    ## Verify similarity
    if (!is.logical(similarity))
        stop("similarity should be a boolean.")
    
    ## Verify orderGenes
    if (!is.logical(orderGenes))
        stop("orderGenes should be a boolean.")
    
    ## Verify returnPlot
    if (!is.logical(returnPlot))
        stop("returnPlot should be a boolean.")
    
    ## Verify saveHeatmapTable
    if (!is.logical(saveHeatmapTable))
        stop("saveHeatmapTable should be a boolean.")
    
    ## Verify width
    if (!is.numeric(width))
        stop("width should be a numeric.")
    
    ## Verify height
    if (!is.numeric(height))
        stop("height should be a numeric.")  
    
}			



#' plotCellHeatmap
#'
#' @description This function plots heatmap with marker genes on rows and
#'  clustered cells on columns.
#'  
#' @param theObject A scRNAseq object with the cluster similarity matrix got 
#' with calculateClustersSimilarity method.
#' @param fileName Name of the output file
#' @param meanCentered Boolean, should mean centering be applied to the 
#' expression data or not. Default = TRUE.
#' @param colorPalette A vector of colors for clusters. Default = "default",
#'  see details.
#' @param statePalette A vector of colors for states or conditions. See details.
#' @param clusteringMethod Clustering method passed to hclust() function. 
#' See ?hclust for a list of method. Default = "ward.D2"
#' @param orderClusters If True, clusters in the similarity matrix of cells will
#'  be ordered by name. Default = FALSE.
#' @param similarity Boolean, should the heatmap be structured by similarities.
#'  Default = FALSE.
#' @param orderGenes Boolean, should the heatmap be structured by gene.
#'  Default = FALSE.
#' @param returnPlot If TRUE export to pdf, if FALSE export to png.
#' @param saveHeatmapTable Boolean, whether to save the expression matrix 
#' used for heatmap into a .csv file or not. 
#' The file will be saved into 'dataDirectory/output_tables' with the same name
#'  as the .pdf plot.
#' @param width Width of the plot. Default = 10.
#' @param height Height of the plot. Default = 8.5.
#' @param ... Additional parameters to pass to pdf() and pheatmap() functions.
#' 
#' @examples
#' experimentName <- "Bergiers"
#' countMatrix <- as.matrix(read.delim(file.path(
#' "tests/testthat/test_data/test_countMatrix.tsv")))
#' outputDirectory <- "./"
#' columnsMetaData <- read.delim(
#' file.path("extdata/Bergiers_colData_filtered.tsv"))
#' 
#' ## Create the initial object
#' scr <- scRNAseq(experimentName = experimentName, 
#'                 countMatrix     = countMatrix, 
#'                 species         = "mouse",
#'                 outputDirectory = outputDirectory)
#' 
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#' 
#' ## Compute the tSNE coordinates
#' scrTsne <- generateTSNECoordinates(scrNorm, cores=5)
#' 
#' ## Perform the clustering with dbScan
#' scrDbscan <- runDBSCAN(scrTsne, cores=5)
#' 
#' ## Compute the cell similarity matrix
#' scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=10, cores=4)
#' 
#' ## Calculate clusters similarity
#' scrCSM <- calculateClustersSimilarity(scrCCI)
#' 
#' ## Ranking genes
#' scrS4MG <- rankGenes(scrCSM)
#' 
#' ## Retrieve the top 10 markers per cluster
#' scrFinal <- retrieveTopClustersMarkers(scrS4MG)
#' 
#' ## Plot the heatmap with marker genes
#' plotCellHeatmap(scrFinal)
#' 
#' @seealso retrieveTopClustersMarkers
#' @exportMethod 
#' @author 
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
setMethod(
    
    f = "plotCellHeatmap",
    
    signature = "scRNAseq",
    
    definition = function(theObject, fileName = NA, meanCentered=TRUE, 
                          colorPalette="default", statePalette="default", 
                          clusteringMethod="ward.D2", orderClusters=FALSE,
                          similarity=FALSE, orderGenes=FALSE, returnPlot=FALSE,
                          saveHeatmapTable=FALSE, width=10, height=8.5, ...){
        
        
        ## Verify parameters
        validObject(theObject)
        
        sceObjectFiltered <- getSceNorm(theObject)
        if (is.na(fileName)){
            fileName=paste0(
                "clusters",
                length(
                    levels(
                        SummarizedExperiment::colData(
                            sceObjectFiltered)$clusters)),
                "_meanCentered",
                meanCentered,
                "_orderClusters",
                orderClusters,
                "_orderGenes",
                orderGenes,
                # "_top",
                # genesNumber,
                "markersPerCluster")
        }
        
        .checkParamCellHeatmap(theObject, fileName, 
                               meanCentered, orderClusters, similarity,
                               orderGenes, returnPlot, saveHeatmapTable,
                               width, height)
        
        markersClusters <- getClustersMarkers(theObject)
        sceObject       <- getSceNorm(theObject)
        dataDirectory   <- getOutputDirectory(theObject)
        experimentName  <- getExperimentName(theObject)
        .plotCellHeatmap(markersClusters, sceObject, dataDirectory, 
                         experimentName, fileName, meanCentered=meanCentered,
                         colorPalette=colorPalette, statePalette=statePalette,
                         clusteringMethod=clusteringMethod, orderClusters=orderClusters,
                         similarity=FALSE, orderGenes=orderGenes, returnPlot=returnPlot,
                         saveHeatmapTable=saveHeatmapTable, width=width, height=height)
    })


################################################################################
## plotGeneExpression
################################################################################

#' .checkParamPlotGeneExpression
#'
#' @description checks parameters of plotGeneExpression
#' 
#' @param theObject A scRNAseq object with the cluster similarity matrix got 
#' with calculateClustersSimilarity method.
#' @param geneName Name of the gene to highlight on the t-SNE plot.
#' @param graphsDirectory Name of the output subdirectory. 
#' Default is "pictures".
#' @param returnPlot If TRUE export to pdf, if FALSE export to png.
#' @param tSNEpicture Number of picture that you want to use for plotting. 
#' Please check "dataDirectory/tsnes" or
#' "dataDirectory/pictures/tSNE_pictures/clusters" to get the number
#' which corresponds to the number of files, it is usually from 1 to 14.
#' Default = 1
#' @param commentName Comment that you want to specify in the filename.
#' @param savePlot {Boolean, should the function export the plot to pdf or not.
#'  Default = TRUE
#' @param width Width of the plot. Default = 6.
#' @param height Height of the plot. Default = 5.
#'
#' @keywords internal
#' @noRd
.checkParamPlotGeneExpression <- function(theObject, geneName, graphsDirectory,
                                          returnPlot, tSNEpicture, commentName,
                                          savePlot, width, height){
    
    ## Verify the object 
    clustersSimilarityMatrix <- getClustersSimilarityMatrix(theObject)
    clusters <- colData(getSceNorm(theObject))$clusters
    nbrClusters <- length(unique(clusters))
    
    if (!isTRUE((ncol(clustersSimilarityMatrix) == nbrClusters)))
        stop(paste("You have to calculate the cluster similarity matrix",
                   "before plotting."))  
    
    
    ## Verify the geneName
    markers <- getClustersMarkers(theObject)$geneName 
    if(!geneName %in% markers)
        stop(paste("geneName should be a marker founded by ",
              "retrieveTopClustersMarkers method'. Please see the",
              "documentation about retrieveTopClustersMarkers method."))
    
    ## Verify the graphsDirectory
    if(!is.character(graphsDirectory))
        stop(paste("graphsDirectory should be a string,",
                   "the path of the graphs directory"))
    
    ## Verify tSNEpicture
    if (!is.numeric(tSNEpicture))
        stop("tSNEpicture should be a integer")
    
    ## Verify returnPlot
    if (!is.logical(returnPlot))
        stop("returnPlot should be a boolean.")
    
    ## Verify commentName
    if(!is.character(commentName))
        stop("commentName should be a string.")
    
    ## Verify savePlot
    if (!is.logical(savePlot))
        stop("savePlot should be a boolean.")  
    
    ## Verify width
    if (!is.numeric(width))
        stop("width should be a numeric.")
    
    ## Verify height
    if (!is.numeric(height))
        stop("height should be a numeric.")  
    
}			




#' plotGeneExpression
#' 
#' @description The function saves a t-SNE plot colored by expression of a  
#' given gene. Warning: filename with t-SNE results is hardcoded, so please 
#' don't rename the output file.
#'
#' @param theObject A scRNAseq object with the cluster similarity matrix got 
#' with calculateClustersSimilarity method.
#' @param geneName Name of the gene to highlight on the t-SNE plot.
#' @param graphsDirectory Name of the output subdirectory. 
#' Default is "pictures".
#' @param palette Color palette for the legend
#' @param returnPlot Boolean, should the function return a ggplot object or not.
#'  Default = FALSE
#' @param tSNEpicture Number of picture that you want to use for plotting. 
#' Please check "dataDirectory/tsnes" or
#' "dataDirectory/pictures/tSNE_pictures/clusters" to get the number
#' which corresponds to the number of files, it is usually from 1 to 14.
#' Default = 1
#' @param commentName Comment that you want to specify in the filename.
#' @param savePlot {Boolean, should the function export the plot to pdf or not.
#'  Default = TRUE
#' @param alpha Opacity of the points of the plot. Default = 1
#' @param limits Range of the gene expression shown in the legend. See details.
#' @param pointSize Size of the point. Default = 1.
#' @param width Width of the plot. Default = 6.
#' @param height Height of the plot. Default = 5.
#' @param ... Additional parameters to pass to pdf() functions.
#' 
#' @examples
#' experimentName <- "Bergiers"
#' countMatrix <- as.matrix(read.delim(file.path(
#' "tests/testthat/test_data/test_countMatrix.tsv")))
#' outputDirectory <- "./"
#' columnsMetaData <- read.delim(
#' file.path("extdata/Bergiers_colData_filtered.tsv"))
#' 
#' ## Create the initial object
#' scr <- scRNAseq(experimentName = experimentName, 
#'                 countMatrix     = countMatrix, 
#'                 species         = "mouse",
#'                 outputDirectory = outputDirectory)
#' 
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#' 
#' ## Compute the tSNE coordinates
#' scrTsne <- generateTSNECoordinates(scrNorm, cores=5)
#' 
#' ## Perform the clustering with dbScan
#' scrDbscan <- runDBSCAN(scrTsne, cores=5)
#' 
#' ## Compute the cell similarity matrix
#' scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=10, cores=4)
#' 
#' ## Calculate clusters similarity
#' scrCSM <- calculateClustersSimilarity(scrCCI)
#' 
#' ## t-SNE plot colored by expression of a  given gene.
#' plotGeneExpression(scrCSM)
#' 
#' @seealso calculateClustersSimilarity
#' @exportMethod 
#' @author 
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
setMethod(
    
    f = "plotGeneExpression",
    
    signature = "scRNAseq",
    
    definition = function(theObject, geneName, graphsDirectory="pictures", 
                          palette=c("grey","red", "#7a0f09", "black"), returnPlot=FALSE,
                          tSNEpicture=1, commentName="", savePlot=TRUE, alpha=1, limits=NA,
                          pointSize=1, width=6, height=5, ...){
        
        ## Verify parameters
        validObject(theObject)
        .checkParamPlotGeneExpression(theObject=theObject,
                                      geneName=geneName,
                                      graphsDirectory=graphsDirectory,
                                      returnPlot=returnPlot,
                                      tSNEpicture=tSNEpicture,
                                      commentName=commentName,
                                      savePlot=savePlot,
                                      width=width, 
                                      height=height)
        
        sceObject       <- getSceNorm(theObject)
        dataDirectory   <- getOutputDirectory(theObject)
        experimentName  <- getExperimentName(theObject)
        tSNEDirectory <- "tsnes"
        
        ## Plot all precalculated t-SNEs to show your clusters
        
        coldDatDf <- SummarizedExperiment::colData(sceObject)
        clustersNumber <- length(unique(coldDatDf$clusters))
        tmpList <- getTSNEList(theObject) 
        tSNECoords <- as.data.frame(getCoordinates(tmpList[[tSNEpicture]]))
        tSNECoords <- tSNECoords[coldDatDf$cellName, ]
        
        if(!geneName %in% rownames(Biobase::exprs(sceObject)))
            print("Gene is not found in expression matrix")
        
        stopifnot(all(rownames(tSNECoords) == colnames(sceObject)))
        tSNECoords$expression <- Biobase::exprs(sceObject)[geneName, ]
        
        if(isTRUE(all.equal(length(limits), 1)))
            limits <- c(min(tSNECoords$expression), max(tSNECoords$expression))
        
        
        if(savePlot){
            
            fileName <- paste(experimentName, "tSNE", clustersNumber, 
                              "clusters", geneName, commentName, "tSNEpicture", 
                              tSNEpicture, "_alpha", alpha,sep="_")
            fileName <- paste0(fileName, ".pdf")
            
            
            subdir <- file.path(dataDirectory, graphsDirectory)
            if(!file.exists(subdir))
                dir.create(subdir, showWarnings=F, recursive = TRUE)
            
            pdf(file.path(dataDirectory, graphsDirectory, fileName), 
                width=width, height=height, ...)
        }
        
        tmp <- ggplot2::ggplot(tSNECoords, aes(x=tSNECoords[,1], 
                                               y=tSNECoords[,2], color=expression)) + 
            geom_point(size=I(pointSize), alpha=alpha) + theme_bw() +
            scale_colour_gradientn(colours=alpha(colorRampPalette(palette)(100),
                                                 0.8), limits=limits) +
            ggtitle(geneName)
        
        print(tmp)
        
        if(savePlot)
            dev.off()
        
        if(returnPlot)
            return(tmp)
    })


#################
## plotClustersSimilarity
#################

#' .checkParamClustersSimilarity
#'
#' @description checks parameters of plotCellSimilarity
#' 
#' @param theObject A scRNAseq object with the cluster similarity matrix got 
#' with calculateClustersSimilarity method.
#' @param returnPlot If TRUE export to pdf, if FALSE export to png.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to TRUE, and forced to true if file is a pipe.
#' @param width Width of the plot. Default = 7.
#' @param height Height of the plot. Default = 6.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to TRUE, and forced to true if file is a pipe.
#' 
#' @keywords internal
#' @noRd
.checkParamClustersSimilarity <- function(theObject, returnPlot, width, height,
                                          onefile, fontsize){
    
    
    ## Verify the object 
    clustersSimilarityMatrix <- getClustersSimilarityMatrix(theObject)
    clusters <- colData(getSceNorm(theObject))$clusters
    nbrClusters <- length(unique(clusters))
    
    if (!isTRUE((ncol(clustersSimilarityMatrix) == nbrClusters)))
        stop(paste("You have to calculate the cluster similarity matrix",
             "before plotting."))
    
    ## Verify returnPlot
    if (!is.logical(returnPlot))
        stop("returnPlot should be a boolean.")  
    
    ## Verify width
    if (!is.numeric(width))
        stop("width should be a numeric.")
    
    ## Verify height
    if (!is.numeric(height))
        stop("height should be a numeric.")
    
    ## Verify onefile
    if (!is.logical(onefile))
        stop("onefile should be a boolean.")  
    
    ## Verify fontsize
    if (!is.numeric(fontsize))
        stop("fontsize should be a numeric.")
}




#' plotCellSimilarity
#' 
#' @description This function plots similarity matrix as a heatmap, so one can 
#' see similarity between parts of different clusters.
#'
#' @param theObject A scRNAseq object with the cluster similarity matrix got 
#' with calculateClustersSimilarity method.
#' @param colorPalette A vector of colors for clusters. Default = "default",
#'  see details.
#' @param statePalette A vector of colors for states or conditions. See details.
#' @param clusteringMethod Clustering method passed to hclust() function. 
#' See ?hclust for a list of method. Default = "ward.D2"
#' @param returnPlot Boolean indicating if the pHeatmap object should  be
#' returned by the function. Default = FALSE.
#' @param width Width of the plot. Default = 7.
#' @param height Height of the plot. Default = 6.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to TRUE, and forced to true if file is a pipe.
#' @param fontsize pheatmap parameter. Base fontsize for the plot
#' @param ... other parameters of the pdf() function.
#' 
#' @examples
#' experimentName <- "Bergiers"
#' countMatrix <- as.matrix(read.delim(file.path(
#' "tests/testthat/test_data/test_countMatrix.tsv")))
#' outputDirectory <- "./"
#' columnsMetaData <- read.delim(
#' file.path("extdata/Bergiers_colData_filtered.tsv"))
#' 
#' ## Create the initial object
#' scr <- scRNAseq(experimentName = experimentName, 
#'                 countMatrix     = countMatrix, 
#'                 species         = "mouse",
#'                 outputDirectory = outputDirectory)
#' 
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#' 
#' ## Compute the tSNE coordinates
#' scrTsne <- generateTSNECoordinates(scrNorm, cores=5)
#' 
#' ## Perform the clustering with dbScan
#' scrDbscan <- runDBSCAN(scrTsne, cores=5)
#' 
#' ## Compute the cell similarity matrix
#' scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=10, cores=4)
#' 
#' ## Calculate clusters similarity
#' scrCSM <- calculateClustersSimilarity(scrCCI)
#' 
#' ## Plot similarity matrix as a heatmap
#' plotClustersSimilarity(scrCSM)
#' 
#' 
#' @exportMethod 
#' @importFrom SummarizedExperiment colData
#' @importFrom stats hclust
#' @importFrom pheatmap pheatmap
#' 
#' @seealso calculateClustersSimilarity
#' @author 
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
setMethod(
    
    f = "plotClustersSimilarity",
    
    signature = "scRNAseq",
    
    definition = function(theObject, colorPalette="default", 
                          statePalette="default", clusteringMethod="ward.D2", 
                          returnPlot=FALSE, width=7, height=5.5, onefile=FALSE, fontsize=7.5,
                          ...){
        
        ## Verify parameters
        validObject(theObject)
        .checkParamClustersSimilarity(theObject, returnPlot, width, height,
                                      onefile, fontsize)
        
        clustersSimilarityMatrix <- getClustersSimilarityMatrix(theObject)
        sceObject       <- getSceNorm(theObject)
        dataDirectory   <- getOutputDirectory(theObject)
        experimentName  <- getExperimentName(theObject)
        
        coldDataDf <- SummarizedExperiment::colData(sceObject)
        clusters <- coldDataDf$clusters
        clustersNumber <- length(unique(clusters))
        clustersNames <- levels(clusters)
        dataDirectory <- dataDirectory
        experimentName <- experimentName
        graphsDirectory <- "pictures"
        
        distanceMatrix <- as.dist(sqrt((1-clustersSimilarityMatrix)/2))
        clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
        colDataSimilarity <- data.frame(clusters=clustersNames)
        rownames(colDataSimilarity) <- colDataSimilarity$clusters
        
        annotationColors <- .generateAnnotationColors(coldDataDf, colorPalette,
                                                      statePalette)
        
        tmpNCol <- ncol(clustersSimilarityMatrix)
        tmpNRow <- nrow(clustersSimilarityMatrix)
        main <- paste0("Clusters similarity matrix ", tmpNCol, " columns, ", 
                       tmpNRow, " rows.")
        
        pheatmapObject <- pheatmap::pheatmap(clustersSimilarityMatrix,
                                             annotation_col=colDataSimilarity,
                                             annotation_colors=annotationColors,
                                             cluster_cols=clusteringTree,
                                             cluster_rows=clusteringTree,
                                             fontsize=fontsize, main=main)
        
        message("\nSaving a heatmap with the cluster similarity matrix.")
        tmpFileP <- paste(experimentName,"clusters_similarity", clustersNumber,
                          "clusters.pdf", sep="_")
        
        
        subdir <- file.path(dataDirectory, graphsDirectory)
        if(!file.exists(subdir))
            dir.create(subdir, showWarnings=F, recursive = TRUE)
        
        pdf(file.path(dataDirectory, graphsDirectory, tmpFileP), width=width, 
            height=height, onefile=onefile, ...)
        grid::grid.newpage()
        grid::grid.draw(pheatmapObject$gtable)
        dev.off()
        
        if(returnPlot)
            return(pheatmapObject)			
    })
