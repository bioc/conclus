################################################################################
## plotCellSimilarity
################################################################################

#' .saveCellSim
#'
#'
#' @description
#' Save the cell similarity heatmap in pdf or png format in the directory
#' defined in theObject (?getOutputDirectory) and in the sub-directory pictures.
#'
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized (see ?normaliseCountMatrix), tSNE were calculated (see
#' ?generateTSNECoordinates), dbScan was run (see ?runDBSCAN), cells were
#' clustered (see ?clusterCellsInternal), as clusters themselves (see
#' ?calculateClustersSimilarity).
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param colDf A data frame with information about cells.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 7.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 6.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param plotPDF If TRUE export heatmap in pdf format, if FALSE export it in
#' png format. Default=TRUE.
#' @param pheatmapObject Object containing the pheatmap
#'
#' @keywords internal
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @noRd
.saveCellSim <- function(theObject, onefile, colDf, width, height, widthPNG,
        heightPNG, plotPDF, pheatmapObject){

    message("\nSaving a heatmap with the cell similarity matrix.")

    experimentName <- getExperimentName(theObject)
    dataDirectory  <- getOutputDirectory(theObject)
    graphsDirectory <- "pictures"
    subdir <- file.path(dataDirectory, graphsDirectory)
    clustersNumber <- length(unique(colDf$clusters))
    filePath <- file.path(subdir, paste(experimentName,
                    "cells_correlation", clustersNumber, "clusters",
                    sep="_"))

    if(!file.exists(subdir))
        dir.create(subdir, showWarnings=FALSE, recursive = TRUE)

    if(plotPDF)
        pdf(file=paste0(filePath, ".pdf"), width=width, height=height,
                onefile=onefile)
    else
        png(filename=paste0(filePath, ".png"), width=widthPNG, height=heightPNG,
                type="cairo")

    grid::grid.newpage()
    grid::grid.draw(pheatmapObject$gtable)
    dev.off()
}


#' .generateAnnotationColors
#'
#' @description
#' Generates annotation_colors used by pheatmap(). It is a list that specifies
#' annotation_row and annotation_col track colors manually.
#'
#' @param colDf A data frame with information about cells.
#' @param colorPalette A vector of colors for clusters.
#' See details of ?plotCellSimilarity.
#' @param statePalette A vector of colors for states or conditions.
#' See details of ?plotCellSimilarity.
#'
#' @keywords internal
#' @return A list of colors.
#' @include sharedInternals.R
#' @noRd
.generateAnnotationColors <- function(colDf, colorPalette, statePalette){

    clusters <- levels(colDf$clusters)
    states <- unique(colDf$state)
    clusterNumber <- length(unique(colDf$clusters))

    colorAnnotationClusters <- .choosePalette(colorPalette, clusterNumber)
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
#' @param colDf A data frame with information about cells.
#' @param mtx Expression count matrix
#' @param clusteringMethod Clustering method passed to hclust() function.
#' See ?hclust for a list of method. Default = "ward.D2"
#'
#' @keywords internal
#' @return Cells ordered by clustering results.
#'
#' @importFrom stats hclust
#' @importFrom stats dist
#' @noRd
.orderCellsInCluster <- function(cluster, colDf, mtx,
        clusteringMethod="ward.D2"){

    cells <- colDf[colDf$clusters == cluster, ]$cellName

    if(length(cells) > 2){
        tree <- hclust(dist(t(mtx[, cells])), method=clusteringMethod)
        return(cells[tree$order])
    }else
        return(cells)
}


#' .checkParamCellSimilaritySub
#'
#' @description Check parameters of plotCellSimilarity
#' 
#' @param showColnames pheatmap parameter. Boolean specifying if column names
#' are displayed. Default=FALSE.
#' @param fontsize pheatmap parameter. Base fontsize for the plot. Default=7.5.
#' @param fontsizeRow pheatmap parameter. Fontsize for rownames. Default=0.03.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param silentPlot If TRUE, does not plot the pheatmap. Default=FALSE.
#' @param returnPlot Boolean indicating if the pHeatmap object should  be
#' returned by the function. Default = TRUE.
#' @param savePlot If TRUE, the heatmap is saved in the directory defined in
#' theObject (?getOutputDirectory) and in the sub-directory "pictures".
#'
#' @keywords internal
#' @noRd
.checkParamCellSimilaritySub <- function(showColnames, fontsize, fontsizeRow,
        widthPNG, heightPNG, silentPlot, returnPlot, savePlot){
    
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
    
    ## Verify silentPlot
    if (!is.logical(silentPlot))
        stop("silentPlot should be a boolean.")
    
    if(silentPlot && !returnPlot && !savePlot)
        stop("You do not plot, neither save the heatmap or return the object.",
                " Nothing will happen. You should either plot the results, ",
                "return the object or save the heatmap.")
}

#' .checkParamCellSimilarity
#'
#' @description Check parameters of plotCellSimilarity
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized (see ?normaliseCountMatrix), tSNE were calculated (see
#' ?generateTSNECoordinates), dbScan was run (see ?runDBSCAN), cells were
#' clustered (see ?clusterCellsInternal), as clusters themselves (see
#' ?calculateClustersSimilarity).
#' @param orderClusters If TRUE, clusters in the similarity matrix of cells will
#' be ordered by name. Default = FALSE.
#' @param savePlot If TRUE, the heatmap is saved in the directory defined in
#' theObject (?getOutputDirectory) and in the sub-directory "pictures".
#' @param plotPDF If TRUE export heatmap in pdf format, if FALSE export it in
#' png format. Default=TRUE.
#' @param returnPlot Boolean indicating if the pHeatmap object should  be
#' returned by the function. Default = TRUE.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 7.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 6.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param showRowNames pheatmap parameter. Boolean specifying if row names
#' are displayed. Default=FALSE.
#' @param showColnames pheatmap parameter. Boolean specifying if column names
#' are displayed. Default=FALSE.
#' @param fontsize pheatmap parameter. Base fontsize for the plot. Default=7.5.
#' @param fontsizeRow pheatmap parameter. Fontsize for rownames. Default=0.03.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param silentPlot If TRUE, does not plot the pheatmap. Default=FALSE.
#'
#' @keywords internal
#' @noRd
.checkParamCellSimilarity <-function(theObject, orderClusters, savePlot,
        plotPDF, returnPlot, width, height, onefile, showRowNames, showColnames,
        fontsize, fontsizeRow, widthPNG, heightPNG, silentPlot){

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

    ## Verify savePlot
    if(!is.logical(savePlot))
        stop("savePlot should be a boolean.")

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

    .checkParamCellSimilaritySub(showColnames, fontsize, fontsizeRow, widthPNG,
            heightPNG, silentPlot, returnPlot, savePlot)
}




#' .plotCellSimilarityOrderClusters
#'
#' @description
#' Orders the clusters to generate a cell similarity pheatmap.
#' 
#' @param orderClusters If TRUE, clusters in the similarity matrix of cells will
#' be ordered by name. Default = FALSE.
#' @param colDf A data frame with information about cells.
#' @param expressionMatrix The count matrix retrieved with ?Biobase::exprs.
#' @param clusteringMethod Clustering method passed to hclust() function.
#' See ?hclust for a list of method. Default = "ward.D2"
#' @param cellsSimilarityMatrix Matrix retrieved with ?getCellsSimilarityMatrix.
#' @param clusterCols Boolean indicating if the clustering should be performed
#' on columns.
#' @param clusterRows Boolean indicating if the clustering should be performed
#' on rows.
#'
#' @keywords internal
#' @importFrom stats as.dist
#' @importFrom stats hclust
#' @noRd
#' @return The update cellsSimilarityMatrix (if orderClusters is TRUE) and the
#' booleans clusterCols and clusterRows.
.plotCellSimilarityOrderClusters <- function(orderClusters, colDf, 
        expressionMatrix, clusteringMethod, cellsSimilarityMatrix, clusterCols,
        clusterRows){
    
    if(orderClusters){
        # Ordering expressionMatrixrix
        newOrder <- lapply(levels(colDf$clusters), function(cluster){
                    .orderCellsInCluster(cluster, colDf, expressionMatrix,
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
    
    return(list(cellsSimilarityMatrix, clusterCols, clusterRows))
}


#' .plotCellSimilarityHeatmap
#'
#' @description
#' Generate a cell similarity pheatmap.
#' 
#' @param colDf A data frame with information about cells.
#' @param colorPalette A vector of colors for clusters. Default = "default",
#' see details.
#' @param statePalette A vector of colors for states or conditions. See details.
#' @param cellsSimilarityMatrix Matrix retrieved with ?getCellsSimilarityMatrix.
#' @param showColnames pheatmap parameter. Boolean specifying if column names
#' are displayed. Default=FALSE.
#' @param showRowNames pheatmap parameter. Boolean specifying if row names
#' are displayed. Default=FALSE.
#' @param fontsizeRow pheatmap parameter. Fontsize for rownames. Default=0.03.
#' @param clusterCols Boolean indicating if the clustering should be performed
#' on columns.
#' @param clusterRows Boolean indicating if the clustering should be performed
#' on rows.
#' @param fontsize pheatmap parameter. Base fontsize for the plot. Default=7.5.
#' @param silentPlot If TRUE, does not plot the heatmap. Default=FALSE.
#' @param savePlot If TRUE, the heatmap is saved in the directory defined in
#' theObject (?getOutputDirectory) and in the sub-directory "pictures".
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized (see ?normaliseCountMatrix), tSNE were calculated (see
#' ?generateTSNECoordinates), dbScan was run (see ?runDBSCAN), cells were
#' clustered (see ?clusterCellsInternal), as clusters themselves (see
#' ?calculateClustersSimilarity).
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 7.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 6.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param plotPDF If TRUE export heatmap in pdf format, if FALSE export it in
#' png format. Default=TRUE.
#'
#' @importFrom pheatmap pheatmap
#' @keywords internal
#' @noRd
#' @return A pheatmap object.
.plotCellSimilarityHeatmap <- function(colDf, colorPalette, statePalette, 
        cellsSimilarityMatrix, showColnames, showRowNames, fontsizeRow, 
        clusterCols, clusterRows, fontsize, silentPlot, savePlot, theObject, 
        onefile, width, height, widthPNG, heightPNG, plotPDF){
    
    annotationColors <- .generateAnnotationColors(colDf, colorPalette,
            statePalette)
    
    columnsToPlot <- switch(is.null(colDf$state) + 1, c("clusters",
                    "state"), c("clusters"))
    annotationCol <- as.data.frame(colDf[columnsToPlot])
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
            main=mainTitle,
            silent=silentPlot)
    
    if(savePlot)
        .saveCellSim(theObject, onefile, colDf, width, height, widthPNG,
                heightPNG, plotPDF, pheatmapObject)
    
    return(pheatmapObject)
}


#' plotCellSimilarity
#'
#' @description This function plots similarity matrix as a heatmap, so one can
#' see similarity between parts of different clusters.
#'
#' @usage
#' plotCellSimilarity(theObject, colorPalette="default",
#'             statePalette="default", clusteringMethod="ward.D2",
#'             orderClusters=FALSE, savePlot=FALSE, plotPDF=TRUE,
#'             returnPlot=FALSE, width=7, height=6, onefile=FALSE,
#'             showRowNames=FALSE, showColnames=FALSE, fontsize=7.5,
#'             fontsizeRow=0.03, widthPNG=800, heightPNG=750, silentPlot=FALSE)
#'
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized (see ?normaliseCountMatrix), tSNE were calculated (see
#' ?generateTSNECoordinates), dbScan was run (see ?runDBSCAN), cells were
#' clustered (see ?clusterCellsInternal), as clusters themselves (see
#' ?calculateClustersSimilarity).
#' @param colorPalette A vector of colors for clusters. Default = "default",
#' see details.
#' @param statePalette A vector of colors for states or conditions. See details.
#' @param clusteringMethod Clustering method passed to hclust() function.
#' See ?hclust for a list of method. Default = "ward.D2"
#' @param orderClusters If TRUE, clusters in the similarity matrix of cells will
#' be ordered by name. Default = FALSE.
#' @param savePlot If TRUE, the heatmap is saved in the directory defined in
#' theObject (?getOutputDirectory) and in the sub-directory "pictures".
#' @param plotPDF If TRUE export heatmap in pdf format, if FALSE export it in
#' png format. Default=TRUE.
#' @param returnPlot Boolean indicating if the pHeatmap object should  be
#' returned by the function. Default = TRUE.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 7.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 6.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param showRowNames pheatmap parameter. Boolean specifying if row names
#' are displayed. Default=FALSE.
#' @param showColnames pheatmap parameter. Boolean specifying if column names
#' are displayed. Default=FALSE.
#' @param fontsize pheatmap parameter. Base fontsize for the plot. Default=7.5.
#' @param fontsizeRow pheatmap parameter. Fontsize for rownames. Default=0.03.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param silentPlot If TRUE, does not plot the heatmap. Default=FALSE.
#'
#' @details
#' colorPalette/statePalette -- A vector of colors for clusters/states or
#' 'default' value. If 'default' is selected, the number of clusters is limited
#' to 16. If an error message is thrown, re-run the function with your own color
#' vector.
#'
#' @aliases plotCellSimilarity
#' @rdname plotCellSimilarity-scRNAseq
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
#' columnID="cell_ID")
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
#' ## Plot the heatmap of the similarity matrix
#' plotCellSimilarity(scrCSM)
#'
#' @return A pheatmap object of the similarity heatmap if returnPlot is TRUE.
#'
#' @seealso calculateClustersSimilarity  plotClusteredTSNE plotCellHeatmap
#' plotGeneExpression plotClustersSimilarity
#'
#' @exportMethod plotCellSimilarity
#' @importFrom SummarizedExperiment colData
#' @importFrom pheatmap pheatmap
#' @importFrom Biobase exprs
#' @importFrom stats as.dist
#' @importFrom stats hclust
#' @importFrom methods validObject
#'
#' @author
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
setMethod(

    f = "plotCellSimilarity",

    signature = "scRNAseq",

    definition = function(theObject, colorPalette="default",
            statePalette="default", clusteringMethod="ward.D2",
            orderClusters=FALSE, savePlot=FALSE, plotPDF=TRUE, returnPlot=FALSE,
            width=7, height=6, onefile=FALSE, showRowNames=FALSE,
            showColnames=FALSE, fontsize=7.5, fontsizeRow=0.03, widthPNG=800,
            heightPNG=750, silentPlot=FALSE){

        ## Verify parameters
        validObject(theObject)
        .checkParamCellSimilarity(theObject, orderClusters, savePlot, plotPDF,
                returnPlot, width, height, onefile, showRowNames, showColnames,
                fontsize, fontsizeRow, widthPNG, heightPNG, silentPlot)

        sceObject <- getSceNorm(theObject)
        expressionMatrix <-  Biobase::exprs(sceObject)
        cellsSimilarityMatrix <- getCellsSimilarityMatrix(theObject)
        colDf <- SummarizedExperiment::colData(sceObject)

        ## Ordering clusters
        result <- .plotCellSimilarityOrderClusters(orderClusters, colDf, 
                expressionMatrix, clusteringMethod, cellsSimilarityMatrix, 
                clusterCols, clusterRows)
        cellsSimilarityMatrix <- result[[1]]
        clusterCols <- result[[2]]
        clusterRows <- result[[3]]
        
        ## Plotting the cell similarity heatmap
        pheatmapObject <- .plotCellSimilarityHeatmap(colDf, colorPalette, 
                statePalette, cellsSimilarityMatrix, showColnames, showRowNames,
                fontsizeRow, clusterCols, clusterRows, fontsize, silentPlot, 
                savePlot, theObject, onefile, width, height, widthPNG, 
                heightPNG, plotPDF)
        
        if(returnPlot)
            return(pheatmapObject)
    })


################################################################################
## plotClusteredTSNE
################################################################################


#' .saveTSNEPlot
#'
#' @description Write the tSNE plots to the output directory obtained with
#' ?getOutputDirectory in the sub-directories 'pictures/graphsTSNEDirectory'.
#'
#' @param tSNEList List of tSNE objects obtained with ?getTSNEList.
#' @param tSNEplots List of ggplot objects representing the tSNEs.
#' @param plotPDF If TRUE export heatmap in pdf format, if FALSE export it in
#' png format. Default=TRUE.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 6.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 5.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE, and forced to true if file is a pipe.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param outputDir Path to the output directory.
#'
#' @keywords internal
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @noRd
.saveTSNEPlot <- function(tSNEList, tSNEplots, plotPDF, width, height, onefile,
        widthPNG, heightPNG, outputDir){

    tSNENameList <- lapply(tSNEList, function(currentTsne)
                return(getName(currentTsne)))


    invisible(mapply(function(currentTsne, currentName, plotPDF, width, height,
                            onefile, widthPNG, heightPNG, outputDir){

                        filePath <- file.path(outputDir, currentName)

                        if(plotPDF)
                            pdf(file=paste0(filePath, ".pdf"), width=width,
                                    height=height,onefile=onefile)
                        else
                            png(filename=paste0(filePath,".png"),
                                    width=widthPNG, height=heightPNG,
                                    type="cairo")

                        print(currentTsne)
                        dev.off()

                    }, tSNEplots, tSNENameList, MoreArgs=list(plotPDF, width,
                            height, onefile, widthPNG, heightPNG, outputDir)))
}




#' .computePlotList
#'
#' @description Generates a list of ggplot objects for each tSNE.
#'
#' @param tSNEList List of tSNE objects obtained with ?getTSNEList.
#' @param sceObject SingleCellExperiment object obtained with ?getSceNorm and
#' containing the normalized count matrix.
#' @param columnName Name of the column to color the t-SNE with.
#' Possible values are clusters, noColor, or state.
#' @param colorPalette A vector of colors for clusters. Default = "default",
#' see details of plotClusteredTSNE below.
#'
#' @return A list of ggplot objects.
#' @keywords internal
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 scale_color_manual
#' @noRd
.computePlotList <- function(tSNEList, sceObject, columnName, colorPalette){

    resultList <- lapply(tSNEList, function(currentTsne, sceObj, columnName,
                    colorPalette){

                tSNEres <- as.data.frame(getCoordinates(currentTsne))
                coordinatesName <- getName(currentTsne)
                colDatDf <- SummarizedExperiment::colData(sceObj)

                ## Get coordinates of filtered cells
                cellVec <- colDatDf$cellName
                tSNEres <- tSNEres[rownames(tSNEres) %in% cellVec, ]

                ## Add new columns
                if(!isTRUE(all.equal(columnName, "noColor")))
                    tSNEres[columnName] <- factor(colDatDf[, columnName])

                if(isTRUE(all.equal(columnName, "noColor")))
                    tmp <- ggplot2::ggplot(tSNEres,
                                    aes_string(x=names(tSNEres)[1],
                                            y=names(tSNEres)[2])) +
                            geom_point(size=I(1)) + theme_bw()

                else
                    tmp <- ggplot2::ggplot(tSNEres, aes_string(
                                            x=names(tSNEres)[1],
                                            y=names(tSNEres)[2],
                                            color=columnName)) +
                            geom_point(size=I(1)) +
                            scale_color_manual(values=colorPalette) +
                            theme_bw()

                return(tmp)

            }, sceObject, columnName, colorPalette)

    return(resultList)
}




#' .checkParamPlotTSNESub
#'
#' @description Check parameters of plotClusteredTSNE
#'
#' @param plotPDF If TRUE export heatmap in pdf format, if FALSE export it in
#' png format. Default=TRUE.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' "pictures/tSNE_pictures".
#' @param silentPlot If TRUE, does not plot the pheatmap. Default=FALSE.
#' @param returnPlot Boolean indicating if the pHeatmap object should  be
#' returned by the function. Default = FALSE.
#' @param savePlot If TRUE, the heatmap is saved in the directory defined in
#' theObject (?getOutputDirectory) and in the sub-directory
#' @param tSNENb Give the number of the tSNE to plot. If NA, all tSNE solutions
#' are plotted (14 tSNE by default). Default=NA.
#' @param PCs A vector of first principal components. For example, to take
#' ranges 1:5 and 1:10 write c(5, 10). Default = c(4, 6, 8, 10, 20, 40, 50).
#' See ?generateTSNECoordinates for details.
#' @param perplexities Numeric scalar defining the perplexity parameter,
#' see ?Rtsne and ?generateTSNECoordinates for more details. Default = c(30, 40)
#'
#' @keywords internal
#' @noRd
.checkParamPlotTSNESub <- function(plotPDF, widthPNG, heightPNG, silentPlot, 
        returnPlot, savePlot, tSNENb, PCs, perplexities){
    
    ## Verify plotPDF
    if(!is.logical(plotPDF))
        stop("plotPDF should be a boolean.")
    
    ## Verify widthPNG
    if (!is.numeric(widthPNG))
        stop("widthPNG should be a numeric.")
    
    ## Verify heightPNG
    if (!is.numeric(heightPNG))
        stop("heightPNG should be a numeric.")
    
    ## Verify silentPlot
    if (!is.logical(silentPlot))
        stop("silentPlot should be a boolean.")
    
    
    if(silentPlot && !returnPlot && !savePlot)
        stop("You do not plot, neither save the heatmap or return the object.",
                " Nothing will happen. You should either plot the results, ",
                "return the object or save the heatmap.")
    
    if(!is.na(tSNENb) && !is.numeric(tSNENb))
        stop("tSNENb should be a numeric.")
    
    if(!is.na(tSNENb) && (tSNENb > (length(PCs)*length(perplexities))))
        stop("The chosen tSNENb should be smaller than PCs x perplexities.")
}



#' .checkParamPlotTSNE
#'
#' @description Check parameters of plotClusteredTSNE
#'
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized (see ?normaliseCountMatrix), tSNE were calculated (see
#' ?generateTSNECoordinates), dbScan was run (see ?runDBSCAN), cells were
#' clustered (see ?clusterCellsInternal), as clusters themselves (see
#' ?calculateClustersSimilarity).
#' @param PCs A vector of first principal components. For example, to take
#' ranges 1:5 and 1:10 write c(5, 10). Default = c(4, 6, 8, 10, 20, 40, 50).
#' See ?generateTSNECoordinates for details.
#' @param perplexities Numeric scalar defining the perplexity parameter,
#' see ?Rtsne and ?generateTSNECoordinates for more details. Default = c(30, 40)
#' @param columnName Name of the column to color the t-SNE with.
#' Possible values are clusters, noColor, or state.
#' @param returnPlot Boolean indicating if the pHeatmap object should  be
#' returned by the function. Default = FALSE.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 6.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 5.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param silentPlot If TRUE, the plots are not displayed on the current device.
#' Default=FALSE.
#' @param savePlot If TRUE, the heatmap is saved in the directory defined in
#' theObject (?getOutputDirectory) and in the sub-directory
#' @param plotPDF If TRUE export heatmap in pdf format, if FALSE export it in
#' png format. Default=TRUE.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' "pictures/tSNE_pictures".
#' @param silentPlot If TRUE, does not plot the pheatmap. Default=FALSE.
#'
#' @keywords internal
#' @noRd
.checkParamPlotTSNE <- function(theObject, PCs, perplexities, columnName,
                                returnPlot, width, height, onefile, silentPlot,
                                savePlot, plotPDF, widthPNG, heightPNG, tSNENb){

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
        !isTRUE(all.equal(columnName, "state")))
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

    ## Verify savePlot
    if (!is.logical(savePlot))
        stop("savePlot should be a boolean.")

    .checkParamPlotTSNESub(plotPDF, widthPNG, heightPNG, silentPlot, 
            returnPlot, savePlot, tSNENb, PCs, perplexities)
}



#' .createTSNEDir
#'
#' @description Create subfolder for TSNE plots and returns the colorPalette.
#'
#'
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized (see ?normaliseCountMatrix), tSNE were calculated (see
#' ?generateTSNECoordinates), dbScan was run (see ?runDBSCAN), cells were
#' clustered (see ?clusterCellsInternal), as clusters themselves (see
#' ?calculateClustersSimilarity).
#' @param columnName Name of the column to color the t-SNE with.
#' Possible values are clusters, noColor, or state.
#' @param sceObject A SingleCellExperiment object with your experiment.
#' @param colorPalette A vector of colors for clusters. Default = "default",
#' see details of plotClusteredTSNE below.
#'
#' @return List containing path to the subdirectory
#' 'pictures/graphsTSNEDirectory' and the color palette
#' @importFrom SummarizedExperiment colData
#' @include sharedInternals.R
#' @keywords internal
#' @noRd
.createTSNEDir <- function(theObject, columnName, sceObject, colorPalette){

    dataDirectory  <- getOutputDirectory(theObject)
    graphsDirectory     <- "pictures"
    graphsTSNEDirectory <- "tSNE_pictures"

    if(isTRUE(all.equal(columnName, "noColor")))
        numberElements <- NULL
    else{
        coldata <- SummarizedExperiment::colData(sceObject)
        nb <- unique(coldata[, columnName])
        numberElements <- length(nb)
        colorPalette <- .choosePalette(colorPalette, numberElements)
    }


    outputDir <- file.path(dataDirectory, graphsDirectory,
                            graphsTSNEDirectory, paste("tSNE", numberElements,
                                    columnName, sep="_"))

    if(!file.exists(outputDir))
        dir.create(outputDir, showWarnings=FALSE, recursive = TRUE)

    return(list(outputDir, colorPalette))
}


#' .plotAndSaveTSNE
#'
#' @description 
#' Plots and/or save the tSNE.
#'
#' @param silentPlot If TRUE, the plots are not displayed on the current device.
#' Default=FALSE.
#' @param tSNENb Give the number of the tSNE to plot. If NA, all tSNE solutions
#' are plotted (14 tSNE by default). Default=NA.
#' @param tSNEplots Result of the .computePlotList function.
#' @param savePlot If TRUE, the heatmap is saved in the directory defined in
#' theObject (?getOutputDirectory) and in the sub-directory
#' "pictures/tSNE_pictures".
#' @param tSNEList List of tSNE objects obtained with ?getTSNEList.
#' @param plotPDF If TRUE export heatmap in pdf format, if FALSE export it in
#' png format. Default=TRUE.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 6.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 5.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param outputdir Directory to write the tSNE to and retrieved from 
#' .createTSNEDir.
#' 
#' @importFrom grDevices dev.new
#' @noRd
.plotAndSaveTSNE <- function(silentPlot, tSNENb, tSNEplots, savePlot, tSNEList, 
        plotPDF, width, height, onefile, widthPNG, heightPNG, outputdir){
    
    ## Plotting tSNE
    if(!silentPlot){
        if(is.na(tSNENb))
            ## !!! Cette partie de code crée des fichiers Rplots.pdf !!!
            ## !!! en dehors du dossier YourOutputDirectory !!!
            ## !!!  alors que savePlot = False !!!
            invisible(lapply(tSNEplots, function(currentTSNE){
                                print(currentTSNE)
                                dev.new()}))
        else
            print(tSNEplots[[tSNENb]])
    }
    
    ## Saving tSNEs
    if(savePlot)
        .saveTSNEPlot(tSNEList, tSNEplots, plotPDF, width, height, onefile,
                widthPNG, heightPNG, outputdir)
}

#' plotClusteredTSNE
#'
#' @description Plot t-SNE generated with different PCs and perplexities.
#' It can also use a coloring scheme by clusters or states.
#'
#' @usage plotClusteredTSNE(theObject, colorPalette="default",
#'             PCs=c(4, 6, 8, 10, 20, 40, 50), perplexities=c(30, 40),
#'             columnName="clusters", savePlot=FALSE, plotPDF=TRUE,
#'             returnPlot=FALSE, width=6, height=5, onefile=FALSE, widthPNG=800,
#'             heightPNG=750, silentPlot=FALSE, tSNENb=NA)
#'
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized (see ?normaliseCountMatrix), tSNE were calculated (see
#' ?generateTSNECoordinates), dbScan was run (see ?runDBSCAN), cells were
#' clustered (see ?clusterCellsInternal), as clusters themselves (see
#' ?calculateClustersSimilarity).
#' @param colorPalette A vector of colors for clusters. Default = "default",
#' see details.
#' @param PCs A vector of first principal components. For example, to take
#' ranges 1:5 and 1:10 write c(5, 10). Default = c(4, 6, 8, 10, 20, 40, 50).
#' See ?generateTSNECoordinates for details.
#' @param perplexities Numeric scalar defining the perplexity parameter,
#' see ?Rtsne and ?generateTSNECoordinates for more details. Default = c(30, 40)
#' @param columnName Name of the column to color the t-SNE with.
#' Possible values are clusters, noColor, or state.
#' @param savePlot If TRUE, the heatmap is saved in the directory defined in
#' theObject (?getOutputDirectory) and in the sub-directory
#' "pictures/tSNE_pictures".
#' @param plotPDF If TRUE export heatmap in pdf format, if FALSE export it in
#' png format. Default=TRUE.
#' @param returnPlot Boolean indicating if the pHeatmap object should  be
#' returned by the function. Default = FALSE.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 6.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 5.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param silentPlot If TRUE, the plots are not displayed on the current device.
#' Default=FALSE.
#' @param tSNENb Give the number of the tSNE to plot. If NA, all tSNE solutions
#' are plotted (14 tSNE by default). Default=NA.
#'
#' @aliases plotClusteredTSNE
#' @rdname plotClusteredTSNE-scRNAseq
#'
#' @details
#' colorPalette -- A vector of colors for clusters/states or
#' 'default' value. If 'default' is selected, the number of clusters is limited
#' to 16. If an error message is thrown, re-run the function with your own color
#' vector.
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
#' columnID="cell_ID")
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
#' ## Plot the heatmap of the similarity matrix
#' plotClusteredTSNE(scrCSM)
#'
#' @return A list of ggplot objects if returnPlot is TRUE.
#'
#' @seealso calculateClustersSimilarity plotCellSimilarity plotCellHeatmap
#' plotGeneExpression plotClustersSimilarity
#'
#' @exportMethod plotClusteredTSNE
#' @importFrom grDevices dev.new
#' @importFrom methods validObject
#'
#' @author
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.

setMethod(

    f = "plotClusteredTSNE",

    signature = "scRNAseq",

    definition = function(theObject, colorPalette="default",
            PCs=c(4, 6, 8, 10, 20, 40, 50), perplexities=c(30, 40),
            columnName="clusters", savePlot=FALSE, plotPDF=TRUE,
            returnPlot=FALSE, width=6, height=5, onefile=FALSE, widthPNG=800,
            heightPNG=750, silentPlot=FALSE, tSNENb=NA){

        ## Verify parameters
        validObject(theObject)
        .checkParamPlotTSNE(theObject, PCs, perplexities, columnName,
                returnPlot, width, height, onefile, silentPlot,
                savePlot, plotPDF, widthPNG, heightPNG, tSNENb)

        ## Creating output folder
        sceObject <- getSceNorm(theObject)
        resultDir <- .createTSNEDir(theObject, columnName, sceObject,
                                    colorPalette)
        outputdir <- resultDir[[1]]
        colorPalette <- resultDir[[2]]

        # plots picture based on t-SNE coordinates from
        # generateTSNECoordinates() and clustering results
        # from clusterCellsInternal() or runClustering()

        ### Plot all precalculated tSNEs to show your clusters
        PCA <- rep(PCs, length(perplexities))
        perp <- rep(perplexities, each=length(PCs))
        tSNEList <- getTSNEList(theObject)
        totalLength <- length(PCs)*length(perplexities)

        if(!isTRUE(all.equal(length(tSNEList), totalLength)))
            stop("The number of elements of TSNEList is not equal to ",
                "PCs x perplexities. Contact the developper.")

        ## Generating the list of clustered tSNEs
        tSNEplots <- .computePlotList(tSNEList, sceObject, columnName,
                colorPalette)

        .plotAndSaveTSNE(silentPlot, tSNENb, tSNEplots, savePlot, tSNEList, 
                plotPDF, width, height, onefile, widthPNG, heightPNG, outputdir)
        
        if(returnPlot)
            return(tSNEplots)

    })


#################
## plotCellHeatmap
#################


#' .saveHeatmap
#'
#' @description Save the heatmap in pdf or png format in the directory defined
#' in theObject (?getOutputDirectory) and in the sub-directory 'pictures'.
#'
#' @param theObject A scRNAseq object with the cluster similarity matrix
#' obtained with ?calculateClustersSimilarity method and the top markers
#' obtained with ?retrieveTopClustersMarkers.
#' @param plotPDF If TRUE, the heatmap is saved in pdf format and in png
#' otherwise. Default = TRUE.
#' @param fileName Name of the output file to which the heatmap is saved.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 10.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 8.5.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param pheatmapObject Object containing the cell heatmap.
#'
#' @keywords internal
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @noRd
.saveHeatmap <- function(theObject, plotPDF, fileName, width, height, onefile,
        widthPNG, heightPNG, pheatmapObject){

    subdir <- file.path(getOutputDirectory(theObject), "pictures")
    if(!file.exists(subdir))
        dir.create(subdir, showWarnings=FALSE, recursive = TRUE)

    if(plotPDF)
        pdf(file=file.path(subdir, paste0(fileName, ".pdf")),
                width=width, height=height, onefile=onefile)
    else
        png(filename=file.path(subdir, paste0(fileName, ".png")),
                width=widthPNG, height=heightPNG, type="cairo")
    grid::grid.newpage()
    grid::grid.draw(pheatmapObject$gtable)
    dev.off()
}


#' .callOrderGenes
#'
#' @description Call the function .orderGenesInCluster on each cluster.
#'
#' @param colDf A data frame with information about cells.
#' @param markersClusters Data frame containing two columns geneName and
#' clusters.
#' @param expressionMatrix The count matrix retrieved with ?Biobase::exprs.
#' @param clusteringMethod Clustering method passed to hclust() function.
#' See ?hclust for a list of method. Default = "ward.D2".
#'
#' @keywords internal
#' @return A vector of ordered cells names
#'
#' @noRd
.callOrderGenes <- function(colDf, markersClusters, expressionMatrix,
        clusteringMethod){
    mat <- unname(unlist(lapply(levels(colDf$clusters), function(cluster)
        .orderGenesInCluster(cluster,
                            markersClusters,
                            expressionMatrix,
                            clusteringMethod))))
    suppressWarnings(mat <- matrix(mat, ncol = length(levels(colDf$clusters))))
    return(mat)
}


#' .callOrderCells
#'
#' @description Call the function .orderCellsInCluster on each cluster.
#'
#' @param colDf A data frame with information about cells.
#' @param expressionMatrix The count matrix retrieved with ?Biobase::exprs.
#' @param clusteringMethod Clustering method passed to hclust() function.
#' See ?hclust for a list of method. Default = "ward.D2".
#'
#' @keywords internal
#' @return A vector of ordered cells names
#'
#' @noRd
.callOrderCells <- function(colDf, expressionMatrix, clusteringMethod){
    
    vec <-  unname(unlist(lapply(levels(colDf$clusters), function(cluster)
        .orderCellsInCluster(cluster, colDf,
                            expressionMatrix,
                            clusteringMethod))))
    return(vec)
}


#' .orderGenesInCluster
#'
#' @description Order cells according to clustering results.
#' Uses for ordering matrix to further plot it with pheatmap()
#'
#' @param cluster Label of the current cluster considered in the sapply of
#' callOrderGenes.
#' @param markersClusters Data frame containing two columns geneName and
#' clusters.
#' @param mtx Expression matrix.
#' @param clusteringMethod Clustering method passed to hclust() function.
#' See ?hclust for a list of method. Default = "ward.D2".
#'
#' @keywords internal
#' @return The markers clusters data frame ordered by genes.
#'
#' @importFrom stats hclust
#' @importFrom stats dist
#' @noRd
.orderGenesInCluster <- function(cluster, markersClusters, mtx,
                                clusteringMethod="ward.D2"){

    genes <- markersClusters[markersClusters$clusters == cluster, ]$geneName

    if(length(genes) > 2){
        tree <- hclust(dist(mtx[as.vector(genes), ]), method=clusteringMethod)
        return(genes[tree$order])
    } else {
        return(genes)
    }
}



#' .checkParamSubFunction
#'
#' @description checks parameters of plotCellHeatmap.
#'
#' @param clusterCols If TRUE, the columns representing the clusters are also
#' taken into account in the hierarchical clustering. Default=FALSE.
#' @param showColnames Shoud the names of the columns (clusters) be indicated on
#' the heatmap. Default = FALSE.
#' @param plotPDF If TRUE, the heatmap is saved in pdf format and in png
#' otherwise. Default = TRUE.
#' @param fontsize base fontsize for the plot. Default = 7.5.
#' @param fontsizeRow fontsize for rownames. Default = 8.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param silentPlot If TRUE, does not plot the pheatmap. Default=FALSE.
#' @param returnPlot If TRUE returns a pheatmap object. Default=FALSE.
#' @param savePlot If TRUE and plotPDF=TRUE, save the heatmap in pdf format.
#' The heatmap is saved in the output directory defined in theObject
#' (?getOutputDirectory) and in the sub-directory 'pictures'.
#'
#' @keywords internal
#' @noRd
.checkParamSubFunction <- function(clusterCols, showColnames, plotPDF, 
        fontsize, fontsizeRow, widthPNG, heightPNG, silentPlot, returnPlot, 
        savePlot){
    
    ## Verify clusterCols
    if(!is.logical(clusterCols))
        stop("clusterCols should be a boolean.")
    
    ## Verify showColnames
    if(!is.logical(showColnames))
        stop("showColnames should be a boolean.")
    
    ## Verify plotPDF
    if(!is.logical(plotPDF))
        stop("plotPDF should be a boolean.")
    
    ## Verify fontsize
    if(!is.numeric(fontsize))
        stop("fontsize should be a numeric.")
    
    ## Verify fontsizeRow
    if(!is.numeric(fontsizeRow))
        stop("fontsizeRow should be a numeric.")
    
    ## Verify widthPNG
    if(!is.numeric(widthPNG))
        stop("widthPNG should be a numeric.")
    
    ## Verify heightPNG
    if(!is.numeric(heightPNG))
        stop("heightPNG should be a numeric.")
    
    ## Verify silentPlot
    if (!is.logical(silentPlot))
        stop("silentPlot should be a boolean.")
    
    if(silentPlot && !returnPlot && !savePlot)
        stop("You do not plot, neither save the heatmap or return the object.",
                " Nothing will happen. You should either plot the results, ",
                "return the object or save the heatmap.")
}

#' .checkParamCellHeatmap
#'
#' @description checks parameters of plotCellHeatmap
#'
#' @param fileName Name of the output file to which the heatmap is saved.
#' @param meanCentered Boolean indicating if mean centering should be applied
#' to the expression matrix. Default = TRUE.
#' @param orderClusters If True, clusters in the similarity matrix of cells will
#' be ordered by name. Default = FALSE.
#' @param orderGenes Boolean, should the heatmap be structured by gene.
#' Default = FALSE.
#' @param returnPlot If TRUE returns a pheatmap object. Default=FALSE.
#' @param savePlot If TRUE and plotPDF=TRUE, save the heatmap in pdf format.
#' The heatmap is saved in the output directory defined in theObject
#' (?getOutputDirectory) and in the sub-directory 'pictures'.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 10.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 8.5.
#' @param markersClusters Data frame containing two columns geneName and
#' clusters.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param clusterCols If TRUE, the columns representing the clusters are also
#' taken into account in the hierarchical clustering. Default=FALSE.
#' @param showColnames Shoud the names of the columns (clusters) be indicated on
#' the heatmap. Default = FALSE.
#' @param fontsize base fontsize for the plot. Default = 7.5.
#' @param fontsizeRow fontsize for rownames. Default = 8.
#' @param plotPDF If TRUE, the heatmap is saved in pdf format and in png
#' otherwise. Default = TRUE.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param silentPlot If TRUE, does not plot the pheatmap. Default=FALSE.
#'
#' @keywords internal
#' @noRd
.checkParamCellHeatmap <- function(fileName, meanCentered,
        orderClusters, orderGenes, returnPlot, savePlot, width,
        height, markersClusters, onefile, clusterCols, showColnames,
        fontsize, fontsizeRow, plotPDF, widthPNG, heightPNG, silentPlot){

    if(!isTRUE((nrow(markersClusters) > 1)))
        stop(paste("You have to calculate the cluster markers before plotting.",
                    "Please see retrieveTopClustersMarkers method."))

    ## Verify fileName
    if(!is.character(fileName) || grepl("/", fileName , fixed = TRUE))
        stop("fileName should be a string, no path.")

    ## Verify meanCentered
    if(!is.logical(meanCentered))
        stop("meanCentered should be a boolean.")

    ## Verify orderClusters
    if(!is.logical(orderClusters))
        stop("orderClusters should be a boolean.")

    ## Verify orderGenes
    if(!is.logical(orderGenes))
        stop("orderGenes should be a boolean.")

    ## Verify returnPlot
    if(!is.logical(returnPlot))
        stop("returnPlot should be a boolean.")

    ## Verify savePlot
    if(!is.logical(savePlot))
        stop("savePlot should be a boolean.")

    ## Verify width
    if(!is.numeric(width))
        stop("width should be a numeric.")

    ## Verify height
    if(!is.numeric(height))
        stop("height should be a numeric.")

    ## Verify onefile
    if(!is.logical(onefile))
        stop("onefile should be a boolean.")

    .checkParamSubFunction(clusterCols, showColnames, plotPDF, fontsize, 
            fontsizeRow, widthPNG, heightPNG, silentPlot, returnPlot, savePlot)
}



#' .orderClustersForHeatmap
#'
#' @description 
#' Centers the values of the matrix according to the mean if meanCentered is
#' TRUE and orders the cells using a hierarchical clustering.
#'
#' @param sceObject SingleCellExperiment object retrieved with ?getSceNorm. See
#'?SingleCellExperiment::SingleCellExperiment for more details.
#' @param markersClusters Data frame containing two columns geneName and
#' clusters.
#' @param meanCentered If TRUE, centers the values according to the mean.
#' @param orderClusters If True, clusters in the similarity matrix of cells will
#' be ordered by name. Default = FALSE.
#' @param colDF Data frame of Metadata of the scRNA-Seq object.
#' @param clusteringMethod Clustering method passed to hclust() function.
#' See ?hclust for a list of method. Default = "ward.D2".
#' @param orderGenes Boolean, should the heatmap be structured by gene.
#' Default = FALSE.
#' 
#' @importFrom Biobase exprs
#' @importFrom stats hclust
#' @keywords internal
#' @noRd
.orderClustersForHeatmap <- function(sceObject, markersClusters, 
        meanCentered, orderClusters, colDf, clusteringMethod, orderGenes){
    
    exprsTmp <- Biobase::exprs(sceObject)
    rowTmp <- rownames(exprsTmp)
    expressionMatrix <- exprsTmp[rowTmp %in% markersClusters$geneName, ]
    
    if(meanCentered){
        meanRows <- rowSums(expressionMatrix) / ncol(expressionMatrix)
        expressionMatrix <- expressionMatrix - meanRows
    }
    
    
    if(orderClusters){
        
        # Ordering expressionMatrixrix
        newOrder <- .callOrderCells(colDf, expressionMatrix,
                clusteringMethod)
        expressionMatrix <- expressionMatrix[, newOrder]
        clusterCols <- FALSE
        
        if(orderGenes){
            
            newOrder <- .callOrderGenes(colDf, markersClusters,
                    expressionMatrix, clusteringMethod)
            expressionMatrix <- expressionMatrix[newOrder, ]
            clusterRows <- FALSE
        }
        
    }else{
        
        distanceMatrix <- dist(t(expressionMatrix))
        clusterCols <- hclust(distanceMatrix, method="ward.D2")
    }
    
    if(!orderGenes)
        clusterRows <- hclust(dist(expressionMatrix), method="ward.D2")
    
    return(list(expressionMatrix, clusterCols, clusterRows))
}


#' .plotCellH
#'
#' @description 
#' This function plots heatmap with marker genes on rows and clustered cells
#' on columns.
#'
#' @param colDF Data frame of Metadata of the scRNA-Seq object.
#' @param colorPalette A vector of colors for clusters. Default = "default",
#' See details.
#' @param statePalette A vector of colors for states or conditions. See details.
#' @param expressionMatrix 
#' @param showColnames Shoud the names of the columns (clusters) be indicated on
#' the heatmap. Default = FALSE.
#' @param fontsizeRow fontsize for rownames. Default = 8.
#' @param clusterCols If TRUE, the columns representing the clusters are also
#' taken into account in the hierarchical clustering. Default=FALSE.
#' @param clusterRows If TRUE, the rows are ordered in the hierarchical
#' clustering.
#' @param fontsize base fontsize for the plot. Default = 7.5.
#' @param silentPlot If TRUE, does not plot the pheatmap. Default=FALSE.
#' @param plotPDF If TRUE, the heatmap is saved in pdf format and in png
#' otherwise. Default = TRUE.
#' @param fileName Name of the output file to which the heatmap is saved.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 10.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 8.5.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param theObject A scRNAseq object with the cluster similarity matrix
#' obtained with ?calculateClustersSimilarity method and the top markers
#' obtained with ?retrieveTopClustersMarkers.
#' 
#' @keywords internal
#' @noRd
#' @return The pheatmap object of the clustering.
.plotCellH <- function(colDf, colorPalette, statePalette, expressionMatrix,
        showColnames, fontsizeRow, clusterCols, clusterRows, fontsize, 
        silentPlot, savePlot, plotPDF, fileName, width, height, onefile,
        widthPNG, heightPNG, theObject){
    
    annotationColors <- .generateAnnotationColors(colDf, colorPalette,
            statePalette)
    
    columnsToPlot <- switch(is.null(colDf$state) + 1, c("clusters",
                    "state"), c("clusters"))
    
    if(is.null(colDf$clusters)){
        
        annCol <- switch(is.null(colDf$state) + 1,
                as.data.frame(colDf["state"]), NA)
        annColors <- switch(is.null(colDf$state) + 1,
                annotationColors[names(annotationColors) == "state"], NA)
    }else{
        annCol <- as.data.frame(colDf[columnsToPlot])
        annColors <- annotationColors
    }
    
    color <- colorRampPalette(c("#023b84", "#4b97fc", "#c9d9ef", "#FEE395",
                    "#F4794E", "#D73027", "#a31008", "#7a0f09"))(100)
    pheatmapObject <- pheatmap::pheatmap(expressionMatrix,
            show_colnames=showColnames, annotation_col=annCol,
            annotation_colors=annColors, fontsize_row=fontsizeRow,
            cluster_cols=clusterCols, cluster_rows=clusterRows,
            color=color, fontsize=fontsize, silent=silentPlot)
    
    
    if(savePlot)
        .saveHeatmap(theObject, plotPDF, fileName, width, height, onefile,
                widthPNG, heightPNG, pheatmapObject)
    
    return(pheatmapObject)
    
}


#' plotCellHeatmap
#'
#' @description This function plots heatmap with marker genes on rows and
#' clustered cells on columns.
#'
#' @usage
#' plotCellHeatmap(theObject, fileName = NA, meanCentered=TRUE,
#'                 colorPalette="default", statePalette="default",
#'                 clusteringMethod="ward.D2", orderClusters=FALSE,
#'                 orderGenes=FALSE, returnPlot=FALSE, savePlot=FALSE, width=10,
#'                 height=8.5, onefile=FALSE, clusterCols=FALSE,
#'                 showColnames=FALSE, fontsize=7.5,  fontsizeRow=8,
#'                 plotPDF=TRUE, widthPNG=800, heightPNG=750, silentPlot)
#'
#' @param theObject A scRNAseq object with the cluster similarity matrix
#' obtained with ?calculateClustersSimilarity method and the top markers
#' obtained with ?retrieveTopClustersMarkers.
#' @param fileName Name of the output file to which the heatmap is saved.
#' @param meanCentered Boolean indicating if mean centering should be applied
#' to the expression matrix. Default = TRUE.
#' @param colorPalette A vector of colors for clusters. Default = "default",
#' See details.
#' @param statePalette A vector of colors for states or conditions. See details.
#' @param clusteringMethod Clustering method passed to hclust() function.
#' See ?hclust for a list of method. Default = "ward.D2".
#' @param orderClusters If True, clusters in the similarity matrix of cells will
#' be ordered by name. Default = FALSE.
#' @param orderGenes Boolean, should the heatmap be structured by gene.
#' Default = FALSE.
#' @param returnPlot If TRUE returns a pheatmap object. Default=FALSE.
#' @param savePlot If TRUE and plotPDF=TRUE, save the heatmap in pdf format.
#' The heatmap is saved in the output directory defined in theObject
#' (?getOutputDirectory) and in the sub-directory 'pictures'.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 10.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 8.5.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param clusterCols If TRUE, the columns representing the clusters are also
#' taken into account in the hierarchical clustering. Default=FALSE.
#' @param showColnames Shoud the names of the columns (clusters) be indicated on
#' the heatmap. Default = FALSE.
#' @param fontsize base fontsize for the plot. Default = 7.5.
#' @param fontsizeRow fontsize for rownames. Default = 8.
#' @param plotPDF If TRUE, the heatmap is saved in pdf format and in png
#' otherwise. Default = TRUE.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param silentPlot If TRUE, does not plot the pheatmap. Default=FALSE.
#'
#' @details
#' colorPalette/statePalette -- A vector of colors for clusters/states or
#' 'default' value. If 'default' is selected, the number of clusters is limited
#' to 16. If an error message is thrown, re-run the function with your own color
#' vector.
#'
#' @return A pheatmap object of the heatmap if returnPlot is TRUE.
#'
#' @aliases plotCellHeatmap
#' @rdname plotCellHeatmap
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
#' columnID="cell_ID")
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
#' ## Retrieve the top 10 markers per cluster
#' scrFinal <- retrieveTopClustersMarkers(scrS4MG)
#'
#' ## Plot the heatmap with marker genes
#' plotCellHeatmap(scrFinal)
#'
#' @seealso calculateClustersSimilarity  plotClusteredTSNE plotCellSimilarity
#' plotGeneExpression plotClustersSimilarity
#'
#' @exportMethod plotCellHeatmap
#' @importFrom SummarizedExperiment colData
#' @importFrom Biobase exprs
#' @importFrom pheatmap pheatmap
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom grDevices colorRampPalette
#' @importFrom methods validObject
#'
#' @author
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
setMethod(

    f = "plotCellHeatmap",

    signature = "scRNAseq",

    definition = function(theObject, fileName = NA, meanCentered=TRUE,
            colorPalette="default", statePalette="default",
            clusteringMethod="ward.D2", orderClusters=FALSE, orderGenes=FALSE,
            returnPlot=FALSE, savePlot=FALSE, width=10, height=8.5,
            onefile=FALSE, clusterCols=FALSE, showColnames=FALSE, fontsize=7.5,
            fontsizeRow=8, plotPDF=TRUE, widthPNG=800, heightPNG=750,
            silentPlot=FALSE){

        ## Verify parameters
        validObject(theObject)

        sceObject <- getSceNorm(theObject)
        colDf <- SummarizedExperiment::colData(sceObject)
        markersClusters <- getClustersMarkers(theObject)

        if(is.na(fileName))
            fileName <- paste0(getExperimentName(theObject), "_", "clusters",
                    length(levels(colDf$clusters)), "_meanCentered",
                    meanCentered, "_orderClusters", orderClusters,
                    "_orderGenes", orderGenes, "markersPerCluster")

        .checkParamCellHeatmap(fileName, meanCentered, orderClusters,
                orderGenes, returnPlot, savePlot, width,
                height, markersClusters, onefile, clusterCols, showColnames,
                fontsize, fontsizeRow, plotPDF, widthPNG, heightPNG, silentPlot)


        # plots correlation between clusters
        results <- .orderClustersForHeatmap(sceObject, markersClusters, 
                meanCentered, orderClusters, colDf, clusteringMethod, 
                orderGenes)
        expressionMatrix <- results[[1]]
        clusterCols <- results[[2]]
        clusterRows <- results[[3]]
        
        pheatmapObject <- .plotCellH(colDf, colorPalette, statePalette, 
                expressionMatrix, showColnames, fontsizeRow, clusterCols, 
                clusterRows, fontsize, silentPlot, savePlot, plotPDF, fileName,
                width, height, onefile, widthPNG, heightPNG, theObject)

        if(returnPlot)
            return(pheatmapObject)
})


################################################################################
## plotGeneExpression
################################################################################

#' .checkParamPlotGeneExpression
#'
#' @description checks parameters of plotGeneExpression
#'
#' @param theObject A scRNAseq object with the top markers retrieved. See
#' ?retrieveTopClustersMarkers.
#' @param geneName Name of the gene to highlight on the t-SNE plot.
#' @param returnPlot If TRUE, returns a ggplot object of the tSNE.
#' Default = FALSE.
#' @param savePlot If TRUE, save the tSNE in pdf or png format. Default=FALSE.
#' @param width Width of the plot. Default = 6.
#' @param height Height of the plot. Default = 5.
#' @param expMat Count matrix retrieved with ?Biobase::exprs.
#' @param tSNECoords Coordinates of the tSNE plot retrieved with
#' ?getCoordinates.
#' @param silentPlot If TRUE, the plots are not displayed on the current device.
#' Default=FALSE.
#' @param plotPDF If TRUE export tSNE in pdf format, if FALSE export it in
#' png format. Default=TRUE.
#'
#' @keywords internal
#' @noRd
.checkParamPlotGeneExpression <- function(theObject, geneName, returnPlot,
        savePlot, width, height, expMat, tSNECoords, silentPlot, plotPDF,
        sceObject){


    if(!geneName %in% rownames(expMat))
        stop("Gene is not found in expression matrix.")

    if(!isTRUE(all.equal(rownames(tSNECoords), colnames(sceObject))))
        stop("The row names of the tSNE coordinates matrix should be equal",
                " to the colnames of the expression matrix.")

    ## Verify the geneName
    markers <- getClustersMarkers(theObject)$geneName
    if(!geneName %in% markers)
        stop(paste("geneName should be a marker founded by ",
                "retrieveTopClustersMarkers method'. Please see the",
                "documentation about retrieveTopClustersMarkers method."))

    ## Verify returnPlot
    if (!is.logical(returnPlot))
        stop("returnPlot should be a boolean.")

    ## Verify savePlot
    if (!is.logical(savePlot))
        stop("savePlot should be a boolean.")

    ## Verify width
    if (!is.numeric(width))
        stop("width should be a numeric.")

    ## Verify height
    if (!is.numeric(height))
        stop("height should be a numeric.")

    ## Verify silentPlot
    if (!is.logical(silentPlot))
        stop("silentPlot should be a boolean.")

    ## Verify plotPDF
    if (!is.logical(plotPDF))
        stop("plotPDF should be a boolean.")

    if(silentPlot && !returnPlot && !savePlot)
        stop("You do not plot, neither save the heatmap or return the object.",
                " Nothing will happen. You should either plot the results, ",
                "return the object or save the heatmap.")
}


#' .saveAndPlotGeneExpression
#'
#' @description
#' Saves and/or plots the tSNE colored by expression.
#'
#' @param theObject A scRNAseq object with the top markers retrieved. See
#' ?retrieveTopClustersMarkers.
#' @param tSNECoords Coordinates of the tSNE.
#' @param pointSize Size of the points on the tSNE. Default = 1.
#' @param alpha Opacity of the points of the plot. Default = 1.
#' @param palette Color palette for the expression levels.
#' @param limits Range of the gene expression shown in the legend. Default = NA.
#' See details.
#' @param geneName Name of the gene to highlight on the t-SNE plot.
#' @param savePlot If TRUE, save the tSNE in pdf or png format. Default=FALSE.
#' @param List of unique cluster numbers.
#' @param tSNEpicture Number of the tSNE picture that you want to use for
#' plotting the gene expression. Default = 1.
#' @param plotPDF If TRUE export tSNE in pdf format, if FALSE export it in
#' png format. Default=TRUE.
#' @param width Width of the plot. Default = 6.
#' @param height Height of the plot. Default = 5.
#' @param silentPlot If TRUE, the plots are not displayed on the current device.
#' Default=FALSE.
#'
#' @keywords internal
#' @noRd
.saveAndPlotGeneExpression <- function(theObject, tSNECoords, pointSize, alpha,
        palette, limits, geneName, savePlot, clustersNumber, tSNEpicture, 
        plotPDF, width, height, silentPlot){
    
    ggres <- ggplot2::ggplot(tSNECoords, aes(x=tSNECoords[,1],
                            y=tSNECoords[,2], color=expression)) +
            geom_point(size=I(pointSize), alpha=alpha) + theme_bw() +
            scale_colour_gradientn(
                    colours=alpha(colorRampPalette(palette)(100), 0.8),
                    limits=limits) + ggtitle(geneName)
    
    if(savePlot){
        
        dataDirectory  <- getOutputDirectory(theObject)
        experimentName <- getExperimentName(theObject)
        subdir <- file.path(dataDirectory, "pictures")
        fileName <- paste(experimentName, "tSNE", clustersNumber,
                "clusters", geneName, "tSNEpicture",
                tSNEpicture, "_alpha", alpha,sep="_")
        
        if(!file.exists(subdir))
            dir.create(subdir, showWarnings=FALSE, recursive = TRUE)
        
        ggsave(filename= paste0(fileName, if(plotPDF) ".pdf" else ".png"),
                plot=ggres, device= if(plotPDF) "pdf" else "png",
                path=subdir, width= width, height = height)
    }
    
    if(!silentPlot)
        print(ggres)
    
    return(ggres)
}

#' plotGeneExpression
#'
#' @description The function saves a t-SNE plot colored by expression of a
#' given gene.
#'
#' @usage
#' plotGeneExpression(theObject, geneName,
#' palette=c("grey","red", "#7a0f09", "black"), returnPlot=FALSE,
#'             tSNEpicture=1, savePlot=FALSE, alpha=1, limits=NA,
#'             pointSize=1, width=6, height=5, plotPDF=TRUE, silentPlot=FALSE)
#'
#' @param theObject A scRNAseq object with the top markers retrieved. See
#' ?retrieveTopClustersMarkers.
#' @param geneName Name of the gene to highlight on the t-SNE plot.
#' @param palette Color palette for the expression levels.
#' @param returnPlot If TRUE, returns a ggplot object of the tSNE.
#' Default = FALSE.
#' @param tSNEpicture Number of the tSNE picture that you want to use for
#' plotting the gene expression. Default = 1.
#' @param savePlot If TRUE, save the tSNE in pdf or png format. Default=FALSE.
#' @param alpha Opacity of the points of the plot. Default = 1.
#' @param limits Range of the gene expression shown in the legend. Default = NA.
#' See details.
#' @param pointSize Size of the points on the tSNE. Default = 1.
#' @param width Width of the plot. Default = 6.
#' @param height Height of the plot. Default = 5.
#' @param plotPDF If TRUE export tSNE in pdf format, if FALSE export it in
#' png format. Default=TRUE.
#' @param silentPlot If TRUE, the plots are not displayed on the current device.
#' Default=FALSE.
#'
#' @aliases plotGeneExpression
#' @rdname plotGeneExpression-scRNAseq
#'
#' @details
#' limits -- This option allows generating t-SNE plots with equal color scale to
#' compare the expression of different genes. By default, limits are the range
#' of expression of a selected gene.
#'
#' @return A ggplot object of the gene expression colored tSNE.
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
#' columnID="cell_ID")
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
#' ## Retrieve top clusters markers
#' scrFinal <- retrieveTopClustersMarkers(scrS4MG, removeDuplicates=FALSE)
#'
#' ## t-SNE plot colored by expression of a  given gene.
#' plotGeneExpression(scrFinal, getClustersMarkers(scrFinal)[1,1])
#'
#' @seealso retrieveTopClustersMarkers plotCellSimilarity plotCellHeatmap
#' plotClusteredTSNE plotClustersSimilarity
#'
#' @exportMethod plotGeneExpression
#' @importFrom SummarizedExperiment colData
#' @importFrom Biobase exprs
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_colour_gradientn
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 ggsave
#' @importFrom grDevices colorRampPalette
#' @importFrom methods validObject
#' @importFrom scales alpha
#'
#' @author
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
setMethod(

    f = "plotGeneExpression",

    signature = "scRNAseq",

    definition = function(theObject, geneName,
            palette=c("grey","red", "#7a0f09", "black"), returnPlot=FALSE,
            tSNEpicture=1, savePlot=FALSE, alpha=1, limits=NA,
            pointSize=1, width=6, height=5, plotPDF=TRUE, silentPlot=FALSE){

        ## Verify parameters
        validObject(theObject)
        sceObject <- getSceNorm(theObject)
        colDf <- SummarizedExperiment::colData(sceObject)
        clustersNumber <- length(unique(colDf$clusters))
        tmpList <- getTSNEList(theObject)
        tSNECoords <- as.data.frame(getCoordinates(tmpList[[tSNEpicture]]))
        tSNECoords <- tSNECoords[colDf$cellName, ]
        expMat <- Biobase::exprs(sceObject)

        .checkParamPlotGeneExpression(theObject, geneName, returnPlot,
                savePlot, width, height, expMat, tSNECoords, silentPlot,
                plotPDF, sceObject)

        ## Plot all precalculated t-SNEs to show your clusters
        tSNECoords$expression <- expMat[geneName, ]

        if(isTRUE(all.equal(length(limits), 1)))
            limits <- c(min(tSNECoords$expression), max(tSNECoords$expression))

        ggres <- .saveAndPlotGeneExpression(theObject, tSNECoords, pointSize, 
                alpha, palette, limits, geneName, savePlot, clustersNumber, 
                tSNEpicture, plotPDF, width, height, silentPlot)
        
        if(returnPlot)
            return(ggres)
    })


#################
## plotClustersSimilarity
#################

#' .checkParamClustersSimilarity
#'
#' @description checks parameters of plotCellSimilarity
#'
#' @param theObject A scRNAseq object with the cluster similarity matrix
#' obtained with ?calculateClustersSimilarity.
#' @param returnPlot Boolean indicating if the pHeatmap object should  be
#' returned by the function. Default = FALSE.
#' @param savePlot If TRUE and plotPDF=TRUE, save the heatmap in pdf format.
#' The heatmap is saved in the output directory defined in theObject
#' (?getOutputDirectory) and in the sub-directory 'pictures'.
#' @param plotPDF If TRUE, the heatmap is saved in pdf format and in png
#' otherwise. Default = TRUE.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 7.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 5.5.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param fontsize pheatmap parameter. Base fontsize for the plot. Default=7.5.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param silentPlot If TRUE, does not plot the pheatmap. Default=FALSE.
#'
#' @keywords internal
#' @noRd
.checkParamClustersSimilarity <- function(theObject, returnPlot, savePlot,
        plotPDF, width, height, onefile, fontsize, widthPNG, heightPNG,
        silentPlot){

    ## Verify the object
    clustersSimilarityMatrix <- getClustersSimilarityMatrix(theObject)
    clusters <- colData(getSceNorm(theObject))$clusters
    nbrClusters <- length(unique(clusters))

    if(!isTRUE((ncol(clustersSimilarityMatrix) == nbrClusters)))
        stop(paste("You have to calculate the cluster similarity matrix",
            "before plotting."))

    ## Verify returnPlot
    if(!is.logical(returnPlot)) stop("returnPlot should be a boolean.")

    ## Verify savePlot
    if(!is.logical(savePlot)) stop("savePlot should be a boolean.")

    ## Verify plotPDF
    if(!is.logical(plotPDF)) stop("plotPDF should be a boolean.")

    ## Verify width
    if(!is.numeric(width)) stop("width should be a numeric.")

    ## Verify height
    if(!is.numeric(height)) stop("height should be a numeric.")

    ## Verify onefile
    if(!is.logical(onefile)) stop("onefile should be a boolean.")

    ## Verify fontsize
    if(!is.numeric(fontsize)) stop("fontsize should be a numeric.")

    ## Verify widthPNG
    if(!is.numeric(widthPNG)) stop("widthPNG should be a numeric.")

    ## Verify heightPNG
    if(!is.numeric(heightPNG)) stop("heightPNG should be a numeric.")

    ## Verify silentPlot
    if(!is.logical(silentPlot)) stop("silentPlot should be a boolean.")

    if(silentPlot && !returnPlot && !savePlot)
        stop("You do not plot, neither save the heatmap or return the object.",
                " Nothing will happen. You should either plot the results, ",
                "return the object or save the heatmap.")
}


#' .pheatmapClusterSim
#'
#' @description Generates the heatmap of the clusters similarities.
#'
#' @param theObject A scRNAseq object with the cluster similarity matrix
#' obtained with ?calculateClustersSimilarity.
#' @param clusteringMethod Clustering method passed to hclust() function.
#' See ?hclust for a list of method. Default = "ward.D2"
#' @param colorPalette A vector of colors for clusters. Default = "default",
#' see details.
#' @param statePalette A vector of colors for states or conditions. See details.
#' @param fontsize pheatmap parameter. Base fontsize for the plot. Default=7.5.
#' @param silentPlot If TRUE, does not plot the pheatmap. Default=FALSE.
#'
#' @keywords internal
#' @importFrom SummarizedExperiment colData
#' @importFrom pheatmap pheatmap
#' @importFrom stats as.dist
#' @importFrom stats hclust
#' @noRd
#' @return The clusters numbers and the pheatmap object.
.pheatmapClusterSim  <- function(theObject, clusteringMethod, colorPalette, 
        statePalette, fontsize, silentPlot){
    
    clustersSimilarityMatrix <- getClustersSimilarityMatrix(theObject)
    colDf <- SummarizedExperiment::colData(getSceNorm(theObject))
    clusters <- colDf$clusters
    clustersNames <- levels(clusters)
    distanceMatrix <- as.dist(sqrt((1-clustersSimilarityMatrix)/2))
    clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
    colDataSimilarity <- data.frame(clusters=clustersNames)
    rownames(colDataSimilarity) <- colDataSimilarity$clusters
    
    annotationColors <- .generateAnnotationColors(colDf, colorPalette,
            statePalette)
    
    pheatmapObject <- pheatmap::pheatmap(clustersSimilarityMatrix,
            annotation_col=colDataSimilarity,
            annotation_colors=annotationColors,
            cluster_cols=clusteringTree,
            cluster_rows=clusteringTree,
            fontsize=fontsize,
            main="Clusters similarity matrix",
            silent=silentPlot)
    
    return(list(clusters, pheatmapObject))
}


#' .savePlotClustersSim
#'
#' @description Write the heatmap to a specified output.
#'
#' @param theObject A scRNAseq object with the cluster similarity matrix
#' obtained with ?calculateClustersSimilarity.
#' @param savePlot If TRUE and plotPDF=TRUE, save the heatmap in pdf format.
#' The heatmap is saved in the output directory defined in theObject
#' @param clusters Cluster numbers contained in the columns meta-data.
#' @param plotPDF If TRUE, the heatmap is saved in pdf format and in png
#' otherwise. Default = TRUE.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 7.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 5.5.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param pheatmapObject Object returned by the pheatmap function in 
#' .pheatmapClusterSim.
#' 
#' @importFrom grDevices pdf
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @keywords internal
#' @noRd
.savePlotClustersSim <- function(savePlot, theObject, clusters, plotPDF, width,
        height, onefile, widthPNG, heightPNG, pheatmapObject){
    
    if(savePlot){
        
        dataDirectory   <- getOutputDirectory(theObject)
        experimentName  <- getExperimentName(theObject)
        clustersNumber <- length(unique(clusters))
        subdir <- file.path(dataDirectory, "pictures")
        
        if(!file.exists(subdir))
            dir.create(subdir, showWarnings=FALSE, recursive = TRUE)
        
        fileName <- paste(experimentName,"clusters_similarity",
                clustersNumber, "clusters", sep="_")
        filePath <- file.path(subdir, fileName)
        
        if(plotPDF)
            pdf(file=paste0(filePath, ".pdf"), width=width, height=height,
                    onefile=onefile)
        else
            png(filename=paste0(filePath, ".png"), width=widthPNG,
                    height=heightPNG, type="cairo")
        grid::grid.newpage()
        grid::grid.draw(pheatmapObject$gtable)
        dev.off()
    }
}

#' plotClustersSimilarity
#'
#' @description This function plots the clusters similarity matrix as a heatmap.
#'
#' @usage
#' plotClustersSimilarity(theObject, colorPalette="default",
#' statePalette="default", clusteringMethod="ward.D2", returnPlot=FALSE,
#' savePlot=FALSE, plotPDF=TRUE, width=7, height=5.5, onefile=FALSE,
#' fontsize=7.5, widthPNG=800, heightPNG=750, silentPlot=FALSE)
#'
#' @param theObject A scRNAseq object with the cluster similarity matrix
#' obtained with ?calculateClustersSimilarity.
#' @param colorPalette A vector of colors for clusters. Default = "default",
#' see details.
#' @param statePalette A vector of colors for states or conditions. See details.
#' @param clusteringMethod Clustering method passed to hclust() function.
#' See ?hclust for a list of method. Default = "ward.D2"
#' @param returnPlot Boolean indicating if the pHeatmap object should  be
#' returned by the function. Default = FALSE.
#' @param savePlot If TRUE and plotPDF=TRUE, save the heatmap in pdf format.
#' The heatmap is saved in the output directory defined in theObject
#' (?getOutputDirectory) and in the sub-directory 'pictures'.
#' @param plotPDF If TRUE, the heatmap is saved in pdf format and in png
#' otherwise. Default = TRUE.
#' @param width Width of the plot in the pdf file. See ?pdf for more details.
#' Default = 7.
#' @param height Height of the plot in the pdf file. See ?pdf for more details.
#' Default = 5.5.
#' @param onefile Logical: if TRUE allow multiple figures in one file. If FALSE,
#' generate a file with name containing the page number for each page.
#' Defaults to FALSE.
#' @param fontsize pheatmap parameter. Base fontsize for the plot. Default=7.5.
#' @param widthPNG Width of the png. See ?png for details. Default=800.
#' @param heightPNG Height of the png. See ?png for details. Default=750.
#' @param silentPlot If TRUE, does not plot the pheatmap. Default=FALSE.
#'
#' @details
#' colorPalette/statePalette -- A vector of colors for clusters/states or
#' 'default' value. If 'default' is selected, the number of clusters is limited
#' to 16. If an error message is thrown, re-run the function with your own color
#' vector.
#'
#' @aliases plotClustersSimilarity
#' @rdname plotClustersSimilarity
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
#' columnID="cell_ID")
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
#' ## Plot similarity matrix as a heatmap
#' plotClustersSimilarity(scrCSM)
#'
#' @return A pheatmap object of the clusters similarity matrix.
#'
#' @seealso alculateClustersSimilarity  plotClusteredTSNE plotCellHeatmap
#' plotGeneExpression plotCellSimilarity
#'
#' @exportMethod plotClustersSimilarity
#' @importFrom SummarizedExperiment colData
#' @importFrom pheatmap pheatmap
#' @importFrom stats as.dist
#' @importFrom stats hclust
#' @importFrom grDevices pdf
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom methods validObject
#'
#' @author
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
setMethod(

    f = "plotClustersSimilarity",

    signature = "scRNAseq",

    definition = function(theObject, colorPalette="default",
                            statePalette="default", clusteringMethod="ward.D2",
                            returnPlot=FALSE, savePlot=FALSE, plotPDF=TRUE,
                            width=7, height=5.5, onefile=FALSE, fontsize=7.5,
                            widthPNG=800, heightPNG=750, silentPlot=FALSE){

        ## Verify parameters
        validObject(theObject)
        .checkParamClustersSimilarity(theObject, returnPlot, savePlot, plotPDF,
                width, height, onefile, fontsize, widthPNG, heightPNG,
                silentPlot)

        ## Generate heatmap
        result <- .pheatmapClusterSim(theObject, clusteringMethod, colorPalette,
                statePalette, fontsize, silentPlot)
        clusters <- result[[1]]
        pheatmapObject <- result[[2]]                        
        
        ## Save plot
        .savePlotClustersSim(savePlot, theObject, clusters, plotPDF, width, 
                height, onefile, widthPNG, heightPNG, pheatmapObject)
        
        if(returnPlot)
            return(pheatmapObject)
    })
