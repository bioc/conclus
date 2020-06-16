#################
## plotCellSimilarity
#################


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


.orderCellsInCluster <- function(cluster, colData, mtx, 
		clusteringMethod="ward.D2"){
	
	# Order cells according to clustering results
	# Uses for ordering matrix to further plot it with pheatmap()
	
	cells <- colData[colData$clusters == cluster, ]$cellName
	if(length(cells) > 2){
		tree <- hclust(dist(t(mtx[, cells])), method=clusteringMethod)
		return(cells[tree$order])
	}else
		return(cells)
}


setMethod(
		
    f = "plotCellSimilarity",
	
    signature = "scRNAseq",
	
    definition = function(theObject, colorPalette="default", 
			statePalette="default", clusteringMethod="ward.D2", 
			orderClusters=FALSE, plotPDF=TRUE, returnPlot=FALSE, width=7, 
			height=6, onefile=FALSE, showRowNames=FALSE, showColnames=FALSE,
			fontsize=7.5, fontsizeRow=0.03, widthPNG=800, heightPNG=750){
        
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


#################
## plotClusteredTSNE
#################


.checkParamPlotTSNE <- function(theObject, columnName){
	
	## Verify experiment name 
	experimentName <- getExperimentName(theObject)
	
	if(!isTRUE(all.equal(columnName, "clusters")) && 
			!isTRUE(all.equal(columnName, "noColor")) &&
			!isTRUE(all.equal(columnName,"state")))
		stop("columnName should be: clusters, noColor, or state.")
}			

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
		dir.create(outputDir, showWarnings=F)
	
	return(list(outputDir, colorPalette))
}			


setMethod(
		
    f = "plotClusteredTSNE",
	
    signature = "scRNAseq",
	
    definition = function(theObject, colorPalette="default",
                        PCs=c(4, 6, 8, 10, 20, 40, 50), perplexities=c(30, 40),
                        columnName="clusters", returnPlot=FALSE, width=6,
                        height=5, onefile=FALSE, ...){
					
		## Verify parameters
		validObject(theObject)
		.checkParamPlotTSNE(theObject, columnName)
			
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
					"PCs*perplexities. Contact the developper.")
		
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


.orderGenesInCluster <- function(cluster, markersClusters, mtx, 
		clusteringMethod="ward.D2"){
	
    # Order cells according to clustering results
    # Uses for ordering matrix to further plot it with pheatmap()
    
    genes <- markersClusters[markersClusters$clusters == cluster, ]$geneName
    if(length(genes) > 2){
        tree <- hclust(dist(mtx[genes, ]), method=clusteringMethod)
        return(genes[tree$order])
    } else {
        return(genes)
    }
}

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



setMethod(
		
    f = "plotCellHeatmap",
	
    signature = "scRNAseq",
	
    definition = function(theObject, fileName, meanCentered=TRUE, 
			colorPalette="default", statePalette="default", 
			clusteringMethod="ward.D2", orderClusters=FALSE, similarity=FALSE,
			orderGenes=FALSE, returnPlot=FALSE, saveHeatmapTable=FALSE, 
			width=10, height=8.5){
        
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
    

#################
## plotGeneExpression
#################


setMethod(
		
    f = "plotGeneExpression",
	
    signature = "scRNAseq",
	
    definition = function(theObject, geneName, graphsDirectory="pictures", 
			palette=c("grey","red", "#7a0f09", "black"), returnPlot=FALSE,
			tSNEpicture=1, commentName="", savePlot=TRUE, alpha=1, limits=NA,
			pointSize=1, width=6, height=5, ...){
        
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

setMethod(
		
    f = "plotClustersSimilarity",
	
    signature = "scRNAseq",
	
    definition = function(theObject, colorPalette="default", 
			statePalette="default", clusteringMethod="ward.D2", 
			returnPlot=FALSE, width=7, height=5.5, onefile=FALSE, fontsize=7.5,
			...){
        
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
					
		message("\nSaving a heatmap with the clusters similarity matrix.")
		tmpFileP <- paste(experimentName,"clusters_similarity", clustersNumber,
				"clusters.pdf", sep="_")
		pdf(file.path(dataDirectory, graphsDirectory, tmpFileP), width=width, 
				height=height, onefile=onefile, ...)
		grid::grid.newpage()
		grid::grid.draw(pheatmapObject$gtable)
		dev.off()
					
		if(returnPlot)
			return(pheatmapObject)			
})
    