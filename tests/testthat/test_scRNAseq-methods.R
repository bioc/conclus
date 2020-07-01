### Test the methods of scRNAseq class.


if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::valid()  
BiocManager::install("biomaRt")
# BiocManager::install("grimbough/biomaRt")

if (!require("pacman")) {
    install.packages("pacman",  repos = "https://cran.biotools.fr/")
}

pacman::p_load(prettyunits, rlist, foreach, ggplot2, pheatmap, zoo,
               dynamicTreeCut, factoextra,
               digest, RColorBrewer,devtools, BiocParallel, scran, scater,
               monocle, SingleCellExperiment , KEGGREST, AnnotationDbi,
               Cairo, rvest, curl,  Matrix, dbscan, fpc, matrixStats,
               dplyr, org.Mm.eg.db, grDevices, S4Vectors, Biobase,
               DataCombine, zoo, rvest, DataCombine, doParallel,
               testthat, tidyr, biomaRt)

# source("R/AllGenerics.R")
# source("R/AllClasses.R")
# source("R/sharedInternals.R")
# source("R/checkFunctions.R")
# source("R/getters.R")
# source("R/setters.R")
# source("R/methods-normalization.R")
# source("R/methods-tsne.R")
# source("R/methods-dbscan.R")
# source("R/methods-clustering.R")
# source("R/methods-plot.R")
# source("R/methods-export.R")
# source("R/methods-markers.R")


library(conclus)
## Data

outputDirectory <- "YourOutputDirectory"
experimentName <- "Bergiers"
columnsMetaData <- read.delim(
    file.path("inst/extdata/Bergiers_colData_filtered.tsv"))

## Creation of the count Matrix

countMatrix <- as.matrix(read.delim(
    file.path("tests/testthat/test_data/test_countMatrix.tsv")))

smallMatrix <- countMatrix[,1:50]

## Load expected results

load(file = "tests/testthat/test_data/scrLight.Rdat")
load(file = "tests/testthat/test_data/expected_normalizedMatrix.Rdat")


## Construction of the object

scr <- scRNAseq(experimentName = experimentName, 
                  countMatrix     = countMatrix, 
                  species         = "mouse",
                  outputDirectory = outputDirectory)

## Performing the normalization

scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
sceNorm <- getSceNorm(scrNorm)


## Performing tSNE

scrTsne <- generateTSNECoordinates(scrNorm, cores=5)
tsneList <- getTSNEList(scrTsne)
tsneListWrong <- tsneList
setCoordinates(tsneListWrong[[1]]) <- getCoordinates(tsneList[[1]])[1:10,]

## Running DbScan

scrDbscan <- runDBSCAN(scrTsne, cores=5)
dbscanList <- getDbscanList(scrDbscan)
clusteringList <- lapply(dbscanList, getClustering)
dbscanListWrong <- dbscanList
setClustering(dbscanListWrong[[1]]) <- 
		as.matrix(getClustering(dbscanList[[1]])[,1:10])

## Running cluster cells internal 

scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber = 10,
		deepSplit = 4, cores = 4, clusteringMethod = "ward.D2")
cci <- getCellsSimilarityMatrix(scrCCI)
scrCCiwrong <- scrCCI
setCellsSimilarityMatrix(scrCCiwrong) <-  matrix(nrow=0,ncol=0)

## Calculate clusters similarity

scrCSM <- calculateClustersSimilarity(scrCCI)
csm <- getClustersSimilarityMatrix(scrCSM)
orderedCLusters <- getClustersSimiliratyOrdered(scrCSM)

## Ranking genes

scrS4MG <- rankGenes(scrCSM)
markers <- getMarkerGenesList(scrS4MG)


## Getting marker genes

scrFinal <- retrieveTopClustersMarkers(scrS4MG, removeDuplicates = F)
scrFinalWrong <- scrFinal
setTSNEList(scrFinalWrong) <- list(new("Tsne"))

## Getting genes info

scrInfos <- retrieveGenesInfo(scrFinal, species = "mouse", cores=5)
wrongInfo <- data.frame(uniprot_gn_symbol=c("symbol1", "symbol2"), 
		clusters=c("1", "3"), external_gene_name=c("gene1", "gene2"), 
		go_id=c("GO1,GO2", "GO1,GO3"), 
		mgi_description=c("description1", "description2"), 
		entrezgene_description=c("description1", "description2"),
		gene_biotype=c("coding", "coding"), chromosome_name=c("1", "2"), 
		Symbol=c("symbol1", "symbol2"), ensembl_gene_id=c("ENS1","ENS2"), 
		mgi_id=c("MGI1", "MGI2"), entrezgene_id=c("1", "2"),
		uniprot_gn_id=c("ID1", "ID2"))


####################  Construction of the object  ####################
		  
test_that("scr is created properly", {
			 
			 expect_identical(getExperimentName(scrLight), 
					 getExperimentName(scr))
			 
			 expect_identical(getCountMatrix(scrLight), getCountMatrix(scr))
			 
			 expect_identical(getSpecies(scrLight), getSpecies(scr))
			 
			 expect_identical(getOutputDirectory(scrLight), 
					 getOutputDirectory(scr))
		 })
		  
		  
test_that("Errors are thrown when creating scr", {
			
			expM <- "'experimentName' slot is empty. Please fill it."
			expect_error(scRNAseq(experimentName = "", 
                  countMatrix     = countMatrix, 
                  species         = "mouse",
                  outputDirectory = outputDirectory), regexp = expM)
			
			expM <- paste0("Experiment name should contain a single string ",
					"describing the experiment, 'My experiment' is not ",
					"correct.")
			expect_error(scRNAseq(experimentName  = "My experiment", 
							countMatrix     = countMatrix, 
							species         = "mouse",
							outputDirectory = outputDirectory), regexp = expM)
			
			expM <- paste0("'countMatrix' slot is empty. It should be a matrix",
					" containing at leat 100 cells.")
			expect_error(scRNAseq(experimentName  = experimentName, 
							countMatrix     = matrix(), 
							species         = "mouse",
							outputDirectory = outputDirectory), regexp = expM)
			
			expM <- paste0("Not enough cells in the count matrix. There Should",
					" be at leat 100 cells. The current count matrix contains",
					" 50 cells.")
			expect_error(scRNAseq(experimentName  = experimentName, 
							countMatrix     = smallMatrix, 
							species         = "mouse",
							outputDirectory = outputDirectory), regexp = expM)
						
			expM <- paste0("species should be 'mouse' or 'human'. '' is ",
					"currently not supported.")
			expect_error(scRNAseq(experimentName  = experimentName, 
							countMatrix     = countMatrix, 
							species         = "",
							outputDirectory = outputDirectory), regexp = expM)
			
			expM <- "'outputDirectory' slot is empty. Please fill it."
			expect_error(scRNAseq(experimentName = experimentName, 
					countMatrix     = countMatrix, 
					species         = "mouse",
					outputDirectory = ""), regexp = expM)
			
			expM <- paste0("'outputDirectory' should be a conform folder path:",
					"'toto tata' is not.")
			expect_error(scRNAseq(experimentName  = experimentName, 
							countMatrix     = countMatrix, 
							species         = "human",
							outputDirectory = "toto tata"), regexp = expM)
			
			expM <- "tSNEList is empty. This should be a list of tSNE objects."
			expect_error(scRNAseq(experimentName = experimentName, 
					countMatrix     = countMatrix, 
					species         = "mouse",
					outputDirectory = outputDirectory,
					tSNEList = list()), regexp = expM)
			
			expM <- paste0("The elements in TsneList slot don't have the same ",
					"number of cells or the same class")
			expect_error(scRNAseq(experimentName = experimentName,
							countMatrix     = countMatrix,
							species         = "mouse",
							outputDirectory = outputDirectory,
							tSNEList = tsneListWrong), regexp = expM)
		
			expM <- "Coordinates should be a matrix with two columns X and Y."
			expect_error(scRNAseq(experimentName = experimentName,
							countMatrix     = countMatrix,
							species         = "mouse",
							outputDirectory = outputDirectory,
							tSNEList = list(Tsne(name = "test", pc = 30,
											perplexity = 4,
											coordinates = matrix(seq_len(9), 
													ncol=3)))), regexp = expM)
			
			expM <- paste0("dbscanList is empty. This should be a list of ",
					"dbScan objects.")
			expect_error(scRNAseq(experimentName = experimentName, 
					countMatrix     = countMatrix, 
					species         = "mouse",
					outputDirectory = outputDirectory,
					dbscanList = list()), regexp = expM)
			
			expM <- paste0("The elements in DbscanList slot don't have the ",
					"same number of cells or the same class")
			expect_error(scRNAseq(experimentName = experimentName, 
					countMatrix     = countMatrix, 
					species         = "mouse",
					outputDirectory = outputDirectory,
					dbscanList = dbscanListWrong), regexp = expM)
	
			expM <- paste0("'cellsSimilarityMatrix' slot should contain a ",
					"square matrix.")
			expect_error(scRNAseq(experimentName = experimentName, 
			countMatrix     = countMatrix, 
			species         = "mouse",
			outputDirectory = outputDirectory,
			cellsSimilarityMatrix = csm[1:2,]), regexp = expM)

			expM <- paste0("'clustersSimilarityMatrix' slot should contain a ",
					"square matrix.")
			expect_error(scRNAseq(experimentName = experimentName, 
				countMatrix     = countMatrix, 
				species         = "mouse",
				outputDirectory = outputDirectory,
				clustersSimilarityMatrix = csm[1:2,]), regexp = expM)

			expM <- paste0("'clustersSimiliratyOrdered' slot should contain ",
					"the same clusters as 'clustersSimilarityMatrix'.")
			expect_error(scRNAseq(experimentName = experimentName,
							countMatrix     = countMatrix,
							species         = "mouse",
							outputDirectory = outputDirectory,
							clustersSimilarityMatrix = csm,
							clustersSimiliratyOrdered = factor(c(15,16,17))), 
					regexp = expM)
			
			expM <- paste0("markerGenesList is empty. This should be a list ",
					"of dataframe")
			expect_error(scRNAseq(experimentName = experimentName, 
					countMatrix     = countMatrix, 
					species         = "mouse",
					outputDirectory = outputDirectory,
					markerGenesList = list()), regexp = expM)
			
			expM <- paste0("'markerGenesList' should contain as many ",
					"dataframes as clusters found. Number of dataframes :9 ",
					"and the number of cluters found is :10.") 
			expect_error(scRNAseq(experimentName = experimentName, 
					countMatrix     = countMatrix, 
					species         = "mouse",
					outputDirectory = outputDirectory,
					clustersSimiliratyOrdered = orderedCLusters,
					markerGenesList = markers[-1]), regexp = expM)
	
			expM <- paste0("clusterMarkers should have the same number of ",
					"clusters than the number of clusters found. Nb clusters ",
					"for markers: 2. Nb of clusters: 10")
			expect_error(setClustersMarkers(scrFinal) <- data.frame(
							geneName= c("gene1", "gene2"), clusters=c(1,2)), 
					expM)
			
			expM <- paste0("The clusterMarkers data frame should have the ",
					"columns 'geneName' and 'clusters'")
			expect_error(setClustersMarkers(scrFinal) <- 
							data.frame(geneNam= rep("gene1", 10), 
									clust=seq_len(10)), expM)
			
			expM <- "clusterMarkers is empty. This should be a dataframe"
			expect_error(scRNAseq(experimentName = experimentName, 
					countMatrix     = countMatrix, 
					species         = "mouse",
					outputDirectory = outputDirectory,
					clustersMarkers = data.frame()), expM)
	
			expM <- "genesInfos is empty. This should be a dataframe"
			expect_error(scRNAseq(experimentName = experimentName,
							countMatrix     = countMatrix,
							species         = "mouse",
							outputDirectory = outputDirectory,
							genesInfos=data.frame()), expM)
			
			expM <- paste0("The genesInfos data frame should have the columns:", 
					" uniprot_gn_symbol;clusters;external_gene_name;go_id;",
					"entrezgene_description;gene_biotype;chromosome_name;",
					"Symbol;ensembl_gene_id;entrezgene_id;uniprot_gn_id;",
					"mgi_description;mgi_id")
			expect_error(scRNAseq(experimentName = experimentName,
							countMatrix     = countMatrix,
							species         = "mouse",
							outputDirectory = outputDirectory,
							genesInfos=data.frame(test="test")), expM)
			
			expM <- paste0("genesInfos should have the same number of clusters",
					" than the number of clusters found. Nb clusters for ",
					"genesInfos: 2. Nb of clusters: 10")
			expect_error(setGenesInfos(scrInfos) <- wrongInfo, expM)
	
		})
		  
		  
		  

###########################  Normalization  ###################################

test_that("Normalization works properly", {
    
			expect_match(class(sceNorm), class(expectedNormalizedMatrix))
			
			expect_equal(Biobase::exprs(sceNorm),
					Biobase::exprs(expectedNormalizedMatrix))
			
			expM <- "'sizes' parameter should be a vector of numeric values."
			expect_error(normaliseCountMatrix(scr, sizes="test"), regexp=expM)
			
			expM <- "'rowdata' parameter should be a data frame or be NULL."
			expect_error(normaliseCountMatrix(scr, rowdata="test"), regexp=expM)
			
			expM <- "'coldata' parameter should be a data frame or be NULL."
			expect_error(normaliseCountMatrix(scr, coldata="test"), regexp=expM)
			
			expM <- "'alreadyCellFiltered' parameter should be a boolean."
			expect_error(normaliseCountMatrix(scr, alreadyCellFiltered="test"), 
					regexp=expM)
			
			expM <- "'runQuickCluster' parameter should be a boolean."
			expect_error(normaliseCountMatrix(scr, runQuickCluster="test"), 
					regexp=expM)
			
			expM <- paste0("species should be 'mouse' or 'human'. ",
					"'melanogaster' is currently not supported.")
			expect_error(normaliseCountMatrix(scRNAseq(
									experimentName = experimentName, 
									countMatrix     = countMatrix, 
									species         = "melanogaster")), 
					regexp=expM)
})


################################  TSNE  #######################################


test_that("Tsne works properly", {
   
			expect_match(class(tsneList), "list")
			expect_match(unique(sapply(tsneList, class)), "Tsne")
			
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'generateTSNECoordinates' function doesn't have its ",
					"'sceNorm' slot updated. Please use 'normaliseCountMatrix'",
					" on the object before.") 
			expect_error(generateTSNECoordinates(scr), regexp=expM)
			
			expM <- "'randomSeed' parameter should be an integer."
			expect_error(generateTSNECoordinates(scrNorm, randomSeed="string"),
					regexp=expM)
			
			expM <- "'cores' parameter should be an integer"	 
			expect_error(generateTSNECoordinates(scrNorm, cores="string"),
					regexp=expM)	 
			
			expM <- "'PCs' parameter should be a vector of numeric." 
			expect_error(generateTSNECoordinates(scrNorm, PCs=c("str1", "str2"),
							regexp=expM))
			
			expM <- "'perplexities' parameter should be a vector of numeric."
			expect_error(generateTSNECoordinates(scrNorm,
							perplexities=c("str1", "str2")),
					regexp=expM)
			
			expM <- "'writeOutput' parameter should be a boolean."
			expect_error(generateTSNECoordinates(scrNorm, writeOutput="str"),
					regexp=expM)
			
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'generateTSNECoordinates' function doesn't have its ",
					"'sceNorm' slot updated. Please use 'normaliseCountMatrix'",
					" on the object before.")
			expect_error(generateTSNECoordinates(scr, cores=5),
			             regexp=expM)
})



#################################  Dbscan  #####################################

test_that("Dbscan works properly", {
			
			# Test class of the output
			expect_match(class(dbscanList), "list")
			
			# Test the class of the first element
			expect_match(class(dbscanList[[1]]), "Dbscan")
			
			# Test the class of the last element
			expect_match(class(dbscanList[[length(dbscanList)]]), "Dbscan")
			
			# Test with a empty class scRNAseq
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'runDBSCAN' function doesn't have its 'sceNorm' slot ",
					"updated. Please use 'normaliseCountMatrix' on the object ",
					"before.")
			expect_error(runDBSCAN(scr), regexp=expM)
			
			## Test with incorrect cores
			expM <- "'cores' parameter should be an integer"
			expect_error(runDBSCAN(scrTsne, cores="1"), regexp=expM)
			
			## Test with incorrect epsilon
			expM <- "'epsilon' parameter should be a vector of numeric"
			epsVec <- c("str1", "str2")
			expect_error(runDBSCAN(scrTsne, epsilon=epsVec, regexp=expM))
			
			## Test with incorrect minPoints
			expM <- "'minPoints' parameter should be a vector of numeric"
			minPvec <- c("str1", "str2")
			expect_error(runDBSCAN(scrTsne,minPoints=minPvec), regexp=expM)
			
			## Test with incorrect writeOutput
			expM <- "'writeOutput' parameter should be a boolean"
			expect_error(runDBSCAN(scrTsne, writeOutput="str"), regexp=expM)
			
			expM <- paste("The 'scRNAseq' object that you're using with",
			"'runDBSCAN' function doesn't have its 'sceNorm'",
			"slot updated. Please use 'normaliseCountMatrix'",
			"on the object before.")
			expect_error(runDBSCAN(scr, cores=5),
			             regexp=expM)
			
			expM <- paste("The 'scRNAseq' object that you're using with",
			              "'runDBSCAN' function doesn't have its 'tSNEList'",
			              "slot updated. Please use 'generateTSNECoordinates'",
			              "on the object before.")
			expect_error(runDBSCAN(scrNorm, cores=5),
			             regexp=expM)
})



#########################  clusterCellsInternal  ###############################

test_that("clusterCellsInternal works properly", {
			
			## Test class of the output
			expect_equal(class(cci), c("matrix", "array"))

			## Test with a empty sceNorm slot
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'clusterCellsInternal' function doesn't have its ",
					"'sceNorm' slot updated. Please use 'normaliseCountMatrix'",
					" on the object before.")
			expect_error(clusterCellsInternal(scr), regexp=expM)
    
			## Test with a empty dbscan slot
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'clusterCellsInternal' function doesn't have its ",
					"'dbscanList' slot updated. Please use 'runDBSCAN' on the ",
					"object before.")
			expect_error(clusterCellsInternal(scrNorm), regexp=expM)
   		
			## Test with incorrect clusterNumber
			expM <- "'clusterNumber' parameter should be a numeric."
			expect_error(clusterCellsInternal(scrDbscan, clusterNumber="str1"),
					regexp=expM)
					
			## Test with incorrect deepSplit
			expM <- "'deepSplit' parameter should be a numeric."
			expect_error(clusterCellsInternal(scrDbscan, deepSplit="str1"),
					regexp=expM)
			
			## Test with incorrect cores
			expM <- "'cores' parameter should be a numeric."
			expect_error(clusterCellsInternal(scrDbscan, cores="str1"),
					regexp=expM)
			
			## Test with incorrect clusteringMethod
			expM <- paste0("'clusteringMethod' should be one of: ward.D; ",
					"ward.D2; single; complete; average; mcquitty; median; ",
					"centroid")
			expect_error(clusterCellsInternal(scrDbscan, 
							clusteringMethod="str1"), regexp=expM)
})



#########################  calculateClustersSimilarity  ########################

test_that("calculateClustersSimilarity works properly", {
		
	## Test with incorrect clusteringMethod
	expM <- paste0("'clusteringMethod' should be one of: ward.D; ",
			"ward.D2; single; complete; average; mcquitty; median; ",
			"centroid")
    expect_error(calculateClustersSimilarity(scrDbscan, 
					clusteringMethod="str1"), regexp=expM)
    
	## Test with unnormalized object 
	expM <- paste0("The 'scRNAseq' object that you're using with ",
			"'calculateClustersSimilarity' function doesn't have its ",
			"'sceNorm' slot updated. Please use 'normaliseCountMatrix' on",
			" the object before.")
	expect_error(calculateClustersSimilarity(scr), regexp=expM)
	
    ## Test with default normalizeCountMatrix 
    expM <- paste0("The 'scRNAseq' object that you're using with ",
			"'calculateClustersSimilarity' function doesn't have a correct ",
			"'sceNorm' slot. This slot should be a 'SingleCellExperiment' ",
			"object containing 'clusters' column in its colData. Please check ",
			"if you correctly used 'clusterCellsInternal' on the object.")
    expect_error(calculateClustersSimilarity(scrDbscan), regexp=expM)
    
	## Test with  
	expM <- paste0("The 'scRNAseq' object that you're using with ",
			"'calculateClustersSimilarity' function doesn't have its ",
			"'cellsSimilarityMatrix' slot updated by clusterCellsInternal. ",
			"Please use 'clusterCellsInternal' on the object before.")
	expect_error(calculateClustersSimilarity(scrCCiwrong), regexp=expM)
	
})



################################## Plotting ####################################

test_that("plotting methods work properly", {
			
			expM <- "columnName should be: clusters, noColor, or state."
			expect_error(plotClusteredTSNE(scrFinal, columnName="toto"), expM)
			
			expM <- paste0("The number of elements of TSNEList is not equal ",
					"to PCs\\*perplexities. Contact the developper.")
			expect_error(plotClusteredTSNE(scrFinalWrong), expM)

		})


test_that("plotCellSimilarity work properly", {
    ## Test with object doesn't have consensus clusters
    expM <- paste0("You have to calculate the cells similarity", 
                   " matrix before plotting.")
    expect_error(plotCellSimilarity(scr), expM)
    
    ## Test with incorrect colorPalette
    expM <- paste0("The number of clusters is greater than the number of",
                   " given colors.")
    expect_error(plotCellSimilarity(scrFinal, colorPalette="str1" ), expM)
    
    ## Test with incorrect statePalette
    expM <- paste0("The number of clusters is greater than the number of",
                   " given colors.")
    expect_error(plotCellSimilarity(scrFinal, statePalette="str1" ), expM)
    
    ## Test with incorrect clusteringMethod
    expM <- paste0("invalid clustering method")
    expect_error(plotCellSimilarity(scrFinal, clusteringMethod="str1" ), expM)
    
    ## Test with incorrect orderClusters
    expM <- paste0("orderClusters should be a boolean.")
    expect_error(plotCellSimilarity(scrFinal, orderClusters="str1" ), expM)

    ## Test with incorrect plotPDF
    expM <- paste0("plotPDF should be a boolean.")
    expect_error(plotCellSimilarity(scrFinal, plotPDF="str1" ), expM)
    
    ## Test with incorrect returnPlot
    expM <- paste0("returnPlot should be a boolean.")
    expect_error(plotCellSimilarity(scrFinal, returnPlot="str1" ), expM)
    
    ## Test with incorrect width
    expM <- paste0("width should be a numeric.")
    expect_error(plotCellSimilarity(scrFinal, width="str1" ), expM)
    
    ## Test with incorrect height
    expM <- paste0("height should be a numeric.")
    expect_error(plotCellSimilarity(scrFinal, height="str1" ), expM)

    ## Test with incorrect onefile
    expM <- paste0("onefile should be a boolean.")
    expect_error(plotCellSimilarity(scrFinal, onefile="str1" ), expM)
    
    ## Test with incorrect showRowNames
    expM <- paste0("showRowNames should be a boolean.")
    expect_error(plotCellSimilarity(scrFinal, showRowNames="str1" ), expM)
    
    ## Test with incorrect showColnames
    expM <- paste0("showColnames should be a boolean.")
    expect_error(plotCellSimilarity(scrFinal, showColnames="str1" ), expM)
    
    ## Test with incorrect fontsize
    expM <- paste0("fontsize should be a numeric.")
    expect_error(plotCellSimilarity(scrFinal, fontsize="str1" ), expM)
    
    ## Test with incorrect fontsizeRow
    expM <- paste0("fontsizeRow should be a numeric.")
    expect_error(plotCellSimilarity(scrFinal, fontsizeRow="str1" ), expM)
    
    ## Test with incorrect widthPNG
    expM <- paste0("widthPNG should be a numeric.")
    expect_error(plotCellSimilarity(scrFinal, widthPNG="str1" ), expM)
    
    ## Test with incorrect heightPNG
    expM <- paste0("heightPNG should be a numeric.")
    expect_error(plotCellSimilarity(scrFinal, heightPNG="str1" ), expM)
   
})


test_that("plotClusteredTSNE work properly", {
    ## Test with no consensus clusters
    expM <- paste0("You have to calculate the cells similarity", 
                   " matrix before plotting.")
    expect_error(plotClusteredTSNE(scr), expM)
    
    ## Test with incorrect colorPalette
    expM <- paste0("The number of clusters is greater than the number of",
                   " given colors.")
    expect_error(plotClusteredTSNE(scrFinal, colorPalette="str1" ), expM)
    
    ## Test with incorrect PCs
    expM <- "'PCs' parameter should be a vector of numeric." 
    expect_error(plotClusteredTSNE(scrFinal, PCs=c("str1", "str2"),
                                   regexp=expM))
    
    ## Test with incorrect perplexities
    expM <- "'perplexities' parameter should be a vector of numeric."
    expect_error(plotClusteredTSNE(scrFinal, perplexities=c("str1", "str2")),
                 regexp=expM)
    
    ## Test with incorrect columnName
    expM <- "columnName should be: clusters, noColor, or state."
    expect_error(plotClusteredTSNE(scrFinal, columnName="toto"), expM)
    
    ## Test with incorrect returnPlot
    expM <- paste0("returnPlot should be a boolean.")
    expect_error(plotClusteredTSNE(scrFinal, returnPlot="str1" ), expM)
    
    ## Test with incorrect width
    expM <- paste0("width should be a numeric.")
    expect_error(plotClusteredTSNE(scrFinal, width="str1" ), expM)
    
    ## Test with incorrect height
    expM <- paste0("height should be a numeric.")
    expect_error(plotClusteredTSNE(scrFinal, height="str1" ), expM)
    
    ## Test with incorrect onefile
    expM <- paste0("onefile should be a boolean.")
    expect_error(plotClusteredTSNE(scrFinal, onefile="str1" ), expM)
    
})


test_that("plotCellHeatmap work properly", {
    ## Test with object doesn't have consensus clusters
    expM <- paste0("You have to calculate the cells similarity", 
                   " matrix before plotting.")
    expect_error(plotCellHeatmap(scr), expM)
    
    ## Test with incorrect fileName
    expM <- paste0("fileName should be a simple string.")
    expect_error(plotCellHeatmap(scrFinal, fileName=TRUE), expM)
    
    expM <- paste0("fileName should be a simple string.")
    expect_error(plotCellHeatmap(scrFinal, fileName="dir/file"), expM)
    
    ## Test with incorrect meanCentered
    expM <- paste0("meanCentered should be a boolean.")
    expect_error(plotCellHeatmap(scrFinal, fileName="str1", 
                                 meanCentered="str2"), expM)
    
    ## Test with incorrect orderClusters
    expM <- paste0("orderClusters should be a boolean.")
    expect_error(plotCellHeatmap(scrFinal, fileName="str1", 
                                 orderClusters="str2"), expM)

    ## Test with incorrect similarity
    expM <- paste0("similarity should be a boolean.")
    expect_error(plotCellHeatmap(scrFinal, fileName="str1", 
                                 similarity="str2"), expM)
    
    ## Test with incorrect orderGenes
    expM <- paste0("orderGenes should be a boolean.")
    expect_error(plotCellHeatmap(scrFinal, fileName="str1", 
                                 orderGenes="str2"), expM)
    
    ## Test with incorrect returnPlot
    expM <- paste0("returnPlot should be a boolean.")
    expect_error(plotCellHeatmap(scrFinal, fileName="str1", 
                                 returnPlot="str2"), expM)
    
    ## Test with incorrect saveHeatmapTable
    expM <- paste0("saveHeatmapTable should be a boolean.")
    expect_error(plotCellHeatmap(scrFinal, fileName="str1", 
                                 saveHeatmapTable="str2"), expM)
    
    ## Test with incorrect width
    expM <- paste0("width should be a numeric.")
    expect_error(plotCellHeatmap(scrFinal, fileName="str1", 
                                 width="str2"), expM)
    
    ## Test with incorrect height
    expM <- paste0("height should be a numeric.")
    expect_error(plotCellHeatmap(scrFinal, fileName="str1", 
                                 height="str2"), expM)
})


test_that("plotGeneExpression work properly", {
    ## Test with object doesn't have consensus clusters
    expM <- paste0("You have to calculate the cells similarity", 
                   " matrix before plotting.")
    expect_error(plotGeneExpression(scr), expM)
    
    ## Test with incorrect geneName
    expM <- paste0("The gene should be one of the normalized count matrix.")
    expect_error(plotGeneExpression(scrFinal, geneName = "gene1"), expM)
    
    ## Test with incorrect graphsDirectory
    expM <- paste0("graphsDirectory should be a string.")
    geneName <- as.character(getClustersMarkers(scrFinal)[1,1])
    expect_error(plotGeneExpression(scrFinal, geneName = geneName, 
                                    graphsDirectory = TRUE), expM)
    
    ## Test with incorrect returnPlot
    expM <- paste0("returnPlot should be a boolean.")
    geneName <- as.character(getClustersMarkers(scrFinal)[1,1])
    expect_error(plotGeneExpression(scrFinal, geneName = geneName, 
                                    returnPlot = "str1"), expM)
    
    ## Test with incorrect tSNEpicture
    expM <- paste0("tSNEpicture should be a integer")
    geneName <- as.character(getClustersMarkers(scrFinal)[1,1])
    expect_error(plotGeneExpression(scrFinal, geneName = geneName, 
                                    tSNEpicture = "str1"), expM)
    
    ## Test with incorrect commentName
    expM <- paste0("commentName should be a string.")
    geneName <- as.character(getClustersMarkers(scrFinal)[1,1])
    expect_error(plotGeneExpression(scrFinal, geneName = geneName, 
                                    commentName = TRUE), expM)
    
    ## Test with incorrect savePlot
    expM <- paste0("savePlot should be a boolean.")
    geneName <- as.character(getClustersMarkers(scrFinal)[1,1])
    expect_error(plotGeneExpression(scrFinal, geneName = geneName, 
                                    savePlot = "str1"), expM)
    
    ## Test with incorrect width
    expM <- paste0("width should be a numeric.")
    geneName <- as.character(getClustersMarkers(scrFinal)[1,1])
    expect_error(plotGeneExpression(scrFinal, geneName = geneName, 
                                    width = "str1"), expM)
    
    ## Test with incorrect height
    expM <- paste0("height should be a numeric.")
    geneName <- as.character(getClustersMarkers(scrFinal)[1,1])
    expect_error(plotGeneExpression(scrFinal, geneName = geneName, 
                                    height = "str1"), expM)
})


test_that("plotClustersSimilarity work properly", {
    ## Test with object doesn't have consensus clusters
    expM <- paste0("You have to calculate the cells similarity", 
                   " matrix before plotting.")
    expect_error(plotClustersSimilarity(scr), expM)
    
    ## Test with incorrect returnPlot
    expM <- paste0("returnPlot should be a boolean.")
    expect_error(plotClustersSimilarity(scrFinal, returnPlot="str1" ), expM)
    
    ## Test with incorrect width
    expM <- paste0("width should be a numeric.")
    expect_error(plotClustersSimilarity(scrFinal, width="str1" ), expM)
    
    ## Test with incorrect height
    expM <- paste0("height should be a numeric.")
    expect_error(plotClustersSimilarity(scrFinal, height="str1" ), expM)
    
    ## Test with incorrect onefile
    expM <- paste0("onefile should be a boolean.")
    expect_error(plotClustersSimilarity(scrFinal, onefile="str1" ), expM)
    
    ## Test with incorrect fontsize
    expM <- paste0("fontsize should be a numeric.")
    expect_error(plotCellSimilarity(scrFinal, fontsize="str1" ), expM)
})


################################  markers  ###################################

test_that("rankGenes method works properly", {
			
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'rankGenes' function doesn't have its 'SceNorm' slot ",
					"updated. Please use 'normaliseCountMatrix' on the object",
					" before.")
			expect_error(rankGenes(scr), expM)
			
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'rankGenes' function doesn't have a correct 'SceNorm' ",
					"slot. This slot should be a 'SingleCellExperiment' object",
					" containing 'clusters' column in its colData. Please ",
					"check if you correctly used 'clusterCellsInternal' on the",
					" object.")
			expect_error(rankGenes(scrDbscan), expM)
			
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'rankGenes' function doesn't have its ",
					"'clustersSimilarityMatrix' slot updated. Please use ",
					"'clusterCellsInternal' on the object before.")
			expect_error(rankGenes(scrCCiwrong), expM)
		})

test_that("retrieveGenesInfo method works properly", {
			
			expM <- "orderGenes should be 'initial' or 'alphabetical'."
			expect_error(retrieveGenesInfo(scr, species="mouse", 
							orderGenes="test"), expM)
			
			expM <- "species should be: mouse or human"
			expect_error(retrieveGenesInfo(scrFinal,  species = "droso"), expM)
			
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'retrieveGenesInfo' function doesn't have its 'SceNorm' ",
					"slot updated. Please use 'normaliseCountMatrix' on the ",
					"object before.")
			expect_error(retrieveGenesInfo(scr, species="mouse"), expM)
			
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'retrieveGenesInfo' function doesn't have a correct ",
					"'SceNorm' slot. This slot should be a ",
					"'SingleCellExperiment' object containing 'clusters' ",
					"column in its colData. Please check if you correctly used",
					" 'clusterCellsInternal' on the object.")
			expect_error(retrieveGenesInfo(scrDbscan, species="mouse"), expM)
			
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'retrieveGenesInfo' function doesn't have a similarity ",
					"matrix, Please use 'calculateClustersSimilarity' on the ",
					"object before.")
			expect_error(retrieveGenesInfo(scrCCI, species="mouse"), expM)
			
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'retrieveGenesInfo' does not have marker genes. Please ",
					"use 'bestClustersMarkers' before.")
			expect_error(retrieveGenesInfo(scrCSM, species="mouse"), expM)
			expect_error(retrieveGenesInfo(scrS4MG, species="mouse"), expM)
		})
			

test_that("saveGenesInfo method works properly", {
			
			expM <- paste0("Your object does not contain genes information. ",
					"Please run 'retrieveGenesInfo' before.")
			expect_error(saveGenesInfo(scr), expM)
			expect_error(saveGenesInfo(scrNorm), expM)
			expect_error(saveGenesInfo(scrTsne), expM)
			expect_error(saveGenesInfo(scrDbscan), expM)
			expect_error(saveGenesInfo(scrCCI), expM)
			expect_error(saveGenesInfo(scrCSM), expM)
			expect_error(saveGenesInfo(scrS4MG), expM)
			expect_error(saveGenesInfo(scrFinal), expM)
		})

test_that("retrieveTopClustersMarkers method works properly", {
			
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'retrieveTopClustersMarkers' function doesn't have its ",
					"'sceNorm' slot updated. Please use 'normaliseCountMatrix'",
					" on the object before.")
			expect_error(retrieveTopClustersMarkers(scr), expM)
			
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'retrieveTopClustersMarkers'function doesn't have a ",
					"correct 'sceNorm' slot. This slot should be a ",
					"'SingleCellExperiment' object containing 'clusters' ",
					"column in its colData. Please check if you correctly used",
					" 'clusterCellsInternal' on the object.")
			expect_error(retrieveTopClustersMarkers(scrNorm), expM)
			expect_error(retrieveTopClustersMarkers(scrTsne), expM)
			expect_error(retrieveTopClustersMarkers(scrDbscan),expM)
			
			expM <- paste0("Something wrong with number of clusters. It is ",
					"supposed to be equal to : 10. Current number: 1. Did you",
					" use 'calculateClustersSimilarity' and 'rankGenes'?")
			expect_error(retrieveTopClustersMarkers(scrCCI), expM)
			expect_error(retrieveTopClustersMarkers(scrCSM), expM)
			
		})			
			
			
		


##################################  Export  ####################################

test_that("exportResults works properly", {
    
    ## Test with object doesn't have consensus clusters
    expM <- paste0("You have to follow all the steps until to have the", 
                   " clusters similarity matrix to use the exportResults",
                   " function.")
    expect_error(exportResults(scr), expM)
    
    ## Test with incorrect saveCellsSimilarityMatrix
    expM <- paste0("saveCellsSimilarityMatrix should be a boolean.")
    expect_error(exportResults(scrFinal, saveCellsSimilarityMatrix="str1" ),
                 expM)
    
    ## Test with incorrect saveClustersSimilarityMatrix
    expM <- paste0("saveClustersSimilarityMatrix should be a boolean.")
    expect_error(exportResults(scrFinal, saveClustersSimilarityMatrix="str1" ),
                 expM)
    
    ## Test with incorrect saveNormalizedMatrix
    expM <- paste0("saveNormalizedMatrix should be a boolean.")
    expect_error(exportResults(scrFinal, saveNormalizedMatrix="str1" ), expM)
    
    ## Test with incorrect saveColData
    expM <- paste0("saveColData should be a boolean.")
    expect_error(exportResults(scrFinal, saveColData="str1" ), expM)
    
    ## Test with incorrect saveRowData
    expM <- paste0("saveRowData should be a boolean.")
    expect_error(exportResults(scrFinal, saveRowData="str1" ), expM)
    
    ## Test with incorrect saveWorkspace
    expM <- paste0("saveWorkspace should be a boolean.")
    expect_error(exportResults(scrFinal, saveWorkspace="str1" ), expM)
    
    ## Test with incorrect saveClusteringResults
    expM <- paste0("saveClusteringResults should be a boolean.")
    expect_error(exportResults(scrFinal, saveClusteringResults="str1" ), expM)

})

##################################  testClustering  ###########################
	
test_that("testClustering works properly", {
			
			expM <- paste0("The 'scRNAseq' object that you're using with ",
					"'testClustering' method doesn't have its 'sceNorm' slot ",
					"updated. Please use 'normaliseCountMatrix' on the object",
					" before.")
			expect_error(testClustering(scr), expM)
			
			expM <- "'dbscanEpsilon' parameter should be an integer."
			expect_error(testClustering(scrNorm, dbscanEpsilon="test"), 
					expM)
			
			expM <- "'minPts' parameter should be an integer"
			expect_error(testClustering(scrNorm, dbscanEpsilon=1.4, 
							minPts="test"), expM)
			
			expM <- paste0("'perplexities' parameter should be a vector of ",
					"numeric.")
			expect_error(testClustering(scrNorm, dbscanEpsilon=1.4, 
							minPts=5, perplexities="test"), expM)
			
			expM <- "'PCs' parameter should be a vector of numeric."
			expect_error(testClustering(scrNorm, dbscanEpsilon=1.4, 
							minPts=5, perplexities=30, PCs="test"), expM)
			
			expM <- "'randomSeed' parameter should be an integer."
			expect_error(testClustering(scrNorm, dbscanEpsilon=1.4, 
							minPts=5, perplexities=30, PCs=4, 
							randomSeed="test"), expM)
			
			expM <- paste0("dbscanEpsilon, minPts, perplexities, PCs, and ",
					"randomSeed should be a single value.")
			expect_error(testClustering(scrNorm, dbscanEpsilon=c(1,21)), expM)
		})
				
	
	
	
	
	
	
	