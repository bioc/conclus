## Data

outputDirectory <- "YourOutputDirectory"
experimentName <- "Bergiers"
columnsMetaData <- read.delim(
    system.file("extdata/test_colData_filtered.tsv", 
                package="conclus"))

## Creation of the count Matrix

countMatrix <- as.matrix(read.delim(
    system.file("extdata/test_countMatrix.tsv", 
                package="conclus")))

smallMatrix <- countMatrix[, 1:50]

## Retrieve the clustering to add
clustAddTab <- read.delim(
    system.file("extdata/Bergiers_clusters_table.tsv", package="conclus"))
clustAddTabColThree <- cbind(clustAddTab, mock=rep(1, nrow(clustAddTab)))
clustWrongName <- clustAddTab
colnames(clustWrongName) <- c("test", "test")
clustWrongcells <- clustAddTab
clustWrongcells$cells <- paste0("test", clustWrongcells$cells)

## Load expected results
## !!! Move to a data folder
load(file = system.file("extdata/scrLight.Rdat", package="conclus"))
load(file = system.file("extdata/expected_normalizedMatrix.Rdat", 
                        package="conclus"))


## Construction of the object

scr <- singlecellRNAseq(experimentName = experimentName, 
                        countMatrix     = countMatrix ,
                        species         = "mouse",
                        outputDirectory = outputDirectory)


############################ scRNAseq setters ##################################


test_that("setExperimentName works properly", {
    
    ## Setting with correct value
    setExperimentName(scr) <- "newName"
    expect_equal("newName", getExperimentName(scr))
})


test_that("setCountMatrix works properly", {
    
    ## Setting with correct matrix
    newCountMatrix <- countMatrix[1:100, 1:100]
    setCountMatrix(scr) <- newCountMatrix
    expect_equal(newCountMatrix, getCountMatrix(scr))

    ## Setting with matrix of NAs
    wrongCountMatrix <- matrix(rep(NA,1000), ncol = 100)
    expM <- paste0("The count matrix is empty or does not contain ",
            "whole numbers. Please check your count matrix.\n")
    expect_error(setCountMatrix(scr) <- wrongCountMatrix, regexp=expM)
    
    ## Setting with matrix of characters
    wrongCountMatrix <- matrix(rep("toto", 1000), ncol = 100)
    expM <- paste0("The count matrix is empty or does not contain ",
            "whole numbers. Please check your count matrix.\n")
    expect_error(setCountMatrix(scr) <- wrongCountMatrix, regexp=expM)
    
    ## Setting with matrix without rownames
    wrongCountMatrix <- matrix(seq(10000), ncol = 100)
    expM <- paste0("The name of the lines should be character class. ",
                    "Please check your count matrix.")
    expect_error(setCountMatrix(scr) <- wrongCountMatrix, regexp=expM)
    
    ## Setting with matrix without colnames
    wrongCountMatrix <- matrix(seq(10000), ncol = 100)
    rownames(wrongCountMatrix) <- paste0(rep("id_", 100) , seq(100))
    expM <- paste0("The name of the columns should be character class. ",
        "Please check your count matrix.")
    expect_error(setCountMatrix(scr) <- wrongCountMatrix, regexp=expM)
    
    ## Setting with too small matrix
    wrongCountMatrix <- countMatrix[1:50, 1:50]
    expM <- paste0("Not enough cells in the count matrix. There Should be at",
            " leat 100 cells. The current count matrix contains 50 cells.\n")
    expect_error(setCountMatrix(scr) <- wrongCountMatrix, regexp=expM)
})


test_that("setSpecies works properly", {
    ## Setting with correct species
    newSpecies <- "human"
    setSpecies(scr) <- newSpecies
    expect_equal(newSpecies, getSpecies(scr))

    ## Setting with wrong species
    wrongSpecies <- "toto"
    expM <- paste0("species should be 'mouse' or 'human'. 'toto' is currently",
            " not supported.")
    expect_error(setSpecies(scr) <- wrongSpecies, regexp=expM)
})


test_that("setOutputDirectory works properly", {
    ## Setting with correct path
    newOutputDir <- "dir"
    setOutputDirectory(scr) <- newOutputDir
    expect_equal(newOutputDir, getOutputDirectory(scr))
    
    ## Setting with wrong path
    wrongOutputDir <- "Single Cell Experiment"
    expM <- paste0("'outputDirectory' should be a conform folder path:",
                    "'Single Cell Experiment' is not.")
    expect_error(setOutputDirectory(scr) <- wrongOutputDir, regexp=expM)
    
    ## Setting with empty value
    wrongOutputDir <- character()
    expM <- "'outputDirectory' slot is empty. Please fill it."
    expect_error(setOutputDirectory(scr) <- wrongOutputDir, regexp=expM)
})


test_that("setSceNorm works properly", {
    ## Setting with correct object
    newSceNorm <- SingleCellExperiment::SingleCellExperiment()
    setSceNorm(scr) <- newSceNorm
    expect_equal(newSceNorm, getSceNorm(scr))
})


test_that("setTSNEList works properly", {
    ## Setting with correct value
    newList <- list(new("Tsne"))
    setTSNEList(scr) <- newList
    expect_equal(newList, getTSNEList(scr))
    
    ## Setting with empty list
    newList <- list()
    expM <- "tSNEList is empty. This should be a list of tSNE objects.\n"
    expect_error(setTSNEList(scr) <- newList, regexp=expM)
    
    ## Setting with wrong list
    newList <- list(1, 2, 3)
    expM <- "tSNEList should be a list of Tsne objects."
    expect_error(setTSNEList(scr) <- newList, regexp=expM)
})


test_that("setDbscanList works properly", {
    ## Setting with correct value
    newList <- list(new("Dbscan"))
    setDbscanList(scr) <- newList
    expect_equal(newList, getDbscanList(scr))
    
    ## Setting with empty list
    newList <- list()
    expM <- "dbscanList is empty. This should be a list of dbScan objects.\n"
    expect_error(setDbscanList(scr) <- newList, regexp=expM)
    
    ## Setting with wrong list
    newList <- list(1, 2, 3)
    expM <- "dbscanList should be a list of Dbscan objects."
    expect_error(setDbscanList(scr) <- newList, regexp=expM)
})


test_that("cellsSimilarityMatrix works properly", {
    ## Setting with correct matrix
    newCCI <- matrix(ncol=3, nrow=3, data=seq(9))
    rownames(newCCI) <- c("c1", "c2", "c3")
    colnames(newCCI) <- c("c1", "c2", "c3")
    setCellsSimilarityMatrix(scr) <- newCCI
    expect_equal(newCCI, getCellsSimilarityMatrix(scr))
    
    ## Setting with wrong value
    wrongCCI <- matrix(ncol=3, nrow=2, data=seq(6))
    expM <- paste0("'cellsSimilarityMatrix' should have column and row names",
                " corresponding to cell names.")
    expect_error(setCellsSimilarityMatrix(scr) <- wrongCCI, regexp=expM)

    ## Setting with matrix non identical names in row and columns
    wrongCCI <- matrix(ncol=3, nrow=3, data=seq(9))
    rownames(wrongCCI) <- c("c1", "c2", "c3")
    colnames(wrongCCI) <- c("c3", "c1", "c2")
    expM <- paste0("'cellsSimilarityMatrix' should be a square matrix with",
            " identical names in rows and columns.")
    expect_error(setCellsSimilarityMatrix(scr) <- wrongCCI, regexp=expM)
    
    
    ## Setting with matrix of non numeric values
    wrongCCI <- matrix(rep("toto", 9), ncol=3, nrow=3)
    rownames(wrongCCI) <- c("c1", "c2", "c3")
    colnames(wrongCCI) <-  c("c1", "c2", "c3")
    
    expM <- "'cellsSimilarityMatrix' should contain only numeric values."
    expect_error(setCellsSimilarityMatrix(scr) <- wrongCCI, regexp=expM)
})


test_that("clustersSimilarityMatrix works properly", {
    ## Setting with correct matrix
    newCSM<- matrix(ncol=3, nrow=3, data = seq(9))
    rownames(newCSM) <- c("1", "2", "3")
    colnames(newCSM) <- c("1", "2", "3")
    setClustersSimilarityMatrix(scr) <- newCSM
    expect_equal(newCSM, getClustersSimilarityMatrix(scr))
    
    ## Setting with wrong value
    wrongCSM <- matrix(ncol=3, nrow=2, data=seq(6))
    expM <- paste0("'clustersSimilarityMatrix' should have column and row ",
            "names corresponding to cluster names.")
    expect_error(setClustersSimilarityMatrix(scr) <- wrongCSM, regexp=expM)
    
    ## Setting with matrix non identical names in row and columns
    wrongCSM <- matrix(ncol=3, nrow=3, data=seq(9))
    rownames(wrongCSM) <- c("1", "2", "3")
    colnames(wrongCSM) <- c("3", "2", "1")
    expM <- paste0("'clustersSimilarityMatrix' should be a square matrix with",
        " identical names in rows and colums.")
    expect_error(setClustersSimilarityMatrix(scr) <- wrongCSM, regexp=expM)
    
    ## Setting with matrix of non numeric values
    wrongCSM <- matrix(rep("toto", 9), ncol=3, nrow=3)
    rownames(wrongCSM) <- c("1", "2", "3")
    colnames(wrongCSM) <- c("1", "2", "3")
    expM <- paste0("'clustersSimilarityMatrix' should contain only numeric",
        " values.")
    expect_error(setClustersSimilarityMatrix(scr) <- wrongCSM, regexp=expM)
})


test_that("markerGenesList works properly", {
    
    ## Setting with correct list
    newList <- list(data.frame(Gene = c("gene1", "gene2"),
                                mean_log10_fdr = c(1, 1),
                                n_05 = c(1,1), score = c(1,2)))
    setMarkerGenesList(scr) <- newList
    expect_equal(newList, getMarkerGenesList(scr))
    
    ## Setting with wrong list
    wrongList <- list(1, 2, 3)
    expM <- paste0("markerGenesList' slot should contain a list of dataframes",
                    " with at least following columns : 'Gene'," ,
                    " 'mean_log10_fdr', 'n_05', 'score")
    expect_error(setMarkerGenesList(scr) <- wrongList, regexp=expM)
    
})


test_that("clustersMarkers works properly", {
    
    ## Setting with correct df
    newDF <- data.frame(geneName="gene1", clusters=1)
    setClustersMarkers(scr) <- newDF
    expect_equal(newDF, getClustersMarkers(scr))
    
    ## Setting with empty df
    wrongDF <- data.frame()
    expM <- "clusterMarkers is empty. This should be a dataframe"
    expect_error(setClustersMarkers(scr) <- wrongDF, regexp=expM)
    
    ## Setting df with wrong colnames
    wrongDF <- data.frame(3)
    expM <- paste0("The clusterMarkers data frame should have the columns", 
                    " 'geneName' and 'clusters'")
    expect_error(setClustersMarkers(scr) <- wrongDF, regexp=expM)
})


test_that("genesInfos works properly", {
    
    ## Setting with correct df
    newDF <- data.frame(uniprot_gn_symbol=c("S100a6"), 
            clusters=c("3"), external_gene_name=c("S100a6"),
            go_id=c("GO:0016020, GO:0046872"), mgi_description=c("descrip1"), 
            entrezgene_description=c("20200"), gene_biotype=c("protein_coding"),
            chromosome_name=c("3"), Symbol=c("S100a6"),
            ensembl_gene_id=c("ENSMUSG00000001025"), mgi_id=c("MGI:1339467"), 
            entrezgene_id= c("20200") , uniprot_gn_id=c("P14069"))
    setGenesInfos(scr) <- newDF
    expect_equal(newDF, getGenesInfos(scr))
    
    ## Setting with incorrect number of clusters
    wrongDF <- data.frame(uniprot_gn_symbol=c("S100a6", "S100a6"), 
                clusters=c("1", "2"), external_gene_name=c("S100a6", "S100a6"),
                go_id=c("GO:0016020", "GO:0046872"),
                mgi_description=c("descrip1", "descrip2"), 
                entrezgene_description=c("20200", "19202"),
                gene_biotype=c("protein_coding", "protein_coding"),
                chromosome_name=c("3", "4"), Symbol=c("S100a6", "S100a6"),
                ensembl_gene_id=c("ENS1", "ENS2"), mgi_id=c("MGI:1", "MGI:2"), 
                entrezgene_id= c("20200", "19202"), 
                uniprot_gn_id=c("P14069", "P14069"))
    expM <- paste0("genesInfos should have the same number of clusters than ",
            "the number of clusters found. Nb clusters for genesInfos: 2. ",
            "Nb of clusters: 1")
    expect_error(setGenesInfos(scr) <- wrongDF, regexp=expM)
    
    ## Setting with empty df
    wrongDF <- data.frame()
    expM <- paste0("genesInfos is empty. This should be a dataframe")
    expect_error(setGenesInfos(scr) <- wrongDF, regexp=expM)
    
    ## Setting df with wrong colnames
    wrongDF <- data.frame(uniprot_gn_symbol=c("S100a6"), 
            clusters=c("3"), external_gene_name=c("S100a6"),
            go_id=c("GO:0016020, GO:0046872"), mgi_description=c("descrip1"), 
            entrezgene_description=c("20200"), gene_biotype=c("protein_coding"),
            chromosome_name=c("3"), Symbol=c("S100a6"),
            ensembl_gene_id=c("ENSMUSG00000001025"), mgi_id=c("MGI:1339467"), 
            entrezgene_id= c("20200"))
    expM <- paste0("The genesInfos data frame should have the columns: ",
                    "uniprot_gn_symbol;clusters;external_gene_name;go_id;",
                    "entrezgene_description;gene_biotype;chromosome_name;",
                    "Symbol;ensembl_gene_id;entrezgene_id;uniprot_gn_id;",
                    "mgi_description;mgi_id")
    expect_error(setGenesInfos(scr) <- wrongDF, regexp=expM)
})



############################### Tsne setters ###################################

tsne <- new("Tsne", coordinates=matrix(data=c(1,2),
                                        dimnames=list(NA,c("X", "Y")), ncol=2))

test_that("setCoordinates works properly", {
    
    ## Seting a correct matrix
    newMatrix <- matrix(data=seq(3), dimnames=list(seq(3), c("X", "Y")),
                                                    ncol=2, nrow=3)
    setCoordinates(tsne) <- newMatrix
    expect_equal(newMatrix, getCoordinates(tsne))

    ## Setting matrix with no column name
    wrongMatrix <- matrix(data=c(1,2), dimnames=list(NA, c("c1", "c2")), ncol=2)
    expM <- "Coordinates should be a matrix with two columns X and Y."
    expect_error(setCoordinates(tsne) <- wrongMatrix, regexp=expM)
    
    ## Setting matrix 
    wrongMatrix <- matrix(data=c("c1", "c2"),
                    dimnames=list(NA, c("X", "Y")), ncol=2)
    expM <- "Coordinates should be a matrix of numeric values."
    expect_error(setCoordinates(tsne) <- wrongMatrix, regexp=expM)
})



############################### Dbscan setters #################################

dbscan <- new("Dbscan", clustering=matrix(data=seq(4),
                                        dimnames=list(c("clust.1", "clust.2"),
                                                    c("c1", "c2")),
                                        ncol=2))

test_that("setClustering works properly", {
    
    ## Setting a correct clustering
    newClustering <- matrix(data=seq(4), 
                        dimnames=list(c("clust.3", "clust.4"), c("c5", "c12")),
                        ncol=2)
    setClustering(dbscan) <- newClustering
    expect_equal(newClustering, getClustering(dbscan))
    
    ## Setting a non conform clustering
    wrongClustering <- matrix(data=c("c1", "c2"),
        dimnames=list(NA, c("X", "Y")), ncol=2)
    expM <- "'Clustering' slot should be a matrix of integer values."
    expect_error(setClustering(dbscan) <- wrongClustering, regexp=expM)
    
})

