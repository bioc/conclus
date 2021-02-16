## Data

outputDirectory <- "YourOutputDirectory"
experimentName <- "Bergiers"

## Load the coldata
coldataPath <- system.file("extdata/colData.tsv", package="conclus")
columnsMetaData <- loadDataOrMatrix(file=coldataPath, type="coldata",
                                    columnID="cell_type")

## Load the count Matrix
countMatrixPath <- file.path(system.file("extdata", package = "conclus"),
                            "countMatrix.tsv")

countMatrix <- loadDataOrMatrix(file=countMatrixPath, type="countMatrix",
                                ignoreCellNumber=TRUE)


## Load expected results
load(file = system.file("extdata/scrLight.Rdat", package="conclus"))
load(file = system.file("extdata/expected_normalizedMatrix.Rdat",
                        package="conclus"))


## Construction of the object

scr <- singlecellRNAseq(experimentName = experimentName,
    countMatrix     = countMatrix ,
    species         = "mouse",
    outputDirectory = outputDirectory)


############################ scRNAseq getters ##################################

test_that("getExperimentName works properly", {
    expect_equal(getExperimentName(scr), "Bergiers")
})

test_that("getCountMatrix works properly", {
    expect_equivalent(getCountMatrix(scr), countMatrix)
})

test_that("getSpecies works properly", {
    expect_equal(getSpecies(scr), "mouse")
})

test_that("getOutputDirectory works properly", {
    expect_equal(getOutputDirectory(scr), outputDirectory)
})

test_that("getSceNorm works properly", {
    expect_equal(getSceNorm(scr), SingleCellExperiment())
})

test_that("getTSNEList works properly", {
    expect_equal(getTSNEList(scr), list(new("Tsne")))
})

test_that("getDbscanList works properly", {
    expect_equal(getDbscanList(scr), list(new("Dbscan")))
})

test_that("getCellsSimilarityMatrix works properly", {
    mat <- matrix(nrow = 1, ncol = 1, dimnames = list("c1", "c1"), data = 1)
    expect_equal(getCellsSimilarityMatrix(scr), mat)
})

test_that("getClustersSimilarityMatrix works properly", {
    mat <- matrix(nrow = 1, ncol = 1, dimnames = list("1", "1"), data = 1)
    expect_equal(getClustersSimilarityMatrix(scr), mat)
})

test_that("getClustersSimilarityOrdered works properly", {
    expect_equal(getClustersSimilarityOrdered(scr), factor(1))
})

test_that("getMarkerGenesList works properly", {
    l <- list(data.frame(Gene = c("gene1"), mean_log10_fdr = c(NA),
                            n_05 = c(NA), score = c(NA)))
    expect_equal(getMarkerGenesList(scr), l)
})

test_that("getClustersMarkers works properly", {
    df <- data.frame(geneName="gene1", clusters=NA)
    expect_equal(getClustersMarkers(scr), df)
})

test_that("getGenesInfos works properly", {
    df <- data.frame(uniprot_gn_symbol=c("symbol"), clusters="1",
        external_gene_name="gene", go_id="GO1,GO2",
        mgi_description="description", entrezgene_description="descr",
        gene_biotype="gene", chromosome_name="1", Symbol="symbol",
        ensembl_gene_id="ENS", mgi_id="MGI", entrezgene_id="1",
        uniprot_gn_id="ID")
    expect_equal(getGenesInfos(scr), df)
})

############################### Tsne getters ###################################

name="Bergiers_tsne_coordinates_1_4PCs_30perp"
pc=4
perplexity=30
coordinates=matrix(data=c(1,2), dimnames=list(NA,c("X", "Y")), ncol=2)

tsne <- new("Tsne", name=name, pc=pc, perplexity=perplexity,
        coordinates=coordinates)

test_that("getName works properly", {
    expect_equal(getName(tsne), name)
})

test_that("getPerplexity works properly", {
    expect_equal(getPerplexity(tsne), perplexity)
})

test_that("getDbscanList works properly", {
    expect_equal(getPC(tsne), pc)
})

test_that("getCoordinates works properly", {
    expect_equivalent(getCoordinates(tsne), coordinates)
})


# ############################### Dbscan getters #################################

name <- "Clustering_1"
minPoints <- 1.3
epsilon <- 3
clustering <- matrix(data=seq(4), dimnames=list(c("clust.1", "clust.2"),
                                                c("c1", "c2")), ncol=2)
dbscan <- new("Dbscan", name=name, minPoints=minPoints, epsilon=epsilon,
        clustering=clustering)

test_that("getName works properly", {
    expect_equal(getName(dbscan), name)
})

test_that("getMinPoints works properly", {
    expect_equal(getMinPoints(dbscan), minPoints)
})

test_that("getEpsilon works properly", {
    expect_equal(getEpsilon(dbscan), epsilon)
})

test_that("getClustering works properly", {
    expect_equivalent(getClustering(dbscan), clustering)
})
