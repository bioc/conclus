source("R/AllClasses.R")
source("R/AllGenerics.R")
source("R/setters.R")
library(testthat)

scr = scRNAseq()
tsne = Tsne()
dbscan = Dbscan()

############################ scRNAseq setters ##################################

setExperimentName(scr) <- "Bergiers"
test_that("setExperimentName works properly", {
    expect_equal(, getExperimentName(scr))
})

test_that("setCountMatrix works properly", {
    expect_equivalent(setCountMatrix(scr), matrix(nrow = 0, ncol = 0))
})

test_that("setNormalizedCountMatrix works properly", {
    expect_equal(setNormalizedCountMatrix(scr), SingleCellExperiment())
})

test_that("setSpecies works properly", {
    expect_equal(setSpecies(scr), character(0))
})

test_that("setOutputDirectory works properly", {
    expect_equal(setOutputDirectory(scr), character(0))
})

test_that("setTSNEList works properly", {
    expect_equal(setTSNEList(scr), list(new("Tsne")))
})

test_that("setDbscanList works properly", {
    expect_equal(setDbscanList(scr), list(new("Dbscan")))
})


############################### Tsne setters ###################################

test_that("setName works properly", {
    expect_equal(setName(tsne), character(0))
})

test_that("setPerplexity works properly", {
    expect_equal(setPerplexity(tsne), numeric(0))
})

test_that("setDbscanList works properly", {
    expect_equal(setPC(tsne), numeric(0))
})

test_that("setCoordinates works properly", {
    expect_equivalent(setCoordinates(tsne),  matrix(nrow = 0, ncol = 0))
})



############################### Dbscan setters #################################

test_that("setName works properly", {
    expect_equal(setName(dbscan), character(0))
})

test_that("setMinPoints works properly", {
    expect_equal(setMinPoints(dbscan), numeric(0))
})


test_that("setClustering works properly", {
    expect_equivalent(setClustering(dbscan),  matrix(nrow = 0, ncol = 0))
})

