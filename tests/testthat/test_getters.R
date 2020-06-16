source("R/AllClasses.R")
source("R/AllGenerics.R")
source("R/getters.R")
library(testthat)

scr = scRNAseq()
tsne = Tsne()
dbscan = Dbscan()

############################ scRNAseq getters ##################################

test_that("getExperimentName works properly", {
    expect_equal(getExperimentName(scr), character(0))
})

test_that("getCountMatrix works properly", {
    expect_equivalent(getCountMatrix(scr), matrix(nrow = 0, ncol = 0))
})

test_that("getSceNorm works properly", {
    expect_equal(getSceNorm(scr), SingleCellExperiment())
})

test_that("getSpecies works properly", {
    expect_equal(getSpecies(scr), character(0))
})

test_that("getOutputDirectory works properly", {
    expect_equal(getOutputDirectory(scr), character(0))
})

test_that("getTSNEList works properly", {
    expect_equal(getTSNEList(scr), list(new("Tsne")))
})

test_that("getDbscanList works properly", {
    expect_equal(getDbscanList(scr), list(new("Dbscan")))
})


############################### Tsne getters ###################################

test_that("getName works properly", {
    expect_equal(getName(tsne), character(0))
})

test_that("getPerplexity works properly", {
    expect_equal(getPerplexity(tsne), numeric(0))
})

test_that("getDbscanList works properly", {
    expect_equal(getPC(tsne), numeric(0))
})

test_that("getCoordinates works properly", {
    expect_equivalent(getCoordinates(tsne),  matrix(nrow = 0, ncol = 0))
})



############################### Dbscan getters #################################

test_that("getName works properly", {
    expect_equal(getName(dbscan), character(0))
})

test_that("getMinPoints works properly", {
    expect_equal(getMinPoints(dbscan), numeric(0))
})


test_that("getClustering works properly", {
    expect_equivalent(getClustering(dbscan),  matrix(nrow = 0, ncol = 0))
})


