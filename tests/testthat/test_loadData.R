coldataPath <- file.path(system.file("extdata", package = "conclus"),
        "test_colData_filtered.tsv")
load(file = system.file("extdata/expected_colData.Rdat", 
                package="conclus"))

wrongColdataNA <- expectedColData
wrongColdataNA$cellName[1] <- NA
colnames(wrongColdataNA)[1] <- "cell_ID"

#########################  Test of loadColdata  ################################

test_that("loadColdata works properly", {
    
    expect_equal(loadColdata(file=coldataPath, columnCell="cell_ID",
                            header=TRUE, dec=".", sep='\t'),
                expectedColData)
    
    ## Test with a NA value in the cell ID column
    expM <- paste0("There are some NA values in the column you choose. ",
            "Please fill theses values or choose another column.")
    expect_error(loadColdata(file=wrongColdataNA, columnCell="cell_ID",
                    header=TRUE, dec=".", sep='\t'), expM)

    ## Test with duplicates in IDs
    wrongColdataPath <- file.path(system.file("extdata", package = "conclus"),
                                    "test_wrong_colData_duplicate_id.tsv")
    expM <- ("Cell IDs should be unique. Please check the selected column.")
    expect_error(loadColdata(file=wrongColdataPath, columnCell="cell_ID",
                                header=TRUE, dec=".", sep='\t'), expM)
    
    
    expM <- "'file' parameter should be a path of existing file."
    expect_error(loadColdata(file="coldata", columnCell=1,
                            header=TRUE, sep='\t', dec="."), regexp = expM)
    
    expM <- "There is no column 'cellName' in the submitted coldata"
    expect_error(loadColdata(file=coldataPath, columnCell="cellName",
                            header=TRUE, sep='\t', dec="."), regexp = expM)

    
    expM <- paste0("'columnCells' parameter should be the name of", 
                    " the column containing cell names/id")
    expect_error(loadColdata(file=coldataPath, columnCell=10,
                            header=TRUE, sep='\t', dec="."), regexp = expM)
    
    
    expM <- paste("'header' parameter should be a boolean. Set TRUE if the",
                "first row of the table corresponds to the column names,",
                "and FALSE if it doesn't.")
    expect_error(loadColdata(file=coldataPath, columnCell="cell",
                            header="yes", sep='\t', dec="."), regexp = expM)


    expM <- paste("'dec' parameter should be the character used in the table",
                    "for decimal points, usually '.' or ',' .")
    expect_error(loadColdata(file=coldataPath, columnCell="cell",
                            header=TRUE, sep='\t', dec=" "), regexp = expM)


    expM <- paste("'sep' parameter should be the character used in the table",
                "to separate the fields. Usually it's ' ' ';' ',' or '\t'.")
    expect_error(loadColdata(file=coldataPath, columnCell="cell",
                            header=TRUE, sep='+', dec="."), regexp = expM)

})



#########################  Test of loadRowdata  ################################

test_that("loadRowdata works properly", {

    rowdataPath <- file.path(system.file("extdata", package = "conclus"),
                                    "test_rowData_filtered.tsv")
    
    load(file = system.file("extdata/expected_rowData.Rdat", 
                package="conclus"))
    expect_equal(loadRowdata(file=rowdataPath, columnGene="gene_ID",
                            header=TRUE, dec=".",
                                 sep='\t'),
                    expectedRowData)
    

    ## Test with a NA value in the gene ID column
    
    wrongRowdataPath <- file.path(system.file("extdata", package = "conclus"),
                                    "test_wrong_rowData_NA.tsv")

    expM <- paste0("There are some NA values in the column you choose. ",
                    "Please fill theses values or choose another column.")
    expect_error(loadRowdata(file=wrongRowdataPath, columnGene="gene_ID",
                                header=TRUE, dec=".", sep='\t'), expM)

    ## Test with duplicate in IDs

    wrongRowdataPath <- file.path(system.file("extdata", package = "conclus"),
                                    "test_wrong_rowData_duplicate_id.tsv")
    expM <- ("Gene IDs should be unique. Please check the selected column.")
    expect_error(loadRowdata(file=wrongRowdataPath, columnGene="gene_ID",
                                header=TRUE, dec=".", sep='\t'), expM)
    
    expM <- "'file' parameter should be a path of existing file."
    expect_error(loadRowdata(file="file", header=TRUE, dec=".",
                sep='\t'), regexp = expM)

    expM <- paste("'header' parameter should be a boolean. Set TRUE if the",
                "first row of the table corresponds to the column names,",
                "and FALSE if it doesn't.")
    expect_error(loadRowdata(file=rowdataPath, header="yes", dec=".",
                sep='\t'), regexp = expM)

    expM <- paste("'dec' parameter should be the character used in the table",
                    "for decimal points, usually '.' or ',' .")
    expect_error(loadRowdata(file=rowdataPath, header=TRUE, dec=" ",
                sep='\t'), regexp = expM)
    
    expM <- paste("'sep' parameter should be the character used in the table",
                "to separate the fields. Usually it's ' ' ';' ',' or '\t'.")
    expect_error(loadRowdata(file=rowdataPath, header=TRUE, sep='+',
                                dec="."), regexp = expM)

})


test_that("loadCountMatrix works properly", {

    countmatrixPath <- file.path(system.file("extdata", package = "conclus"),
                                    "test_countMatrix.tsv")

    load(file = system.file("extdata/expected_countMatrix.Rdat", 
                package="conclus"))
    expect_equal(loadCountMatrix(file=countmatrixPath, header=TRUE, dec=".",
                                 sep='\t'), 
                expectedCountMatrix)

    expM <- "'file' parameter should be a path of existing file."
    expect_error(loadCountMatrix(file="countMatrix", header=TRUE, dec=".",
                sep='\t'), regexp = expM)

    expM <- paste("'header' parameter should be a boolean. Set TRUE if the",
                "first row of the table corresponds to the column names,",
                "and FALSE if it doesn't.")
    expect_error(loadCountMatrix(file=countmatrixPath, header="yes", dec=".",
                sep='\t'), regexp = expM)

    expM <- paste("'dec' parameter should be the character used in the table",
                    "for decimal points, usually '.' or ',' .")
    expect_error(loadCountMatrix(file=countmatrixPath, header=TRUE, dec=" ",
                sep='\t'), regexp = expM)
    
    expM <- paste("'sep' parameter should be the character used in the table",
                "to separate the fields. Usually it's ' ' ';' ',' or '\t'.")
    expect_error(loadCountMatrix(file=countmatrixPath, header=TRUE, sep='+',
                                dec="."), regexp = expM)

})