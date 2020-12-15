## Prepare col data

coldataPath <- system.file("extdata/test_colData_filtered.tsv", 
        package="conclus")
columnsMetaData <- read.delim(file=coldataPath, header=TRUE)
wrongColdataNA <- columnsMetaData
wrongColdataNA$cellName[1] <- NA
colnames(wrongColdataNA)[1] <- "cell_ID"
wrongColdataDup <- wrongColdataNA
wrongColdataDup$cellName[1] <- "c2"
wrongColDataName <- columnsMetaData
colnames(wrongColDataName)[1] <- "cell_ID"

## Prepare row data
rowdataPath <- file.path(system.file("extdata", package = "conclus"),
        "test_rowData_filtered.tsv")
rowMetaData <- read.delim(file=rowdataPath, header=TRUE)
load(file = system.file("extdata/expected_rowData.Rdat", package="conclus"))
wrongRowMetaDataNA <- rowMetaData
wrongRowMetaDataNA$gene_ID[476] <- NA
wrongRowMetaDataDup <- rowMetaData
wrongRowMetaDataDup$gene_ID[1] <- "Tspan32"

## Prepare count matrix
countmatrixPath <- file.path(system.file("extdata", package = "conclus"),
        "test_countMatrix.tsv")
expectedCountMatrix <- as.matrix(read.delim(countmatrixPath, header=TRUE, 
                row.names=1))

####################  Test loadDataOrMatrix with col data ####################

test_that("loadDataOrMatrix works properly with colData", {
    
    expect_equal(loadDataOrMatrix(file=coldataPath, type="coldata", 
                    columnID="cell_ID"), columnsMetaData)
    
    ## Test with a NA value in the cell ID column
    expM <- paste0("There are some NA values in the column you choose. ",
            "Please fill theses values or choose another column.")
    expect_error(loadDataOrMatrix(file=wrongColdataNA, type="coldata", 
                    columnID="cell_ID"), expM)

    ## Test with duplicates in IDs
    expM <- ("IDs should be unique. Please check the selected column.")
    expect_error(loadDataOrMatrix(file=wrongColdataDup, type="coldata", 
                    columnID="cell_ID"), expM)
        
    expM <- "'file' parameter should be a path of existing file."
    expect_error(loadDataOrMatrix(file="coldata", type="coldata", columnID=1),
            regexp = expM)
    
    expM <- "There is no column 'cellName' in the submitted data"
    expect_error(loadDataOrMatrix(file=wrongColDataName, type="coldata", 
                    columnID="cellName"), regexp = expM)

    expM <- paste0("'columnCells' parameter should be the name of", 
                    " the column containing cell names/id")
    expect_error(loadDataOrMatrix(file=columnsMetaData, type="coldata", 
                    columnID=10), regexp = expM)
    
    expM <- paste("'header' parameter should be a boolean. Set TRUE if the",
                "first row of the table corresponds to the column names,",
                "and FALSE if it doesn't.")
    expect_error(loadDataOrMatrix(file=columnsMetaData, type="coldata", 
                    columnID="cell_ID", header="yes"), regexp = expM)

    expM <- paste("'dec' parameter should be the character used in the table",
                    "for decimal points, usually '.' or ',' .")
    expect_error(loadDataOrMatrix(file=columnsMetaData, type="coldata", 
                    columnID="cell_ID", dec=" "), regexp = expM)

    expM <- paste("'sep' parameter should be the character used in the table",
                "to separate the fields. Usually it's ' ' ';' ',' or '\t'.")
    expect_error(loadDataOrMatrix(file=columnsMetaData, type="coldata", 
                    columnID="cell_ID", sep='+'), regexp = expM)
})



####################  Test loadDataOrMatrix with row data ####################

test_that("loadDataOrMatrix works properly with row data", {
    
    expect_equal(loadDataOrMatrix(file=rowMetaData, type="rowdata", 
                    columnID="gene_ID"), expectedRowData)
    
    ## Test with a NA value in the gene ID column
    expM <- paste0("There are some NA values in the column you choose. ",
            "Please fill theses values or choose another column.")
    expect_error(loadDataOrMatrix(file=wrongRowMetaDataNA, type="rowdata", 
                    columnID="gene_ID"), expM)

    ## Test with duplicate in IDs
    expM <- ("IDs should be unique. Please check the selected column.")
    expect_error(loadDataOrMatrix(file=wrongRowMetaDataDup, type="rowdata", 
                    columnID="gene_ID"), expM)
})


##################  Test loadDataOrMatrix with count matrix ##################

test_that("loadDataOrMatrix works properly with a count matrix", {

    expect_equal(loadDataOrMatrix(file=countmatrixPath, type="countMatrix"), 
                expectedCountMatrix)
})
