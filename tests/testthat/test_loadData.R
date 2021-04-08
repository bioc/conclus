## Prepare col data

coldataPath <- system.file("extdata/colData.tsv", package="conclus")
expectedColData <- read.delim(file=coldataPath, header=TRUE)
colnames(expectedColData)[colnames(expectedColData) == "cell_ID"] <- "cellName"

## Creation of a wrong col data
wrongColdataNA <- expectedColData
colnames(wrongColdataNA)[1] <- "cell_ID"
wrongColdataNA$cell_ID[1] <- NA


wrongColdataDup <- wrongColdataNA
wrongColdataDup$cell_ID[1] <- "c2"

wrongColDataName <- expectedColData
colnames(wrongColDataName)[1] <- "name"

## Prepare row data
rowdataPath <- file.path(system.file("extdata", package="conclus"), 
                         "rowData.tsv")
expectedRowData <- read.delim(file=rowdataPath, header=TRUE)
colnames(expectedRowData)[colnames(expectedRowData) == "gene_ID"] <- 
                                                            "nameInCountMatrix"

## Creation of a wrong row data with NA value
wrongRowMetaDataNA <- expectedRowData
colnames(wrongRowMetaDataNA)[1] <- "gene_ID"
wrongRowMetaDataNA$gene_ID[1] <- NA

## Wrong data with duplicated ID
wrongRowMetaDataDup <- wrongRowMetaDataNA
wrongRowMetaDataDup$gene_ID[1] <- "Mcts1"

## Prepare count matrix
countmatrixPath <- file.path(system.file("extdata", package = "conclus"),
                                "countMatrix.tsv")
expectedCountMatrix <- as.matrix(read.delim(countmatrixPath, header=TRUE,
                                            row.names=1))



####################  Test loadDataOrMatrix with col data ####################

test_that("loadDataOrMatrix works properly with colData", {

    expect_equal(loadDataOrMatrix(file=coldataPath, type="coldata",
                    columnID="cell_ID"), expectedColData)

    ## Test with a NA value in the cell ID column
    expM <- paste0("There are some NA values in the column you choose. ",
            "Please fill theses values or choose another column.")
    expect_error(loadDataOrMatrix(file=wrongColdataNA, type="coldata",
                                  columnID="cell_ID"), regexp=expM)

    ## Test with duplicates in IDs
    expM <- ("IDs should be unique. Please check the selected column.")
    expect_error(loadDataOrMatrix(file=wrongColdataDup, type="coldata",
                                   columnID="cell_ID"), regexp=expM)

    expM <- "'file' parameter should be a path of existing file."
    expect_error(loadDataOrMatrix(file="coldata", type="coldata", columnID=1),
           regexp=expM)

    expM <- "There is no column 'cell_ID' in the submitted data"
    expect_error(loadDataOrMatrix(file=wrongColDataName, type="coldata",
                    columnID="cell_ID"), regexp=expM)

    expM <- paste0("'columnCells' parameter should be the name of",
                    " the column containing cell names/id")
    expect_error(loadDataOrMatrix(file=expectedColData, type="coldata",
                    columnID=10), regexp=expM)

    expM <- paste("'header' parameter should be a boolean. Set TRUE if the",
                "first row of the table corresponds to the column names,",
                "and FALSE if it doesn't.")
    expect_error(loadDataOrMatrix(file=expectedColData, type="coldata",
                    columnID="cell_ID", header="yes"), regexp=expM)

    expM <- paste("'dec' parameter should be the character used in the table",
                    "for decimal points, usually '.' or ',' .")
    expect_error(loadDataOrMatrix(file=expectedColData, type="coldata",
                    columnID="cell_ID", dec=" "), regexp=expM)

    expM <- paste("'sep' parameter should be the character used in the table",
                "to separate the fields. Usually it's ' ' ';' ',' or '\t'.")
    expect_error(loadDataOrMatrix(file=expectedColData, type="coldata",
                    columnID="cell_ID", sep='+'), regexp=expM)
})



####################  Test loadDataOrMatrix with row data ####################

test_that("loadDataOrMatrix works properly with row data", {

    expect_equal(loadDataOrMatrix(file=rowdataPath, type="rowdata",
                    columnID="gene_ID"), expectedRowData)

    colnames(expectedRowData)
    colnames(loadDataOrMatrix(file=rowdataPath, type="rowdata",
                    columnID="gene_ID"))
    
    ## Test with a NA value in the gene ID column
    expM <- paste0("There are some NA values in the column you choose. ",
                    "Please fill theses values or choose another column.")
    expect_error(loadDataOrMatrix(file=wrongRowMetaDataNA, type="rowdata",
                    columnID="gene_ID"), expM)

    ## Test with duplicate in IDs
    expM <- ("IDs should be unique. Please check the selected column.")
    expect_error(loadDataOrMatrix(file=wrongRowMetaDataDup, type="rowdata",
                    columnID="gene_ID"), expM)
    
    expM <- paste("'columnGenes' parameter should be the name of the column", 
                    "containing gene SYMBOLS.")
    expect_error(loadDataOrMatrix(file=expectedRowData, type="rowdata",
                    columnID=10), regexp=expM)

    expM <- paste("'header' parameter should be a boolean. Set TRUE if the",
                "first row of the table corresponds to the column names,",
                "and FALSE if it doesn't.")
    expect_error(loadDataOrMatrix(file=expectedColData, type="rowdata",
                    columnID="gene_ID", header="yes"), regexp=expM)

    expM <- paste("'dec' parameter should be the character used in the table",
                    "for decimal points, usually '.' or ',' .")
    expect_error(loadDataOrMatrix(file=expectedColData, type="rowdata",
                    columnID="gene_ID", dec=" "), regexp=expM)

    expM <- paste("'sep' parameter should be the character used in the table",
                "to separate the fields. Usually it's ' ' ';' ',' or '\t'.")
    expect_error(loadDataOrMatrix(file=expectedColData, type="rowdata",
                    columnID="gene_ID", sep='+'), regexp=expM)
})



##################  Test loadDataOrMatrix with count matrix ##################

test_that("loadDataOrMatrix works properly with a count matrix", {

    expect_equal(loadDataOrMatrix(file=countmatrixPath, type="countMatrix", 
                                    ignoreCellNumber=TRUE), expectedCountMatrix)
    
    ## Test with small matrix 
    expM <- paste0("Not enough cells in the count matrix. There ",
                   "should be at leat 100 cells. The current count matrix ",
                   "contains 20 cells.\n")
    expect_error(loadDataOrMatrix(file=countmatrixPath, type="countMatrix", 
                 ignoreCellNumber=FALSE), expM)


})
