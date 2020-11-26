test_that("loadColdata works properly", {
            
            coldataPath <- file.path(system.file("extdata", package = "conclus"),
                    "test_colData_filtered.tsv")
            
            expect_equal(loadColdata(file=coldataPath, columnCell="cellName",
                            header=TRUE, dec=".", sep='\t'),
                    read.delim(coldataPath))
            
            expM <- "'file' parameter should be a path of existing file."
            expect_error(loadColdata(file="coldata", columnCell=1,
                            header=TRUE, sep='\t', dec="."), regexp = expM)
            
            wrongColdataPath <- file.path(system.file("extdata", package = "conclus"),
                    "coldataWithoutCellName.tsv")
            expM <- "There is no column '10' in the table."
            expect_error(loadColdata(file=wrongColdataPath, columnCell=10,
                            header=TRUE, sep='\t', dec="."), regexp = expM)
            
            expM <- "There is no column 'cell_ID' in the table."
            expect_error(loadColdata(file=wrongColdataPath, columnCell="cell_ID",
                            header=TRUE, sep='\t', dec="."), regexp = expM)
            
            expM <- paste("'header' parameter should be a boolean. Set TRUE if the",
                    "first row of the table corresponds to the column names,",
                    "and FALSE if it doesn't.")
            expect_error(loadColdata(file=coldataPath, columnCell=10,
                            header="yes", sep='\t', dec="."), regexp = expM)
            
            expM <- paste("'dec' parameter should be the character used in the table",
                    "for decimal points, usually '.' or ',' .")
            expect_error(loadColdata(file=coldataPath, columnCell="cellName",
                            header=TRUE, sep='\t', dec=" "), regexp = expM)
            
            expM <- paste("'sep' parameter should be the character used in the table",
                    "to separate the fields. Usually it's ' ' ';' ',' or '\t'.")
            expect_error(loadColdata(file=coldataPath, columnCell="cellName",
                            header=TRUE, sep='+', dec="."), regexp = expM)
            
        })



test_that("loadCountMatrix works properly", {
            
            countmatrixPath <- file.path(system.file("extdata", package = "conclus"),
                    "Bergiers_counts_matrix_filtered.tsv")
            
            expect_equal(loadCountMatrix(file=countmatrixPath,header=TRUE, dec=".",
                            sep='\t'), 
                    as.matrix(read.delim(countmatrixPath)))
            
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
