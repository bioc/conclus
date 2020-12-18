.checkLoadedData <- function(file, column, header, sep, dec, type){
    
    if(!is.data.frame(file) && !file.exists(file))
        stop("'file' parameter should be a path of existing file.")
    
    if (!is.logical(header))
        stop("'header' parameter should be a boolean. Set TRUE if the first ", 
                "row of the table corresponds to the column names, and FALSE ",
                "if it doesn't.")
    
    if (!is.character(sep) || !sep %in% c(" ", ",", "\t", ";"))
        stop("'sep' parameter should be the character used in the table ",
                "to separate the fields. Usually it's ' ' ';' ',' or '\t'.")
    
    if (!is.character(dec) || !dec %in% c(".", ","))
        stop("'dec' parameter should be the character used in the table ",
                "for decimal points, usually '.' or ',' .")
    
    if (!type %in% c("coldata", "countMatrix", "rowdata"))
        stop("'type' parameter should be the string 'coldata' 'rowdata' ",
                "or 'countMatrix'. ")
    
    if(!is.character(column) && isTRUE(type == "coldata"))
        stop("'columnCells' parameter should be the name of the column",
            " containing cell names/id")
    
    if(!is.character(column) && isTRUE(type == "rowdata"))
        stop("'columnGenes' parameter should be the name of the column",
                " containing gene SYMBOLS. ")
}


#' loadDataOrMatrix
#' 
#' @description 
#' This function allows to import the coldata, rowData or the countMatrix. It 
#' formats each type of data to follow the requirements of CONCLUS.
#' 
#' @usage 
#' loadDataOrMatrix(file, type, columnID=NULL, header=TRUE, sep='\t', dec=".")
#'             
#' @param file Path to the rowData, colData or Matrix.
#' @param type Values should be "coldata", "rowdata", or "countMatrix".
#' @param columnID For row and col data, column name containing cells/genes 
#' names/id. Should not be used when inporting a matrix. Default=NULL.
#' @param header Set TRUE if the first row of the table corresponds to the
#' column names, and FALSE if it doesn't. Default=TRUE.
#' @param sep Character used in the table to separate the fields. 
#' Usually it's ' ' ';' ',' or '\\t'. Default='\\t'.
#' @param dec Character used in the table for decimal points.
#' Usually '.' or ',' . Default=".".
#'
#' @return The formatted row, col data or the matrix.
#' @export loadDataOrMatrix
#'
#' @examples
#' 
#' ## ColData
#' coldataPath <- system.file("extdata/test_colData_filtered.tsv", 
#'                             package = "conclus")                 
#' loadDataOrMatrix(file=coldataPath, type="coldata", columnID="cell_ID")
#' 
#' ## RowData
#' rowdataPath <- system.file("extdata/test_rowData_filtered.tsv",
#'                             package="conclus")
#' loadDataOrMatrix(file=rowdataPath, type="rowdata", columnID="gene_ID")
#' 
#' ## CountMatrix
#' countmatrixPath <- file.path(system.file("extdata", package = "conclus"),
#'                                 "test_countMatrix.tsv")
#' loadDataOrMatrix(file=countmatrixPath, type="countMatrix")
#'
#' @author Ilyess RACHEDI and Nicolas DESCOSTES
#' @importFrom utils read.delim
#'
loadDataOrMatrix <- function(file, type, columnID=NULL, header=TRUE, sep='\t', 
        dec="."){
    
    .checkLoadedData(file, columnID, header, sep, dec, type)
    
    if(!is.data.frame(file))
        if(isTRUE(all.equal(type, "countMatrix")))
            df <- read.delim(file=file, header=header, sep=sep, dec=dec,
                    na.strings=c("", "NA",  "<NA>"), stringsAsFactors=FALSE,
                    row.names=1)
        else
            df <- read.delim(file=file, header=header, sep=sep, dec=dec,
                    na.strings=c("", "NA",  "<NA>"), stringsAsFactors=FALSE)
            
    else
        df <- file
    
    if(isTRUE(all.equal(type, "coldata")))
        refColName <- "cellName"
    else if(isTRUE(all.equal(type, "rowdata")))
        refColName <- "nameInCountMatrix"
    else
        return(as.matrix(df))
    
    if(isFALSE(refColName %in% colnames(df))){
        
        if(isTRUE(columnID %in% colnames(df))){
            
            if(isFALSE(any(is.na(df[, columnID])))){
                ## Change the column name
                names(df)[names(df) == columnID] <- refColName
                ## Re-order the columns
                df <- df[, c(refColName, 
                                colnames(df)[!colnames(df) %in% refColName])]
            }else 
                stop("There are some NA values in the column you choose. ",
                        "Please fill theses values or choose another column.")
            
        }else 
            stop("There is no column '", columnID, "' in the submitted ",
                    "data")
    }
    
    if(length(unique(df[, refColName])) != nrow(df))
        stop("IDs should be unique. Please check the selected column.")
    else
        rownames(df) <- df[[refColName]]
    
    return(df)
}
