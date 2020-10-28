.checkLoadMatrix <- function(file, columnCell, header, sep, dec, type){
    
    if(!file.exists(file))
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
    
    if (!type %in% c("coldata", "countMatrix"))
        stop("'type' parameter should be the string 'coldata' ",
                "or 'countMatrix'. ")
    
    if(is.na(columnCell) && isTRUE(columnCell == "coldata"))
        stop("'columnCell' parameter should be the column name ",
                "or the indice of the column containing cell names/id. ")
}


#' loadColData
#' 
#' @description 
#' This function import the coldata. All parameters are mandatory and no
#' default values, so the user has to know how is the coldata before the 
#' import. This function also replaces the name of the column containing cell
#' names and informed with 'columnCell' parameter by "cellName" for the good 
#' walk of CONCLUS.
#' 
#' @usage 
#' loadColData(file=coldataPath, columnCell="cellName", header=TRUE, 
#'             dec=".", sep='\t')
#'             
#' @param file Path of the coldata to load.
#' @param columnCell Column name or indice of the column containing 
#' cell names/id.
#' @param header Set TRUE if the first row of the table corresponds to the
#' column names, and FALSE if it doesn't.
#' @param sep Character used in the table to separate the fields. 
#' Usually it's ' ' ';' ',' or '\t'.
#' @param dec Character used in the table for decimal points.
#' Usually '.' or ',' .
#'
#' @return The coldata with cellName column
#' @export
#'
#' @examples
#' coldataPath <- file.path(system.file("extdata", package = "conclus"),
#'                         "test_colData_filtered.tsv")
#'                         
#' loadColData(file=coldataPath, columnCell="cellName", header=TRUE, 
#'             dec=".", sep='\t')
#'
#'
loadColData<- function(file, columnCell, header, sep, dec){

    .checkLoadMatrix(file, columnCell, header, sep, dec, type="coldata")

    df <- read.delim(file=file,
                        header=header,
                        sep=sep,
                        dec=dec,
                        na.strings=c("", "NA",  "<NA>"),
                        stringsAsFactors=FALSE)
    
    ## The submitted coldata should have 'cellName' column to allow the merge
    ## this one with the coldata created by conclus in method-normalisation.R
    if(isFALSE("cellName" %in% colnames(df))){
        
        if(isTRUE(columnCell %in% colnames(df)) || columnCell <= ncol(df)){
            names(df)[names(df) == columnCell] <- "cellName"
            df <- df[, c("cellName", 
                        colnames(df)[!colnames(df) %in% "cellName"])]
            
        }else 
            stop("There is no column '", columnCell, "' in the table.")
    }
    
    if(isFALSE("state" %in% colnames(df))){
         
        exp <- grep("state", colnames(df), ignore.case = TRUE, 
                    value = TRUE)
        colnames(df)[colnames(df) == exp] <- "state"
        df$cellName <- coldata$cellName
     }
   
    rownames(df) <- df$cellName
    
    return(df)
}


loadCountMatrix<- function(file, header, sep, dec){

    .checkLoadMatrix(file, columnCell=NA, header, sep, 
                    dec, type="countMatrix")

    df <- read.delim(file=file,
                    header=header,
                    sep=sep,
                    dec=dec,
                    na.strings=c("", "NA",  "<NA>"),
                    stringsAsFactors=FALSE,
                    row.names = 1)
    
    mat <- as.matrix(df)
    
    return(mat)
}

