.checkLoadMatrix <- function(file, column, header, sep, dec, type){
    
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
    
    if (!type %in% c("coldata", "countMatrix", "rowdata"))
        stop("'type' parameter should be the string 'coldata' 'rowdata' ",
                "or 'countMatrix'. ")
    
    if(is.na(column) && isTRUE(type == "coldata"))
        stop("'columnCells' parameter should be the column name ",
                "or the indice of the column containing cell names/id. ")
    
    if(is.na(column) && isTRUE(type == "rowdata"))
        stop("'columnGenes' parameter should be the column name ",
                "or the indice of the column containing gene SYMBOLS. ")
}


#' loadColdata
#' 
#' @description 
#' This function import the coldata. All parameters are mandatory and no
#' default values, so the user has to know how is the coldata before the 
#' import. This function also replaces the name of the column containing cell
#' names and informed with 'columnCells' parameter by "cellName" for the good 
#' walk of CONCLUS.
#' 
#' @usage 
#' loadColdata<- function(file, columnCells, header, sep, dec)
#'             
#' @param file Path of the coldata to load.
#' @param columnCells Column name or indice of the column containing 
#' cell names/id.
#' @param header Set TRUE if the first row of the table corresponds to the
#' column names, and FALSE if it doesn't.
#' @param sep Character used in the table to separate the fields. 
#' Usually it's ' ' ';' ',' or '\\t'.
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
#' loadColdata(file=coldataPath, columnCells="cellName", header=TRUE, 
#'             dec=".", sep="\t")
#'
#'
loadColdata<- function(file, columnCells, header, sep, dec){
    
    .checkLoadMatrix(file, columnCells, header, sep, dec, type="coldata")
    
    df <- read.delim(file=file,
            header=header,
            sep=sep,
            dec=dec,
            na.strings=c("", "NA",  "<NA>"),
            stringsAsFactors=FALSE)
    
    ## The submitted coldata should have 'cellName' column to allow the merge
    ## this one with the coldata created by conclus in method-normalisation.R
    if(isFALSE("cellName" %in% colnames(df))){
        
        if(isTRUE(columnCells %in% colnames(df)) || columnCells <= ncol(df)){
            ## Change the column name of the cell column
            names(df)[names(df) == columnCells] <- "cellName"
            ## Re-order the columns
            df <- df[, c("cellName", 
                            colnames(df)[!colnames(df) %in% "cellName"])]
            
        }else 
            stop("There is no column '", columnCells, "' in the submitted ",
                    "coldata")
    }
    
    rownames(df) <- df$cellName
    
    return(df)
}



#' loadRowdata
#' 
#' @description 
#' This function import the coldata. All parameters are mandatory and no
#' default values, so the user has to know how is the coldata before the 
#' import. This function also replaces the name of the column containing cell
#' names and informed with 'columnGenes' parameter by "cellName" for the good 
#' walk of CONCLUS.
#' 
#' @usage 
#' loadRowdata(file, columnGenes, header, sep, dec)
#'             
#' @param file Path of the coldata to load.
#' @param columnGenes Column name or indice of the column containing 
#' genes SYMBOL.
#' @param header Set TRUE if the first row of the table corresponds to the
#' column names, and FALSE if it doesn't.
#' @param sep Character used in the table to separate the fields. 
#' Usually it's ' ' ';' ',' or '\\t'.
#' @param dec Character used in the table for decimal points.
#' Usually '.' or ',' .
#'
#' @return The rowdata with nameInCountMatrix column
#' @export
#'
#' @examples
#' rowdataPath <- file.path(system.file("extdata", package = "conclus"),
#'                         "test_colData_filtered.tsv")
#'                         
#' loadColdata(file=coldataPath, columnGenes="nameInCountMatrix", header=TRUE, 
#'             dec=".", sep="\t")
#'
#'
loadRowdata <- function(file, columnGenes, header, sep, dec){
    
    .checkLoadMatrix(file, columnGenes, header, sep, dec, type="rowdata")
    
    df <- read.delim(file=file,
            header=header,
            sep=sep,
            dec=dec,
            na.strings=c("", "NA",  "<NA>"),
            stringsAsFactors=FALSE)
    
    ## The submitted coldata should have 'nameInCountMatrix' column to allow 
    ## to merge the one with the coldata created by conclus in 
    ## method-normalisation.R
    if(isFALSE("nameInCountMatrix" %in% colnames(df))){
        
        if(isTRUE(columnGenes %in% colnames(df)) || columnGenes <= ncol(df)){
            ## Change the column name of the gene column
            names(df)[names(df) == columnGenes] <- "nameInCountMatrix"
            ## Re-order the columns
            df <- df[, c("nameInCountMatrix", 
                            colnames(df)[!colnames(df) %in% 
                                            "nameInCountMatrix"])]
            
        }else 
            stop("There is no column '", columnGenes, "' in the submitted ",
                    "rowdata")
    }
    
    rownames(df) <- df$nameInCountMatrix
    
    return(df)
}


loadCountMatrix<- function(file, header, sep, dec){
    
    .checkLoadMatrix(file, column=NA, header, sep, 
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
