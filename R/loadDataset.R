.checkLoadMatrix <- function(file, column, header, sep, dec, type){
    
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


#' loadColdata
#' 
#' @description 
#' This function allows to import the coldata. All parameters are mandatory 
#' and no default values, so it is necessary to know how is the coldata before
#' importing. This function also replaces via 'columnCells' parameter the 
#' name of the column containing cell names by 'cellName' for the 
#' good walk of CONCLUS.
#' 
#' @usage 
#' loadColdata(file, columnCells="cell_ID", header=TRUE, sep='\t', dec=".")
#'             
#' @param file Path or a data.frame of the coldata to load.
#' @param columnCells Column name containing cell names/id. Default="cell_ID".
#' @param header Set TRUE if the first row of the table corresponds to the
#' column names, and FALSE if it doesn't. Default=TRUE.
#' @param sep Character used in the table to separate the fields. 
#' Usually it's ' ' ';' ',' or '\\t'. Default='\\t'.
#' @param dec Character used in the table for decimal points.
#' Usually '.' or ',' . Default=".".
#'
#' @return The coldata with cellName column
#' @export
#'
#' @examples
#' coldataPath <- system.file("extdata/test_colData_filtered.tsv", 
#'                             package = "conclus")
#'                         
#' loadColdata(file=coldataPath, columnCells="cell_ID", header=TRUE, 
#'             sep='\t', dec=".")
#'
#' @importFrom utils read.delim
#'
loadColdata<- function(file, columnCells="cell_ID", header=TRUE, sep='\t', 
        dec="."){

    .checkLoadMatrix(file, columnCells, header, sep, dec, type="coldata")

    if(!is.data.frame(file))
        df <- read.delim(file=file, header=header, sep=sep, dec=dec,
                na.strings=c("", "NA",  "<NA>"), stringsAsFactors=FALSE)
    else
        df <- file
    
    ## The submitted coldata should have 'cellName' column to allow the merge
    ## of this one with the coldata created by conclus in 
    ## method-normalisation.R
    
    if(isFALSE("cellName" %in% colnames(df))){
        
       if(isTRUE(columnCells %in% colnames(df))){
            
            if(isFALSE(any(is.na(df[, columnCells])))){
                ## Change the column name of the cell column
                names(df)[names(df) == columnCells] <- "cellName"
                ## Re-order the columns
                df <- df[, c("cellName", 
                            colnames(df)[!colnames(df) %in% "cellName"])]
            }else 
                stop("There are some NA values in the column you choose. ",
                    "Please fill theses values or choose another column.")
            
        }else 
            stop("There is no column '", columnCells, "' in the submitted ",
                    "coldata")
    }
   
    if(length(unique(df[, "cellName"])) != nrow(df))
        stop("Cell IDs should be unique. Please check the selected column.")
    else
        rownames(df) <- df$cellName
    
    return(df)
}



#' loadRowdata
#' 
#' @description 
#' This function allows to import the rowdata. All parameters are mandatory 
#' and no default values, so it is necessary to know how is the rowdata before 
#' importing. This function also replaces via 'columnGenes' parameter the 
#' name of the column containing gene names by 'nameInCountMatrix' for the 
#' good walk of CONCLUS.
#' 
#' @usage 
#' loadRowdata(file, columnGenes, header, sep, dec)
#'             
#' @param file Path of the rowdata. to load.
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
#' rowdataPath <- system.file("extdata/test_rowData_filtered.tsv",
#'                             package="conclus")
#'                         
#' loadRowdata(file=rowdataPath, columnGenes="gene_ID", header=TRUE, 
#'             sep='\t', dec=".")
#'
#' @importFrom utils read.delim
#' 
loadRowdata <- function(file, columnGenes, header, sep, dec){

    .checkLoadMatrix(file, columnGenes, header, sep, dec, type="rowdata")

    df <- read.delim(file=file,
                        header=header,
                        sep=sep,
                        dec=dec,
                        na.strings=c("", "NA", "<NA>"),
                        stringsAsFactors=FALSE)
    
    ## The submitted rowdata should have 'nameInCountMatrix' column to allow 
    ## the merge of this one with the rowdata created by conclus in 
    ## method-normalisation.R
    
    if(isFALSE("nameInCountMatrix" %in% colnames(df))){

        if(isTRUE(columnGenes %in% colnames(df))){
            
            if(isFALSE(any(is.na(df[, columnGenes])))){
                ## Change the column name of the gene column
                names(df)[names(df) == columnGenes] <- "nameInCountMatrix"
                ## Re-order the columns
                df <- df[, c("nameInCountMatrix", 
                            colnames(df)[!colnames(df) %in% "nameInCountMatrix"])]
            }else 
                stop("There are some NA values in the column you choose. ",
                    "Please fill theses values or choose another column.")
            
        }else 
            stop("There is no column '", columnGenes, "' in the submitted ",
                    "rowdata")
    }
   
    if(length(unique(df[, "nameInCountMatrix"])) != nrow(df))
        stop("Gene IDs should be unique. Please check the selected column.")
    else
        rownames(df) <- df$nameInCountMatrix
    
    return(df)
}


#' loadCountMatrix
#' 
#' @description 
#' This function allows to import the count matrix. All parameters are 
#' mandatory and no default values, so it is necessary to know how is your 
#' count matrix before importing. 
#' 
#' @usage 
#' loadCountMatrix(file, header, sep, dec)
#' 
#' @param file Path of the count matrix to load.
#' @param header Set TRUE if the first row of the table corresponds to the
#' column names, and FALSE if it doesn't.
#' @param sep Character used in the table to separate the fields. 
#' Usually it's ' ' ';' ',' or '\\t'.
#' @param dec Character used in the table for decimal points.
#' Usually '.' or ',' .
#'
#' @return The count matrix.
#' @export loadCountMatrix
#'
#' @examples
#' 
#' countmatrixPath <- system.file("extdata/Bergiers_counts_matrix_filtered.tsv",
#'                                 package = "conclus")
#'
#' loadCountMatrix(file=countmatrixPath, header=TRUE, sep='\t', dec=".")
#'
#' @importFrom utils read.delim
#' 
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

