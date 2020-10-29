#' retrieveFromGEO
#' 
#' @description
#' This function retrieves the count matrix and columns meta-data from GEO.
#' They are formatted to be suitable inputs for conclus.
#' 
#' @usage
#' retrieveFromGEO(matrixURL, countMatrixPath, seriesMatrixName,
#'        species, convertToSymbols=TRUE)
#' 
#' @param matrixURL URL of the count matrix. The matrix must be un-normalized.
#' @param countMatrixPath Path to the file to which the downloaded count 
#' matrix will be saved.
#' @param seriesMatrixName Name of the columns meta-data file hosted on GEO.
#' This name can usually be found in the 'Series Matrix File(s)' section.
#' @param species Values should be 'mouse' or 'human'.
#' @param convertToSymbols Boolean indicating if the genes IDs 
#' contained in the row names of the matrix should be converted to official 
#' genes symbols. Default: TRUE. To choose the type of IDs contained in the 
#' count matrix, see the annoType parameter just below.
#' @param annoType Type of the genes annotations contained in the row names of 
#' the count Matrix. Default: "ENSEMBL". See details.
#' 
#' @return
#' A list. The first element contains the count matrix and the second element 
#' contains the columns meta-data.
#' 
#' @aliases retrieveFromGEO
#' @rdname retrieveFromGEO
#' 
#' @details
#' The conversion (TRUE by default) of the row genes IDs (ENSEMBL by default) 
#' to official genes symbols is done with the function 'bitr' of the 
#' 'clusterProfiler' package. To see a list of all possible values to pass to 
#' the annoType parameter use 'keytypes' method on "org.Mm.eg.db" (for mouse) 
#' or "org.Hs.eg.db" (for human). For example, copy/paste in a R terminal: 
#' library(org.Mm.eg.db);keytypes(org.Mm.eg.db) 
#' 
#' @examples
#' outputDirectory <- "./YourOutputDirectory"
#' dir.create(outputDirectory, showWarnings=FALSE)
#' species <- "mouse"
#'
#' countMatrixPath <- file.path(outputDirectory, "countmatrix.txt")
#' matrixURL <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=",
#' "GSE96982&format=file&file=GSE96982%5FcountMatrix%2Etxt%2Egz") 
#' seriesMatrix <- "GSE96982-GPL19057_series_matrix.txt.gz"
#' 
#' result <- retrieveFromGEO(matrixURL, countMatrixPath, seriesMatrix, species)
#' countMatrix <- result[[1]]
#' columnsMetaData <- result[[2]]
#' 
#' @author
#' Nicolas DESCOSTES
#' 
#' @importFrom GEOquery getGEO
#' @importFrom clusterProfiler bitr
#' @importFrom stringr str_extract
#' @export retrieveFromGEO
retrieveFromGEO <- function(matrixURL, countMatrixPath, seriesMatrixName,
        species, convertToSymbols=TRUE, annoType="ENSEMBL"){
    
    ## Retrieving the count matrix
    message("Downloading the count matrix.")
    download.file(matrixURL, countMatrixPath)
    countMatrix <- as.matrix(read.table(countMatrixPath, header=TRUE, 
                    row.names=1, stringsAsFactors = FALSE))
    
    ## Retrieving the columns meta-data
    message("Downloading the columns meta-data.")
    GEOnb <- stringr::str_extract(seriesMatrixName, "GSE[0-9]+")
    gpl <- GEOquery::getGEO(GEOnb, GSEMatrix=TRUE)
    
    if(!any(names(gpl) == seriesMatrixName))
        stop("The series matrix was not found. Contact the developper with ",
                "a reproducible example.")
    
    gpl <- gpl[[seriesMatrixName]]
    gpl <- as(gpl, "data.frame")
    columnsMetaData <- data.frame(state=gpl$sampletype.ch1, 
            cellBarcode=gpl$wellbarcode.ch1)
    
    ## Removing columns of the matrix not having meta-data
    message("Formating data.")
    idx <- match(colnames(countMatrix), columnsMetaData$cellBarcode)
    if(isTRUE(all.equal(length(which(!is.na(idx))), 0)))
        stop("The cell barcodes were not found in the matrix. Are you sure ",
                "that the count matrix and the series matrix correspond?")
    
    idxNA <- which(is.na(idx))
    if(!isTRUE(all.equal(length(idxNA), 0))){
        countMatrix <- countMatrix[,-idxNA]
        idx <- idx[-idxNA]
    }
    
    ## Reordering columns meta-data
    columnsMetaData <- columnsMetaData[idx,]
    
    ## Adding cells names
    cellsnames <- paste0("c", seq_len(ncol(countMatrix)))
    colnames(countMatrix) <- cellsnames
    columnsMetaData <- cbind(cellName=cellsnames, columnsMetaData)
    rownames(columnsMetaData) <- cellsnames
    
    ## Converting ensembl IDs to gene symbols
    if(convertToSymbols){
       
        message("Converting ENSEMBL IDs to symbols.")
        if(isTRUE(all.equal(species, "mouse")))
            organismDatabase <- "org.Mm.eg.db"
        else if(isTRUE(all.equal(species, "human")))
            organismDatabase <- "org.Hs.eg.db"
        else
            stop("species should be mouse or human.")
        
        
        matrixSym <- rownames(countMatrix)
        symbolsVec <- clusterProfiler::bitr(matrixSym, fromType=annoType, 
                toType=c("SYMBOL"), OrgDb= organismDatabase)
        
        if(isTRUE(is(symbolsVec, "try-error")))
            stop("The matrix does not contains ENSEMBL IDs. Please verify ",
                    "the row names of the count matrix.")
        
        idxSymRemove <- which(duplicated(symbolsVec$ENSEMBL))
        
        if(!isTRUE(all.equal(length(idxSymRemove), 0)))
            symbolsVec <- symbolsVec[-idxSymRemove,]
        
        idxSym <- match(symbolsVec$ENSEMBL, matrixSym)
        rownames(countMatrix)[idxSym] <- symbolsVec$SYMBOL
    }
    
    idxDup <- duplicated(rownames(countMatrix))
    ldup <- length(which(idxDup))
    
    if(!isTRUE(all.equal(ldup, 0))){
        countMatrix <- countMatrix[!idxDup,]
        warning("Nb of lines removed due to duplication of row names: ", ldup)
    }
        
    
    return(list(countMatrix, columnsMetaData))
}
