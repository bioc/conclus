retrieveFromGEO <- function(matrixURL, countMatrixPath, seriesMatrixName,
        species, convertToSymbols=TRUE, annoType="ENSEMBL"){
    
    ## Retrieving the count matrix
    message("Downloading the count matrix.")
    download.file(matrixURL, countMatrixPath)
    countMatrix <- as.matrix(read.table(countMatrixPath, header=TRUE, 
                    row.names=1, stringsAsFactors = FALSE))
    
    ## Retrieving the columns meta-data
    message("Downloading the columns meta-data.")
    gpl <- GEOquery::getGEO(filename=seriesMatrixName)
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
    if(convertEnsemblToSymbols){
       
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
    
    
    return(list(countMatrix, columnsMetaData))
}

