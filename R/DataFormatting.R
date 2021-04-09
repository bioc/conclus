#' conclusCacheClear
#'
#' @description
#' This function deletes the cache of conclus.
#'
#' @usage
#' conclusCacheClear()
#' 
#' @return Nothing, it deletes the cache of conclus.
#'
#' @aliases conclusCacheClear
#' @rdname conclusCacheClear
#' 
#' 
#' @details
#' This function don't return anything. It deletes the current contents
#' of the cache.
#' 
#' @examples 
#' NULL
#'
#' @import BiocFileCache
#' @importFrom  tools R_user_dir
#' 
#' @author Ilyess RACHEDI & Nicolas DESCOSTES
#'
#' @export conclusCacheClear
conclusCacheClear <- function() {
    cache <- Sys.getenv(x = "CONCLUS_CACHE", 
                        unset = tools::R_user_dir("conclus", which="cache"))
    bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
    removebfc(bfc, ask = FALSE)
}



#' .checkCache
#'
#' @description
#' Check before download if the data is one the cache. If it is, the cached 
#' file is copied the destination where it should have been downloaded.
#'
#' @usage
#' .checkCache(bfc, name, dataPath=NULL
#' 
#' @param bfc BiocFileCache object corresponding to the cache of CONCLUS
#' @param name ID (string) of the data extract from the url or series.
#' @param dataPath Path where the file must be downloaded
#' 
#' @return Path of the file in the cache if it exists, instead NULL.
#' @keywords internal
#' @noRd
#' @import BiocFileCache
.checkCache <- function(bfc, name, dataPath=NULL){

    rname <- bfcquery(bfc, name, field="rname")
    
    ## The cache exists and the path of the data is into it.
    if(isTRUE(as.logical(nrow(rname)))){
        cacheDataPath <- bfcrpath(bfc, rnames=name)
        
        if(file.exists(cacheDataPath)){
            if(!is.null(dataPath))
                ## The data is copied from the cache to the destination
                file.copy(cacheDataPath, dataPath, overwrite=TRUE) 
            return(cacheDataPath)
        }
    }
    
    return(NULL)
}



#' .retrieveMatrix
#'
#' @description
#' Download and read a count matrix from a URL
#'
#' @usage
#' .retrieveMatrix(matrixURL, countMatrixPath)
#'
#' @param matrixURL URL of the count matrix. The matrix must be un-normalized.
#' @param countMatrixPath Path to the file to which the downloaded count
#' matrix will be saved.
#' @param bfc BiocFileCache object corresponding to the cache of CONCLUS
#' 
#' @return The count matrix as a matrix.
#' @keywords internal
#' @noRd
#' @importFrom utils download.file
#' @import BiocFileCache
.retrieveMatrix <- function(bfc, matrixURL, countMatrixPath){

    name <- gsub("http.+file=", "", matrixURL)
    cacheDataPath <- .checkCache(bfc, name, countMatrixPath)
    
    if(is.null(cacheDataPath)){
        ## The count matrix isn't in the cache
        message("Downloading the count matrix.")
        download.file(matrixURL, countMatrixPath, timeout=120)
        add <- bfcadd(bfc, name, countMatrixPath, action="copy")
        
        countMatrix <- as.matrix(read.table(countMatrixPath, header=TRUE, 
                                            row.names=1, 
                                            stringsAsFactors=FALSE))
    }else 
        countMatrix <- as.matrix(read.table(cacheDataPath, header=TRUE, 
                                        row.names=1, stringsAsFactors=FALSE))

    return(countMatrix)

}


#' .retrieveColMetaDataFromURL
#'
#' @description
#' Download and read the columns meta-data from a URL
#'
#' @usage
#' .retrieveColMetaDataFromURL(bfc, colMetaDataURL, metaDataPath)
#' 
#' @param bfc BiocFileCache object corresponding to the cache of CONCLUS
#' @param colMetaDataURL URL of the meta-data. The columns should be cell ID,
#' cell barcode and state.
#' @param metaDataPath Path to the file to which the columns meta-data will be
#' saved.
#' 
#' @return The columns meta-data as a dataframe.
#' @keywords internal
#' @noRd
#' @importFrom utils download.file
#' @import BiocFileCache
.retrieveColMetaDataFromURL <- function(bfc, colMetaDataURL, metaDataPath){

    name <- gsub("http.+file=", "", colMetaDataURL)
    cacheDataPath <- .checkCache(bfc, name, metaDataPath)

    if(is.null(cacheDataPath)){
        message("Downloading the columns meta-data.")
        download.file(colMetaDataURL, metaDataPath, timeout=120)
        add <- bfcadd(bfc, name, metaDataPath, action="copy")
        metadata <- read.delim(metaDataPath, header=TRUE, 
                                stringsAsFactors=FALSE)

    }else 
        metadata <- read.table(cacheDataPath, header=TRUE, 
                                stringsAsFactors=FALSE, sep="\t")

    if(!isTRUE(all.equal(colnames(metadata), c("cellName", "state",
                        "cellBarcode"))))
        warning("The columns of the cells meta-data should be: cellName, ",
                "state, and cellBarcode. Please correct the dataframe.")   

    return(metadata)
}


#' .retrieveColMetaDataFromSeries
#'
#' @description
#' Download the columns meta-data from GEO.
#'
#' @usage
#' .retrieveColMetaDataFromSeries(seriesMatrixName)
#'
#' @param seriesMatrixName Name of the columns meta-data file hosted on GEO.
#' This name can usually be found in the 'Series Matrix File(s)' section.
#' @param bfc BiocFileCache object corresponding to the cache of CONCLUS
#'
#' @details
#' The method GEOquery::getGEO downloads all the series matrices of the GEO
#' record but only the one of interest is kept.
#'
#' @return A data.frame of the columns meta-data.
#' @keywords internal
#' @noRd
#'
#' @importFrom stringr str_extract
#' @importFrom GEOquery getGEO
#' @importFrom methods as
#' @importFrom  tools R_user_dir
#' @import BiocFileCache
.retrieveColMetaDataFromSeries <- function(bfc, seriesMatrixName){

    name <- tools::file_path_sans_ext(seriesMatrixName)
    metaDataCachePath <- .checkCache(bfc, name)
    
    if(is.null(metaDataCachePath)){
            
        message("Downloading the columns meta-data.")
        GEOnb <- stringr::str_extract(seriesMatrixName, "GSE[0-9]+")
        gpl <- GEOquery::getGEO(GEOnb, GSEMatrix=TRUE)
    
        if(!any(names(gpl) == seriesMatrixName))
            stop("The series matrix was not found. If the columns meta-data ",
                    "are provided as supplementary file, please use the ",
                    "colMetaDataURL parameter.")
    
        gpl <- gpl[[seriesMatrixName]]
        gpl <- methods::as(gpl, "data.frame")
        columnsMetaData <- data.frame(state=gpl$sampletype.ch1,
                cellBarcode=gpl$wellbarcode.ch1)
        
        ## Save on the cache
        savepath <- bfcnew(bfc, name, ext=".RData")
        save(columnsMetaData, file=savepath)
        
    }else 
            load(metaDataCachePath)
    
    return(columnsMetaData)
}

#' .filteringAndOrdering
#'
#' @description
#' Filter and order the count matrix and columns meta-data.
#'
#' @usage
#' .filteringAndOrdering(countMatrix, columnsMetaData)
#'
#' @param countMatrix Matrix containing single-cell gene expression retrieved
#' from GEO.
#' @param columnsMetaData Corresponding meta-data of the columns/cells of the
#' count matrix.
#'
#' @details
#' The columns of the matrix not having meta-data are removed. The meta-data
#' are then re-ordered accordingly to have a match between the matrix columns
#' and the rows of the meta-data data.frame. The cells names being "c[0-9]+"
#' are added to the matrix and meta-data.
#'
#' @return A list with the first element being the count matrix, and the second
#' the columns meta-data.
#' @keywords internal
#' @noRd
.filteringAndOrdering <- function(countMatrix, columnsMetaData){

    ## Removing columns of the matrix not having meta-data
    message("Formating data.")
    idx <- match(colnames(countMatrix), columnsMetaData$cellBarcode)
    if(isTRUE(all.equal(length(which(!is.na(idx))), 0)))
        warning("The cell barcodes were not found in the matrix. The columns ",
                "of the count matrix and the rows of the meta-data will ",
                "not be re-ordered. Are you sure that the count matrix and ",
                "the meta-data correspond?")
    else{

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
    }
    return(list(countMatrix, columnsMetaData))
}


#' .convertToSymbols
#'
#' @description
#' Convert the ENSEMBL IDs of the count matrix row names to symbols.
#'
#' @usage
#' .convertToSymbols(species, countMatrix, annoType)
#'
#' @param species Values should be 'mouse' or 'human'. Other organisms can be
#' added on demand.
#' @param countMatrix Matrix countaining single-cell gene expression retrieved
#' from GEO.
#' @param annoType Type of the genes annotations contained in the row names of
#' the count Matrix. Default: "ENSEMBL".
#'
#' @details
#' The conversion (TRUE by default) of the row genes IDs (ENSEMBL by default)
#' to official genes symbols is done with the function 'bitr' of the
#' 'clusterProfiler' package. To see a list of all possible values to pass to
#' the annoType parameter use 'keytypes' method on "org.Mm.eg.db" (for mouse)
#' or "org.Hs.eg.db" (for human). For example, copy/paste in a R terminal:
#' library(org.Mm.eg.db);keytypes(org.Mm.eg.db)
#'
#' @return The count matrix with row names converted to gene symbols.
#' @keywords internal
#' @noRd
#' @importFrom clusterProfiler bitr
.convertToSymbols <- function(species, countMatrix, annoType){

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

    return(countMatrix)
}



#' retrieveFromGEO
#'
#' @description
#' This function retrieves the count matrix and columns meta-data from GEO.
#' They are formatted to be suitable inputs for conclus.
#'
#' @usage
#' retrieveFromGEO(matrixURL, countMatrixPath, species,
#' seriesMatrixName=NA, metaDataPath=NA, colMetaDataURL=NA,
#' convertToSymbols=TRUE, annoType="ENSEMBL")
#'
#' @param matrixURL URL of the count matrix. The matrix must be un-normalized.
#' @param countMatrixPath Path to the file to which the downloaded count
#' matrix will be saved.
#' @param species Values should be 'mouse' or 'human'. Other organisms can be
#' added on demand.
#' @param seriesMatrixName Name of the columns meta-data file hosted on GEO.
#' This name can usually be found in the 'Series Matrix File(s)' section.
#' Should not be used if colMetaDataURL is defined. Default=NA.
#' @param metaDataPath If colMetaDataURL is used, defines the path to the file
#' to which the downloaded meta-data will be saved.
#' @param colMetaDataURL URL of the columns meta-data file hosted on GEO.
#' This file can be found in 'supplementary file'.
#' Should not be used if seriesMatrixName is defined. Default=NA.
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
#' result <- retrieveFromGEO(matrixURL, countMatrixPath, species,
#' seriesMatrixName=seriesMatrix)
#'
#' countMatrix <- result[[1]]
#' columnsMetaData <- result[[2]]
#' 
#' @author
#' Nicolas DESCOSTES & Ilyess RACHEDI
#' 
#' @importFrom  tools R_user_dir
#' @export retrieveFromGEO
retrieveFromGEO <- function(matrixURL, countMatrixPath, species,
        seriesMatrixName=NA, metaDataPath=NA, colMetaDataURL=NA,
        convertToSymbols=TRUE, annoType="ENSEMBL"){

    if(is.na(seriesMatrixName) && is.na(colMetaDataURL))
        stop("You should define at least seriesMatrixName or colMetaDataURL")

    cachePath <- tools::R_user_dir("conclus", which="cache")
    bfc <- BiocFileCache::BiocFileCache(cachePath, ask=FALSE)
    
    ## Retrieving the count matrix and columns meta-data
    countMatrix <- .retrieveMatrix(bfc, matrixURL, countMatrixPath)

    if(!is.na(seriesMatrixName))
        columnsMetaData <- .retrieveColMetaDataFromSeries(bfc, seriesMatrixName)
    else
        columnsMetaData <- .retrieveColMetaDataFromURL(bfc, colMetaDataURL,
                metaDataPath)

    ## Filtering, ordering and adding cells names
    result <- .filteringAndOrdering(countMatrix, columnsMetaData)
    countMatrix <- result[[1]]
    columnsMetaData <- result[[2]]


    ## Converting ensembl IDs to gene symbols
    if(convertToSymbols)
        countMatrix <- .convertToSymbols(species, countMatrix, annoType)

    idxDup <- duplicated(rownames(countMatrix))
    ldup <- length(which(idxDup))

    if(!isTRUE(all.equal(ldup, 0))){
        countMatrix <- countMatrix[!idxDup,]
        warning("Nb of lines removed due to duplication of row names: ", ldup)
    }

    return(list(countMatrix, columnsMetaData))
}
