#' .removeNoSymbol
#'
#' @description
#' Remove genes with no SYMBOL after the normalization.
#'
#' @param sceNorm SingleCellExperiment object
#'
#' @keywords internal
#' @importFrom SingleCellExperiment logcounts counts rowData colData
#' @importFrom rlang .data
#' @return SingleCellExperiment object
#' @noRd
.removeNoSymbol <- function(sceNorm){

    norm_mat <- SingleCellExperiment::logcounts(sceNorm)
    cm <- SingleCellExperiment::counts(sceNorm)
    rd <- SingleCellExperiment::rowData(sceNorm)
    cd <- SingleCellExperiment::colData(sceNorm)
    
    ## Remove SYMBOL genes on the count matrix
    genesToRemove <- rownames(subset(rd[, c("ENSEMBL", "SYMBOL")], 
                                is.na(rd$SYMBOL)))
    
    norm_mat <- norm_mat[!(row.names(norm_mat) %in% genesToRemove), ] 
    rowDataFiltered <- rd[!(row.names(rd) %in% genesToRemove), ] 
    countMatrixSCEFiltered <- cm[!(row.names(cm) %in% genesToRemove), ] 
    
    newSceNorm <- SingleCellExperiment(list(counts=countMatrixSCEFiltered,
                                            logcounts=norm_mat),
                                        rowData=rowDataFiltered,
                                        colData=cd)
    message("Genes with no SYMBOL IDs removed.")
    
    rownames(newSceNorm) <- rowData(newSceNorm)$SYMBOL
    return(newSceNorm)
}


#' .normalizeCon
#'
#' @description
#' This function is inspired from the function normalizeSCE of scater v1.12.2
#' and normalize the count matrix.
#'
#' @param object Result of scran::computeSumFactors filtering on cells having
#' a size factor higher than zero.
#' '
#' @keywords internal
#' @importFrom SummarizedExperiment assay
#' @return A SingleCellExperiment object containing the normalized matrix in
#' the logcounts slot.
#' @noRd
.normalizeCon <- function(object) {

    cur_exprs <- SummarizedExperiment::assay(object, i = "counts")
    col_list <-  split(cur_exprs, rep(seq_len(ncol(cur_exprs)),
                    each = nrow(cur_exprs)))
    norm_exprs <- mapply(function(exprsCol, sizeFact){

                result <- log2((exprsCol/sizeFact)+1)

            }, col_list, SingleCellExperiment::sizeFactors(object),
            SIMPLIFY=FALSE)
    norm_exprs <- do.call("cbind", norm_exprs)

    SummarizedExperiment::assay(object, "logcounts",
            withDimnames=FALSE) <- norm_exprs

    return(object)
}


#' .defineMartVar
#'
#' @description
#' This function retrieves the correct biomaRt parameters according to the
#' species.
#'
#' @param species The studied species.
#' '
#' @keywords internal
#' @import org.Mm.eg.db org.Hs.eg.db
#' @importFrom biomaRt useMart
#' @return Returns the genome annotation name, the pattern to use to retrieve
#' ensembl IDs and the biomaRt database.
#' @noRd
.defineMartVar <- function(species){

    if(isTRUE(all.equal(species, "human"))){

        # suppressMessages(library(org.Hs.eg.db, warn.conflicts=FALSE))
        genomeAnnot <- org.Hs.eg.db::org.Hs.eg.db
        ensemblPattern <- "ENSG"
        dataset <- "hsapiens_gene_ensembl"

    }else if(isTRUE(all.equal(species, "mouse"))){

        # suppressMessages(library(org.Mm.eg.db, warn.conflicts=FALSE))
        genomeAnnot <- org.Mm.eg.db::org.Mm.eg.db
        ensemblPattern <- "ENSMUSG"
        dataset <-"mmusculus_gene_ensembl"

    }else
        stop("Species should be human or mouse: ", species, "is not ",
                "currently supported.")

    ensembl <- .tryUseMart(biomart="ensembl", dataset)

    return(list(genomeAnnot, ensemblPattern, ensembl))
}


#' .testMattrixNames
#'
#' @description
#' This function tests that the count matrix has at least ensembl IDs or gene
#' symbols.
#'
#' @param countMatrix The raw count matrix
#' @param ensemblPattern Expression used to retrieve the ensembl IDs.
#' @param genomeAnnot Bioconductor package name to be used with biomaRt.
#'
#' @keywords internal
#' @importFrom AnnotationDbi keys
#' @return Returns the number of lines of count matrix having a gene symbol.
#' @noRd
.testMattrixNames <- function(countMatrix, ensemblPattern, genomeAnnot){

    allnames <- rownames(countMatrix)
    lengthEnsemblGenes <- length(grep(ensemblPattern, allnames))
    lengthSymbols <- length(intersect(AnnotationDbi::keys(genomeAnnot,
                            keytype = "SYMBOL"), allnames))

    if(isTRUE(all.equal(lengthEnsemblGenes, 0)) &&
            isTRUE(all.equal(lengthSymbols, 0)))
        stop("The row names of the matrix neither contain ensembl genes or",
                "official gene symbols.")

    return(lengthSymbols)
}


#' .annotateEnsembl
#'
#' @description
#' This function selects the symbol and the gene name from the ensembl IDs.
#'
#' @param ensemblGenes Row names of the counts matrix having an ensembl ID.
#' @param ensemblPattern Expression used to retrieve the ensembl IDs.
#' @param genomeAnnot Bioconductor package name to be used with biomaRt.
#'
#' @keywords internal
#' @importFrom AnnotationDbi keys select
#' @return Returns the rowData with ensembl IDs completed with symbols and gene
#' name.
#' @noRd
.annotateEnsembl <- function(ensemblGenes, ensemblPattern, genomeAnnot){

    message("Annotating ",length(ensemblGenes), " genes containing ",
            ensemblPattern, " pattern.")

    ## If none of ENSEMBL genes are in database, set to NA
    lengthEnsInDb <- length(intersect(AnnotationDbi::keys(genomeAnnot,
                            keytype="ENSEMBL"), ensemblGenes))

    if(!isTRUE(all.equal(lengthEnsInDb, 0))){

        rowdataEnsembl <- AnnotationDbi::select(genomeAnnot,
                keys = ensemblGenes, keytype = "ENSEMBL",
                columns = c("SYMBOL", "GENENAME"), multiVals="first")
        rowdataEnsembl <-
                rowdataEnsembl[!duplicated(rowdataEnsembl$ENSEMBL),]
    }else
        rowdataEnsembl <- data.frame(ENSEMBL=ensemblGenes, SYMBOL = NA,
                GENENAME=NA)

    rowdataEnsembl$nameInCountMatrix <- ensemblGenes

    return(rowdataEnsembl)
}


#' .annotateSymbols
#'
#' @description
#' This function selects the ensembl ID and the gene name from the symbol.
#'
#' @param genomeAnnot Bioconductor package name to be used with biomaRt.
#' @param symbolGenes Names of the rows of the count matrix having a gene symbol
#' as name instead of an ensembl ID.
#'
#' @keywords internal
#' @importFrom AnnotationDbi select
#' @return Returns the rowData with symbols completed with ensembl IDs and gene
#' name.
#' @noRd
.annotateSymbols <- function(symbolGenes, genomeAnnot){

    message("Annotating ", length(symbolGenes),
            " genes considering them as SYMBOLs.")

    rowdataSymbol <- AnnotationDbi::select(genomeAnnot,
            keys=symbolGenes,
            keytype="SYMBOL",
            columns=c("ENSEMBL",
                    "GENENAME"),
            multiVals="first")
    rowdataSymbol <- rowdataSymbol[!duplicated(rowdataSymbol$SYMBOL),]
    rowdataSymbol$nameInCountMatrix <- symbolGenes

    return(rowdataSymbol)
}


#' .annotateRowData
#'
#' @description
#' This function calls the annotation methods on ensembl IDs and gene symbols.
#'
#' @param ensemblGenes Row names of the counts matrix having an ensembl ID.
#' @param ensemblPattern Expression used to retrieve the ensembl IDs.
#' @param genomeAnnot Bioconductor package name to be used with biomaRt.
#' @param lengthSymbols Number of rows in the count matrix having a gene symbol
#' instead of an ensembl ID.
#' @param symbolGenes Names of the rows of the count matrix having a gene symbol
#' as name instead of an ensembl ID.
#'
#' @keywords internal
#' @return Returns the rowData annotated with the sub-functions considering
#' ensembl IDs or symbols.
#' @noRd
.annotateRowData <- function(ensemblGenes, ensemblPattern, genomeAnnot,
        lengthSymbols, symbolGenes){

    if(!isTRUE(all.equal(length(ensemblGenes), 0)))
        rowdataEnsembl <- .annotateEnsembl(ensemblGenes, ensemblPattern,
                genomeAnnot)
    else
        rowdataEnsembl <- NULL

    if(!isTRUE(all.equal(lengthSymbols, 0)))
        rowdataSymbol <- .annotateSymbols(symbolGenes, genomeAnnot)
    else
        rowdataSymbol <- NULL

    rowdata <- base::rbind(rowdataSymbol, rowdataEnsembl)
    return(rowdata)
}


#' .retrieveGenesInfoBiomart
#'
#' @description
#' This function retrieves the ensembl_gene_id, go_id, name_1006,
#' chromosome_name, and gene_biotype for each gene retrieved from the count
#' matrix.
#'
#' @param ensembl Value returned by the method UseMart called in the function
#' .defineMartVar.
#' @param rowdata rowData of class data.frame, it contains gene names of the
#' count matrix.
#'
#' @keywords internal
#' @importFrom biomaRt getBM
#' @return Returns the rowData filled with the bioMart annotations.
#' @noRd
.retrieveGenesInfoBiomart <- function(ensembl, rowdata){

    attributes <- c("ensembl_gene_id", "go_id", "name_1006", "chromosome_name",
                    "gene_biotype")

    genes <- rowdata$ENSEMBL
    res <- .tryGetBM(attributes, ensembl, values=genes, "ensembl_gene_id")

    tmp <- res[!duplicated(res$ensembl_gene_id),]

    rowdata <- merge(rowdata, tmp[c("ensembl_gene_id", "chromosome_name",
                            "gene_biotype")], by.x = "ENSEMBL",
            by.y = "ensembl_gene_id", all.x = TRUE, all.y = FALSE, sort = FALSE)

    rowdataGO <- merge(rowdata, res[c("ensembl_gene_id",
                            "go_id", "name_1006")], by.x = "ENSEMBL",
            by.y = "ensembl_gene_id", all.x = TRUE, all.y = FALSE, sort = FALSE)
    rowdataGO <- rowdataGO[!is.na(rowdataGO$name_1006) &
                    ((rowdataGO$name_1006 == "cell surface") |
                        (rowdataGO$name_1006 == strwrap("cell surface
                                            receptor signaling
                                            pathway"))), ]
    rowdataGO$name_1006[duplicated(rowdataGO$ENSEMBL)] <-
            "cell surface receptor signaling pathway"
    rowdataGO <- rowdataGO[!duplicated(rowdataGO$ENSEMBL),]
    rowdata <- merge(rowdata, rowdataGO[c("nameInCountMatrix", "go_id",
                            "name_1006")], by.x = "nameInCountMatrix",
            by.y = "nameInCountMatrix", all.x = TRUE, all.y = TRUE,
            sort = FALSE)
    return(rowdata)
}


#' .mergeRowDataDf
#'
#' @description
#' This function merges the rowdataDF to the rowData if not null.
#'
#' @param rowdataDF rowData of class data.frame, it contains gene names of the
#' count matrix.
#' @param Data frame containing annotations retrieved from biomaRt.
#'
#' @keywords internal
#'
#' @return Returns the rowData filled with the rowdataDF annotations.
#' @noRd
.mergeRowDataDf <- function(rowdataDF, rowdata){

    if("nameInCountMatrix" %in% colnames(rowdata) &&
            "nameInCountMatrix" %in% colnames(rowdataDF)){
        rowdataMerged <- merge(rowdataDF, rowdata, all.x = TRUE, all.y = TRUE,
                sort = FALSE)

        return(rowdataMerged)

        }else
            stop("No 'nameInCountMatrix' column in the rowdata or rowdataDF",
                    " to perform the merge." )
}


#' .mergeColDataDf
#'
#' @description
#' This function merge the coldata created by CONCLUS with  the coldata
#' provided by the user if he has entered it.
#'
#' @param coldataDF Data frame provided by the user containing information
#' about cells.
#' @param coldata Data frame created by CONCLUS containing information
#' about cells.
#'
#' @keywords internal
#'
#' @return Returns the final coldata
#' @noRd
.mergeColDataDf <- function(coldataDF, coldata){

    if("cellName" %in% colnames(coldata) &&
        "cellName" %in% colnames(coldataDF)){
            coldataMerged <- merge(coldataDF, coldata, all.x = TRUE,
                                    all.y = TRUE, sort = FALSE)

        return(coldataMerged)

    }else
        stop("No 'cellName' column in the coldata or coldataDF",
                " to perform the merge." )

}

#' .annotateGenes
#'
#' @description
#' This function returns a rowData with the same rownames as in countMatrix
#' with annotations. Only genes having an ENSEMBL IDs or SYMBOL ID
#' will receive the annotation. This function use \code{\link{bioMart}} and
#' \code{\link{AnnotationDbi}} to retrieve annotations.
#'
#' @param countMatrix The raw count matrix
#' @param species The studied species. Currently limited to mouse and human.
#' Other organisms can be added on demand.
#' @param rowdataDF rowData of class data.frame, it contains gene names of the
#' count matrix.
#' @param info Logical. If TRUE, additional annotations like ensembl_gene_id,
#' go_id, name_1006, chromosome_name and gene_biotype are added to the 
#' row data, for all the genes from the count matrix with ENSEMBL IDs or 
#' SYMBOL ID. Default: TRUE.
#' 
#' @keywords internal
#'
#' @return Returns the rowData filled with annotations
#' @noRd
.annotateGenes <- function(countMatrix, species, rowdataDF, info=TRUE){

    martResult <- .defineMartVar(species)
    genomeAnnot <- martResult[[1]]
    ensemblPattern <- martResult[[2]]
    ensembl <- martResult[[3]]
    lengthSymbols <- .testMattrixNames(countMatrix, ensemblPattern, genomeAnnot)

    ensemblGenes <- rownames(countMatrix)[grep(ensemblPattern,
                    rownames(countMatrix))]
    symbolGenes <- rownames(countMatrix)[!grepl(ensemblPattern,
                    rownames(countMatrix))]

    rowdata <- .annotateRowData(ensemblGenes, ensemblPattern, genomeAnnot,
            lengthSymbols, symbolGenes)

    ## Filtering duplicated symbols
    (multSym <- rowdata$SYMBOL[!is.na(rowdata$SYMBOL) &
                                duplicated(rowdata$SYMBOL)])
    ## We don't combine such ensembls, but leave them unique with
    ## "symbol_ensembl"
    (rowdata$SYMBOL[rowdata$SYMBOL %in% multSym] <-
                paste0(rowdata$SYMBOL[rowdata$SYMBOL %in% multSym],
                        rowdata$ENSEMBL[rowdata$SYMBOL %in% multSym]))
    
    if(isTRUE(info))
        rowdata <- .retrieveGenesInfoBiomart(ensembl, rowdata)

    if(!is.null(rowdataDF))
        rowdata <- .mergeRowDataDf(rowdataDF, rowdata)

    rownames(rowdata) <- rowdata$nameInCountMatrix
    rowdata <- rowdata[rownames(countMatrix),]
    rowdata$SYMBOL[(S4Vectors::isEmpty(rowdata$SYMBOL)) |
                        (rowdata$SYMBOL == "")] <- NA
    stopifnot(all(rownames(rowdata) == rownames(countMatrix)))

    return(rowdata)
}



#' .createReportTable
#'
#' @description
#' This function creates a boolean vector that is used to filter the cells
#' (see .filterCells below).
#'
#' @details
#' Rather 'columnName' is a positive or a negative feature,
#' the function applies a different threshold on the data.
#' This threshold is used to create a boolean vector that will be used to
#' filter the cells. This function is used by .filterCells to create the report
#' table.
#'
#' @param columnName Name of the column to calculate quantiles
#' @param colData Data frame containing informations about cells
#' @param mb Vector of positive features reflecting the quality of cells
#' @param mw Vector of negative features reflecting the lack of quality of cell
#'
#' @importFrom methods is
#' @importFrom stats quantile
#' @import doParallel
#' @keywords internal
#'
#' @return A boolean vector
#' @noRd
.createReportTable <- function(columnName, colData, mb, mw){

    quan <- quantile(colData[, colnames(colData) == columnName])
    
    if (columnName %in% mb) {
        threshold <- 2.5 * quan[2] - 1.5 * quan[4]
        if (threshold < 0){
            threshold <- (quan[1] + quan[2]) / 2
            }
        vec <- as.numeric(colData[ ,columnName] >=  as.numeric(threshold))
        
    } else if (columnName %in% mw) {
        threshold <- 2.5 * quan[4] - 1.5 * quan[2]
        if (threshold > quan[5]){
            threshold <- (quan[3] + quan[4]) / 2
        }
        vec <- as.numeric(colData[, columnName] <=  as.numeric(threshold))
    }

    return(vec)
}



#' .filterCells
#'
#' @description
#' This function removes not significant cells from the analysis. It filters
#' the cells by creating a report table that is used to define which cell is
#' informative or not.
#'
#' @param countMatrix Matrix containing the transcripts counts.
#' @param colData Data frame containing information about cells
#' @param genesSumThr Threshold of the minimum nb of s that should be
#' significant
#' @param MoreBetter Vector with the name of positive features in the colData
#' reflecting the quality of cells. Default: c("genesNum", "sumCodPer",
#' "genesSum")
#' @param MoreWorse Vector with the name of negative features in the colData
#' reflecting the quality of cells. Default: "sumMtPer"
#' @importFrom dplyr mutate
#' @keywords internal
#'
#' @return Returns a rlist containing the count matrix and the colData filtered
#' @noRd
.filterCells <- function(countMatrix, colData, genesSumThr=100,
                        MoreBetter=c("genesNum", "sumCodPer", "genesSum"),
                        MoreWorse=c("sumMtPer")){
    
    message("Running filterCells.")
    countMatrix <- countMatrix[, colSums(countMatrix) >= genesSumThr]
    if (isTRUE(all.equal(ncol(countMatrix), 0)))
        stop("None of your cells has at least 100 genes expressed. Since the ",
            "filtering keeps only those cells, nothing will be kept. ",
            "Please check the count matrix.")

    colData <- colData[colnames(countMatrix), ]
    mb <- MoreBetter
    mw <- MoreWorse
    columnNames <- c(mb, mw)

    ## Create the report table
    reportTable <- unlist(lapply(columnNames, .createReportTable, colData,
                                    mb, mw))
    reportTable <- matrix(reportTable, ncol=length(columnNames),
                        dimnames=list(colData$cellName, columnNames))

    ## Add cell names column to the report table
    reportTable <- data.frame(cellName=colData$cellName, reportTable)

    ### add columns with filtering score and verdict ###
    reportTable <- dplyr::mutate(reportTable, score=NA)
    reportTable$score <- rowSums(reportTable[, colnames(reportTable) %in%
                                                c(mb,mw)])
    reportTable <- dplyr::mutate(reportTable, filterPassed=NA)
    reportTable$filterPassed[reportTable$score >= length(mb) + length(mw)] <- 1
    reportTable$filterPassed[reportTable$score <  length(mb) + length(mw)] <- 0

    ### add filtering verdict to colData ###
    colData <- dplyr::mutate(colData, filterPassed=NA)
    colData$filterPassed[colData$cellName %in%
                    reportTable$cellName[reportTable$filterPassed == 1]] <- 1
    colData$filterPassed[colData$cellName %in%
                    reportTable$cellName[reportTable$filterPassed == 0]] <- 0

    reportTable <- reportTable[order(reportTable$score, decreasing = FALSE), ]

    rownames(colData) <- colData$cellName
    colData <- colData[colnames(countMatrix), ]
    stopifnot(all(rownames(colData) == colnames(countMatrix)))
    countMatrix <- countMatrix[, colData$filterPassed == 1]
    colData <- colData[colData$filterPassed == 1, ]

    return(list(countMatrix, colData))
}


#' .fillGenesNumColumn
#'
#' @description
#' This function counts the number of genes having read counts > 1 in a cell
#'
#' @param countMatrix The count matrix.
#' @param cell Name of the cell considered in the sapply (see.
#' @param coldata Data frame containing some informations about cells.
#'
#' @keywords internal
#' @noRd
.fillGenesNumColumn <- function(cell, coldata, countMatrix){

    readsCell <- countMatrix[, cell]
    coldata$genesNum[coldata$cellName == cell] <- length(
            readsCell[readsCell > 0])
}



#' .fillOneUmmiColumn
#'
#' @description
#' This function fills the OneUMiColumn of the report table for one cell.
#'
#' @param countMatrix The count matrix.
#' @param cell Cell to count
#' @param coldata  Data frame containing informations about cells
#'
#' @keywords internal
#' @noRd
.fillOneUmmiColumn <- function(cell, coldata, countMatrix){

    readsCell <- countMatrix[, cell]
    coldata$oneUMI[coldata$cellName == cell] <- length(
            readsCell[readsCell == 1])
}

#' .addGenesInfoCell
#'
#' @description
#' This function fills the genesNum, genesSum, and oneUMI that represents
#' columns meta-data. This information will be used to filter out cells.
#'
#' @param coldata  Data frame containing informations about cells
#' @param countMatrix The count matrix.
#'
#' @importFrom dplyr mutate
#' @keywords internal
#' @noRd
#' @return The updated columns meta-data.
.addGenesInfoCell <- function(coldata, countMatrix){

    coldata <- dplyr::mutate(coldata, genesNum=NA, genesSum=NA, oneUMI=NA)
    coldata$genesSum <- colSums(countMatrix)
    coldata$genesNum <- vapply(colnames(countMatrix), .fillGenesNumColumn,
            coldata, countMatrix, FUN.VALUE=integer(1))
    coldata$oneUMI <- vapply(colnames(countMatrix), .fillOneUmmiColumn, coldata,
            countMatrix, FUN.VALUE=integer(1))
    coldata <- dplyr::mutate(coldata,
            oneUMIper =100 * coldata$oneUMI / coldata$genesNum)
    return(coldata)
}

#' .addCellsInfo
#'
#' @description
#' This function creates colData or add columns mtGenes, genesNum, codGenes,
#' genesSum, codSum, mtPer, codPer, sumMtPer, sumCodPer to the existing colData.
#'
#' @param countMatrix The count matrix.
#' @param rowdataDF  Data frame containing informations about genes.
#' @param coldataDF  Optionnal. Data frame containing information about cells.
#' @importFrom  dplyr mutate
#' @keywords internal
#'
#' @return Returns the colData
#' @noRd
.addCellsInfo <- function(countMatrix, rowdataDF, coldataDF = NULL){

    message("Adding cell info for cells filtering.")
    coldata <- data.frame(cellName = colnames(countMatrix),
            stringsAsFactors = FALSE)

    ### add info about all genes in a cell
    coldata <- .addGenesInfoCell(coldata, countMatrix)
    ### add info about mitochondrial and protein-coding genes
    coldata <- dplyr::mutate(coldata, mtGenes=NA, mtSum=NA, codGenes=NA, 
                                codSum=NA)

    for(i in seq_len(ncol(countMatrix))){

        mt <- countMatrix[rownames(countMatrix) %in%
                rowdataDF$nameInCountMatrix[
                        rowdataDF$chromosome_name == "MT"], i]

        coldata$mtGenes[coldata$cellName == colnames(countMatrix)[i]] <-
            length(mt[mt > 0])
        coldata$mtSum[coldata$cellName == colnames(countMatrix)[i]] <- sum(mt)
        cod <- countMatrix[rownames(countMatrix) %in%
                        rowdataDF$nameInCountMatrix[rowdataDF$gene_biotype ==
                                                        "protein_coding"], i]
        coldata$codGenes[coldata$cellName == colnames(countMatrix)[i]] <-
            length(cod[cod > 0])
        coldata$codSum[coldata$cellName == colnames(countMatrix)[i]] <- sum(cod)
    }

    coldata <- dplyr::mutate(coldata,
            mtPer=100 * coldata$mtGenes / coldata$genesNum,
            codPer=100 * coldata$codGenes / coldata$genesNum,
            sumMtPer=100 * coldata$mtSum / coldata$genesSum,
            sumCodPer=100 * coldata$codSum / coldata$genesSum)

    if (!is.null(coldataDF))
        coldata <- .mergeColDataDf(coldataDF, coldata)

    rownames(coldata) <- coldata$cellName
    coldata <- coldata[colnames(countMatrix), ]
    stopifnot(all(rownames(coldata) == colnames(countMatrix)))

    return(coldata)
}



#' .filterGenes
#'
#' @description
#' This function filters the genes of the rowData and the count matrix.
#' If a gene has less than ten counts among all cells, it is removed.
#'
#' @param rowData Data frame containing informations about genes.
#' @param countMatrix Class matrix. It's the count matrix.
#'
#' @keywords internal
#'
#' @return Returns a rlist containing the count matrix and the colData
#' @noRd
.filterGenes <- function(countMatrix, rowData){

    ## internal function, filters genes which are more than
    ## in 10 cells and less than (all-10) cells
    message("Running filterGenes.")
    selRows <- ((rowSums(countMatrix[, ] >= 1)) > 10)
    countMatrix <- countMatrix[selRows, ]
    rowData <- rowData[rowData$nameInCountMatrix %in% rownames(countMatrix), ]

    return(list(countMatrix, rowData))
}


#' .checkParamNorm
#'
#' @description
#' This function checks the parameters of the method normaliseCountMatrix.
#'
#' @param sizes  Vector of size factors from scran::computeSumFactors()
#' function.
#' @param rowdata Data frame containing genes informations. Default is NULL.
#' @param coldata Data frame containing cells informations. Default is NULL.
#' @param alreadyCellFiltered Logical. If TRUE, quality check and
#' filtering will not be applied.
#' @param runQuickCluster Logical. If TRUE scran::quickCluster() function will
#' be applied. It usually improves the normalization for medium-size count
#' matrices. However, it is not recommended for datasets with less than 200
#' cells and may take too long for datasets with more than 10000 cells.
#' @param info Logical. If TRUE, additional annotations like ensembl_gene_id,
#' go_id, name_1006, chromosome_name and gene_biotype are added to the 
#' row data, for all the genes from the count matrix with ENSEMBL IDs or 
#' SYMBOL ID. Default: TRUE.
#' @keywords internal
#' @noRd
.checkParamNorm <- function(sizes, rowdata, coldata, alreadyCellFiltered,
        runQuickCluster, info){

    if(!is.numeric(sizes))
        stop("'sizes' parameter should be a vector of numeric values.")

    if(!is.data.frame(rowdata) && !is.null(rowdata))
        stop("'rowdata' parameter should be a data frame or be NULL.")

    if(!is.data.frame(coldata) && !is.null(coldata))
        stop("'coldata' parameter should be a data frame or be NULL.")

    if (!is.logical(alreadyCellFiltered))
        stop("'alreadyCellFiltered' parameter should be a boolean.")

    if (!is.logical(runQuickCluster))
        stop("'runQuickCluster' parameter should be a boolean.")
    
    if (!is.logical(info))
        stop("'info' parameter should be a boolean.")
}


#' .filterSCE
#'
#' @description
#' This function calls sub-functions that filters non informative cells and
#' genes.
#'
#' @param alreadyCellFiltered Logical. If TRUE, quality check and
#' filtering will not be applied.
#' @param countMatrix Class matrix representing the genes expression.
#' @param rowdata Data frame containing genes informations.
#' @param coldata Data frame containing cells informations.
#'
#' @keywords internal
#' @noRd
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @return A single cell experiment object
.filterSCE <- function(alreadyCellFiltered, countMatrix, coldata, rowdata){


    if(!alreadyCellFiltered){

        nbCellsBefore <- nrow(coldata)
        filterCellsResult <- .filterCells(countMatrix, coldata)
        countMatrix <- filterCellsResult[[1]]
        coldata <- filterCellsResult[[2]]
        nbCellsAfter <- nrow(coldata)
        message("Number of cells removed after filtering: ",
                nbCellsBefore - nbCellsAfter)
        message("Number of cells remaining after filtering: ", nbCellsAfter)
    }

    if(isTRUE(nrow(coldata) == 0))
        stop("There are no more cells after filtering. Maybe the ",
            "count matrix does not contain enough informative cells. ",
            "Please check the count matrix.")

    nbGenesBefore <- nrow(rowdata)
    filterGenesResult <- .filterGenes(countMatrix, rowdata)
    countMatrix <- filterGenesResult[[1]]
    rowdata <- filterGenesResult[[2]]
    nbGenesAfter <- nrow(rowdata)

    message("Number of genes removed after filtering: ",
            nbGenesBefore - nbGenesAfter)
    message("Number of genes remaining after filtering: ", nbGenesAfter)

    if(isTRUE(nrow(rowdata) == 0))
        stop("There are no more genes after filtering. Maybe the count matrix ",
            "contains only genes which are less than in 10 cells or more than ",
            "all-10 cells. Please check the count matrix.")

    stopifnot(all(rownames(countMatrix) == rownames(rowdata)))
    stopifnot(all(colnames(countMatrix) == rownames(coldata)))
    sce <- SingleCellExperiment::SingleCellExperiment(assays=list(
                    counts=as.matrix(countMatrix)),
            colData=coldata,
            rowData=rowdata)

    return(sce)

}


#' .checkRowAndColdata
#'
#' @description
#' This function verifies that the meta-data have correct dimensions.
#'
#' @param countMatrix Class matrix representing the genes expression.
#' @param rowdata Data frame containing genes informations. Default is NULL.
#' @param coldata Data frame containing cells informations. Default is NULL.
#'
#' @keywords internal
#' @noRd
.checkRowAndColdata <- function(countMatrix, rowdata, coldata){

    if(!is.null(rowdata) && !isTRUE(all.equal(nrow(rowdata),
                                                nrow(countMatrix))))
        stop("The provided row metadata should contain the same number ",
                "of rows than the matrix.")

    if(!is.null(coldata) && !isTRUE(all.equal(nrow(coldata),
                                            ncol(countMatrix))))
        stop("The provided col metadata should contain the same number ",
                "of rows than the matrix number of columns.")
}


#' .quickClusterScran
#'
#' @description
#' This function uses quickCluster of scran package to assign clusters to
#' cells. When the data is too small, the PC of internal function of
#' quickCluster can be too big. In this case, the function is restarted
#' with 'd=NA' paramater. In this way, no dimensionality reduction is
#' performed and the gene expression values are directly used in clustering.
#'
#' @param sce Single Cell experiment object
#'
#' @keywords internal
#' @noRd
.quickClusterScran <- function(sce){

    cl <- tryCatch(
        scran::quickCluster(sce),

        warning = function(w){

                    m <- paste("You're computing too large a percentage of",
                        "total singular values, use a standard svd instead.")

                    if(grepl(m, w$message) && ncol(sce) < 200 )
                        scran::quickCluster(sce, d=NA)
                },

        error=function(e) NULL)

    return(cl)
}

#' normaliseCountMatrix
#'
#' @description
#' This function uses coldata (cells informations) and rowdata (genes
#' informations) to filter the count matrix. It also normalizes by using
#' deconvolution with size factors.
#'
#' @usage
#' normaliseCountMatrix(theObject, sizes=c(20,40,60,80,100), rowdata=NULL,
#'                     coldata=NULL, alreadyCellFiltered=FALSE,
#'                     runQuickCluster=TRUE, info=TRUE)
#'
#'
#' @param theObject A scRNAseq object
#' @param sizes  Vector of size factors used by scran::computeSumFactors().
#' It is a numeric vector of pool sizes, i.e., number of cells per pool. See
#' ?scran::computeSumFactors for more details.
#' @param coldata Data frame containing cells informations. Default is NULL.
#' @param rowdata Data frame containing genes informations. Default is NULL.
#' @param alreadyCellFiltered Logical. If TRUE, quality check and
#' filtering will not be applied.
#' @param runQuickCluster Logical. If TRUE scran::quickCluster() function will
#' be applied. It usually improves the normalization for medium-size count
#' matrices. However, it is not recommended for datasets with less than 200
#' cells and may take too long for datasets with more than 10000 cells.
#' @param info Logical. If TRUE, additional annotations like ensembl_gene_id,
#' go_id, name_1006, chromosome_name and gene_biotype are added to the 
#' row data, for all the genes from the count matrix with ENSEMBL IDs or 
#' SYMBOL ID. Default: TRUE.
#' @param removeNoSymbol Logical. If TRUE, genes with no SYMBOL are removed
#' after the normalization
#' 
#' @details
#' This function uses the normalization method of the scater package. For more
#' details about the normalization used see ?scater::normalize. The size factors
#' used in the normalization are calculated with scran::computeSumFactors.
#'
#' Beforehand, the function will annotate genes creating rowData and add
#' statistics about cells into columnsMetaData. If you already have
#' columnsMetaData and rowData, you can give it to the function (see manual).
#' It will keep your columns and add new ones at the end. If you do not want
#' to lose any cell after quality metrics check, select
#' alreadyCellFiltered = TRUE, by default it is FALSE. Before scater
#' normalization, the function will call scran::quickCluster (see manual for
#' details). If you want to skip this step, set runQuickCluster = FALSE, by
#' default it is TRUE. We advice to first try the analysis with this option
#' and to set it to FALSE if no rare populations are found.
#'
#'
#' @rdname normaliseCountMatrix-scRNAseq
#' @aliases normaliseCountMatrix
#'
#' @importFrom scran quickCluster computeSumFactors
#' @importFrom scater logNormCounts
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SingleCellExperiment sizeFactors
#'
#' @return Returns a scRNASeq object with its sceNorm slot updated. This slot
#' contains a SingleCellExperiment object having the normalized count matrix,
#' the colData (table with cells informations), and the rowData (table with the
#' genes informations). See ?SingleCellExperiment for more details.
#'
#' @examples
#' ## Load the count matrix
#' countmatrixPath <- system.file("extdata/countMatrix.tsv", package="conclus")
#' countMatrix <- loadDataOrMatrix(file=countmatrixPath, type="countMatrix",
#'                                 ignoreCellNumber=TRUE)
#'
#' ## Load the coldata
#' coldataPath <- system.file("extdata/colData.tsv", package="conclus")
#' columnsMetaData <- loadDataOrMatrix(file=coldataPath, type="coldata",
#' columnID="cell_ID")
#'
#' ## Create the initial object
#' scr <- singlecellRNAseq(experimentName = "Bergiers",
#'                 countMatrix     = countMatrix,
#'                 species         = "mouse",
#'                 outputDirectory = "YourOutputDirectory")
#'
#' ## Normalize and filter the raw counts matrix
#' scr <- normaliseCountMatrix(scr, coldata = columnsMetaData, info=FALSE)
#'
#' @exportMethod normaliseCountMatrix
#' @importFrom methods validObject
#'
#' @author Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas
#' DESCOSTES.

setMethod(

    f = "normaliseCountMatrix",

    signature = "scRNAseq",

    definition = function(theObject, sizes=c(20,40,60,80,100), rowdata=NULL,
                            coldata=NULL, alreadyCellFiltered=FALSE,
                            runQuickCluster=TRUE, info=TRUE, 
                            removeNoSymbol=FALSE){
        
        validObject(theObject)
        .checkParamNorm(sizes, rowdata, coldata, alreadyCellFiltered,
                runQuickCluster, info)
        countMatrix <- getCountMatrix(theObject)
        species <- getSpecies(theObject)

        .checkRowAndColdata(countMatrix, rowdata, coldata)

        rowdata <- .annotateGenes(countMatrix, species=species,
                    rowdataDF=rowdata, info=info)

        coldata <- .addCellsInfo(countMatrix, rowdataDF=rowdata, 
                                coldataDF=coldata)
        sce <- .filterSCE(alreadyCellFiltered, countMatrix, coldata, rowdata)

        # normalization
        message("Running normalization. It can take a while depending on the",
                " number of cells.")
        
        cl <- NULL
        if (runQuickCluster) cl <- .quickClusterScran(sce)

        # Compute sizeFactors which will be used for normalization
        sceNorm <- suppressMessages(
            scran::computeSumFactors(sce, sizes=sizes, clusters=cl, 
                                    positive=FALSE))
        message("summary(sizeFactors(sceObject)):")
        print(summary(SingleCellExperiment::sizeFactors(sceNorm)))

        nbrNegSFCells <- length(SingleCellExperiment::sizeFactors(sceNorm)[
                        SingleCellExperiment::sizeFactors(sceNorm) <= 0])

        if (nbrNegSFCells > 0){
            message(nbrNegSFCells, " cells with negative sizeFactors will be ",
                    "deleted before the downstream analysis.")
            sceNorm <- sceNorm[, SingleCellExperiment::sizeFactors(sceNorm) > 0]
        }
        
        ## Put SYMBOL names in rownames
        sceNorm <- .normalizeCon(sceNorm)
        if(removeNoSymbol == TRUE) sceNorm <- .removeNoSymbol(sceNorm)
        
        setSceNorm(theObject) <- sceNorm
        
        return(theObject)
    })
