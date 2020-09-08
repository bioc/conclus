#' .annotateGenes
#' Filling the rowData
#'
#' This function returns a rowData with the same rownames as in countMatrix
#' with annotations. Only genes having an ENSEMBL IDs or SYMBOL ID
#' will receive the annotation. This function use \code{\link{bioMart}} and
#' \code{\link{AnnotationDbi}} to retrieve annotations.
#'
#' @param countMatrix The raw count matrix
#' @param species The studied species
#' @param rowdataDF rowData of class data.frame, it contains gene names of the
#' count matrix.
#'
#' @keywords internal
#'
#' @importFrom biomaRt useMart getBM
#' @importFrom AnnotationDbi keys select
#' @return Returns the rowData filled with annotations
#' @noRd
.annotateGenes <- function(countMatrix, species, rowdataDF){


    if(isTRUE(all.equal(species, "human"))){

        # suppressMessages(library(org.Hs.eg.db, warn.conflicts=F))
        genomeAnnot <- org.Hs.eg.db::org.Hs.eg.db
        ensemblPattern <- "ENSG"
        ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

    }else if(isTRUE(all.equal(species, "mouse"))){

        # suppressMessages(library(org.Mm.eg.db, warn.conflicts=F))
        genomeAnnot <- org.Mm.eg.db::org.Mm.eg.db
        ensemblPattern <- "ENSMUSG"
        ensembl <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")

    }else
        stop("Species should be human or mouse: ", species, "is not ",
                "currently supported.")



    allnames <- rownames(countMatrix)
    lengthEnsemblGenes <- length(grep(ensemblPattern, allnames))
    lengthSymbols <- length(intersect(AnnotationDbi::keys(genomeAnnot,
                            keytype = "SYMBOL"), allnames))

    if(isTRUE(all.equal(lengthEnsemblGenes, 0)) &&
            isTRUE(all.equal(lengthSymbols, 0)))
        stop("The row names of the matrix neither contain ensembl genes or",
                "official gene symbols.")



    ensemblGenes <- rownames(countMatrix)[grep(ensemblPattern,
                    rownames(countMatrix))]
    symbolGenes <- rownames(countMatrix)[!grepl(ensemblPattern,
                    rownames(countMatrix))]

    if(!isTRUE(all.equal(length(ensemblGenes), 0))){

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
    }

    if(!isTRUE(all.equal(lengthSymbols, 0))){

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
    }


    rowdata <- base::rbind(rowdataSymbol, rowdataEnsembl)
    ## Filtering duplicated symbols
    (multSym <- rowdata$SYMBOL[!is.na(rowdata$SYMBOL) &
                                duplicated(rowdata$SYMBOL)])

    ## We don't combine such ensembls, but leave them unique with
    ## "symbol_ensembl"
    (rowdata$SYMBOL[rowdata$SYMBOL %in% multSym] <-
                paste0(rowdata$SYMBOL[rowdata$SYMBOL %in% multSym],
                        rowdata$ENSEMBL[rowdata$SYMBOL %in% multSym]))



    message("Retrieving information about genes from biomaRt.")

    res <- getBM(attributes=c("ensembl_gene_id", "go_id", "name_1006",
                    "chromosome_name", "gene_biotype"), mart=ensembl)
    tmp <- res[!duplicated(res$ensembl_gene_id),]

    rowdata <- merge(rowdata, tmp[c("ensembl_gene_id", "chromosome_name",
                            "gene_biotype")], by.x = "ENSEMBL",
            by.y = "ensembl_gene_id", all.x = TRUE, all.y = FALSE,
            sort = FALSE)

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


    if (!is.null(rowdataDF)){

        rowdataDF$nameInCountMatrix <- rownames(rowdataDF)
        rowdata <- merge(rowdataDF, rowdata, by.x = "nameInCountMatrix",
                by.y = "nameInCountMatrix", all.x = TRUE, all.y = TRUE,
                sort = FALSE)
    }

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
#' This function creates a boolean vector that is used to filer the cells
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
        threshold <- 2.5*quan[2] - 1.5*quan[4]
        if (threshold < 0){
            threshold <- (quan[1]+quan[2]) / 2
            }
        vec <- as.numeric(colData[ ,columnName] >=  as.numeric(threshold))
    } else if (columnName %in% mw) {
        threshold <- 2.5*quan[4] - 1.5*quan[2]
        if (threshold > quan[5]){
            threshold <- (quan[3]+quan[4]) / 2
        }
        vec <- as.numeric(colData[ ,columnName] <=  as.numeric(threshold))
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
    countMatrix <- countMatrix[,colSums(countMatrix) > genesSumThr]
    colData <- colData[colnames(countMatrix),]
    mb <- MoreBetter
    mw <- MoreWorse
    columnNames <- c(mb,mw)

    ## Create the report table
    reportTable <- unlist(lapply(columnNames, .createReportTable, colData,
                                    mb, mw))
    reportTable <- matrix(reportTable, ncol=length(columnNames), )
    colnames(reportTable) <- columnNames

    ## Add rownames to the report table
    rownames(reportTable) <- colData$cellName

    ## Add cell names column to the report table
    reportTable <- data.frame(cellName=colData$cellName, reportTable)

    ### add columns with filtering score and verdict ###
    reportTable <- dplyr::mutate(reportTable, score=NA)
    reportTable$score <- rowSums(reportTable[,colnames(reportTable) %in%
                                                c(mb,mw)])
    reportTable <- dplyr::mutate(reportTable, filterPassed=NA)
    reportTable$filterPassed[reportTable$score >= length(mb)+length(mw)] <- 1
    reportTable$filterPassed[reportTable$score < length(mb)+length(mw)] <- 0

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

    countMatrix <- countMatrix[,colData$filterPassed == 1]
    colData <- colData[colData$filterPassed == 1,]

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
    coldata <- dplyr::mutate(coldata, genesNum=NA, genesSum=NA, oneUMI=NA)
    coldata$genesSum <- colSums(countMatrix)
    coldata$genesNum <- vapply(colnames(countMatrix), .fillGenesNumColumn,
            coldata, countMatrix, FUN.VALUE=integer(1))
    coldata$oneUMI <- vapply(colnames(countMatrix), .fillOneUmmiColumn, coldata,
            countMatrix, FUN.VALUE=integer(1))
    coldata <- dplyr::mutate(coldata,
                            oneUMIper =100 * coldata$oneUMI / coldata$genesNum)

    ### add info about mitochondrial and protein-coding genes
    coldata <- dplyr::mutate(coldata, mtGenes = NA, mtSum = NA, codGenes = NA,
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
            mtPer = 100*coldata$mtGenes/coldata$genesNum,
            codPer = 100*coldata$codGenes/coldata$genesNum,
            sumMtPer = 100*coldata$mtSum/coldata$genesSum,
            sumCodPer = 100*coldata$codSum/coldata$genesSum)

    if (!is.null(coldataDF)){

        coldataDF$cellName <- rownames(coldataDF)
        coldata <- merge(coldataDF, coldata, by.x = "cellName",
                by.y = "cellName", all.x = FALSE, all.y = TRUE, sort = FALSE)
    }

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

    selRows <- ((rowSums(countMatrix[, ] >= 1)) > 10)
    countMatrix <- countMatrix[selRows, ]
    rowData <- rowData[rowData$nameInCountMatrix %in% rownames(countMatrix), ]

    return(list(countMatrix, rowData))
}

.checkParamNorm <- function(sizes, rowdata, coldata, alreadyCellFiltered,
        runQuickCluster){

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
#'                     runQuickCluster=TRUE)
#' 
#' 
#' @param theObject A scRNAseq object
#' @param sizes  Vector of size factors from scran::computeSumFactors()
#' function.
#' @param coldata Data frame containing cells informations. Default is NULL.
#' @param rowdata Data frame containing genes informations. Default is NULL.
#' @param alreadyCellFiltered Logical. If TRUE, quality check and
#' filtering will not be applied.
#' @param runQuickCluster Logical. If TRUE scran::quickCluster() function will
#' be applied. It usually improves the normalization for medium-size count
#' matrices. However, it is not recommended for datasets with less than 200
#' cells and may take too long for datasets with more than 10000 cells.
#'
#' @details
#' This function uses the normalization method of the scater package. For more
#' details about the normalization used see ?scater::normalize. The size factors
#' used in the normalization are calculated with scran::computeSumFactors.
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
#' @return Returns a scRNASeq object with its SCEnorm slot updated. This slot
#' contains a SingleCellExperiment object having the normalized count matrix,
#' the colData (table with cells informations), and the rowData (table with the
#' genes informations). See ?SingleCellExperiment for more details.
#'
#' @examples
#' experimentName <- "Bergiers"
#' countMatrix <- as.matrix(read.delim(system.file(
#' "extdata/test_countMatrix.tsv", package="conclus")))
#' outputDirectory <- "./"
#' columnsMetaData <- read.delim(
#' system.file("extdata/Bergiers_colData_filtered.tsv", package="conclus"))
#'
#' scr <- singlecellRNAseq(experimentName = experimentName,
#'                 countMatrix     = countMatrix,
#'                 species         = "mouse",
#'                 outputDirectory = outputDirectory)
#'
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
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
                            runQuickCluster=TRUE){

        validObject(theObject)

        .checkParamNorm(sizes, rowdata, coldata, alreadyCellFiltered,
                runQuickCluster)

        countMatrix <- getCountMatrix(theObject)
        species <- getSpecies(theObject)
        rowdata <- .annotateGenes(countMatrix, species = species,
                rowdataDF = rowdata)
        coldata <- .addCellsInfo(countMatrix, rowdataDF = rowdata,
                coldataDF = coldata)


        if (!alreadyCellFiltered){

            filterCellsResult <- .filterCells(countMatrix, coldata)
            countMatrix <- filterCellsResult[[1]]
            coldata <- filterCellsResult[[2]]
        }

        filterGenesResult <- .filterGenes(countMatrix, rowdata)
        countMatrix <- filterGenesResult[[1]]
        rowdata <- filterGenesResult[[2]]

        stopifnot(all(rownames(countMatrix) == rownames(rowdata)))
        stopifnot(all(colnames(countMatrix) == rownames(coldata)))

        sce <-
            SingleCellExperiment::SingleCellExperiment(assays=list(
                                                counts=as.matrix(countMatrix)),
                                                colData=coldata,
                                                rowData=rowdata)

        # normalization
        message("Running normalization. It can take a while depending on the",
                " number of cells.")

        if(runQuickCluster)
            cl <- tryCatch(scran::quickCluster(sce), error=function(e) NULL)
        else
            cl <- NULL

        # Compute sizeFactors which will be used for normalization
        sceNorm <- scran::computeSumFactors(sce, sizes=sizes, clusters=cl)

        message("summary(sizeFactors(sceObject)):")
        print(summary(SingleCellExperiment::sizeFactors(sceNorm)))

        if (length(
            SingleCellExperiment::sizeFactors(
                sceNorm)[SingleCellExperiment::sizeFactors(sceNorm) <= 0]) > 0)
            message("Cells with negative sizeFactors will be deleted before the
                        downstream analysis.")

        sceNorm <- sceNorm[, SingleCellExperiment::sizeFactors(sceNorm) > 0]
        sceNorm <- scater::logNormCounts(sceNorm)
        setSceNorm(theObject) <- sceNorm
        return(theObject)
    })






