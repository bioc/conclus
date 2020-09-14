#' @description
#' Constructor of the class scRNAseq (might be extended in the future).
#'
#' @rdname constructors
#' @name constructors
#' @title constructors
NULL

#' singlecellRNAseq
#'
#' @usage
#' singlecellRNAseq(experimentName, countMatrix, species, outputDirectory)
#'
#' @param experimentName character string representing the name of the
#' experiment. Many output of scRNAseq will use this name.
#' @param countMatrix An integer matrix representing the raw count matrix with
#' reads or unique molecular identifiers (UMIs).
#' @param species character string representing the species of interest.
#' Currently limited to "mouse" and "human".
#' @param outputDirectory A character string of the path to the root output
#' folder.
#' 
#' @return Object of class scRNAseq
#'
#' @rdname constructors
#' 
#' @examples
#' experimentName <- "Bergiers"
#' countMatrix <- as.matrix(read.delim(system.file(
#' "extdata/test_countMatrix.tsv", package="conclus")))
#' outputDirectory <- "YourOutputDirectory"
#' columnsMetaData <- read.delim(
#' system.file("extdata/Bergiers_colData_filtered.tsv", package="conclus"))
#'
#' ## Create the initial object
#' scr <- singlecellRNAseq(experimentName = experimentName,
#'                 countMatrix     = countMatrix,
#'                 species         = "mouse",
#'                 outputDirectory = outputDirectory)
#'                 
#' @aliases singlecellRNAseq
#' @seealso scRNAseq-class
#' @export singlecellRNAseq
singlecellRNAseq <- function(experimentName, countMatrix, species,
        outputDirectory, tSNElist=list(new("Tsne"))){

    new("scRNAseq",
            experimentName=experimentName,
            countMatrix=countMatrix,
            species=species,
            outputDirectory=outputDirectory,
            tSNEList=tSNElist)
}
