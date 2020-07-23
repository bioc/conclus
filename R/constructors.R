#' @description
#' Constructor of the class scRNAseq (might be extended in the future).
#'   
#' @rdname constructors
#' @name constructors
#' @title constructors
NULL

#' scRNAseq
#' 
#' @usage
#' scRNAseq(experimentName, countMatrix, species, outputDirectory)
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
#' @rdname constructors
#' @aliases scRNAseq
#' @seealso scRNAseq-class
#' @export scRNAseq
scRNAseq <- function(experimentName, countMatrix, species, outputDirectory){
	
	new("scRNAseq",
			experimentName=experimentName, 
			countMatrix=countMatrix, 
			species=species, 
			outputDirectory=outputDirectory)
}



