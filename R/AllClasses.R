## All classes defined for the scRNAseq package


################################################################################
############################### Tsne class  ####################################
################################################################################

#' The Tsne class
#'
#' S4 class containing the features to plot tSNEs. This constructor is internal
#' and is used by the method generateTSNECoordinates.
#'
#' @name Tsne-class
#' @rdname Tsne-class
#' @aliases getTsneName getPC getPerplexity getCoordinates
#' @aliases setTsneName setPC setPerplexity setCoordinates
#' @exportClass Tsne
#' 
#' @section Details:
#'  \describe{
#'	  Tsne is a vector of principal commponents (PC) and perplexity that are
#'    the parameters necessary to reduce the dimensionality of the data in the
#'   the form of a t-distributed stochastic neighbor embedding (t-SNE). For 
#'   details about perplexities parameter see ‘?Rtsne’. This information is 
<<<<<<< HEAD
#'   stored in three components:
#'
#'    \item{\code{name}:}{A \code{"character"} string representing the name of 
#'    the tSNE coordinates.}
#'    \item{\code{pc}:}{A \code{"numeric"} vector representing the number
#'    of principal components used by CONCLUS to perfom a PCA before performing
#'    a tSNE. Default: c(8, 10, 20, 40, 50, 80, 100).}
#'    \item{\code{perplexity}:}{A \code{"numeric"} vector. Default: c(30, 40) }
#'    \item{\code{coordinates}:}{A \code{"numeric"} \code{"matrix"} that 
#'    contains the coordinates of one tSNE solution.}
#'    }
#'
#'
#' @section Constructor:
#' \describe{
#'    Tsne(name = "character", pc = "numeric", perplexity = "numeric", 
#'         coordinates = "matrix")
#' 
#'    \item{\code{name}:}{Empty character string or the name of the tSNE.}
#'    \item{\code{pc}:}{Empty \code{"numeric"} or vector of the number
#'    of PCs.}
#'    \item{\code{perplexity}:}{Empty \code{"numeric"} or vector of perplexity 
#'     values.}
#'    \item{\code{coordinates}:}{Empty \code{"numeric"} \code{"matrix"} or 
#'    matrix of coordinates.}}
#'
#'  
#' @section Accessors:
#'   \describe{
#'     In the following snippets, x is a Tsne object.
#' 
#'     \item{\code{getTsneName(x)}: Get the name of the tSNE.}
#'     \item{\code{getPC(x)}: Get all PC parameters.} 
#'     \item{\code{getPerplexity(x)}: Get all perplexity parameters.}
#'     \item{\code{getCoordinates(x)}: Get the matrix of tSNE coordinates.}}
#' 
#' 
#' @section Subsetting:
#'   \describe{
#'     In the following snippets, x is a Tsne object.
#' 
#'     \item{\code{setTsneName(x) <- value}: Set the name of the tSNE.}
#'     \item{\code{setPC(x) <- value}: Set all PC parameters.} 
#'     \item{\code{setPerplexity(x) <- value}: Set all perplexity parameters.}
#'     \item{\code{setCoordinates(x) <- value}: Set the matrix of tSNE coordinates.}}
#' 
#' @section Authors:
#' \describe{Ilyess Rachedi}
#' 
#' #' @section See also:
#' \describe{
#'   \code{\link{generateTSNECoordinates}}
#' }  
#' }

=======
#'   stored in four components:
#'
#'    \item{\code{name}:}{A \code{"character"} string representing the name of 
#'    the tSNE coordinates.}
#'    \item{\code{pc}:}{A \code{"numeric"} value representing the number
#'    of principal components used by CONCLUS to perfom a PCA before performing
#'    a tSNE. 
#'    \item{\code{perplexity}:}{A \code{"numeric"}. Default: c(30, 40) }
#'    \item{\code{coordinates}:}{A \code{"numeric"} \code{"matrix"} that 
#'    contains the coordinates of one tSNE solution.}
#'    }}
#'
#'
#' @section Constructor:
#' \describe{
#'    Tsne(name = "character", pc = "numeric", perplexity = "numeric", 
#'         coordinates = "matrix")
#' 
#'    \item{\code{name}:}{Empty character string or the name of the tSNE.}
#'    \item{\code{pc}:}{Empty \code{"numeric"} number of PCs.}
#'    \item{\code{perplexity}:}{Empty \code{"numeric"} perplexity values.}
#'    \item{\code{coordinates}:}{Empty \code{"numeric"} \code{"matrix"} or 
#'    matrix of coordinates.}}
#'
#'  
#' @section Accessors:
#'   \describe{
#'     In the following snippets, x is a Tsne object.
#' 
#'     \item{\code{getTsneName(x)}: Get the name of the tSNE.}
#'     \item{\code{getPC(x)}: Get the PC used} 
#'     \item{\code{getPerplexity(x)}: Get the perplexity used}
#'     \item{\code{getCoordinates(x)}: Get the matrix of tSNE coordinates.}}
#' 
#' 
#' @section Subsetting:
#'   \describe{
#'     In the following snippets, x is a Tsne object.
#' 
#'     \item{\code{setTsneName(x) <- value}: Set the name of the tSNE.}
#'     \item{\code{setPC(x) <- value}: Set the PC parameter.} 
#'     \item{\code{setPerplexity(x) <- value}: Set the perplexity parameter.}
#'     \item{\code{setCoordinates(x) <- value}: Set the matrix of tSNE 
#'     coordinates.}}
#' 
#' @author Ilyess Rachedi
#' @seealso \code{\link{generateTSNECoordinates}}
>>>>>>> branch 'roxygen' of https://git.embl.de/rachedi/concluspoo.git
Tsne <- setClass(
    "Tsne",
    slots = c(
        name = "character",
        pc = "numeric",
        perplexity = "numeric",
        coordinates = "matrix"
    ),
	
    validity = function(object){
		
		coordinates <- getCoordinates(object)
		
        if(!isTRUE(all.equal(ncol(coordinates),2)))
			stop("Coordinates should be a matrix with two columns X and Y.")
            
    })



################################################################################
############################### Dbscan class ###################################
################################################################################

#' The Dbscan class
#'
#' S4 class containing the features to plot DBSCAN. This constructor is internal
#' and is used by the method runDBSCAN.
#'
#' @name Dbscan-class
#' @rdname Dbscan-class
#' @aliases getDbscanName getEpsilon getMinPoints getClustering
#' @aliases setDbscanName setEpsilon setMinPoints setClustering
#' @exportClass Dbscan
#' 
#' @section Details:
#'  \describe{
#'    \item{\code{name}:}{A \code{"character"} string representing the name of 
#'    the Dbscan clustering.}
#'    \item{\code{epsilon}:}{\code{"numeric"} value. The epsilon
#'    is the distance to consider two points belonging to the same cluster.
#'    Default epsilon= c(1.3, 1.4, 1.5)}
#'    \item{\code{minPoints}:}{\code{"numeric"} value. The minPoints
#'    is the minimum number of points to construct a cluster.}
#'    \item{\code{clustering}:}{Variable of class \code{"matrix"}, it contains
#'    the result of one DBSCAN clustering solution.}}
#'
#'
#' @section Constructor:
#' \describe{
#'    Dbscan(name = "character", epsilon = "numeric", minPoints = "numeric", 
#'         clustering = "matrix")
#' 
#'    \item{\code{name}:}{Empty character string or the name of the tSNE.}
#'    \item{\code{epsilon}:}{Empty \code{"numeric"} representing the epsilon.}
#'    \item{\code{minPoints}:}{Empty \code{"numeric"} representing the
#'     minPoints value.}
#'    \item{\code{clustering}:}{Empty \code{"numeric"} \code{"matrix"} or 
#'    matrix of clustering}}
#'
#'  
#' @section Accessors:
#'   \describe{
#'     In the following snippets, x is a Tsne object.
#' 
#'     \item{\code{getDbscanName(x)}: Get the name of the Dbscan}
#'     \item{\code{getEpsilon(x)}: Get the epsilon used.} 
#'     \item{\code{getMinPoints(x)}: Get the MinPoint used.}
#'     \item{\code{getClustering(x)}: Get the matrix of clustering DBSCAN}}
#' 
#' 
#' @section Subsetting:
#'   \describe{
#'     In the following snippets, x is a Tsne object.
#' 
#'     \item{\code{setDbscanName(x) <- value}: Set the name of the Dbscan}
#'     \item{\code{setEpsilon(x) <- value}: Set the epsilon used.} 
#'     \item{\code{setMinPoints(x) <- value}: Set the minPoints used}
#'     \item{\code{setClustering(x) <- value}: Set the matrix of Dbscan
#'     clustering..}}
#' 
#' @author Ilyess Rachedi
#' 
#' #' @section See also:
#' \describe{
#'   \code{\link{runDBSCAN}}
#' }  

Dbscan <- setClass(
    "Dbscan",
    slots = c(
        name = "character",
        epsilon    = "numeric",
        minPoints  = "numeric",
        clustering = "matrix"
    ))



################################################################################
############################## scRNAseq class ##################################
################################################################################

#' The scRNAseq class
#'
#' S4 class and the main class used by CONCLUS containing the different steps
#' to analyse rare cell populations.
#'
#' @name scRNAseq-class
#' @rdname scRNAseq-class
#' @aliases getExperimentName getCountMatrix getSceNorm getSpecies
#' getOutputDirectory getTSNEList getDbscanList getCellsSimilarityMatrix
#' getClustersSimilarityMatrix getClustersSimiliratyOrdered
#' getMarkerGenesList getClustersMarkers getGenesInfos
#' @aliases setExperimentName setCountMatrix setSceNorm setSpecies
#' setOutputDirectory setTSNEList setDbscanList setCellsSimilarityMatrix
#' setClustersSimilarityMatrix setClustersSimiliratyOrdered
#' setMarkerGenesList setClustersMarkers setGenesInfos
#' @exportClass scRNAseq
#' @import SingleCellExperiment
#' @importFrom methods new
#'
#' @section Details:
#' \describe{
#'    \item{\code{experimentName}:}{\code{"character"} string representing 
#'    the name of the experiment. Many output of scRNAseq will use
#'    this name.}
#'    \item{\code{countMatrix}:}{\code{"integer"}{\code{"matrix"} representing,
#'    the raw matrix with reads or unique molecular identifiers (UMIs).}}
#'    \item{\code{sceNorm}:}{Object of class SingleCellExperiment.
#'    Contain the colData giving informations about cell;
#'    the rowData giving informations about genes and the normalized
#'    count matrix. Fill it with \code{\link{normaliseCountMatrix}.}}
#'    \item{\code{species}:}{\code{"character"} string representing the species
#'     of interest. Actually it's limited to "mouse" or "human".}
#'    \item{\code{outputDirectory}:}{\code{"character"} string representing the
#'     path where the outputs have to go.}
#'    \item{\code{tSNEList}:}{List of \code{"Tsne"} objects representing the
#'    different tSNE coordinates generated by CONCLUS.}
#'    \item{\code{dbscanList}:}{List of \code{"Dbscan"} objects representing the
#'    different Dbscan clustering generated by CONCLUS.}
#'    \item{\code{cellsSimilarityMatrix}:}{ cells * cells similarity
#'     \code{matrix}. Define how many times twos cell have been seen together
#'     across the 84 solutions of clustering }
#'    \item{\code{clustersSimilarityMatrix}:}{\code{"numeric"} \code{"matrix"}
#'    comparing the robustness of the consensus clusters.}
#'    \item{\code{clustersSimiliratyOrdered}:}{\code{factor} representing 
#'    the clusters ordered by similarity.}
#'    \item{\code{markerGenesList}:}{\code{list} of data.frames. Each data frame 
#'    contains the ranked genes of one cluster.}
#'    \item{\code{clustersMarkers}:}{\code{data.frame} containing the top 10 
#'    by default of marker genes of each clusters.}
#'    \item{\code{genesInfos}:}{\code{data.frame} containing informations of
#'     the markers genes each clusters.}
#' }
#'     
#'     
#' @section Constructor:
#' \describe{
#'    \item{\code{experimentName}:}{string of the name of the experiment.}
#'    \item{\code{countMatrix}:}{\code{matrix} containing the counts}
#'    \item{\code{sceNorm}:}{Empty or \code{SingleCellExperiment object }
#'    containing, the count matrix, the colData and the rowData.}
#'    \item{\code{species}:}{\code{"character"} representig the species of
#'    interest.}
#'    \item{\code{outputDirectory}:}{\code{"character"} representing the path
#'    of the output directory.}
#'    \item{\code{tSNEList}:}{Empty or list of \code{Tsne} object.}
#'    \item{\code{dbscanList}:}{Empty or list of \code{Dbscan} object.}
#'    \item{\code{cellsSimilarityMatrix}:}{Empty or cells * cells \code{matrix} 
#'    with similarity scores.}
#'    \item{\code{clustersSimilarityMatrix}:}{Empty or clusters * clusters
#'    \code{matrix} with similarity scores.}
#'    \item{\code{clustersSimiliratyOrdered}:}{Empty or \code{factor}
#'     representing the clusters ordered by similarity.}
#'    \item{\code{markerGenesList}:}{Empty of \code{list} of data.frames. 
#'    Each data frame contains the ranked genes of one cluster.}
#'    \item{\code{clustersMarkers}:}{Empty or \code{data.frame} containing the 
#'    top 10 by default of marker genes of each clusters.}
#'    \item{\code{genesInfos}:}{Empty or \code{data.frame} containing 
#'    informations about the markers genes each clusters.}}
#'    
#'    
#' @section Accessors:
#'   \describe{
#'     In the following snippets, x is a scRNAseq object.
#' 
#'     \item{\code{getExperimentName(x)}: get the name of the experiment}
#'     \item{\code{getCountMatrix(x)}: get the count matrix.} 
#'     \item{\code{getSceNorm(x)}: get the SingleCellExperiment object used.}
#'     \item{\code{getSpecies(x)}: get the species.}  
#'     \item{\code{getOutputDirectory(x)}: get the path of the output directory.} 
#'     \item{\code{getTSNEList(x)}: get the list of Tsne objects.}
#'     \item{\code{getDbscanList(x)}: get the list of Dbscan objects.}
#'     \item{\code{getCellsSimilarityMatrix(x)}: get the cell similarity 
#'     matrix.} 
#'     \item{\code{getClustersSimilarityMatrix(x)}: get the cluster similarity
#'      matrix.}
#'     \item{\code{getClustersSimiliratyOrdered(x)}: get the clusters ordered
#'      by similarity}
#'     \item{\code{getMarkerGenesList(x)}: get the list of marker genes by 
#'     clusters.}
#'     \item{\code{getClustersMarkers(x)}: get the more significant markers by 
#'     clusters into a data.frame}
#'     \item{\code{getGenesInfos(x)}: get a data.frame containing informations
#'      about marker genes}}
#' 
#' @section Subsetting:
#'   \describe{
#'     In the following snippets, x is a scRNAseq object.
#' 
#'     \item{\code{setExperimentName(x)}: set the name of the experiment}
#'     \item{\code{setCountMatrix(x)}: set the count matrix.} 
#'     \item{\code{setSceNorm(x)}: set the SingleCellExperiment object used.}
#'     \item{\code{setSpecies(x)}: set the species.} 
#'     \item{\code{setOutputDirectory(x)}: set the path of the output directory.} 
#'     \item{\code{setTSNEList(x)}: set the list of Tsne objects.}
#'     \item{\code{setDbscanList(x)}: set the list of Dbscan objects.}
#'     \item{\code{setCellsSimilarityMatrix(x)}: set the cell similarity 
#'     matrix.} 
#'     \item{\code{setClustersSimilarityMatrix(x)}: set the cluster similarity
#'      matrix.}
#'     \item{\code{setClustersSimiliratyOrdered(x)}: set the clusters ordered
#'      by similarity}   
#'     \item{\code{setMarkerGenesList(x)}: set the list of marker genes by 
#'     clusters.}
#'     \item{\code{setClustersMarkers(x)}: set the more significant markers by 
#'     clusters into a data.frame}
#'     \item{\code{setGenesInfos(x)}: set a data.frame containing informations
#'      about marker genes}}
#' 
#' @author Ilyess Rachedi
#' @seealso \code{\link{scRNAseq}}

scRNAseq <- setClass(
    "scRNAseq",
    slots = c(
        experimentName = "character",
        countMatrix = "matrix",
        sceNorm = "SingleCellExperiment",
        species = "character",
        outputDirectory = "character",
        tSNEList = "list",
        dbscanList = "list",
        cellsSimilarityMatrix = "matrix",
        clustersSimilarityMatrix = "matrix",
        clustersSimiliratyOrdered = "factor",
        markerGenesList = "list",
		clustersMarkers = "data.frame",
		genesInfos = "data.frame"
    ),
    prototype = list(
        sceNorm = SingleCellExperiment(),
        tSNEList = list(new("Tsne")),
        dbscanList = list(new("Dbscan")),
        clustersSimilarityMatrix = as.matrix(data.frame(1, row.names=1)),  
        clustersSimiliratyOrdered = factor(1),
        markerGenesList = list(data.frame(Gene = c("gene1"),
                                          mean_log10_fdr = c(NA),
                                          n_05 = c(NA), score = c(NA))),
	    clustersMarkers = data.frame(geneName="gene1", clusters=NA),
		
		genesInfos = data.frame(uniprot_gn_symbol=c("symbol"), clusters="1",
				external_gene_name="gene", go_id="GO1,GO2", 
				mgi_description="description", entrezgene_description="descr",
				gene_biotype="gene", chromosome_name="1", Symbol="symbol",
				ensembl_gene_id="ENS", mgi_id="MGI", entrezgene_id="1",
				uniprot_gn_id="ID")
    ),
	
    validity = function(object) {
		     
        ## Test experimentName slot
        experimentName <- getExperimentName(object)
		
        if(isTRUE(all.equal(length(experimentName),0)) ||
				isTRUE(all.equal(experimentName,"")))
			stop("'experimentName' slot is empty. Please fill it.")
            
		if(!is.character(experimentName) | grepl(" ", experimentName))
			stop("Experiment name should contain a single string ",
                "describing the experiment, '", experimentName, 
				"' is not correct.")
            
        ## Test countMatrix slot
        countMatrix <- getCountMatrix(object)
		
		if(isTRUE(all.equal(nrow(countMatrix), 0)) && 
				isTRUE(all.equal(ncol(countMatrix), 0)))
			stop("The count matrix should not be empty. Please fill it.")
		
        if(all(is.na(countMatrix)))
			stop("'countMatrix' slot is empty. It should be a matrix ",
					"containing at leat 100 cells.\n")

        if(ncol(countMatrix) < 100)
			stop("Not enough cells in the count matrix. There ",
                "Should be at leat 100 cells. ", "The current count matrix ",
				"contains ", ncol(countMatrix), " cells.\n")
            
	
		## Test sceNorm slot
		sceNorm <- getSceNorm(object)
			
		if(class(sceNorm) != "SingleCellExperiment")
			stop("Normalized count matrix should be a SingleCellExperiment ",
					"object and not a '", is(sceNorm), "'.")
			
        ## Test species slot
        species <- getSpecies(object)
		
		if(isTRUE(all.equal(length(species), 0)))
			stop("The species is empty. Please fill it.")
		
        if (!is.element(species, c("mouse","human")))
			stop("species should be 'mouse' or 'human'. '", species, 
					"' is currently not supported.\n")
			
        ## Test outputDirectory slot
        outputDirectory <- getOutputDirectory(object)
		
        if (isTRUE(all.equal(length(outputDirectory), 0)) ||
				isTRUE(all.equal(outputDirectory, "")))
			stop("'outputDirectory' slot is empty. Please fill it.")
            		
		if(!is.character(outputDirectory) | grepl(" ", outputDirectory))
			stop("'outputDirectory' should be a conform folder path:",
                "'", outputDirectory, "' is not.")
            
        ## Test tSNEList slot
        tSNEList <- getTSNEList(object)
		
		if(isTRUE(all.equal(length(tSNEList), 0)))
			stop("tSNEList is empty. This should be a list of tSNE objects.\n")
				
		invisible(checkList(tSNEList, getCoordinates, "Tsne"))
		
        
        ## Test dbscan slot
        dbscanList <- getDbscanList(object)
        
		if(isTRUE(all.equal(length(dbscanList), 0)))
			stop("dbscanList is empty. This should be a list of dbScan ",
					"objects.\n")
		
		invisible(checkList(dbscanList, getClustering, "Dbscan"))
        
        
        ## Test cellsSimilarityMatrix slot
        cellsSimilarityMatrix <- getCellsSimilarityMatrix(object)
        
		if(!isTRUE(all.equal(nrow(cellsSimilarityMatrix), 
						ncol(cellsSimilarityMatrix))))
			stop("'cellsSimilarityMatrix' slot should contain a square matrix.")
            
		
        ## Test clustersSimilarityMatrix slot
        clustersSimilarityMatrix <- getClustersSimilarityMatrix(object)
		
        if (!isTRUE(all.equal(nrow(clustersSimilarityMatrix), 
						ncol(clustersSimilarityMatrix))))
			stop("'clustersSimilarityMatrix' slot should contain a square ",
					"matrix. ")
        
        
        ## Test clustersSimiliratyOrdered slot
        clustersSimiliratyOrdered <- getClustersSimiliratyOrdered(object)
		
		if(isTRUE(all.equal(nrow(clustersSimiliratyOrdered), 0)) &&
				isTRUE(all.equal(ncol(clustersSimiliratyOrdered), 0)))
			stop("'clustersSimiliratyOrdered' is empty. It should be a matrix.")
		
        if (all(!clustersSimiliratyOrdered %in% 
						rownames(clustersSimilarityMatrix)))
			stop("'clustersSimiliratyOrdered' slot should contain the same ",
					"clusters as 'clustersSimilarityMatrix'.")
        
        
        ## Test getMarkerGenesList slot
        markerGenesList <- getMarkerGenesList(object)
		
		if(isTRUE(all.equal(length(markerGenesList), 0)))
			stop("markerGenesList is empty. This should be a list of dataframe")
				
        invisible(checkMarkerGenesList(markerGenesList, 
						clustersSimiliratyOrdered))
		
		## Test clustersMarkers
		clusterMarkers <- getClustersMarkers(object)
		
		if(isTRUE(all.equal(length(clusterMarkers), 0)))
			stop("clusterMarkers is empty. This should be a dataframe")
		
		invisible(checkClusterMarkers(clusterMarkers, 
						clustersSimiliratyOrdered))
		
		## Test genesInfos slot
		genesInfos <- getGenesInfos(object)
		
		if(isTRUE(all.equal(length(genesInfos), 0)))
			stop("genesInfos is empty. This should be a dataframe")
		
		invisible(checkGenesInfos(genesInfos, species, 
						clustersSimiliratyOrdered))
    }
)

