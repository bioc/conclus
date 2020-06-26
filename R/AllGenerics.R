################################################################################
########################### scRNAseq class methods #############################
################################################################################

#' @rdname normaliseCountMatrix-scRNAseq
#' @aliases normaliseCountMatrix,normaliseCountMatrix-method
setGeneric(
		
	name = "normaliseCountMatrix",
		
	def = function(theObject, sizes=c(20,40,60,80,100), rowdata=NULL,
			coldata=NULL, alreadyCellFiltered=FALSE, runQuickCluster=TRUE){
			standardGeneric("normaliseCountMatrix")    
    },
    signature="theObject") 



setGeneric(
		
	name = "testClustering",
		
	def = function(theObject, dbscanEpsilon=1.4, minPts=5, 
				perplexities=c(30), PCs=c(4), randomSeed=42, width=7, height=7,
				onefile=FALSE, ...){
			standardGenericic("testClustering")    
    },
    signature="theObject")



setGeneric(
		
	name = "generateTSNECoordinates",
		
	def = function(theObject, randomSeed=42, cores=1, 
			PCs=c(4, 6, 8, 10, 20, 40, 50), perplexities=c(30,40), 
			writeOutput=FALSE){
		standardGeneric("generateTSNECoordinates")    
    },
    signature="theObject")



setGeneric(
		
	name = "runDBSCAN",
		
	def = function(theObject, cores=1, epsilon=c(1.3, 1.4, 1.5), 
			minPoints=c(3, 4), writeOutput=FALSE){
		standardGeneric("runDBSCAN")    
    },
    signature="theObject")



setGeneric(
    
	name = "clusterCellsInternal",
		
	def = function(theObject, clusterNumber=0, deepSplit=4, cores=1,
				clusteringMethod="ward.D2"){
		standardGeneric("clusterCellsInternal")    
    },
    signature="theObject")



setGeneric(
		
	name = "calculateClustersSimilarity",
		
	def = function(theObject, clusteringMethod="ward.D2"){
		standardGeneric("calculateClustersSimilarity")    
    },
    signature="theObject")



setGeneric(
		
	name = "addClusteringManually",
		
	def = function(theObject, fileName, columnName="clusters"){
		standardGeneric("addClusteringManually")    
    },
    signature="theObject")



setGeneric(
		
	name = "plotCellSimilarity",
		
	def = function(theObject, colorPalette="default", 
				statePalette="default", clusteringMethod="ward.D2",
				orderClusters=FALSE, plotPDF=TRUE, returnPlot=FALSE, width=7, 
				height=6, ...){
		standardGeneric("plotCellSimilarity")    
    },
    signature="theObject")



setGeneric(
		
	name = "plotClusteredTSNE",
		
	def = function(theObject, tSNEresExp="", colorPalette="default",
				PCs=c(4, 6, 8, 10, 20, 40, 50), perplexities=c(30, 40),
				columnName="clusters", returnPlot=FALSE, width=6, height=5, 
				onefile=FALSE, ...){
		standardGeneric("plotClusteredTSNE")
			
    },
    signature="theObject")



setGeneric(
		
	name = "plotCellHeatmap",
		
	def = function(theObject, fileName, meanCentered=TRUE, 
				colorPalette="default", statePalette="default", 
				clusteringMethod="ward.D2", orderClusters=FALSE, 
				similarity=FALSE, orderGenes=FALSE, returnPlot=FALSE,
				saveHeatmapTable=FALSE, width=10, height=8.5, ...){
		standardGeneric("plotCellHeatmap")    
    },
    signature="theObject")



setGeneric(
		
	name = "plotGeneExpression",
		
	def = function(theObject, geneName, graphsDirectory="pictures", 
			palette=c("grey","red", "#7a0f09", "black"), returnPlot=FALSE,
			tSNEpicture=1, commentName="", savePlot=TRUE, alpha=1, limits=NA,
			pointSize=1, width=6, height=5, ...){
		standardGeneric("plotGeneExpression")    
    },
    signature="theObject")



setGeneric(
		
	name = "plotClustersSimilarity",
		
	def = function(theObject, colorPalette="default", 
				statePalette="default", clusteringMethod="ward.D2", 
				returnPlot=FALSE, width=7, height=5.5, ...){
		standardGeneric("plotClustersSimilarity")    
    },
    signature="theObject")



setGeneric(
		
	name = "exportResults",
		
	def = function(theObject, saveCellsSimilarityMatrix=TRUE,
                 saveClustersSimilarityMatrix=TRUE, saveNormalizedMatrix=TRUE,
                 saveColData=TRUE, saveRowData=TRUE, saveWorkspace=TRUE,
                 saveClusteringResults=TRUE){
			 standardGeneric("exportResults")
        },
    signature="theObject")



setGeneric(
		
	name = "rankGenes",
		
	def = function(theObject, column="clusters", writeMarkerGenes=FALSE){
		standardGeneric("rankGenes")
    },
    signature="theObject")



setGeneric(
    name="bestClustersMarkers",
    def=function(theObject,
                 genesNumber=10,
                 removeDuplicates = TRUE){
        standardGeneric("bestClustersMarkers")
    },
    signature="theObject")



setGeneric(
    name="getGenesInfos",
    def=function(theObject,
                 groupBy="clusters",
                 orderGenes="initial",
                 getUniprot=TRUE,
                 silent=FALSE,
                 cores=1){
        standardGeneric("getGenesInfos")
    },
    signature="theObject")



setGeneric(
    name = "saveMarkersLists",
    def = function(theObject,
                   pattern = "genes.tsv",
                   Ntop = 100) {
        standardGeneric("saveMarkersLists")
    },
    signature = "theObject"
)



setGeneric(
    name = "saveGenesInfos",
    def = function(theObject,
                   sep = ";",
                   header = TRUE,
                   startFromCluster = 1,
                   groupBy = "clusters",
                   # getGenesInfos params
                   orderGenes = "initial",
                   getUniprot = TRUE,
                   silent = FALSE,
                   cores = 1) {
        standardGeneric("saveGenesInfos")
    },
    signature = "theObject"
)


################################################################################
############################### Getter methods #################################
################################################################################

#' Getter of experimentName slot of scRNAseq class
#'
#' Get the experiment name assigned to the experimentName slot. Experiment name
#' describes the name of th experiment of the Conclus analysis.
#' 
#' @name getExperimentName
#' @rdname getExperimentName-scRNAseq
#' @importFrom methods setMethod
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' 
#' @return a string corresponding to the experiment name
#' @export
setGeneric(
		
	name = "getExperimentName",
	
    def = function(theObject){
        standardGeneric("getExperimentName")    
    },
    signature="theObject")



#' Getter of countMatrix slot of scRNAseq class
#'
#' Get the the count matrix assigned to the countMatrix slot.
#' 
#' @name getCountMatrix
#' @rdname getCountMatrix-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#'
#' @return The matrix assigned to the countMatrix slot.
#' @export
setGeneric(
		
    name = "getCountMatrix",
	
    def = function(theObject){
        standardGeneric("getCountMatrix")    
    },
    signature = "theObject")



#' Getter of sceNorm slot of scRNAseq class
#'
#' Get the SCE object assigned to the sceNorm slot. This contains the colData,
#' the rowData, and the normalized count matrix.
#'
#' @name getSceNorm
#' @rdname getSceNorm-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#' @param theObject scRNAseq object
#'
#' @return The SCE object assigned to the sceNorm slot.
#' @export
setGeneric(
		
    name = "getSceNorm",
	
    def = function(theObject){
        standardGeneric("getSceNorm")    
    },
    signature = "theObject")



#' Getter of getSpecies slot of scRNAseq class
#'
#' Get the species assigne to the slot. This is the studied species.
#'
#' @name getSpecies
#' @rdname getSpecies-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#' 
#' @param theObject scRNAseq object
#'
#' @return The species assigned to the slot.
#' @export
setGeneric(
		
    name = "getSpecies",
	
    def = function(theObject){
        standardGeneric("getSpecies")    
    },
    signature = "theObject")




#' Getter of the outputDirectory slot of scRNAseq class
#'
#' Get the path of the output directory. This is the folder where all the 
#' outputs are.
#'
#' @name getOutputDirectory
#' @rdname getOutputDirectory-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' 
#' @return The path of the output directory
#' @export
setGeneric(
		
    name = "getOutputDirectory",
	
    def = function(theObject){
        standardGeneric("getOutputDirectory")    
    },
    signature = "theObject")


#' Getter of the PCs slot of scRNAseq class
#'
#' Get the vector of PCs. 
#'
#' @name getPCs
#' @rdname getPCs-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' 
#' @return PCs
#' @export
setGeneric(
		
    name = "getPCs",
	
    def = function(theObject){
        standardGeneric("getPCs")    
    },
    signature = "theObject")


#' Getter of the perplexities slot of scRNAseq class
#'
#' Get the vector of perplexities 
#'
#' @name getPerplexities
#' @rdname getPerplexities-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' 
#' @return Vector of perplexities
#' @export
setGeneric(
		
    name = "getPerplexities",
	
    def = function(theObject){
        standardGeneric("getPerplexities")    
    },
    signature = "theObject")


#' Getter of the tSNEList slot of scRNAseq class
#'
#' Get a list of tSNE results 
#'
#' @name getTSNEList
#' @rdname getTSNEList-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' 
#' @return List of tSNE results 
#' @export
setGeneric(
		
    name = "getTSNEList",
	
    def = function(theObject){
        standardGeneric("getTSNEList")    
    },
    signature = "theObject")


#' Getter of the dbscanList slot of scRNAseq class
#'
#' Get a list of DBSCAN results 
#'
#' @name getDbscanList
#' @rdname getDbscanList-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' 
#' @return List of DBSCAN results 
#' @export
setGeneric(
		
    name = "getDbscanList",
	
    def = function(theObject){
        standardGeneric("getDbscanList")    
    },
    signature = "theObject")


#' Getter of the cellsSimilarityMatrix slot of scRNAseq class
#'
#' Get a Cells Similarity Matrix
#'
#' @name getCellsSimilarityMatrix
#' @rdname getCellsSimilarityMatrix-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' 
#' @return Cells Similarity Matrix
#' @export
setGeneric(
		
    name = "getCellsSimilarityMatrix",
	
    def = function(theObject){
        standardGeneric("getCellsSimilarityMatrix")    
    },
    signature = "theObject")


#' Getter of the clusters similarity matrix slot of scRNAseq class
#'
#' Get clusters similarity matrix
#'
#' @name getClustersSimilarityMatrix
#' @rdname getClustersSimilarityMatrix-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' 
#' @return clusters similarity matrix
#' @export
setGeneric(
		
    name = "getClustersSimilarityMatrix",
	
    def = function(theObject){
        standardGeneric("getClustersSimilarityMatrix")    
    },
    signature = "theObject")


#' Getter of the clustersSimiliratyOrdered slot of scRNAseq class
#'
#' Get the order of the consensus clusters
#'
#' @name getClustersSimiliratyOrdered
#' @rdname getClustersSimiliratyOrdered-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' 
#' @return Factor. Order of the consensus clusters.
#' @export
setGeneric(
		
    name = "getClustersSimiliratyOrdered",
	
    def = function(theObject){
        standardGeneric("getClustersSimiliratyOrdered")    
    },
    signature = "theObject")



#' Getter of the markerGenesList slot of scRNAseq class
#'
#' Get the order of the consensus clusters
#'
#' @name getMarkerGenesList
#' @rdname getMarkerGenesList-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' 
#' @return markerGenesList
#' @export
setGeneric(
		
    name = "getMarkerGenesList",
	
    def = function(theObject){
        standardGeneric("getMarkerGenesList")    
    },
    signature = "theObject")



#' Getter of the clustersMarkers slot of scRNAseq class
#'
#' Get the clusters markers
#'
#' @name getClustersMarkers
#' @rdname getClustersMarkers-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' 
#' @return Clusters markers
#' @export
setGeneric(
		
		name = "getClustersMarkers",
		
		def=function(theObject){
			standardGeneric("getClustersMarkers")
		},
		signature = "theObject")


#' Getter of the getGenesInfos slot of scRNAseq class
#'
#' Get the clusters markers from the object
#'
#' @name getGenesInfos
#' @rdname getGenesInfos-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' 
#' @return Data.frame containing informations about markers genes.
#' @export
setGeneric(
    name="getGenesInfos",
    def=function(theObject){
        standardGeneric("getGenesInfos")
    },
    signature = "theObject")



## Tsne Class

#' Getter of the name slot of Tsne class
#'
#' Get the name from the Tsne object
#'
#' @name getTsneName
#' @rdname getTsneName-Tsne
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Tsne object
#' 
#' @return Name of the Tsne object
#' @export
setGeneric(
		
    name = "getTsneName",
	
    def = function(theObject){
        standardGeneric("getTsneName")    
    },
    signature = "theObject")


#' Getter of the PC slot of Tsne class
#'
#' Get the PC from the Tsne object
#'
#' @name getPC
#' @rdname getPC-Tsne
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Tsne object
#' 
#' @return PC of the Tsne object
#' @export
setGeneric(
		
    name = "getPC",
	
    def = function(theObject){
        standardGeneric("getPC")    
    },
    signature = "theObject")


#' Getter of the perplexity slot of Tsne class
#'
#' Get the perplexity from the Tsne object
#'
#' @name getPerplexity
#' @rdname getPerplexity-Tsne
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Tsne object
#' 
#' @return perplexity of the Tsne object
#' @export
setGeneric(
		
    name = "getPerplexity",
	
    def = function(theObject){
        standardGeneric("getPerplexity")    
    },
    signature = "theObject")


#' Getter of the coordinates slot of Tsne class
#'
#' Get the coordinates from the Tsne object
#'
#' @name getCoordinates
#' @rdname getCoordinates-Tsne
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Tsne object
#' 
#' @return coordinates of the Tsne object
#' @export
setGeneric(
		
    name = "getCoordinates",
	
    def = function(theObject){
        standardGeneric("getCoordinates")    
    },
    signature = "theObject")



## Dbscan getters

#' Getter of the name slot of Dbscan class
#'
#' Get the name from the Dbscan object
#'
#' @name getDbscanName
#' @rdname getDbscanName-Dbscan
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Dbscan object
#' 
#' @return Name of the Dbscan object
#' @export
setGeneric(
    
    name = "getDbscanName",
    
    def = function(theObject){
        standardGeneric("getDbscanName")    
    },
    signature = "theObject")


#' Getter of the epsilon slot of Dbscan class
#'
#' Get the name from the Dbscan object
#'
#' @name getEpsilon
#' @rdname getEpsilon-Dbscan
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Dbscan object
#' 
#' @return Name of the Dbscan object
#' @export
setGeneric(
		
    name = "getEpsilon",
	
    def = function(theObject){
        standardGeneric("getEpsilon")    
    },
    signature = "theObject")


#' Getter of the minPoints slot of Dbscan class
#'
#' Get the minPoints from the Dbscan object
#'
#' @name getMinPoints
#' @rdname getMinPoints-Dbscan
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Dbscan object
#' 
#' @return minPoints of the Dbscan object
#' @export
setGeneric(
		
    name = "getMinPoints",
	
    def = function(theObject){
        standardGeneric("getMinPoints")    
    },
    signature = "theObject")


#' Getter of the clustering slot of Dbscan class
#'
#' Get the clustering from the Dbscan object
#'
#' @name getClustering
#' @rdname getClustering-Dbscan
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Dbscan object
#' 
#' @return clustering result of the Dbscan object
#' @export
setGeneric(
		
    name = "getClustering",
	
    def = function(theObject){
        standardGeneric("getClustering")    
    },
    signature = "theObject")



################################################################################
############################### Setter methods #################################
################################################################################

#' Setter of experimentName slot of scRNAseq class
#'
#' Set the experiment name to the experimentName slot. Experiment name
#' describes the name of th experiment of the Conclus analysis.
#' 
#' @name setExperimentName<-
#' @rdname setExperimentName-scRNAseq
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' @param value Name of the experiment
#' @export
setGeneric(
		
    name = "setExperimentName<-",
	
    def = function(theObject, value){
        standardGeneric("setExperimentName<-")    
    },
    signature = "theObject")



#' Setter of countMatrix slot of scRNAseq class
#'
#' Set the the count matrix to the countMatrix slot.
#' 
#' @name setCountMatrix<-
#' @rdname setCountMatrix-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' @param value COunt matrix
#' @export
setGeneric(
		
    name = "setCountMatrix<-",
	
    def = function(theObject, value){
        standardGeneric("setCountMatrix<-")    
    },
    signature = "theObject")



#' Setter of sceNorm slot of scRNAseq class
#'
#' Set the SCE object  to the sceNorm slot. This contains the colData,
#' the rowData, and the normalized count matrix.
#'
#' @name setSceNorm<-
#' @rdname setSceNorm-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#' @param theObject scRNAseq object
#' @param value SingleCellExperiment object with normalized count matrix
#' @export
setGeneric(
		
    name = "setSceNorm<-",
	
    def = function(theObject, value){
        standardGeneric("setSceNorm<-")    
    },
    signature = "theObject")



#' Setter of setSpecies slot of scRNAseq class
#'
#' Set the species  to the slot. This is the studied species.
#'
#' @name setSpecies<-
#' @rdname setSpecies-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#' 
#' @param theObject scRNAseq object
#' @param value Species
#' @export
setGeneric(
		
    name = "setSpecies<-",
	
    def = function(theObject, value){
        standardGeneric("setSpecies<-")    
    },
    signature = "theObject")



#' Setter of the outputDirectory slot of scRNAseq class
#'
#' Set the path of the output directory. This is the folder where all the 
#' outputs are.
#'
#' @name setOutputDirectory<-
#' @rdname setOutputDirectory-scRNAseq
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' @param value Path of the output folder.
#' @export
setGeneric(
		
    name = "setOutputDirectory<-",
	
    def = function(theObject, value){
        standardGeneric("setOutputDirectory<-")    
    },
    signature = "theObject")


#' Setter of PCs slot of scRNAseq class
#'
#' Set the PCs to the PCs slot. a PC is the number of principal components
#' used by CONCLUS to perfom a PCA before to perform the tSNE and create the
#' object. By default, CONCLUS use a range of
#' PCs=c(8, 10, 20, 40, 50, 80, 100).
#' 
#' @name setPCs<-
#' @rdname setPCs-scRNAseq
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' @param value Vector of PCs
#' @export
setGeneric(
		
    name = "setPCs<-",
	
    def = function(theObject, value){
        standardGeneric("setPCs<-")    
    },
    signature = "theObject")



#' Setter of perplexities slot of scRNAseq class
#'
#' Set the perplexities to the perplexities slot.
#' The perplexity is a value used by tSNE algorithm.
#' 
#' @name setPerplexities<-
#' @rdname setPerplexities-scRNAseq
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' @param value Vector of perplexities
#' @export
setGeneric(
		
    name = "setPerplexities<-",
	
    def = function(theObject, value){
        standardGeneric("setPerplexities<-")    
    },
    signature = "theObject")



#' Setter of setTSNEList slot of scRNAseq class
#'
#' Set a list of tSNE solutions to the TSNEList slot.
#' 
#' @name setTSNEList<-
#' @rdname setTSNEList-scRNAseq
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' @param value List of tSNE results
#' @export
setGeneric(
		
    name = "setTSNEList<-",
	
    def = function(theObject, value){
        standardGeneric("setTSNEList<-")    
    },
    signature = "theObject")



#' Setter of dbscanList slot of scRNAseq class
#'
#' Set a list of DBSCAN solutions to the dbscanList slot.
#' 
#' @name setDbscanList<-
#' @rdname setDbscanList-scRNAseq
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' @param value List of DBSCAN results
#' @export
setGeneric(
		
    name = "setDbscanList<-",
	
    def = function(theObject, value){
        standardGeneric("setDbscanList<-")    
    },
    signature = "theObject")



#' Setter of cellsSimilarityMatrix slot of scRNAseq class
#'
#' Set a value to the cellsSimilarityMatrix slot.
#' 
#' @name setCellsSimilarityMatrix<-
#' @rdname setCellsSimilarityMatrix-scRNAseq
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' @param value cells Similarity Matrix
#' @export
setGeneric(
		
    name = "setCellsSimilarityMatrix<-",
	
    def = function(theObject, value){
        standardGeneric("setCellsSimilarityMatrix<-")    
    },
    signature = "theObject")



#' Setter of clustersSimilarityMatrix slot of scRNAseq class
#'
#' Set a value to the clustersSimilarityMatrix slot.
#' 
#' @name setClustersSimilarityMatrix<-
#' @rdname setClustersSimilarityMatrix-scRNAseq
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' @param value clusters Similarity Matrix
#' @export
setGeneric(
		
    name = "setClustersSimilarityMatrix<-",
	
    def = function(theObject, value){
        standardGeneric("setClustersSimilarityMatrix<-")    
    },
    signature = "theObject")


#' Setter of clustersSimiliratyOrdered slot of scRNAseq class
#'
#' Set a value to the clustersSimiliratyOrdered slot.
#' 
#' @name setClustersSimiliratyOrdered<-
#' @rdname setClustersSimiliratyOrdered-scRNAseq
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' @param value Factor. Label of clusters similiraty ordered.
#' @export
setGeneric(

    name = "setClustersSimiliratyOrdered<-",

    def = function(theObject, value){
        standardGeneric("setClustersSimiliratyOrdered<-")    
    },
    signature = "theObject")



#' Setter of markerGenesList slot of scRNAseq class
#'
#' Set a value to the markerGenesList slot.
#' 
#' @name setMarkerGenesList<-
#' @rdname setMarkerGenesList-scRNAseq
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' @param value List of markers genes (data.frames) for each cluster.
#' @export
setGeneric(

    name = "setMarkerGenesList<-",

    def = function(theObject, value){
        standardGeneric("setMarkerGenesList<-")    
    },
    signature = "theObject")



#' Setter of clustersMarkers slot of scRNAseq class
#'
#' Set a value to the clustersMarkers slot.
#' 
#' @name setClustersMarkers<-
#' @rdname setClustersMarkers-scRNAseq
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' @param value Data.frame of best markers genes among all clusters.
#' @export
setGeneric(
    name="setClustersMarkers<-",
    def=function(theObject, value){
        standardGeneric("setClustersMarkers<-")
    },
    signature = "theObject")



#' Setter of genesInfos slot of scRNAseq class
#'
#' Set a value to the genesInfos slot.
#' 
#' @name setGenesInfos<-
#' @rdname setGenesInfos-scRNAseq
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject scRNAseq object
#' @param value Data.frame containing some informations about markers genes.
#' @export
setGeneric(
    name="setGenesInfos<-",
    def=function(theObject, value){
        standardGeneric("setGenesInfos<-")
    },
    signature = "theObject")





##############################  Tsne class setters #############################


#' Setter of name slot of Tsne class
#'
#' Set a name to the name slot.
#' 
#' @name setTsneName<-
#' @rdname setTsneName-Tsne
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Tsne object
#' @param value Name
#' @export
setGeneric(

    name = "setTsneName<-",

    def = function(theObject, value){
        standardGeneric("setTsneName<-")    
    },
    signature = "theObject")



#' Setter of PC slot of Tsne class
#'
#' Set a PC to the PC slot.
#' 
#' @name setPC<-
#' @rdname setPC-Tsne
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Tsne object
#' @param value PC
#' @export
setGeneric(

    name = "setPC<-",

    def = function(theObject, value){
        standardGeneric("setPC<-")    
    },
    signature = "theObject")



#' Setter of perplexity slot of Tsne class
#'
#' Set a perplexity to the perplexity slot.
#' 
#' @name setPerplexity<-
#' @rdname setPerplexity-Tsne
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Tsne object
#' @param value perplexity
#' @export
setGeneric(

    name = "setPerplexity<-",

    def = function(theObject, value){
        standardGeneric("setPerplexity<-")    
    },
    signature = "theObject")



#' Setter of perplexity slot of Tsne class
#'
#' Set a perplexity to the perplexity slot.
#' 
#' @name setCoordinates<-
#' @rdname setCoordinates-Tsne
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Tsne object
#' @param value Coordinates (data.frame)
#' @export
setGeneric(

    name = "setCoordinates<-",

    def = function(theObject, value){
        standardGeneric("setCoordinates<-")    
    },
    signature = "theObject")




#############################  Dbscan class setters ############################

#' Setter of name slot of Dbscan class
#'
#' Set a name to the name slot.
#' 
#' @name setDbscanName<-
#' @rdname setDbscanName-Dbscan
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Dbscan object
#' @param value Name
#' @export
setGeneric(

    name = "setDbscanName<-",

    def = function(theObject, value){
        standardGeneric("setDbscanName<-")    
    },
    signature = "theObject")



#' Setter of epsilon slot of Dbscan class
#'
#' Set a epsilon to the epsilon slot.
#' 
#' @name setEpsilon<-
#' @rdname setEpsilon-Dbscan
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Dbscan object
#' @param value epsilon
#' @export
setGeneric(

    name = "setEpsilon<-",

    def = function(theObject, value){
        standardGeneric("setEpsilon<-")    
    },
    signature = "theObject")



#' Setter of minPoints slot of Dbscan class
#'
#' Set a minPoints to the Dbscan slot.
#' 
#' @name setMinPoints<-
#' @rdname setMinPoints-Dbscan
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Dbscan object
#' @param value epsilon
#' @export
setGeneric(

    name = "setMinPoints<-",

    def = function(theObject, value){
        standardGeneric("setMinPoints<-")    
    },
    signature = "theObject")



#' Setter of clustering slot of Dbscan class
#'
#' Set a clustering to the clustering slot.
#' 
#' @name setClustering<-
#' @rdname setClustering-Dbscan
#' @importFrom methods setReplaceMethod validObject
#' @author Ilyess RACHEDI and Nicolas DESCOSTES, based on code by Polina
#' PAVLOVICH and Artem ADAMOV.
#'  
#' @param theObject Dbscan object
#' @param value clustering
#' @export
setGeneric(
    
    name = "setClustering<-",
    
    def = function(theObject, value){
        standardGeneric("setClustering<-")    
    },
    signature = "theObject")
