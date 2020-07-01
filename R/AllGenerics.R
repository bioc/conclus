################################################################################
########################### scRNAseq class methods #############################
################################################################################

setGeneric(
		
	name = "normaliseCountMatrix",
		
	def = function(theObject, sizes=c(20,40,60,80,100), rowdata=NULL,
			coldata=NULL, alreadyCellFiltered=FALSE, runQuickCluster=TRUE){
			standardGeneric("normaliseCountMatrix")    
    },
    signature="theObject") 



setGeneric(
		
	name = "testClustering",
		
	def = function(theObject, dbscanEpsilon=1.4, minPts=5, perplexities=c(30), 
			PCs=c(4), randomSeed=42, width=7, height=7, onefile=FALSE, ...){
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
		
	def = function(theObject, colorPalette="default", statePalette="default", 
			clusteringMethod="ward.D2", orderClusters=FALSE, plotPDF=TRUE, 
			returnPlot=FALSE, width=7, height=6, ...){
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
		
	def = function(theObject, colorPalette="default", statePalette="default", 
			clusteringMethod="ward.D2", returnPlot=FALSE, width=7, height=5.5, 
			...){
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
		
    name =  "retrieveTopClustersMarkers",
	
    def = function(theObject, nTop=10, removeDuplicates = TRUE, 
			writeMarkerGenes = FALSE){
        standardGeneric("retrieveTopClustersMarkers")
    },
    signature="theObject")



setGeneric(
		
    name = "retrieveGenesInfo",
	
    def = function(theObject, species, groupBy="clusters", orderGenes="initial",
			getUniprot=TRUE, cores=1){
        standardGeneric("retrieveGenesInfo")
    },
    signature="theObject")



setGeneric(
		
    name = "saveGenesInfo",
	
    def = function(theObject) {
        standardGeneric("saveGenesInfo")
    },
    signature = "theObject"
)


################################################################################
############################### Getter methods #################################
################################################################################

setGeneric(
		
	name = "getExperimentName",
	
    def = function(theObject){
        standardGeneric("getExperimentName")    
    },
    signature="theObject")


setGeneric(
		
    name = "getCountMatrix",
	
    def = function(theObject){
        standardGeneric("getCountMatrix")    
    },
    signature = "theObject")



setGeneric(
		
    name = "getSceNorm",
	
    def = function(theObject){
        standardGeneric("getSceNorm")    
    },
    signature = "theObject")



setGeneric(
		
    name = "getSpecies",
	
    def = function(theObject){
        standardGeneric("getSpecies")    
    },
    signature = "theObject")



setGeneric(
		
    name = "getOutputDirectory",
	
    def = function(theObject){
        standardGeneric("getOutputDirectory")    
    },
    signature = "theObject")


setGeneric(
		
    name = "getPCs",
	
    def = function(theObject){
        standardGeneric("getPCs")    
    },
    signature = "theObject")


setGeneric(
		
    name = "getPerplexities",
	
    def = function(theObject){
        standardGeneric("getPerplexities")    
    },
    signature = "theObject")


setGeneric(
		
    name = "getTSNEList",
	
    def = function(theObject){
        standardGeneric("getTSNEList")    
    },
    signature = "theObject")


setGeneric(
		
    name = "getDbscanList",
	
    def = function(theObject){
        standardGeneric("getDbscanList")    
    },
    signature = "theObject")


setGeneric(
		
    name = "getCellsSimilarityMatrix",
	
    def = function(theObject){
        standardGeneric("getCellsSimilarityMatrix")    
    },
    signature = "theObject")


setGeneric(
		
    name = "getClustersSimilarityMatrix",
	
    def = function(theObject){
        standardGeneric("getClustersSimilarityMatrix")    
    },
    signature = "theObject")


setGeneric(
		
    name = "getClustersSimiliratyOrdered",
	
    def = function(theObject){
        standardGeneric("getClustersSimiliratyOrdered")    
    },
    signature = "theObject")



setGeneric(
		
    name = "getMarkerGenesList",
	
    def = function(theObject){
        standardGeneric("getMarkerGenesList")    
    },
    signature = "theObject")



setGeneric(
		
		name = "getClustersMarkers",
		
		def=function(theObject){
			standardGeneric("getClustersMarkers")
		},
		signature = "theObject")


setGeneric(
    name="getGenesInfos",
    def=function(theObject){
        standardGeneric("getGenesInfos")
    },
    signature = "theObject")



## Tsne Class

setGeneric(
		
    name = "getName",
	
    def = function(theObject){
        standardGeneric("getName")    
    },
    signature = "theObject")


setGeneric(
		
    name = "getPC",
	
    def = function(theObject){
        standardGeneric("getPC")    
    },
    signature = "theObject")


setGeneric(
		
    name = "getPerplexity",
	
    def = function(theObject){
        standardGeneric("getPerplexity")    
    },
    signature = "theObject")


setGeneric(
		
    name = "getCoordinates",
	
    def = function(theObject){
        standardGeneric("getCoordinates")    
    },
    signature = "theObject")



## Dbscan getters

setGeneric(
		
    name = "getEpsilon",
	
    def = function(theObject){
        standardGeneric("getEpsilon")    
    },
    signature = "theObject")


setGeneric(
		
    name = "getMinPoints",
	
    def = function(theObject){
        standardGeneric("getMinPoints")    
    },
    signature = "theObject")


setGeneric(
		
    name = "getClustering",
	
    def = function(theObject){
        standardGeneric("getClustering")    
    },
    signature = "theObject")



################################################################################
############################### Setter methods #################################
################################################################################

setGeneric(
		
    name = "setExperimentName<-",
	
    def = function(theObject, value){
        standardGeneric("setExperimentName<-")    
    },
    signature = "theObject")



setGeneric(
		
    name = "setCountMatrix<-",
	
    def = function(theObject, value){
        standardGeneric("setCountMatrix<-")    
    },
    signature = "theObject")



setGeneric(
		
    name = "setSceNorm<-",
	
    def = function(theObject, value){
        standardGeneric("setSceNorm<-")    
    },
    signature = "theObject")



setGeneric(
		
    name = "setSpecies<-",
	
    def = function(theObject, value){
        standardGeneric("setSpecies<-")    
    },
    signature = "theObject")



setGeneric(
		
    name = "setOutputDirectory<-",
	
    def = function(theObject, value){
        standardGeneric("setOutputDirectory<-")    
    },
    signature = "theObject")


setGeneric(
		
    name = "setPCs<-",
	
    def = function(theObject, value){
        standardGeneric("setPCs<-")    
    },
    signature = "theObject")



setGeneric(
		
    name = "setPerplexities<-",
	
    def = function(theObject, value){
        standardGeneric("setPerplexities<-")    
    },
    signature = "theObject")



setGeneric(
		
    name = "setTSNEList<-",
	
    def = function(theObject, value){
        standardGeneric("setTSNEList<-")    
    },
    signature = "theObject")



setGeneric(
		
    name = "setDbscanList<-",
	
    def = function(theObject, value){
        standardGeneric("setDbscanList<-")    
    },
    signature = "theObject")



setGeneric(
		
    name = "setCellsSimilarityMatrix<-",
	
    def = function(theObject, value){
        standardGeneric("setCellsSimilarityMatrix<-")    
    },
    signature = "theObject")



setGeneric(
		
    name = "setClustersSimilarityMatrix<-",
	
    def = function(theObject, value){
        standardGeneric("setClustersSimilarityMatrix<-")    
    },
    signature = "theObject")


setGeneric(

    name = "setClustersSimiliratyOrdered<-",

    def = function(theObject, value){
        standardGeneric("setClustersSimiliratyOrdered<-")    
    },
    signature = "theObject")



setGeneric(

    name = "setMarkerGenesList<-",

    def = function(theObject, value){
        standardGeneric("setMarkerGenesList<-")    
    },
    signature = "theObject")



setGeneric(

    name = "setClustersMarkers<-",
    
    def = function(theObject, value){
        standardGeneric("setClustersMarkers<-")
    },
    signature = "theObject")



setGeneric(

    name="setGenesInfos<-",
    
    def=function(theObject, value){
        standardGeneric("setGenesInfos<-")
    },
    signature = "theObject")





##############################  Tsne class setters #############################


setGeneric(

    name = "setName<-",

    def = function(theObject, value){
        standardGeneric("setName<-")    
    },
    signature = "theObject")



setGeneric(

    name = "setPC<-",

    def = function(theObject, value){
        standardGeneric("setPC<-")    
    },
    signature = "theObject")



setGeneric(

    name = "setPerplexity<-",

    def = function(theObject, value){
        standardGeneric("setPerplexity<-")    
    },
    signature = "theObject")



setGeneric(

    name = "setCoordinates<-",

    def = function(theObject, value){
        standardGeneric("setCoordinates<-")    
    },
    signature = "theObject")




#############################  Dbscan class setters ############################

setGeneric(

    name = "setEpsilon<-",

    def = function(theObject, value){
        standardGeneric("setEpsilon<-")    
    },
    signature = "theObject")



setGeneric(

    name = "setMinPoints<-",

    def = function(theObject, value){
        standardGeneric("setMinPoints<-")    
    },
    signature = "theObject")



setGeneric(
    
    name = "setClustering<-",
    
    def = function(theObject, value){
        standardGeneric("setClustering<-")    
    },
    signature = "theObject")
