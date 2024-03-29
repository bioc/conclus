################################################################################
########################### scRNAseq class methods #############################
################################################################################

setGeneric(

    name = "normaliseCountMatrix",

    def = function(theObject, sizes=c(20,40,60,80,100), rowdata=NULL,
            coldata=NULL, alreadyCellFiltered=FALSE, runQuickCluster=TRUE,
            info=TRUE, removeNoSymbol=FALSE){
            standardGeneric("normaliseCountMatrix")
    },
    signature="theObject")



setGeneric(

    name = "testClustering",

    def = function(theObject, dbscanEpsilon=1.4, minPts=5,
            perplexities=30, PCs=4, randomSeed=42, width=7, height=7,
            cores=2, writeOutput=FALSE, fileTSNE="test_tSNE.pdf",
            fileDist="distance_graph.pdf", fileClust="test_clustering.pdf",
            silent=FALSE, plotKNN=TRUE, ...){
            standardGenericic("testClustering")
    },
    signature="theObject")


setGeneric(

    name = "generateTSNECoordinates",

    def = function(theObject, randomSeed=42, cores=2,
            PCs=c(4, 6, 8, 10, 20, 40, 50), perplexities=c(30,40),
            writeOutput=FALSE){
        standardGeneric("generateTSNECoordinates")
    },
    signature="theObject")



setGeneric(

    name = "runDBSCAN",

    def = function(theObject, cores=2, epsilon=c(1.3, 1.4, 1.5),
            minPoints=c(3, 4), writeOutput=FALSE){
        standardGeneric("runDBSCAN")
    },
    signature="theObject")



setGeneric(

    name = "clusterCellsInternal",

    def = function(theObject, clusterNumber=NULL, deepSplit=4, cores=2,
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

    name = "addClustering",

    def = function(theObject, filePathAdd=NA, clusToAdd=NA){
        standardGeneric("addClustering")
    },
    signature="theObject")



setGeneric(

    name = "plotCellSimilarity",

    def = function(theObject, colorPalette="default", statePalette="default",
            clusteringMethod="ward.D2", orderClusters=FALSE, savePlot=FALSE,
            plotPDF=TRUE, returnPlot=FALSE, width=7, height=6, onefile=FALSE,
            showRowNames=FALSE, showColnames=FALSE, fontsize=7.5,
            fontsizeRow=0.03, widthPNG=800, heightPNG=750, silentPlot=FALSE){
        standardGeneric("plotCellSimilarity")
    },
    signature="theObject")



setGeneric(

    name = "plotClusteredTSNE",

    def = function(theObject, colorPalette="default",
            PCs=c(4, 6, 8, 10, 20, 40, 50), perplexities=c(30, 40),
            columnName="clusters", savePlot=FALSE, plotPDF=TRUE,
            returnPlot=FALSE, width=6, height=5, onefile=FALSE, widthPNG=800,
            heightPNG=750, silentPlot=FALSE, tSNENb=NA){
        standardGeneric("plotClusteredTSNE")

    },
    signature="theObject")



setGeneric(

    name = "plotCellHeatmap",

    def = function(theObject, fileName = NA, meanCentered=TRUE,
            colorPalette="default", statePalette="default",
            clusteringMethod="ward.D2", orderClusters=FALSE, orderGenes=FALSE,
            returnPlot=FALSE, savePlot=FALSE, width=10, height=8.5,
            onefile=FALSE, clusterCols=FALSE, showColnames=FALSE, fontsize=7.5,
            fontsizeRow=8, plotPDF=TRUE, widthPNG=800, heightPNG=750,
            silentPlot=FALSE){
        standardGeneric("plotCellHeatmap")
    },
    signature="theObject")



setGeneric(

    name = "plotGeneExpression",

    def = function(theObject, geneName,
            palette=c("grey","red", "#7a0f09", "black"), returnPlot=FALSE,
            tSNEpicture=1, savePlot=FALSE, alpha=1, limits=NA,
            pointSize=1, width=6, height=5, plotPDF=TRUE, silentPlot=FALSE){
        standardGeneric("plotGeneExpression")
    },
    signature="theObject")



setGeneric(

    name = "plotClustersSimilarity",

    def = function(theObject, colorPalette="default", statePalette="default",
            clusteringMethod="ward.D2", returnPlot=FALSE, savePlot=FALSE,
            plotPDF=TRUE, width=7, height=5.5, onefile=FALSE, fontsize=7.5,
            widthPNG=800, heightPNG=750, silentPlot=FALSE){
        standardGeneric("plotClustersSimilarity")
    },
    signature="theObject")



setGeneric(

    name = "exportResults",

    def = function(theObject, saveClusteringResults=TRUE, saveAll=FALSE,
            saveNormalizedMatrix=FALSE, saveColData=FALSE, saveRowData=FALSE,
            saveTsne=FALSE, saveDBScan=FALSE, saveCellsSimilarityMatrix=FALSE,
            saveClustersSimilarityMatrix=FALSE, saveFullMarkers=FALSE,
            saveTopMarkers=FALSE, saveGenesInfos=FALSE){
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

    def = function(theObject, groupBy="clusters", orderGenes="initial",
            getUniprot=TRUE, cores=2, saveInfos=FALSE){
        standardGeneric("retrieveGenesInfo")
    },
    signature="theObject")


setGeneric(

    name = "retrieveTableClustersCells",

    def= function(theObject){
        standardGeneric("retrieveTableClustersCells")
    },
    signature="theObject")


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

    name = "getSuggestedClustersNumber",

    def = function(theObject){
        standardGeneric("getSuggestedClustersNumber")
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

    name = "getClustersSimilarityOrdered",

    def = function(theObject){
        standardGeneric("getClustersSimilarityOrdered")
    },
    signature = "theObject")



setGeneric(

    name = "getMarkerGenesList",

    def = function(theObject, cluster="all"){
        standardGeneric("getMarkerGenesList")
    },
    signature = "theObject")



setGeneric(

        name = "getTopMarkers",

        def=function(theObject){
            standardGeneric("getTopMarkers")
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
    signature = c("theObject"))


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

    name = "setSuggestedClustersNumber<-",

    def = function(theObject, value){
        standardGeneric("setSuggestedClustersNumber<-")
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

    name = "setTopMarkers<-",

    def = function(theObject, value){
        standardGeneric("setTopMarkers<-")
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
    signature = c("theObject"))



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
