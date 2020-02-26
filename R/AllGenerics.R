################################################################################
########################### scRNAseq class methods #############################
################################################################################

setGeneric(
    name="normaliseCountMatrix",
    def=function(theObject,
                 countMatrix=getCountMatrix(theObject),
                 species=getSpecies(theObject),
                 sizes=c(20,40,60,80,100),
                 rowData=NULL,
                 colData=NULL,
                 alreadyCellFiltered=FALSE, 
                 runQuickCluster=TRUE,
                 databaseDir=TRUE){
        standardGeneric("normaliseCountMatrix")    
    },
    signature="theObject")


setGeneric(
    name="generateTSNECoordinates",
    def=function(theObject,
                 sceObject=getNormalizedCountMatrix(theObject),
                 dataDirectory=getOutputDirectory(theObject),
                 experimentName=getExperimentName(theObject),
                 randomSeed=42,
                 cores=1,
                 PCs=c(4, 6, 8, 10, 20, 40, 50),
                 perplexities=c(30,40)){
        standardGeneric("generateTSNECoordinates")    
    },
    signature="theObject")


################################################################################
############################### Getter methods #################################
################################################################################

setGeneric(
    name="getExperimentName",
    def=function(theObject){
        standardGeneric("getExperimentName")    
    },
    signature="theObject")


setGeneric(
    name="getCountMatrix",
    def=function(theObject){
        standardGeneric("getCountMatrix")    
    },
    signature = "theObject")


setGeneric(
    name="getNormalizedCountMatrix",
    def=function(theObject){
        standardGeneric("getNormalizedCountMatrix")    
    },
    signature = "theObject")


setGeneric(
    name="getColData",
    def=function(theObject){
        standardGeneric("getColData")    
    },
    signature = "theObject")


setGeneric(
    name="getSpecies",
    def=function(theObject){
        standardGeneric("getSpecies")    
    },
    signature = "theObject")


setGeneric(
    name="getOutputDirectory",
    def=function(theObject){
        standardGeneric("getOutputDirectory")    
    },
    signature = "theObject")


setGeneric(
    name="getPCs",
    def=function(theObject){
        standardGeneric("getPCs")    
    },
    signature = "theObject")


setGeneric(
    name="getPerplexities",
    def=function(theObject){
        standardGeneric("getPerplexities")    
    },
    signature = "theObject")


setGeneric(
    name="getTSNEResults",
    def=function(theObject){
        standardGeneric("getTSNEResults")    
    },
    signature = "theObject")


setGeneric(
    name="getClusteringResults",
    def=function(theObject){
        standardGeneric("getClusteringResults")    
    },
    signature = "theObject")


################################################################################
############################### Setter methods #################################
################################################################################

setGeneric(
    name="setExperimentName<-",
    def=function(theObject, value){
        standardGeneric("setExperimentName<-")    
    },
    signature = "theObject")


setGeneric(
    name="setCountMatrix<-",
    def=function(theObject, value){
        standardGeneric("setCountMatrix<-")    
    },
    signature = "theObject")


setGeneric(
    name="setNormalizedCountMatrix<-",
    def=function(theObject, value){
        standardGeneric("setNormalizedCountMatrix<-")    
    },
    signature = "theObject")


setGeneric(
    name="setColData<-",
    def=function(theObject, value){
        standardGeneric("setColData<-")    
    },
    signature = "theObject")


setGeneric(
    name="setSpecies<-",
    def=function(theObject, value){
        standardGeneric("setSpecies<-")    
    },
    signature = "theObject")


setGeneric(
    name="setOutputDirectory<-",
    def=function(theObject, value){
        standardGeneric("setOutputDirectory<-")    
    },
    signature = "theObject")


setGeneric(
    name="setPCs<-",
    def=function(theObject, value){
        standardGeneric("setPCs<-")    
    },
    signature = "theObject")


setGeneric(
    name="setPerplexities<-",
    def=function(theObject, value){
        standardGeneric("setPerplexities<-")    
    },
    signature = "theObject")


setGeneric(
    name="setTSNEResults<-",
    def=function(theObject, value){
        standardGeneric("setTSNEResults<-")    
    },
    signature = "theObject")


setGeneric(
    name="setClusteringResults<-",
    def=function(theObject, value){
        standardGeneric("setClusteringResults<-")    
    },
    signature = "theObject")


