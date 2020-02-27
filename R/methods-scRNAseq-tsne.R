.initialisePath <- function(dataDirectory){
    # creates directories for further writing of results.
    # names of directories are hardcoded.
    # no idea if it is good or bad.
    
    graphsDirectory <- "pictures"
    markerGenesDirectory <- "marker_genes"
    tSNEDirectory <- "tsnes"
    outputDataDirectory <- "output_tables"
    tSNEPicturesDirectory <- "tSNE_pictures"
    
    
    dir.create(dataDirectory, showWarnings=F)
    dir.create(file.path(dataDirectory, graphsDirectory), showWarnings=F)
    dir.create(file.path(dataDirectory, graphsDirectory, tSNEPicturesDirectory),
               showWarnings=F)
    dir.create(file.path(dataDirectory, markerGenesDirectory), showWarnings=F)
    dir.create(file.path(dataDirectory, tSNEDirectory), showWarnings=F)
    dir.create(file.path(dataDirectory, outputDataDirectory), showWarnings=F)
    
}


.getTSNEresults <- function(expressionMatrix = getNormalizedCountMatrix(theObject),
                            cores=1,
                            PCs=c(4, 6, 8, 10, 20, 40, 50),
                            perplexities=c(30, 40),
                            randomSeed=42){
    PCAData <- prcomp(t(expressionMatrix))$x
    myCluster <- parallel::makeCluster(cores, # number of cores to use
                                       type = "PSOCK") # type of cluster
    doParallel::registerDoParallel(myCluster)
    tSNECoordinates <- foreach::foreach(PCA=rep(PCs, length(perplexities)),
                                        perp=rep(perplexities,
                                                 each=length(PCs)),
                                        .combine='cbind') %dopar% {
                                            library(SingleCellExperiment)
                                            tmp <- scater::runTSNE(
                                                SingleCellExperiment(
                                                    assays=list(
                                                    logcounts=t(
                                                        PCAData[, 1:PCA]))),
                                                scale_features=FALSE,
                                                perplexity=perp,
                                                rand_seed=randomSeed,
                                                theme_size=13, 
                                                return_SCESet=FALSE)
                                            
                                            scater::plotTSNE(tmp)
                                        }
    parallel::stopCluster(myCluster)
    message(paste("Calculated", length(PCs)*length(perplexities),
                  "2D-tSNE plots."))
    return(tSNECoordinates)
}


setMethod(
    f="generateTSNECoordinates",
    signature="scRNAseq",
    definition=function(theObject,
                        sceObject=getNormalizedCountMatrix(theObject),
                        dataDirectory=getOutputDirectory(theObject),
                        experimentName=getExperimentName(theObject),
                        randomSeed=42,
                        cores=1,
                        PCs=c(4, 6, 8, 10, 20, 40, 50),
                        perplexities=c(30,40)){
        
        .initialisePath(dataDirectory)
        tSNEDirectory <- "tsnes"
        message(paste0("Running TSNEs using ", cores, " cores."))
        TSNEres <- .getTSNEresults(Biobase::exprs(sceObject),
                                   cores=cores,
                                   PCs=PCs,
                                   perplexities=perplexities,
                                   randomSeed=randomSeed)
        
        PCA <- rep(PCs, length(perplexities))
        perp <- rep(perplexities, each=length(PCs))
        
        outputDir <- file.path(dataDirectory, tSNEDirectory)
        filesList <- list.files(outputDir, pattern="_tsne_coordinates_")
        deletedFiles <- sapply(filesList, function(fileName) 
            file.remove(file.path(outputDir, fileName)))
        
        for (i in 1:(length(PCs)*length(perplexities))){
            write.table(TSNEres[1, i][[1]],
                        file=file.path(dataDirectory, tSNEDirectory,
                                       paste0(experimentName,
                                              '_tsne_coordinates_',
                                              i,
                                              "_" ,
                                              PCA[i],
                                              "PCs_",
                                              perp[i],
                                              "perp.tsv")),
                        quote=FALSE, sep='\t')
        }
        saveRDS(TSNEres, file=file.path(dataDirectory,
                                          "output_tables",
                                          paste0(experimentName,
                                                 "_tSNEResults.rds")))
        rm(tSNEDirectory, PCA, perp)
        setTSNEResults(theObject) <- TSNEres
        return(theObject)
    })
