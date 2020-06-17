################################################################################
############################### Tsne class  ####################################
################################################################################


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
		
		invisible(checkGenesInfos(genesInfos, species))
    }
)

