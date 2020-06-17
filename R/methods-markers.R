#################
## rankGenes
#################

.checkParamsRankGenes <- function(sceObject, simMed){
	
	## Check if the normalized matrix is correct
	if(all(dim(sceObject) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with 'rankGenes' ",
				"function doesn't have its 'SceNorm' slot updated. Please",
				" use 'normaliseCountMatrix' on the object before.")
	
	## Check if the SCE object contain cluster colums in its colData
	if(!("clusters" %in% names(colData(sceObject))))
		stop("The 'scRNAseq' object that you're using with 'rankGenes' ",
				"function doesn't have a correct 'SceNorm' slot. This slot",
				" should be a 'SingleCellExperiment' object containing ",
				"'clusters' column in its colData. Please check if you ",
				"correctly used 'clusterCellsInternal' on the object.")
	
	## Check the cluster similarity matrix
	if(all(dim(simMed) == c(0,0)) || all(dim(simMed) == c(1,1)))
		stop("The 'scRNAseq' object that you're using with 'rankGenes' ",
				"function doesn't have its 'clustersSimilarityMatrix' slot",
				" updated. Please use 'clusterCellsInternal' on the object",
				" before.")
}

.buildTTestPval <- function(otherGroups, tTestPval, colDF, mat, colLabel, 
		currentGroup){
	
	tTestPval <- lapply(otherGroups, function(currentOther){
				
				tTestPval[, paste0("vs_", currentOther)] <- NA
				x <- mat[, colDF[, c(colLabel)] == currentGroup]
				y <- mat[, colDF[, c(colLabel)] == currentOther]
				t <- (rowMeans(x)-rowMeans(y))/sqrt(apply(mat, 1, 
								var)*(1 / ncol(x) + 1 / ncol(y)))
				df <- ncol(x)+ncol(y)-2
				tTestPval[, paste0("vs_", currentOther)] <- pt(t, 
						df, lower.tail=FALSE)
				return(tTestPval)
			})
	
	return(tTestPval)
}

.buildMarkerGenesList <- function(groups, simMedRowList, exprM, column, 
		colNamesVec, colDF, writeMarkerGenes){
	
	result <- mapply(function(currentGroup, currentSimMed, mat, 
					allGroups, colLabel, simMedNames){
				
				message(paste("Working on cluster", currentGroup))
				tTestPval <- data.frame(row.names=rownames(mat))
				tTestFDR <- data.frame(row.names=rownames(mat))
				otherGroups <- allGroups[allGroups!=currentGroup]
				
				## Create a dataframe clustering vs clustering
				tTestPval <- .buildTTestPval(otherGroups, tTestPval, colDF, mat,
						colLabel, currentGroup) 
				tTestPval <- bind_cols(tTestPval)
				tTestPval <- cbind(Gene=rownames(mat), tTestPval)
				
				## Apply the FDR
				tTestFDR <- lapply(otherGroups, function(currentOther){
							
							tTestFDR[, paste0("vs_", currentOther)] <- 
									p.adjust(tTestPval[, paste0("vs_", 
															currentOther)], 
											method="fdr")
							return(tTestFDR)
						})
				tTestFDR <- bind_cols(tTestFDR)
				tTestFDR <- cbind(Gene=as.factor(rownames(mat)), tTestFDR)
				
				submat <- as.matrix(tTestFDR[, 2:(length(otherGroups) + 1 )])
				
				## Add column mean_log10_fdr 
				tTestFDR$mean_log10_fdr <- rowMeans(log10(submat + 1e-300))
				
				## Add column n_05
				tTestFDR$n_05 <- apply(submat, 1, function(x) 
							length(x[!is.na(x) & x < 0.05]))
				
				##Add column score 
				weights <- currentSimMed[match(otherGroups, simMedNames)]
				tTestFDR$score <- apply(submat, 1, function(x) 
							sum(-log10(x+1e-300) * weights) / ncol(submat))
				
				tTestFDR <- tTestFDR[order(tTestFDR$score, decreasing=TRUE), ]
				row.names(tTestFDR) <- NULL
				
				## Write list if option = TRUE
				if(writeMarkerGenes){
					
					outputDir <- file.path(getOutputDirectory(theObject), 
							"marker_genes")
					outputFile <- paste0(getExperimentName(theObject),
							"_cluster_", allGroups[i], "_genes.tsv")
					
					write.table(tTestFDR, file.path(outputDir, outputFile),
							col.names=TRUE, row.names=FALSE, quote=FALSE,
							sep="\t")
				}
				
				return(tTestFDR)
				
			}, groups, simMedRowList, 
			MoreArgs = list(exprM, groups, column, colNamesVec), 
			SIMPLIFY = FALSE)
	
	return(result)
}

setMethod(
		
    f = "rankGenes",
	
    signature = "scRNAseq",
	
    definition = function(theObject, column="clusters", writeMarkerGenes=FALSE){
        
        sceObject  <- getSceNorm(theObject)
		simMed <- getClustersSimilarityMatrix(theObject)
		.checkParamsRankGenes(sceObject, simMed)
        exprM <- Biobase::exprs(sceObject)
        colDF <- SummarizedExperiment::colData(sceObject)
         
        message("Ranking marker genes for each cluster.")
        
        stopifnot(all(colnames(exprM) == rownames(colDF)))
        groups <- unique(colDF[, column])
        simMed <- simMed + 0.05
		simMedRowList <- split(simMed, seq_len(nrow(simMed)))
		colNamesVec <- as.numeric(colnames(simMed))
        		
		markerGenesList <- .buildMarkerGenesList(groups, simMedRowList, exprM,
				column, colNamesVec, colDF, writeMarkerGenes) 		
		
        setMarkerGenesList(theObject) <- markerGenesList
        rm(markerGenesList, simMed, groups)
        
        return(theObject)
    })


###################
## retrieveGenesInfo
###################


.checkParamsretrieveGenesInfo <- function(orderGenes, genes, species,
		databaseDict){
	
	if(!isTRUE(all.equal(orderGenes,"initial")) && 
			!isTRUE(all.equal(orderGenes,"alphabetical")))
		stop("orderGenes should be 'initial' or 'alphabetical'.")
	
	if (!(species %in% names(databaseDict)))
		stop("species should be: ", paste(names(databaseDict), 
						collapse = " or "))
}


.returnDB1 <- function(genes, ensembl){
	
	db1 <- getBM(attributes=c("uniprot_gn_symbol",  # geneName
					"external_gene_name", # Complete name,
					"description",        # Description
					"entrezgene_description",
					"gene_biotype",        # Feature.Type
					"go_id" ),             # Gene Ontology
			values=genes$geneName,
			filters = "uniprot_gn_symbol", mart=ensembl)
	
	return(db1)
}

.returnDB2 <- function(genes, ensembl){
	
	db2 <- getBM(attributes=c("uniprot_gn_symbol",  # geneName
					"external_gene_name", # Complete name,   
					"chromosome_name",    # Chromosome name
					"ensembl_gene_id",    # Ensembl
					"entrezgene_id", #Entrez.Gene.ID
					"uniprot_gn_id"),     # Uniprot.ID
			values=genes$geneName,
			filters = "uniprot_gn_symbol", mart=ensembl)
	
	return(db2)
}

setMethod(
		
		f = "retrieveGenesInfo",
		
		signature = "scRNAseq",
		
		definition = function(theObject, species, groupBy="clusters", 
		orderGenes="initial", getUniprot=TRUE, silent=FALSE, cores=1){
	
	genes <- getClustersMarkers(theObject)
	databaseDict <- c(mouse = "mmusculus_gene_ensembl",
			human = "hsapiens_gene_ensembl")
	.checkParamsretrieveGenesInfo(orderGenes, genes, species, databaseDict)
	
    ## Additional Special database and columns according to species
	if(isTRUE(all.equal(species, "mouse"))){
        speDbID <- "mgi_id"
        speDBdescription <- "mgi_description"
    }else{
        speDbID <- NULL
        speDBdescription <- NULL
        spDB <- NULL
    }
    
    ## Selecting the BioMart database to use
    dataset <- databaseDict[species]
	ensembl <- useEnsembl(biomart="genes", dataset=dataset)
	
    ## Query biomart
	database <- merge(x = .returnDB1(genes, ensembl), 
			y = .returnDB2(genes, ensembl),  
			by = c("uniprot_gn_symbol", "external_gene_name"), all.x = TRUE)
    
    
    ## Group the GO id and the uniprot Id
	options(dplyr.summarise.inform = FALSE)
	database <- database %>% group_by(uniprot_gn_symbol, chromosome_name,
					entrezgene_description) %>% summarise(
					go_id = paste(unique(go_id), collapse=', '),
					uniprot_gn_id = paste(unique(uniprot_gn_id), collapse=', '),
					description = paste(unique(description), collapse=', '),
					external_gene_name = paste(unique(external_gene_name),
                                       collapse=', '),
					gene_biotype = paste(unique(gene_biotype), collapse=', '),
					ensembl_gene_id = paste(unique(ensembl_gene_id), 
							collapse=', '),
					entrezgene_id = paste(unique(entrezgene_id), collapse=', '),
        )
    
    
    ## Group the GO id and the uniprot Id
    if(!(is.null(speDbID) & is.null(speDBdescription)))
        ## Query biomart for specific db according to species
        ## external_gene_name is also searched because one gene can have
        ## several external_gene_name and so several speDbID and 
        ## speDBdescription
        spDB <- getBM(attributes=c("uniprot_gn_symbol",  # geneName
                                   "external_gene_name", speDbID, 
								   speDBdescription),  
                      values=genes$geneName,
                      filters = "uniprot_gn_symbol",
                      mart=ensembl)
        
    
    
    ## Merge the second column, the first already existing in the first table
    if(!is.null(spDB))
        database <- merge(x = database, y = spDB, by = c("uniprot_gn_symbol", 
						"external_gene_name"), all.x = TRUE)
    
	## Add cluster column
    database <- merge(x = database, y = genes, by.x = "uniprot_gn_symbol", 
			by.y = "geneName", all.x = TRUE)
    
    ## Add Symbol column
    database <- cbind(database, Symbol = database$uniprot_gn_symbol)
    
    
    ## Order the data frames colums to be abe to bind them
    colnamesOrder <- c("uniprot_gn_symbol", "clusters", "external_gene_name",
                      "go_id", speDBdescription, "entrezgene_description",
                      "gene_biotype", "chromosome_name", "Symbol", 
					  "ensembl_gene_id", speDbID, "entrezgene_id", 
					  "uniprot_gn_id")
    
    result  <- database[, colnamesOrder]
    
	if(isTRUE(all.equal(orderGenes, "initial"))){
		desiredOrder <- genes$geneName
		result <- merge(list(uniprot_gn_symbol = unique(desiredOrder)),
				result, sort = FALSE)
		rownames(result) <- seq_len(nrow(result))
	}
		 
    rm(database, colnamesOrder)
    biomartCacheClear()
	
	setGenesInfos(theObject) <- result
    return(theObject)
})


###################
## saveMarkersLists
###################


setMethod(
		
		f = "saveMarkersLists",
		
		signature = "scRNAseq",
		
		definition = function(theObject, dataDirectory, outputDir=NA,
                             pattern="genes.tsv", Ntop=100){
						 
						 
	if(is.na(outputDir))
		outputDir <- file.path(dataDirectory, 
				paste0("marker_genes/markers_lists"))
    					 
    if(!file.exists(outputDir))
        dir.create(outputDir, showWarnings=F)
    
    markerGenesList <- getMarkerGenesList(theObject)
    clustersSimiliratyOrdered <- getClustersSimiliratyOrdered(theObject)
    checkMarkerGenesList(markerGenesList, clustersSimiliratyOrdered)
    experimentName <- getExperimentName(theObject)
    clusterIndexes <- seq_len(length(markerGenesList))
	
	invisible(mapply(function(currentMarkerGenes, clustIdx, threshold, expN){
						
						markerGenes <- currentMarkerGenes[seq_len(threshold), ]
						markerGenes <- as.data.frame(markerGenes$Gene)
						outputName  <- paste0(expN,"_cluster_", clustIdx)
						colnames(markerGenes) <- "geneName"
						markerGenes$clusters <- paste0("cluster_", clustIdx)
						fileName <- paste0(outputName, "_markers.csv")
						outputFile <- file.path(outputDir, fileName)
						write.table(markerGenes, outputFile, row.names=FALSE, 
								quote=FALSE, sep=";")
						}, markerGenesList, clusterIndexes, MoreArgs=list(Ntop, 
							experimentName)))
	
    message(paste("You can find all marker genes in :", outputDir, "."))
})



###################
## saveGenesInfo
###################


setMethod(
		
		f = "saveGenesInfo",
		
		signature = "scRNAseq",
		
		definition = function(theObject, outputDir=NA){
			
			if(is.na(outputDir))
				outputDir <- file.path(getOutputDirectory(theObject),
						"marker_genes/saveGenesInfo")
			
			if(!file.exists(outputDir))
					dir.create(outputDir, recursive=TRUE)
			
			species <- getSpecies(theObject)	
			infos <- getGenesInfos(theObject)
			infosList <- split(infos, infos$clusters)
			
			invisible(lapply(infosList, function(clusterDF){
								
								clustNB <- unique(clusterDF$clusters)
								
								if(!isTRUE(all.equal(length(clustNB), 1)))
									stop("Error in saveGenesInfo. Contact the",
											" developper.")
								
								outputFile <- file.path(outputDir, 
										paste0("genes_info_clust", clustNB, 
												".csv"))
								write.csv(clusterDF, file=outputFile, 
										row.names=FALSE)
								
							}))
	})




###################
## bestClustersMarkers
###################

.checkParamsBestClust <- function(sceObject, markerGenesList, 
		numberOfClusters){
	
	## Check if the normalized matrix is correct
	if(all(dim(sceObject) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with 'rankGenes' ",
				"function doesn't have its 'sceNorm' slot updated. Please use",
				" 'normaliseCountMatrix' on the object before.")
	
	## Check if the SCE object contain cluster colums in its colDF
	if(!("clusters" %in% names(colData(sceObject))))
		stop("The 'scRNAseq' object that you're using with 'rankGenes'",
				"function doesn't have a correct 'sceNorm' slot. This slot ",
				"should be a 'SingleCellExperiment' object containing ",
				"'clusters' column in its colData.\n Please check if you ",
				"correctly used 'clusterCellsInternal' on the object.")
	
	## Check if the markerGenesList is correct
	if(!isTRUE(all.equal(length(markerGenesList),numberOfClusters)))
		stop("Something wrong with number of clusters. It is supposed",
				"to be equal to : ", numberOfClusters, ". Current ",
				"number: ", length(markerGenesList))
	
}


setMethod(
		
		f = "bestClustersMarkers",
		
		signature = "scRNAseq",
		
		definition = function(theObject, nTop=10, removeDuplicates = TRUE){
			
			sceObject <- getSceNorm(theObject)
			markerGenesList <-  getMarkerGenesList(theObject)
			clustVec <- SummarizedExperiment::colData(sceObject)$clusters
			clustersVec <- unique(clustVec)
			numberOfClusters <- length(clustersVec)
			
			.checkParamsBestClust(sceObject, markerGenesList, numberOfClusters)
			
			##Retrieve the nTop genes in each element of markerGenesList
			geneName <- unlist(lapply(markerGenesList, 
					function(currentMarkers) 
						currentMarkers$Gene[seq_len(nTop)]))
			clusters <- rep(clustersVec, each=nTop)
			markersClusters <- data.frame(geneName, clusters)
				
			if(removeDuplicates)
				markersClusters <- 
						markersClusters[!duplicated(markersClusters$geneName), ] 
			
			
			setClustersMarkers(theObject) <- markersClusters
			return(theObject)
		})

