.checkParamsRankGenes <- function(sceObject){
	
	## Check if the normalized matrix is correct
	if(all(dim(sceObject) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with 'rankGenes' ",
				"function doesn't have its 'SceNorm' slot updated. Please",
				" use 'normaliseCountMatrix' on the object before.")
	
	## Check if the SCE object contain cluster colums in its colData
	if(!("clusters" %in% names(colData(sceObject))))
		stop("The 'scRNAseq' object that you're using with 'rankGenes' ",
				"function doesn't have a correct 'SceNorm' slot. This slot",
				"should be a 'SingleCellExperiment' object containing ",
				"'clusters' column in its colData. Please check if you ",
				"correctly used 'clusterCellsInternal' on the object.")
	
	## Check the cluster similarity matrix
	if(all(dim(sceObject) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with 'rankGenes' ",
				"function doesn't have its 'clustersSimilarityMatrix' slot",
				" updated. Please use 'clusterCellsInternal' on the object",
				" before.")
}


setMethod(
		
    f = "rankGenes",
	
    signature = "scRNAseq",
	
    definition = function(theObject, column="clusters", writeMarkerGenes=FALSE){
        
        sceObject  <- getSceNorm(theObject)
		.checkParamsRankGenes(sceObject)
		simMed <- getClustersSimilarityMatrix(theObject)
        exprM <- Biobase::exprs(sceObject)
        colDF <- SummarizedExperiment::colData(sceObject)
         
        message("Ranking marker genes for each cluster.")
        
        stopifnot(all(colnames(exprM) == rownames(colDF)))
        groups <- unique(colDF[, column])
        simMed <- simMed + 0.05
        
		markerGenesList <- lapply(groups, function(currentGroup, mat, 
						allGroups, colLabel){
					
					message(paste("Working on cluster", currentGroup))
					tTestPval <- data.frame(row.names=rownames(mat))
					tTestFDR <- data.frame(row.names=rownames(mat))
					otherGroups <- allGroups[allGroups!=currentGroup]
					
					## Create a dataframe clustering vs clustering
					tTestPval <- lapply(otherGroups, function(currentOther){
								
								tTestPval[, paste0("vs_", currentOther)] <- NA
								x <- mat[, colDF[, c(colLabel)] == currentGroup]
								y <- mat[, colDF[, c(colLabel)] == currentOther]
								t <- (rowMeans(x)-rowMeans(y))/sqrt(apply(mat, 1, var)*(1 / ncol(x) + 1 / ncol(y)))
								df <- ncol(x)+ncol(y)-2
								tTestPval[, paste0("vs_", currentOther)] <- pt(t, df, lower.tail=FALSE)
								return(tTestPval)
							})
					
					tTestPval <- bind_cols(tTestPval)
					tTestPval <- cbind(Gene=rownames(mat), tTestPval)
					
					## Apply the FDR
					tTestFDR <- lapply(otherGroups, function(currentOther){
								
								tTestFDR[, paste0("vs_", currentOther)] <- p.adjust(tTestPval[, paste0("vs_", currentOther)], method="fdr")
								return(tTestFDR)
							})
					tTestFDR <- bind_cols(tTestFDR)
					tTestFDR <- cbind(Gene=as.factor(rownames(mat)), tTestFDR)
					
					submat <- as.matrix(tTestFDR[, 2:(length(otherGroups) + 1 )])
					
					## Add column mean_log10_fdr 
					tTestFDR$mean_log10_fdr <- rowMeans(log10(submat + 1e-300))
					
					## Add column n_05
					tTestFDR$n_05 <- apply(submat, 1, function(x) length(x[!is.na(x) & x < 0.05]))
					
					##Add column score 
					weights <- simMed[i, otherGroups]
					tTestFDR$score <- apply(submat, 1, function(x) sum(-log10(x+1e-300) * weights) / ncol(submat))
					
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
					
				}, exprM, groups, column)
							
		
		
							
				
                
                
                
                
                
			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
				
				markerGenesList <- lapply(seq_len(length(groups)), function(i){
							
									
							message(paste("Working on cluster", i))
							tTestPval <- data.frame(row.names=rownames(exprM))
							otherGroups <- groups[groups!=groups[i]]
							
							## Create a dataframe clustering vs clustering
							tTestPval <- lapply(seq_len(length(otherGroups)),
									function(k){
										
												tTestPval[, paste0("vs_", otherGroups[k])] <- NA
										x <- exprM[, colDF[, c(column)] == groups[i]]
										y <- exprM[, colDF[, c(column)] == otherGroups[k]]
										t <- (rowMeans(x) - rowMeans(y)) /
												sqrt(apply(exprM, 1, var) *
																(1 / ncol(x) + 1 / ncol(y)))
										df <- ncol(x) + ncol(y) - 2
										tTestPval[, paste0("vs_", otherGroups[k])] <-
												pt(t, df, lower.tail=FALSE)
										return(tTestPval)})
							
							tTestPval <- bind_cols(tTestPval)
							tTestPval <- cbind(Gene=rownames(exprM), tTestPval)
							
							## Apply the FDR
							tTestFDR <- data.frame(row.names=rownames(exprM))
							
							tTestFDR <- lapply(seq_len(length(otherGroups)), function(l){
												
												tTestFDR[, paste0("vs_", otherGroups[l])] <-
												p.adjust(tTestPval[, 
																paste0("vs_", 
																		otherGroups[l])],
														method="fdr")
										return(tTestFDR)})
							
							tTestFDR <- bind_cols(tTestFDR)
							tTestFDR <- cbind(Gene=as.factor(rownames(exprM)), tTestFDR)
							
							submat <- as.matrix(tTestFDR[, 2:(length(otherGroups) + 1 )])
							
							## Add column mean_log10_fdr 
							tTestFDR$mean_log10_fdr <- rowMeans(log10(submat + 1e-300))
							
							## Add column n_05
							tTestFDR$n_05 <- apply(submat, 1, function(x)
										length(x[!is.na(x) & x < 0.05]))
							
							##Add column score 
							weights <- simMed[i, otherGroups]
							tTestFDR$score <- apply(submat, 1, function(x)
										sum(-log10(x+1e-300) * weights) / ncol(submat))
							
							tTestFDR <- tTestFDR[order(tTestFDR$score, decreasing=TRUE), ]
							row.names(tTestFDR) <- NULL
							
							## Write list if option = TRUE
							if(writeMarkerGenes){
								
								outputDir <- file.path(getOutputDirectory(theObject), 
										"marker_genes")
								outputFile <- paste0(getExperimentName(theObject),
										"_cluster_", groups[i], "_genes.tsv")
								
								write.table(tTestFDR, file.path(outputDir, outputFile),
										col.names=TRUE, row.names=FALSE, quote=FALSE,
										sep="\t")
							}
							
							return(tTestFDR)
						})
				
				
				
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
			
			
        setMarkerGenesList(theObject) <- markerGenesList
        rm(markerGenesList, simMed, groups)
        
        return(theObject)
    })



getGenesInfo <- function(theObject, species, groupBy="clusters", 
		orderGenes="initial", getUniprot=TRUE, silent=FALSE, cores=1){
    
	genes <- getClustersMarkers(theObject)
			
    if (orderGenes != "initial" && orderGenes != "alphabetical")
        stop("orderGenes should be 'initial' or 'alphabetical'.")
    
    if (!("geneName" %in% colnames(genes))){
        stop("genes should be a dataframe containing at least geneName column.")
    }
    
    ## Dictionnary to access to the dataset according the species
    databaseDict     <- c(mouse = "mmusculus_gene_ensembl",
                          human = "hsapiens_gene_ensembl")
    
    if (!(species %in% names(databaseDict))){
        msg <- paste(names(databaseDict), collapse = " or ")
        stop("species should be: ", msg)
    }
    
    ## Additional Special database and columns according to species
    if (species == "mouse"){
        speDbID <- "mgi_id"
        speDBdescription <- "mgi_description"
        
    } else {
        speDbID <- NULL
        speDBdescription <- NULL
        spDB <- NULL
    }
    
    
    ## Selecting the BioMart database to use
    dataset <- databaseDict[species]
	ensembl <- useEnsembl(biomart="genes", dataset=dataset)
	
    ## Query biomart
    db1 <- getBM(attributes=c("uniprot_gn_symbol",  # geneName
                              "external_gene_name", # Complete name,
                              "description",        # Description
                              "entrezgene_description",
                              "gene_biotype",        # Feature.Type
                              "go_id" ),             # Gene Ontology
                 values=genes$geneName,
                 filters = "uniprot_gn_symbol",
                 mart=ensembl)
    db2 <- getBM(attributes=c("uniprot_gn_symbol",  # geneName
                              "external_gene_name", # Complete name,   
                              "chromosome_name",    # Chromosome name
                              "ensembl_gene_id",    # Ensembl
                              "entrezgene_id", #Entrez.Gene.ID
                              "uniprot_gn_id"),     # Uniprot.ID
                 values=genes$geneName,
                 filters = "uniprot_gn_symbol",
                 mart=ensembl)
    
    
    database <- merge(x = db1,
                      y = db2, 
                      by = c("uniprot_gn_symbol", "external_gene_name") ,
                      all.x = TRUE)
    
    
    ## Group the GO id and the uniprot Id
    database <- database %>%
        group_by(
            uniprot_gn_symbol,
            chromosome_name,
            entrezgene_description
        ) %>%
        summarise(
            go_id = paste(unique(go_id),
                          collapse = ', '),
            uniprot_gn_id = paste(unique(uniprot_gn_id),
                                  collapse = ', '),
            description = paste(unique(description),
                                collapse = ', '),
            external_gene_name = paste(unique(external_gene_name),
                                       collapse = ', '),
            gene_biotype = paste(unique(gene_biotype),
                                 collapse = ', '),
            ensembl_gene_id = paste(unique(ensembl_gene_id),
                                    collapse = ', '),
            entrezgene_id = paste(unique(entrezgene_id),
                                  collapse = ', '),
        )
    
    
    ## Group the GO id and the uniprot Id
    if (!(is.null(speDbID) & is.null(speDBdescription))){
        ## Query biomart for specific db according to species
        ## external_gene_name is also searched because one gene can have
        ## several external_gene_name and so several speDbID and 
        ## speDBdescription
        spDB <- getBM(attributes=c("uniprot_gn_symbol",  # geneName
                                   "external_gene_name",
                                   speDbID,
                                   speDBdescription),  
                      values=genes$geneName,
                      filters = "uniprot_gn_symbol",
                      mart=ensembl)
        
    }
    
    ## Merge the second column, the first already existing in the first table
    if (!is.null(spDB)){
        database <- merge(x = database,
                          y = spDB, 
                          by = c("uniprot_gn_symbol", "external_gene_name") ,
                          all.x = TRUE)
    }
    
    ## Add cluster column
    database <- merge(x = database,
                      y = genes, 
                      by.x = "uniprot_gn_symbol", by.y = "geneName",
                      all.x = TRUE)
    
    ## Add Symbol column
    database <- cbind(database, Symbol = database$uniprot_gn_symbol)
    
    
    ## Order the data frames colums to be abe to bind them
    colnamesOrder = c("uniprot_gn_symbol",
                      "clusters", 
                      "external_gene_name",
                      "go_id",
                      speDBdescription,
                      "entrezgene_description",
                      "gene_biotype", 
                      "chromosome_name",
                      "Symbol",
                      "ensembl_gene_id",
                      speDbID,
                      "entrezgene_id",
                      "uniprot_gn_id")
    
    result  <- database[, colnamesOrder]
    
    if (orderGenes == "initial"){
        desiredOrder <- genes$geneName
        result <-
            merge(list(uniprot_gn_symbol = unique(desiredOrder)),
                  result,
                  sort = FALSE)
        rownames(result) <- seq_len(nrow(result)) 
    }
    
    # inserting space for comments
    if (any(colnames(result) %in% groupBy) &
        (orderGenes == "initial") &
        (length(unique(result$clusters)) > 1)){
        resultFinal <- result
        groupingTable <- table(resultFinal[, groupBy])
        
        groupingTable <- groupingTable[unique(resultFinal$clusters)]
        resultFinal <- InsertRow(resultFinal, c("For notes:",
                                                rep("", (ncol(result)) - 1)),
                                 RowNum=1)
        
        RowNum <- groupingTable[1] + 1
        
        for(i in 1:(length(groupingTable)-1)){
            resultFinal <- InsertRow(resultFinal, rep("", ncol(result)),
                                     RowNum=(RowNum + 1))
            resultFinal <- InsertRow(resultFinal, c("For notes:",
                                                    rep("",
                                                        (ncol(result)) - 1)
            ),
            RowNum=(RowNum + 2))
            RowNum <- RowNum + 2 + groupingTable[i +1]
        }
        result <- resultFinal
        rm(resultFinal)
    }
    rm(database, colnamesOrder)
    biomartCacheClear()
    return(result)
}



saveMarkersLists <- function(theObject,
                             dataDirectory,
                             outputDir=file.path(
                                 dataDirectory,
                                 paste0("marker_genes/markers_lists")),
                             pattern="genes.tsv", Ntop=100){
    
    if(!file.exists(outputDir))
        dir.create(outputDir, showWarnings=F)
    
    markerGenesList <- getMarkerGenesList(theObject)
    clustersSimiliratyOrdered <- getClustersSimiliratyOrdered(theObject)
    checkMarkerGenesList(markerGenesList, clustersSimiliratyOrdered)
    experimentName <- getExperimentName(theObject)
    
    tmp <- sapply(seq_len(length(markerGenesList)), function(i){
        markerGenes <- markerGenesList[[i]][1:Ntop, ]
        markerGenes <- as.data.frame(markerGenes$Gene)
        outputName  <- paste0(experimentName,"_cluster_", i)
        colnames(markerGenes) <- "geneName"
        markerGenes$clusters <- paste0("cluster_", i)
        write.table(markerGenes, file.path(outputDir,
                                           paste0(outputName, "_markers.csv")),
                    row.names=FALSE, quote=FALSE, sep=";")
        
    })
    message(paste("You could find all marker genes in :", outputDir, "."))
}



saveGenesInfo <- function(species,
                          dataDirectory="",
                          outputDir="",
                          sep=";", header=TRUE,
                          startFromCluster=1, #outputFormat=c("csv", "xlsx"),
                          groupBy="clusters", # getGenesInfo params
                          orderGenes="initial",
                          getUniprot=TRUE,
                          silent=FALSE, cores=1){
    
    if(dataDirectory != ""){
        inputDir = file.path(dataDirectory, "/marker_genes/markers_lists")
        outputDir = file.path(dataDirectory, "/marker_genes/saveGenesInfo")
        pattern = "markers.csv"
        if(!file.exists(outputDir))
            dir.create(outputDir, showWarnings=F)
    }
    
    filesList <- list.files(inputDir, pattern = pattern)
    
    if (is.null(filesList)){
        stop("There is no file finishing by ", pattern)
    }
    
    filePrefix <- do.call(rbind, strsplit(filesList, "[.]"))[, 1]
    
    
    infos <- sapply(seq(startFromCluster, length(filesList)), 
                    function(i){
                        genes <- read.delim(file.path(inputDir, filesList[i]),
                                            sep = sep, header = header,
                                            stringsAsFactors = FALSE)
                        
                        result <- getGenesInfo(genes, 
                                               species, 
                                               groupBy=groupBy,
                                               orderGenes=orderGenes,
                                               getUniprot=getUniprot,
                                               silent=silent, cores=cores)
                        message("Writing the output file number ", i, "\n")
                        
                        
                        write.table(
                            result,
                            file = file.path(outputDir,
                                             paste0(filePrefix[i],
                                                    "_genesInfo.csv")),
                            quote = FALSE,
                            sep = ";",
                            row.names = FALSE
                        )
                    }
    )
}


setMethod(
		f="bestClustersMarkers",
		signature="scRNAseq",
		definition=function(theObject, genesNumber=10, removeDuplicates = TRUE){
			
			sceObject       <- getSceNorm(theObject)
			dataDirectory   <- getOutputDirectory(theObject)
			experimentName  <- getExperimentName(theObject)
			markerGenesList <-  getMarkerGenesList(theObject)
			
			## Check if the normalized matrix is correct
			if(all(dim(sceObject) == c(0,0)))
				stop("The 'scRNAseq' object that you're using with 'rankGenes'",
						"function doesn't have its 'sceNorm' slot updated.",
						"Please use 'normaliseCountMatrix' on the object ",
						"before.")
				
				
			## Check if the SCE object contain cluster colums in its colDF
			if(!("clusters" %in% names(colData(sceObject))))
				stop("The 'scRNAseq' object that you're using with 'rankGenes'",
						"function doesn't have a correct 'sceNorm' slot. This",
						"slot should be a 'SingleCellExperiment' object",
						"containing 'clusters' column in its colData.\n",
						"Please check if you correctly used ",
						"'clusterCellsInternal' on the object.")
			
			## Check if the markerGenesList is correct
			numberOfClusters <- length(unique(
							SummarizedExperiment::colData(sceObject)$clusters))
			
			if(!isTRUE(all.equal(length(markerGenesList),numberOfClusters)))
				stop("Something wrong with number of clusters. It is supposed",
						"to be equal to : ", numberOfClusters, ". Current ",
						"number: ", length(markerGenesList))
				
			numberOfClusters <- length(unique(
							SummarizedExperiment::colData(sceObject)$clusters))
			
			dir <- file.path(dataDirectory, "marker_genes")
			nTop <- genesNumber
			clusters <- unique(
					SummarizedExperiment::colData(sceObject)$clusters)
			
			markersClusters <- as.data.frame(matrix(ncol = 2,
							nrow = nTop*numberOfClusters))
			
			colnames(markersClusters) = c("geneName", "clusters")
			runUntil = length(markerGenesList)
			
			for(i in seq_len(runUntil)){
				tmpAll <- markerGenesList[[i]]
				
				markersClusters$clusters[(nTop*(i-1)+1):(nTop*i)] <- 
						as.character(clusters[i])
				
				markersClusters$geneName[(nTop*(i-1)+1):(nTop*i)] <- 
						as.character(tmpAll$Gene[1:nTop])
			}
			
			if(removeDuplicates)
				markersClusters <- 
						markersClusters[!duplicated(markersClusters$geneName), ] 
			
			
			setClustersMarkers(theObject) <- markersClusters
			return(theObject)
		})

