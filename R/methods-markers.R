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


#' .buildTTestPval
#'
#' @description 
#' Generate a list of p-values of a cluster comparing it to all the other 
#' clusters.
#'
#' @param otherGroups List of marker genes for the other groups.
#'@param tTestPval An empty data frame to retrieve the results.
#' @param colDF Data frame of Metadata of the scRNA-Seq object.
#' @param mat Normalized expression matrix.
#' @param colLabel Name of the column with a clustering result. 
#' Default="clusters"
#' @param currentGroup Marker genes for the cluster currently considered in the 
#' mapply of .buildMarkerGenesList.
#' 
#' @keywords internal
#' @return A list containing for each cluster the p-value of the t-test.
#' @noRd

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


#' .buildMarkerGenesList
#'
#' @description 
#' Generate genes list for each cluster.
#'
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized (see ?normaliseCountMatrix), tSNE were calculated (see 
#' ?generateTSNECoordinates), dbScan was run (see ?runDBSCAN), cells were
#' clustered (see ?clusterCellsInternal), as clusters themselves (see 
#' ?calculateClustersSimilarity).
#' @param groups Character vector of the different clusters names.
#' @param simMedRowList List of similarity values per cluster.
#' @param exprM Normalized expression matrix.
#' @param column Name of the column with a clustering result. Default="clusters"
#' @param colNamesVec Column names of the cluster similarity matrix.
#' @param colDF Data frame of Metadata of the scRNA-Seq object. 
#' @param writeMarkerGenes If TRUE, output one list of marker genes per cluster 
#' in the output directory defined in theObject and in the sub-directory 
#' 'marker_genes'. Default=FALSE. 
#' 
#' @keywords internal
#' @return A list containing for each cluster the marker genes.
#' @importFrom utils write.table
#' @importFrom stats p.adjust
#' @importFrom stats var
#' @importFrom stats pt
#' @noRd

.buildMarkerGenesList <- function(theObject, groups, simMedRowList, exprM, 
		column, colNamesVec, colDF, writeMarkerGenes){
	
	result <- mapply(function(currentGroup, currentSimMed, mat, 
					allGroups, colLabel, simMedNames){
				
				message(paste("Working on cluster", currentGroup))
				tTestPval <- data.frame(row.names=rownames(mat))
				tTestFDR <- data.frame(row.names=rownames(mat))
				otherGroups <- allGroups[allGroups!=currentGroup]
				
				## Create a dataframe clustering vs clustering
				tTestPval <- .buildTTestPval(otherGroups, tTestPval, colDF, mat,
						colLabel, currentGroup) 
				tTestPval <- dplyr::bind_cols(tTestPval)
				tTestPval <- cbind(Gene=rownames(mat), tTestPval)
				
				## Apply the FDR
				tTestFDR <- lapply(otherGroups, function(currentOther){
							
							tTestFDR[, paste0("vs_", currentOther)] <- 
									p.adjust(tTestPval[, paste0("vs_", 
															currentOther)], 
											method="fdr")
							return(tTestFDR)
						})
				tTestFDR <- dplyr::bind_cols(tTestFDR)
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
							"_cluster_", currentGroup, "_genes.tsv")
					
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


#' rankGenes
#'
#' @description 
#' This function searches marker genes for each cluster.
#'
#' @usage 
#' rankGenes(theObject, column="clusters", writeMarkerGenes=FALSE)
#' 
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized (see ?normaliseCountMatrix), tSNE were calculated (see 
#' ?generateTSNECoordinates), dbScan was run (see ?runDBSCAN), cells were
#' clustered (see ?clusterCellsInternal), as clusters themselves (see 
#' ?calculateClustersSimilarity).
#' @param column Name of the column with a clustering result. 
#' Default="clusters"
#' @param writeMarkerGenes If TRUE, output one list of marker genes per cluster 
#' in the output directory defined in theObject and in the sub-directory 
#' 'marker_genes'. Default=FALSE.
#' 
#' @details 
#'  To understand the nature of the consensus clusters identified by CONCLUS, 
#' it is essential to identify genes which could be classified as marker genes 
#' for each cluster. To this aim, each gene should be "associated" to a 
#' particular cluster. This association is performed by looking at upregulated 
#' genes in a particular cluster compared to the others (multiple comparisons). 
#' 
#' The function *rankGenes* performs multiple comparisons of all genes from 
#' theObject and rank them according to a score reflecting a FDR power.
#' 
#' For each table corresponding to a particular consensus cluster, the first 
#' column is a gene name. The following columns represent adjusted p-values 
#' (FDR) of a one-tailed T-test between the considered cluster and all others. 
#' 
#' Top genes with significant FDR in most of the comparisons can be assumed as 
#' positive markers of a cluster. The column mean_log10_fdr is the mean power of
#'  FDR in all comparisons; the column n_05 is the number of comparisons in 
#' which the gene was significantly upregulated. The score for marker genes is 
#' the average power of FDR among all comparisons for a cluster multiplied to 
#' weights taken from the clustersSimilarityMatrix + 0.05. Taking into account 
#' both FDRs of all comparisons and clustersSimilarityMatrix allows us to keep 
#' the balance between highlighting markers for individual clusters and their 
#' 'families' which makes the final heatmap as informative as possible. 
#' 
#' Note: Adding 0.05 to the clustersSimilarityMatrix in calculating the score 
#' helps avoiding the following problem: in case you have a cluster very 
#' different from all others, it will have the value 1 on the diagonal and 0 
#' similarities to all others groups in the clustersSimilarityMatrix. So all 
#' weights for that cluster will be zeros meaning that the score would also be 
#' zero and genes will be ordered in alphabetical order in the corresponding 
#' marker genes list file.
#' 
#' For a cluster k and a gene G, a scoreG was defined in the following way:
#'  
#' scoreG= ∑i(−log10(fdrk,i + ϵ)∗weightk,i) / nClusters−1, 
#' 
#' Where 
#' 
#' 1. fdrk,i is an adjusted p-value obtained by comparing expression of G in 
#' cluster k versus expression of G in cluster i. 
#' 2. weightk,i is a similarity between these two groups taken from the element 
#' in the clustersSimilarityMatrix 
#' 3. nClusters is a number of consensus clusters given to the rankGenes() 
#' 4. ϵ = 10−300 is a small number which does not influence the ranking and 
#' added to avoid an error when fdr is equal to zero 
#' 5. k = [1,…,nClusters] 
#' 6. I = ([1,…,nClusters]exceptfor[k])
#' 
#' @aliases rankGenes
#'  
#' @author 
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.
#' 
#' @rdname 
#' rankGenes-scRNAseq
#' 
#' @return 
#' An object of class scRNASeq with its markerGenesList slot updated. 
#' 
#' @examples
#' experimentName <- "Bergiers"
#' countMatrix <- as.matrix(read.delim(file.path(
#' "tests/testthat/test_data/test_countMatrix.tsv")))
#' outputDirectory <- "./"
#' columnsMetaData <- read.delim(
#' file.path("extdata/Bergiers_colData_filtered.tsv"))
#' 
#' ## Create the initial object
#' scr <- scRNAseq(experimentName = experimentName, 
#'                 countMatrix     = countMatrix, 
#'                 species         = "mouse",
#'                 outputDirectory = outputDirectory)
#' 
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#' 
#' ## Compute the tSNE coordinates
#' scrTsne <- generateTSNECoordinates(scrNorm, cores=5)
#' 
#' ## Perform the clustering with dbScan
#' scrDbscan <- runDBSCAN(scrTsne, cores=5)
#' 
#' ## Compute the cell similarity matrix
#' scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=10, cores=4)
#' 
#' ## Calculate clusters similarity
#' scrCSM <- calculateClustersSimilarity(scrCCI)
#' 
#' ## Ranking genes
#' scrS4MG <- rankGenes(scrCSM)
#' 
#' @seealso
#' retrieveTopClustersMarkers retrieveGenesInfo 
#' 
#' @exportMethod rankGenes
#' @importFrom SummarizedExperiment colData
#' @importFrom Biobase exprs

setMethod(
		
		f = "rankGenes",
		
		signature = "scRNAseq",
		
		definition = function(theObject, column="clusters", 
				writeMarkerGenes=FALSE){
			
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
			
			markerGenesList <- .buildMarkerGenesList(theObject, groups, 
					simMedRowList, exprM, column, colNamesVec, colDF, 
					writeMarkerGenes) 		
			
			setMarkerGenesList(theObject) <- markerGenesList
			rm(markerGenesList, simMed, groups)
			
			return(theObject)
		})


###################
## retrieveGenesInfo
###################


#' .saveGenesInfo
#' 
#' @param theObject An Object of class scRNASeq for which retrieveGenesInfo was 
#' run. See ?retrieveGenesInfo.
#' 
#' @keywords internal
#' @noRd
#' 
#' @importFrom utils write.table
.saveGenesInfo <- function(theObject){
	
	outputDir <- file.path(getOutputDirectory(theObject),
			"marker_genes/saveGenesInfo")
	
	if(!file.exists(outputDir))
		dir.create(outputDir, recursive=TRUE)
	
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
}


.checkParamsretrieveGenesInfo <- function(theObject, orderGenes, genes, species,
		databaseDict, saveInfos){
	
	if(!isTRUE(all.equal(orderGenes,"initial")) && 
			!isTRUE(all.equal(orderGenes,"alphabetical")))
		stop("orderGenes should be 'initial' or 'alphabetical'.")
			
	## Check if the normalized matrix is correct
	sceObject <- getSceNorm(theObject)
	if(all(dim(sceObject) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with 'retrieveGenesInfo'",
				" function doesn't have its 'SceNorm' slot updated. Please use",
				" 'normaliseCountMatrix' on the object before.")
	
	## Check if the SCE object contain cluster colums in its colData
	if(!("clusters" %in% names(colData(sceObject))))
		stop("The 'scRNAseq' object that you're using with 'retrieveGenesInfo'",
				" function doesn't have a correct 'SceNorm' slot. This slot",
				" should be a 'SingleCellExperiment' object containing ",
				"'clusters' column in its colData. Please check if you ",
				"correctly used 'clusterCellsInternal' on the object.")
	
	## Check the cluster similarity matrix
	simMed <- getClustersSimilarityMatrix(theObject)
	if(all(dim(simMed) == c(0,0)) || all(dim(simMed) == c(1,1)))
		stop("The 'scRNAseq' object that you're using with 'retrieveGenesInfo'",
				" function doesn't have a similarity matrix, Please use ",
				"'calculateClustersSimilarity' on the object before.")
	
	## Check the marker genes
	if(isTRUE(all.equal(dim(genes), c(1,2))) && 
			isTRUE(all.equal(genes$geneName, "gene1")) &&
			is.na(genes$clusters))
		stop("The 'scRNAseq' object that you're using with 'retrieveGenesInfo'",
				" does not have marker genes. Please use ",
				"'retrieveTopClustersMarkers' before.")
	
	## Check saveInfos
	if(!is.logical(saveInfos))
		stop("saveInfos should be a boolean.")
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



#' retrieveGenesInfo
#'
#' @description 
#' This method retrieve information about the marker genes of each cluster 
#' querying the Ensembl database with biomaRt.
#'
#' @usage 
#' retrieveGenesInfo(theObject, groupBy="clusters", 
#' 				orderGenes="initial", getUniprot=TRUE, cores=1, 
#' 				saveInfos=FALSE)
#' 
#' @param theObject An Object of class scRNASeq for which the count matrix was
#' normalized (see ?normaliseCountMatrix), tSNE were calculated (see 
#' ?generateTSNECoordinates), dbScan was run (see ?runDBSCAN), cells were
#' clustered (see ?clusterCellsInternal), as clusters themselves (see 
#' ?calculateClustersSimilarity).
#' @param groupBy A column in the input table used for grouping the genes in 
#' the output tables. This option is useful if a table contains genes from 
#' different clusters. Default = "clusters".
#' @param orderGenes If "initial" then the order of genes will not be changed. 
#' The other option is "alphabetical". Default = "initial".
#' @param getUniprot Boolean, whether to get information from UniProt or not. 
#' Default = TRUE.
#' @param cores Maximum number of jobs that CONCLUS can run in parallel. 
#' Default is 1.
#' @param saveInfos If TRUE, save the genes infos table in the directory 
#' defined in theObject (?getOutputDirectory) and in the sub-directory 
#' 'marker_genes/saveGenesInfo'.
#' 
#' @details 
#' The output dataframe is composed of the following columns:
#'   
#'  - uniprot_gn_symbol: Uniprot gene symbol.
#'  - clusters: The cluster to which the gene is associated.
#'  - external_gene_name: The complete gene name.
#'  - go_id: Gene Ontology (GO) identification number.
#'  - mgi_description: If the species is mouse, description of the gene on MGI.
#'  - entrezgene_description: Description of the gene by Entrez database.
#'  - gene_biotype: protein coding gene, lincRNA gene, miRNA gene, unclassified 
#' non-coding RNA gene, or pseudogene.
#'  - chromosome_name: The chromosome on which the gene is located.
#'  - Symbol: Official gene symbol.
#'  - ensembl_gene_id: ID of the gene on the ensembl database.
#'  - mgi_id: If the species is mouse, ID of the gene on the MGI database.
#'  - entrezgene_id: ID of the gene on the entrez database.
#'  - uniprot_gn_id: ID of the gene on the uniprot database.
#' 
#' @aliases retrieveGenesInfo
#'  
#' @rdname 
#' retrieveGenesInfo-scRNAseq
#' 
#' @return 
#' An object of class scRNASeq with its genesInfos slot updated. 
#' 
#' @examples
#' experimentName <- "Bergiers"
#' countMatrix <- as.matrix(read.delim(file.path(
#' "tests/testthat/test_data/test_countMatrix.tsv")))
#' outputDirectory <- "./"
#' columnsMetaData <- read.delim(
#' file.path("extdata/Bergiers_colData_filtered.tsv"))
#' 
#' ## Create the initial object
#' scr <- scRNAseq(experimentName = experimentName, 
#'                 countMatrix     = countMatrix, 
#'                 species         = "mouse",
#'                 outputDirectory = outputDirectory)
#' 
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#' 
#' ## Compute the tSNE coordinates
#' scrTsne <- generateTSNECoordinates(scrNorm, cores=5)
#' 
#' ## Perform the clustering with dbScan
#' scrDbscan <- runDBSCAN(scrTsne, cores=5)
#' 
#' ## Compute the cell similarity matrix
#' scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=10, cores=4)
#' 
#' ## Calculate clusters similarity
#' scrCSM <- calculateClustersSimilarity(scrCCI)
#' 
#' ## Ranking genes
#' scrS4MG <- rankGenes(scrCSM)
#' 
#' ## Getting marker genes
#' scrFinal <- retrieveTopClustersMarkers(scrS4MG, removeDuplicates = F)
#' 
#' ## Getting genes info
#' scrInfos <- retrieveGenesInfo(scrFinal, species = "mouse", cores=5)
#' 
#' @seealso
#' rankGenes retrieveTopClustersMarkers 
#' 
#' @exportMethod retrieveGenesInfo
#' @importFrom biomaRt useEnsembl biomartCacheClear
#' @importFrom dplyr group_by %>% summarise
#' 
#' @author 
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.


setMethod(
		
		f = "retrieveGenesInfo",
		
		signature = "scRNAseq",
		
		definition = function(theObject, groupBy="clusters", 
				orderGenes="initial", getUniprot=TRUE, cores=1, 
				saveInfos=FALSE){
			
			species <- getSpecies(theObject)
			genes <- getClustersMarkers(theObject)
			databaseDict <- c(mouse = "mmusculus_gene_ensembl",
					human = "hsapiens_gene_ensembl")
			.checkParamsretrieveGenesInfo(theObject, orderGenes, genes, species, 
					databaseDict, saveInfos)
			
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
					by = c("uniprot_gn_symbol", "external_gene_name"), 
					all.x = TRUE)
			
			
			## Group the GO id and the uniprot Id
			options(dplyr.summarise.inform = FALSE)
			database <- database %>% group_by(.data$uniprot_gn_symbol, 
							chromosome_name, entrezgene_description) %>% 
					summarise(go_id = paste(unique(go_id), collapse=', '),
							uniprot_gn_id = paste(unique(uniprot_gn_id), 
									collapse=', '),
							description = paste(unique(description), 
									collapse=', '),
							external_gene_name = paste(
									unique(external_gene_name), collapse=', '),
							gene_biotype = paste(unique(gene_biotype), 
									collapse=', '),
							ensembl_gene_id = paste(unique(ensembl_gene_id), 
									collapse=', '),
							entrezgene_id = paste(unique(entrezgene_id), 
									collapse=', '),
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
			
			
			
			## Merge the second column, the first already existing in the first 
			## table
			if(!is.null(spDB))
				database <- merge(x = database, y = spDB, by = 
								c("uniprot_gn_symbol", "external_gene_name"), 
						all.x = TRUE)
			
			## Add cluster column
			database <- merge(x = database, y = genes, 
					by.x = "uniprot_gn_symbol",  by.y = "geneName", 
					all.x = TRUE)
			
			## Add Symbol column
			database <- cbind(database, Symbol = database$uniprot_gn_symbol)
			
			
			## Order the data frames colums to be abe to bind them
			colnamesOrder <- c("uniprot_gn_symbol", "clusters", 
					"external_gene_name", "go_id", speDBdescription, 
					"entrezgene_description", "gene_biotype", "chromosome_name",
					"Symbol", "ensembl_gene_id", speDbID, "entrezgene_id", 
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
			
			if(saveInfos)
				.saveGenesInfo(theObject)
			
			return(theObject)
		})



###################
## retrieveTopClustersMarkers
###################

.checkParamsTopClust <- function(sceObject, markerGenesList, 
		numberOfClusters){
	
	## Check if the normalized matrix is correct
	if(all(dim(sceObject) == c(0,0)))
		stop("The 'scRNAseq' object that you're using with ",
				"'retrieveTopClustersMarkers' function doesn't have its ",
				"'sceNorm' slot updated. Please use 'normaliseCountMatrix' on",
				" the object before.")
	
	## Check if the SCE object contain cluster colums in its colDF
	if(!("clusters" %in% names(colData(sceObject))))
		stop("The 'scRNAseq' object that you're using with ",
				"'retrieveTopClustersMarkers' function doesn't have a correct",
				" 'sceNorm' slot. This slot should be a 'SingleCellExperiment'",
				" object containing 'clusters' column in its colData. Please ",
				"check if you correctly used 'clusterCellsInternal' on the ",
				"object.")
	
	## Check if the markerGenesList is correct
	if(!isTRUE(all.equal(length(markerGenesList),numberOfClusters)))
		stop("Something wrong with number of clusters. It is supposed",
				" to be equal to : ", numberOfClusters, ". Current ",
				"number: ", length(markerGenesList), ". Did you use ",
				"'calculateClustersSimilarity' and 'rankGenes'?")
}


#' .writeMarkersList
#'
#' @description 
#' Writes one list per cluster in the output folder defined in theObject, and 
#' in the sub-directory marker_genes/markers_lists.
#'
#' @param theObject An Object of class scRNASeq for which rankGenes was 
#' run. See ?rankGenes.
#' @param markersClusters Data frame containing two columns geneName and 
#' clusters.
#' @param nTop Number of marker genes to retrieve per cluster.
#' 
#' @keywords internal
#' @noRd

.writeMarkersList <- function(theObject, markersClusters, nTop){
	
	## Creating the output folder
	dataDirectory <- getOutputDirectory(theObject)
	outputDir <- file.path(dataDirectory, paste0("marker_genes/markers_lists"))
	
	if(!file.exists(outputDir))
		dir.create(outputDir, showWarnings=F)
	
	message("Writing lists of markers to ", outputDir)
	
	## Creating a list of markers per clusters
	markersClustersList <- split(markersClusters, markersClusters$clusters)
	
	## Writing each element of the lists to a corresponding file
	experimentName <- getExperimentName(theObject)
	invisible(mapply(function(currentMarkers, currentClustNb, thres, expN){
						
						outputName  <- paste0(expN,"_cluster_", currentClustNb)
						fileName <- paste0(outputName, "_markers.csv")
						outputFile <- file.path(outputDir, fileName)
						write.table(currentMarkers, file=outputFile, sep="\t",
								quote=FALSE, row.names=FALSE, col.names=TRUE)
						
					}, markersClustersList, names(markersClustersList), 
					MoreArgs=list(nTop, experimentName)))
}

#' retrieveTopClustersMarkers
#'
#' @description 
#' This function retrieves the top N marker genes for each cluster. 
#'
#' @usage 
#' retrieveTopClustersMarkers(theObject, nTop=10, removeDuplicates = TRUE,
#' 				writeMarkerGenes = FALSE)
#' 
#' @param theObject An Object of class scRNASeq for which rankGenes was 
#' run. See ?rankGenes.
#' @param nTop Number of marker genes to retrieve per cluster. Default=10.
#' @param removeDuplicates If TRUE, duplicated markers are removed from the 
#' lists. Default=TRUE.
#' @param writeMarkerGenes If TRUE, writes one list per cluster in the output
#' folder defined in theObject, and in the sub-directory 
#' marker_genes/markers_lists. Default=FALSE.
#' 
#' @aliases retrieveTopClustersMarkers
#' @rdname 
#' retrieveTopClustersMarkers-scRNAseq
#' 
#' @return 
#' Output the list of markers to marker_genes/markers_lists if writeMarkersGenes
#' is TRUE and return a scRNASeq object with its clustersMarkers slot updated.  
#' 
#' @examples
#' experimentName <- "Bergiers"
#' countMatrix <- as.matrix(read.delim(file.path(
#' "tests/testthat/test_data/test_countMatrix.tsv")))
#' outputDirectory <- "./"
#' columnsMetaData <- read.delim(
#' file.path("extdata/Bergiers_colData_filtered.tsv"))
#' 
#' ## Create the initial object
#' scr <- scRNAseq(experimentName = experimentName, 
#'                 countMatrix     = countMatrix, 
#'                 species         = "mouse",
#'                 outputDirectory = outputDirectory)
#' 
#' ## Normalize and filter the raw counts matrix
#' scrNorm <- normaliseCountMatrix(scr, coldata = columnsMetaData)
#' 
#' ## Compute the tSNE coordinates
#' scrTsne <- generateTSNECoordinates(scrNorm, cores=5)
#' 
#' ## Perform the clustering with dbScan
#' scrDbscan <- runDBSCAN(scrTsne, cores=5)
#' 
#' ## Compute the cell similarity matrix
#' scrCCI <- clusterCellsInternal(scrDbscan, clusterNumber=10, cores=4)
#' 
#' ## Calculate clusters similarity
#' scrCSM <- calculateClustersSimilarity(scrCCI)
#' 
#' ## Ranking genes
#' scrS4MG <- rankGenes(scrCSM)
#' 
#' ## Retrieve the top 10 markers per cluster
#' scrFinal <- retrieveTopClustersMarkers(theObject)
#' 
#' @seealso
#' retrieveGenesInfo
#' 
#' @exportMethod retrieveTopClustersMarkers
#' @author 
#' Ilyess RACHEDI, based on code by Polina PAVLOVICH and Nicolas DESCOSTES.

setMethod(
		
		f = "retrieveTopClustersMarkers",
		
		signature = "scRNAseq",
		
		definition = function(theObject, nTop=10, removeDuplicates = TRUE,
				writeMarkerGenes = FALSE){
			
			markerGenesList <-  getMarkerGenesList(theObject)
			sceObject <- getSceNorm(theObject)
			clustVec <- SummarizedExperiment::colData(sceObject)$clusters
			clusterIndexes <- unique(clustVec)
			numberOfClusters <- length(clusterIndexes)
			
			.checkParamsTopClust(sceObject, markerGenesList, numberOfClusters)
			
			##Retrieve the nTop genes in each element of markerGenesList
			geneName <- unlist(lapply(markerGenesList, 
							function(currentMarkers) 
								currentMarkers$Gene[seq_len(nTop)]))
			clusters <- rep(clusterIndexes, each=nTop)
			markersClusters <- data.frame(geneName, clusters)
			
			if(removeDuplicates)
				markersClusters <- 
						markersClusters[!duplicated(markersClusters$geneName), ] 
			
			if(writeMarkerGenes)
				.writeMarkersList(theObject, markersClusters, nTop)
				
			setClustersMarkers(theObject) <- markersClusters
			return(theObject)
		})
