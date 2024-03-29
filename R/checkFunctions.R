.checkObject <- function(obj, val, getElementFunction, className){

    element <- getElementFunction(obj)

    if(isTRUE(all.equal(className, "Dbscan")))
        nbElement <- ncol(element)
    else
        nbElement <- nrow(element)

    if(isTRUE(all.equal(nbElement, val)) && is(obj) == className)
        return(TRUE)
    else
        return(FALSE)
}

checkList <- function(theobjectList, getFunction, className){

    if(!isTRUE(all.equal(className, "Dbscan")) &&
        !isTRUE(all.equal(className, "Tsne")))
        stop("In checkList, className should be Dbscan or Tsne.")

    if(isTRUE(all.equal(className, "Dbscan")))
        val <- ncol(getFunction(theobjectList[[1]]))
    else
        val <- nrow(getFunction(theobjectList[[1]]))

    vec <- vapply(theobjectList, .checkObject, val, getFunction,
            className, FUN.VALUE=logical(1))

    if(!all(vec)){
        theobjectList[!vec]
        names <- vapply(theobjectList , getName, FUN.VALUE=character(1))
        names <- paste(names , collapse="; ")
        stop("The elements in ", className, "List slot don't have the same ",
                "number of cells or the same class")
    }
}


checkMarkerGenesList <- function(markerGeneobjectlist,
        clustersSimiliratyOrdered=NULL){

    expecteColumn <- c("Gene", "mean_log10_fdr", "n_05", "score")

    vec <- vapply(markerGeneobjectlist,

                    function(df, expcol){
                        nameColumn <- colnames(df)
                        if(all(expcol %in% nameColumn))
                            return(TRUE)
                        else
                            return(FALSE)
                        },

                    expecteColumn, FUN.VALUE=logical(1)
                    )

    if(!all(vec))
        stop("'markerGenesList' slot should contain a list of ",
                "dataframes with at least following columns : 'Gene', ",
                "'mean_log10_fdr', ", "'n_05', 'score'.\n")


    if(!is.null(clustersSimiliratyOrdered) &&
            all(!is.na(markerGeneobjectlist[[1]]$score)))
        if(!isTRUE(all.equal(length(markerGeneobjectlist),
                        length(clustersSimiliratyOrdered))))
            stop("'markerGenesList' should contain as many dataframes as ",
                    "clusters found. Number of dataframes :",
                    length(markerGeneobjectlist), " and the number of cluters ",
                    "found is :", length(clustersSimiliratyOrdered), ".")
}


checkTopMarkers <- function(topMarkers, clustersSimiliratyOrdered=NULL){


    expectedColumn <- c("geneName", "clusters")

    if(!all(expectedColumn %in% colnames(topMarkers)))
        stop("The topMarkers data frame should have the columns 'geneName'",
                " and 'clusters'")


    nbClustMark <- length(unique(topMarkers$clusters))
    nbClust <- length(clustersSimiliratyOrdered)

    if(nrow(topMarkers) > 1 && !isTRUE(all.equal(nbClust, nbClustMark)))
        stop("topMarkers should have the same number of clusters than the",
                " number of clusters found. Nb clusters for markers: ",
                nbClustMark, ". Nb of clusters: ", nbClust, ".\n",
                "If you use 'retrieveTopClustersMarkers', try setting ",
                "removeDuplicates = FALSE")
}


checkGenesInfos <- function(genesInfos, species, clustersSimiliratyOrdered){

    expectedColumn <- c("uniprot_gn_symbol", "clusters",
                        "go_id", "entrezgene_description", "gene_biotype",
                        "chromosome_name", "Symbol", "ensembl_gene_id",
                        "entrezgene_id", "uniprot_gn_id")

    if(isTRUE(all.equal(species, "mouse")))
        expectedColumn <- c(expectedColumn, "mgi_description", "mgi_id")

    if(!all(expectedColumn %in% colnames(genesInfos)))
        stop("The genesInfos data frame should have the columns: ",
                paste(expectedColumn, collapse=";"))

    clustersVec <- unique(genesInfos$clusters)
    clustersVec <- clustersVec[clustersVec != ""]
    nbClustInfo <- length(clustersVec)
    nbClust <- length(clustersSimiliratyOrdered)

    if(nrow(genesInfos) > 1 && !isTRUE(all.equal(nbClust, nbClustInfo)))
        stop("genesInfos should have the same number of clusters than the",
                " number of clusters found. Nb clusters for genesInfos: ",
                nbClustInfo, ". Nb of clusters: ", nbClust)
}
