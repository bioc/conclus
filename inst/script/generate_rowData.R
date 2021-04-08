################################################################################
###########################   extdata/rowData.tsv   ############################
################################################################################
##
## Row data used for examples and tests.
## Size : 20 * 5.
##
## 
## This is the code to generate the extdata/rowData.tsv :


library(biomaRt)
genomeAnnot <- org.Mm.eg.db::org.Mm.eg.db

symbolGenes <- c("Sipa1l2", "Mcts1", "Mov10", "Rtca", "Fmr1", "Sp1", "Uck1", 
                 "Gm2a", "Napsa", "Tsc2", "Rab5b", "Efnb2", "Apoe", "S100a4", 
                 "Eef1e1", "Braf", "Usp32", "Rapsn", "Ndufa9", "Ap4e1")


## Get the rowdata containing SYMBOL ID and ENSEMBL_ID equivalent
rowdata <- AnnotationDbi::select(genomeAnnot, 
                                 keys=symbolGenes, 
                                 keytype="SYMBOL",
                                 columns=c("SYMBOL", "ENSEMBL", "GENENAME"), 
                                 multiVals="first")

rowdata <- rowdata[!duplicated(rowdata$SYMBOL), ]

## Use the retrived ensembl id to have other informations
ensemblID <- rowdata$ENSEMBL
ensembl<- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
attributes <- c("ensembl_gene_id", "go_id", "chromosome_name", "gene_biotype")

res <- getBM(attributes=attributes, mart=ensembl, values=ensemblID, 
             filters="ensembl_gene_id")

tmp <- res[!duplicated(res$ensembl_gene_id), ]

rowdata <- merge(rowdata, tmp[c("ensembl_gene_id", "chromosome_name",
                "gene_biotype")], by.x = "ENSEMBL", by.y = "ensembl_gene_id", 
                all.x = TRUE, all.y = FALSE, sort = FALSE)

rownames(rowdata) <- rowdata$SYMBOL
colnames(rowdata)[colnames(rowdata) == "SYMBOL"] <- "gene_ID"
rowdata <- rowdata[, c("gene_ID", 
                        colnames(rowdata)[!colnames(rowdata) == "gene_ID"])]
write.table(rowdata, file="inst/extdata/rowData.tsv", sep="\t", quote=FALSE)
