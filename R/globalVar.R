utils::globalVariables(c(
                ## function .getTSNEresults in foreach
                "PCAGetTSNEresults", "perpGetTSNEresults",

                ## function .mkDbscan in foreach
                "iMkDbscan", "epsMkDbscan", "MinPtsMkDbscan",

                ## function .mkSimMat in foreach
                "iMkSimMat",

                ## function retrieveGenesInfo from dplyr operations
                "uniprot_gn_symbol", "chromosome_name",
                "entrezgene_description", "go_id", "uniprot_gn_id",
                "description", "external_gene_name", "gene_biotype",
                "ensembl_gene_id", "entrezgene_id"))
