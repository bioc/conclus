################################################################################
#######################  extdata/clusters_table.tsv  ###########################
################################################################################
##
## This table gives which cells has been assigned which cluster.
##
## This is the code to generate this table:

library("conclus")
load(system.file("extdata/scrFull.Rdat", package="conclus"))

clustCellsDf <- retrieveTableClustersCells(scr)

write.table(clustCellsDf, file="inst/extdata/clusters_table.tsv", sep="\t", 
            quote=FALSE)
