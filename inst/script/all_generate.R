## Script to generate all files in inst/extdata

source(file.path(system.file("script", package = "conclus"), 
                 "generate_countMatrix.R"))

source(file.path(system.file("script", package = "conclus"), 
                 "generate_colData.R"))

source(file.path(system.file("script", package = "conclus"), 
                 "generate_rowData.R"))

source(file.path(system.file("script", package = "conclus"), 
                 "generate_scrLight.R"))

source(file.path(system.file("script", package = "conclus"), 
                 "generate_scrFull.R"))

source(file.path(system.file("script", package = "conclus"), 
                 "generate_expected_normalizedMatrix.R")) 
