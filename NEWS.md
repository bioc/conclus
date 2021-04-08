Changes in version 0.99.338 (2021-04-08)
----------------------------------------

+ Submission to Bioconductor

+ NAMESPACE
  - Imported BiocFileCache package and R_user_dir() of tools package
  - Exported new function conclusCacheClear()
  
+ DESCRIPTION
  - Removed LazyData: true.
  
+ DataFormatting.R
  - Added a caching system for retrieveFromGEO()
  - Created conclusCacheClear() to delete the cache
  - Updated documentation
  
+ loadDataset.R
  - Updated documentation
  - Simplified nested "if" in loadDataset.R
  
+ methods-normalization.R
  - Modified the use of getBM to retrieve only genes of the count matrix (instead of the all database)
  
+ test_setters.R/test_getters.R
  - Modified documentation for loadDataOrMatrix()
  - Simplified nested "if" in loadDataset.R
  
+ test_loadData.R
  - Adapted the unit tests to the new format of coldata and rowdata
  
+ test_scRNAseq-methods.R
  - Changed the experiment name to "Light_Experience"
  - Used tempdir() for output directory
  
+ inst
  - Added inst/script to generate data on inst/extdata
  - New data generated
  
+ vignette
  - Used tempdir() for output directory
  - Specified other parameters in the first example of runCONCLUS for the very small dataset
  - Replaced old paths by new ones
  
