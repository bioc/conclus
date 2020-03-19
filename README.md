# Conclus with S4 class

All classes and methods are in R repository. <br>
The data are in inst/extdata.

To run the code, run main.R

### 21/02/2020

- Write scRNAseq class
- Write getters & setters methods, and generics for scRNAseq cmethods
- Transform conclus::normaliseCountMatrix to normaliseCountMatrix method in
methods-scRNAseq-normalization.R :
    * Coding style for Bioconductor
    * Modify 2 lines for accessing to data without system.file function
- Write main.R for testing
    * Load all libraires to compensate the missing of NAMESPACE FILE

### 27/02/2020

- Write tsne methods for scRNAseq class with methods:
    * getTSNEresults
    * generateTSNECoordinates
- Write accessors for Tsne class.

### 04/03/2020

- Write dbscan methods for scRNAseq class with :
    * runDBSCAN

### 18/02/2020

- Finish to write dbscan methods for scRNAseq class
- Begin to write similarity method of scRNAseq class