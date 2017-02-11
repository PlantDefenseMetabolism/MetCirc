
MetCirc: Navigating mass spectral similarity in high-resolution MS/MS metabolomics data
===========

One of the main problems in MS/MS metabolomics is the rapid dereplication of 
previously characterized metabolites across a range of biological samples and 
the structural prediction of unknowns from MS/MS data. `MetCirc` aims to 
faciliate these steps by offering functionalities to display, 
(interactively) explore similarities and annotate features of MS/MS metabolomics 
data. The `R` package is especially designed to 
improve the interactive exploration of metabolomics data obtained from 
cross-species/cross-tissues comparative experiments. Notably, `MetCirc`
includes functions to calculate the similarity between individual MS/MS 
spectra based on a normalised dot product calculation taking into account 
shared fragments or main neutral losses.


To install `MetCirc`from this GitHub page enter:

`library(devtools)` 

`install_github(repo = "PlantDefenseMetabolism/MetCirc")`


`MetCirc` is also available via the Bioconductor framework. To install 
`MetCirc` from Bioconductor enter: 

`source("https://bioconductor.org/biocLite.R")`

`biocLite("MetCirc")`
