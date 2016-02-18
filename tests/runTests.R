library(amap)
library(circlize)
data("sd01_outputXCMS", package = "MetabolomicTools")
data("sd02_deconvoluted", package = "MetabolomicTools")
data("binnedMSP", package = "MetabolomicTools")
data("similarityMat", package = "MetabolomicTools")

BiocGenerics:::testPackage("MetabolomicTools")