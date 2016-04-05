library(amap)
library(circlize)
library(scales)
library(shiny)

data("sd01_outputXCMS", package = "MetCirc")
data("sd02_deconvoluted", package = "MetCirc")
data("binnedMSP", package = "MetCirc")
data("similarityMat", package = "MetCirc")

BiocGenerics:::testPackage("MetCirc")