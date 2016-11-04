## 
## no unit test for shinyCircos since it is a shiny application
## 

## create objects which will be used in unit tests
data("idMSMStoMSP", package = "MetCirc")
data("binnedMSP", package = "MetCirc")
## use only a selection
binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
similarityMat <- createSimilarityMatrix(binnedMSP)
## order similarityMat according to mz
simMat <- createOrderedSimMat(similarityMat, order = "mz")
groupname <- rownames(simMat)
linkMat_thr <- createLinkMatrix(simMat, 0.75, 1) 
ind <- 18
linkMatIndsHover <- getLinkMatrixIndices(groupname[ind], linkMat_thr)
## MetCirc:::printInformationHover(groupname = groupname, 
##  msp = NULL, ind = ind, lMatIndHover = linkMatIndsHover, 
##  linkMatrixThreshold = linkMat_thr, highlight = TRUE, 
##  similarityMatrix = simMat)


## START unit test printInformationHover 
test_printInformationHover <- function() {
    checkEquals(MetCirc:::printInformationHover( 
                    groupname = groupname, msp = NULL, ind = ind, 
                    lMatIndHover = linkMatIndsHover, 
                    linkMatrixThreshold = linkMat_thr, 
                    similarityMatrix = orderedSimMat),
        "ANT_0018_1052.939378/1020.95636760976 connects to <br/> ANT_0019_1063.440342/1020.17226498551")
}
## END unit test printInformationHover

## START unit test createOrderedSimMat
test_createOrderedSimMat <- function() {
    checkException(createOrderedSimMat(similarityMat, order = "foo"))
    checkException(createOrderedSimMat(order = "mz"))
    checkEquals(dim(simMat), dim(similarityMat))
    checkEquals(colnames(simMat), rownames(simMat))
    checkTrue(is.matrix(simMat))
    checkTrue(is.numeric(simMat))
}
## END unit test createOrderedSimMat