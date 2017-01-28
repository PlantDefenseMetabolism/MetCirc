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
linkMatInds <- getLinkMatrixIndices(groupname[ind], linkMat_thr)
## MetCirc:::printInformationSelect(groupname = groupname, 
##  msp = NULL, ind = ind, lMatInd = linkMatInds, 
##  linkMatrixThreshold = linkMat_thr, highlight = TRUE, 
##  similarityMatrix = simMat)


## START unit test printInformationSelect 
test_printInformationSelect <- function() {
    checkEquals(MetCirc:::printInformationSelect( 
                    groupname = groupname, msp = NULL, ind = ind, 
                    lMatInd = linkMatInds, 
                    linkMatrixThreshold = linkMat_thr, 
                    similarityMatrix = orderedSimMat),
        "ANT_1052.94/1018.15 connects to <br/> ANT_1063.44/1017.07 <br/>ANT_1398.71/1015.75 <br/>LIM_1398.71/1015.75")
}
## END unit test printInformationSelect

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