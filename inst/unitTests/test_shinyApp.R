
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

## START unit test shinyCircos
test_shinyCircos <- function() {
    checkException(shinyCircos(1:2, NULL))
    checkException(shinyCircos(similarityMat, "a"))
}
## END unit test shinyCircos


## START unit test printInformationSelect 
test_printInformationSelect <- function() {
    checkEquals(MetCirc:::printInformationSelect( 
                    groupname = groupname, msp = NULL, ind = ind, 
                    lMatInd = linkMatInds, 
                    linkMatrixThreshold = linkMat_thr, 
                    similarityMatrix = orderedSimMat),
        "ANT_1052.94/1018.15 connects to <br/> ANT_1063.44/1017.07 <br/>ANT_1398.71/1015.75 <br/>LIM_1398.71/1015.75")
    checkEquals(MetCirc:::printInformationSelect( 
        groupname = rownames(similarityMat), msp = finalMSP, ind = ind, 
        lMatInd = linkMatInds, linkMatrixThreshold = linkMat_thr, 
        similarityMatrix = similarityMat),
        "SPL_966.94/990.15 (Unknown, Unknown, Unknown, Unknown) connects to  <br/>ANT_1063.44/1017.07 (0.01, Unknown, Unknown, Unknown, Unknown)<br/> ANT_1398.71/1015.75 (0.002, Unknown, Unknown, Unknown, Unknown)<br/> LIM_1398.71/1015.75 (0.002, Unknown, Unknown, Unknown, Unknown)<br/> ANT_1052.94/1018.15 (0.043, Unknown, Unknown, Unknown, Unknown)<br/>")
    checkEquals(MetCirc:::printInformationSelect( 
        groupname = rownames(similarityMat), msp = finalMSP, ind = ind, 
        lMatInd = NULL, linkMatrixThreshold = linkMat_thr, 
        similarityMatrix = similarityMat),
        "SPL_966.94/990.15 (Unknown, Unknown, Unknown,Unknown) does not connect to any feature")
}
## END unit test printInformationSelect

## START unit test createOrderedSimMat
test_createOrderedSimMat <- function() {
    checkException(createOrderedSimMat(similarityMat, order = "foo"))
    checkException(createOrderedSimMat(order = "mz"))
    checkException(createOrderedSimMat(order = "retentionTime"))
    checkException(createOrderedSimMat(order = "clustering"))
    checkEquals(dim(simMat), dim(similarityMat))
    checkEquals(colnames(simMat), rownames(simMat))
    checkTrue(is.matrix(simMat))
    checkTrue(is.numeric(simMat))
}
## END unit test createOrderedSimMat