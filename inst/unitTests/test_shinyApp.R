## 
## no unit test for shinyCircos since it is a shiny application
## 

## START unit test createOrderedSimMat
## create objects which will be used in unit tests
data("binnedMSP", package = "MetCirc")
## use only a selection 
binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
similarityMat <- createSimilarityMatrix(binnedMSP)  
namesPrec <- rownames(binnedMSP)

## create dfNameGroup
dfNameGroup <- data.frame(
    group = unlist(lapply(strsplit(namesPrec, "_"), "[[", 1)), name = namesPrec)
dfNameGroupRT <- orderNames(dfNameGroup = dfNameGroup, 
    similarityMatrix = NULL, order = "retentionTime")
dfNameGroupRTMock <- dfNameGroupRT
colnames(dfNameGroupRTMock) <- c("group", "mockname")
orderedSimMat <- createOrderedSimMat(dfNameGroupRT, similarityMat)

test_createOrderedSimMat <- function() {
    checkException(createOrderedSimMat(dfNameGroupRTMock, similarityMat))
    checkException(createOrderedSimMat(dfNameGroupRT))
    ## needs column in dfNameGroup with entries group_number_mz/rt
    checkException(createOrderedSimMat(dfNameGroup, similarityMat))
    checkEquals(dim(orderedSimMat), dim(similarityMat))
    checkEquals(colnames(orderedSimMat), rownames(orderedSimMat))
    checkEquals(colnames(orderedSimMat), dfNameGroupRT$name)
    checkTrue(is.matrix(orderedSimMat))
    checkTrue(is.numeric(orderedSimMat))
}
## END unit test createOrderedSimMat
