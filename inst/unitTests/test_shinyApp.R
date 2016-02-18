## 
## no unit test for shinyCircos since it is a shiny application
## 

## START unit test createOrderedSimMat
namesPrec <- rownames(binnedMSP)
dfNameGroup <- data.frame(group = unlist(lapply(strsplit(namesPrec, "_"), "[[", 1)), name = namesPrec)
dfNameGroupRT <- MetabolomicTools::orderNames(dfNameGroup = dfNameGroup, 
                            similarityMatrix = NULL, order = "retentionTime")
dfNameGroupRTMock <- dfNameGroupRT
colnames(dfNameGroupRTMock) <- c("group", "mockname")
orderedSimMat <- MetabolomicTools::createOrderedSimMat(dfNameGroupRT, similarityMat)

test_createOrderedSimMat <- function() {
    checkException(MetabolomicTools::createOrderedSimMat(dfNameGroupRTMock, similarityMat))
    checkException(MetabolomicTools::createOrderedSimMat(dfNameGroupRT))
    ## needs column in dfNameGroup with entries group_number_mz/rt
    checkException(MetabolomicTools::createOrderedSimMat(dfNameGroup, similarityMat))
    checkEquals(dim(orderedSimMat), dim(similarityMat))
    checkEquals(colnames(orderedSimMat), rownames(orderedSimMat))
    checkEquals(colnames(orderedSimMat), dfNameGroupRT$name)
    checkTrue(is.matrix(orderedSimMat))
    checkTrue(is.numeric(orderedSimMat))
}
## END unit test createOrderedSimMat
