## START unit test NDP
mass <- colnames(binnedMSP)

test_NDP <- function() {
    checkEquals(MetabolomicTools::NDP(binnedMSP[1,], binnedMSP[1,], mass = mass), 1)
    checkEquals(MetabolomicTools::NDP(binnedMSP[1,], binnedMSP[2,], mass = mass), 0.01356902, 
                tolerance = 0.000001)
    checkEquals(MetabolomicTools::NDP(binnedMSP[2,], binnedMSP[1,], mass = mass), 0.01356902, 
                tolerance = 0.000001)    
    checkEquals(MetabolomicTools::NDP(binnedMSP[1,], binnedMSP[2,], mass = mass), 
                MetabolomicTools::NDP(binnedMSP[2,], binnedMSP[1,], mass = mass))
    checkException(MetabolomicTools::NDP(binnedMSP[1,], binnedMSP[2,], mass = c(0:10)))
    checkException(MetabolomicTools::NDP(binnedMSP[1,1:10], binnedMSP[2,], mass = mass))
}
## END unit test NDP 

## START unit test createSimilarityMatrix
simMat <- MetabolomicTools::createSimilarityMatrix(binnedMSP)
test_createSimilarityMatrix <- function() {
    checkEquals(dim(simMat)[1], dim(binnedMSP)[1])
    checkEquals(dim(simMat)[2], dim(binnedMSP)[1])
    checkEquals(colnames(simMat), rownames(simMat))
    checkEquals(rownames(simMat), rownames(binnedMSP))
    checkTrue(is.numeric(simMat))
    checkTrue(is.matrix(simMat))
}
## END unit test createSimilarityMatrix