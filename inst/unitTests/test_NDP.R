## START unit test NDP
## create objects which will be used in unit tests
data("binnedMSP", package = "MetCirc")
## use only a selection 
binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]

## get masses
mass <- colnames(binnedMSP)

test_NDP <- function() {
    checkEquals(NDP(binnedMSP[1,], binnedMSP[1,], mass = mass), 1)
    checkEquals(NDP(binnedMSP[1,], binnedMSP[2,], mass = mass), 0.7980518, 
                tolerance = 0.000001)
    checkEquals(NDP(binnedMSP[2,], binnedMSP[1,], mass = mass), 0.7980518, 
                tolerance = 0.000001)    
    checkEquals(NDP(binnedMSP[1,], binnedMSP[2,], mass = mass), 
                NDP(binnedMSP[2,], binnedMSP[1,], mass = mass))
    checkException(NDP(binnedMSP[1,], binnedMSP[2,], mass = c(0:10)))
    checkException(NDP(binnedMSP[1,1:10], binnedMSP[2,], mass = mass))
}
## END unit test NDP 

## START unit test createSimilarityMatrix
simMat <- createSimilarityMatrix(binnedMSP)
test_createSimilarityMatrix <- function() {
    checkEquals(dim(simMat)[1], dim(binnedMSP)[1])
    checkEquals(dim(simMat)[2], dim(binnedMSP)[1])
    checkEquals(colnames(simMat), rownames(simMat))
    checkEquals(rownames(simMat), sort(rownames(binnedMSP)))
    checkTrue(is.numeric(simMat))
    checkTrue(is.matrix(simMat))
}
## END unit test createSimilarityMatrix