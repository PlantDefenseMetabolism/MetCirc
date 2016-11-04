## START unit test createLink0Matrix
## create objects which will be used in unit tests
data("binnedMSP", package = "MetCirc")
## use only a selection 
binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
similarityMat <- createSimilarityMatrix(binnedMSP)  
groupname <- rownames(similarityMat)
group <- unlist(lapply(strsplit(groupname, "_"), "[", 1))

## create link0Matrix
link0Matrix <- createLink0Matrix(similarityMatrix = similarityMat)
ndps <- as.numeric(link0Matrix[,"NDP"])


test_createLink0Matrix <- function() {
    checkEquals(dim(link0Matrix)[2], 5)
    checkTrue(is.matrix(link0Matrix))
    checkTrue(all(
        colnames(link0Matrix) == c("group1", "name1", "group2", "name2", "NDP")))
    checkTrue(
        all(unique(c(link0Matrix[,"group1"], link0Matrix[,"group2"])) %in% unique(group)))
    checkTrue(
        all(unique(c(link0Matrix[,"name1"], link0Matrix[,"name2"])) %in% unique(groupname)))
    checkTrue(all(0 < ndps & ndps <= 1))
}
## END unit test link0Matrix

## START unit test thresholdLinkMatrix
test_thresholdLinkMatrix <- function() {
    checkEquals(dim(thresholdLinkMatrix(link0Matrix, 0, 1)), dim(link0Matrix))
    checkException(thresholdLinkMatrix(similarityMat, 0, 1))
    checkException(thresholdLinkMatrix(similarityMat, 0.6, 0.5))
    checkException(thresholdLinkMatrix(link0Matrix, 1.05, 1.1))
    checkTrue(
        dim(thresholdLinkMatrix(link0Matrix, 0.2, 1))[1] >= 
            dim(thresholdLinkMatrix(link0Matrix, 0.3, 1))[1])
}
## END unit test thresholdLinkMatrix

## START unit test createLinkMatrix
tLinkMatrix1 <- thresholdLinkMatrix(link0Matrix, 0.9, 1)
tLinkMatrix2 <- createLinkMatrix(similarityMat, 0.9, 1)

test_createLinkMatrix <- function() {
    checkTrue(identical(tLinkMatrix1, tLinkMatrix2))
}
## END unit test createLinkMatrix

## START unit test cutLinkMatrix
cutLMInter <- cutLinkMatrix(tLinkMatrix1, type = "inter")
cutLMIntra <- cutLinkMatrix(tLinkMatrix1, type = "intra")

test_cutLinkMatrix <- function() {
    checkTrue(
        all(dim(cutLinkMatrix(tLinkMatrix1, type = "all")) == dim(tLinkMatrix2))
    )
    checkException(cutLinkMatrix(tLinkMatrix1, type = "foo"))
    checkTrue(
        all(unlist(lapply(1:dim(cutLMInter)[1], 
                          function(x) cutLMInter[x,1] != cutLMInter[x,3])))
    )
    checkTrue(
        all(unlist(lapply(1:dim(cutLMIntra)[1], 
                          function(x) cutLMIntra[x,1] == cutLMIntra[x,3])))
    )
}
## END unit test cutLinkMatrix
