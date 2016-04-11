## create objects which will be used in unit tests
data("binnedMSP", package = "MetCirc")
## use only a selection 
binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
similarityMat <- createSimilarityMatrix(binnedMSP)  
namesPrec <- rownames(binnedMSP)

## create dfNameGroup
dfNameGroup <- data.frame(
    group = unlist(lapply(strsplit(namesPrec, "_"), "[[", 1)), name = namesPrec)
## order according to compartment
dfNameGroup <- dfNameGroup[order(dfNameGroup[,"group"]),] 
dfNameGroupRT <- orderNames(dfNameGroup = dfNameGroup, 
                            similarityMatrix = NULL, order = "retentionTime")
linkMat <- createLinkMatrix(similarityMatrix = similarityMat, 
                threshold=0.95, dfNameGroup = dfNameGroup)

## START unit test for plotCircos
circos.clear()
circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), 
           track.margin = c(0.0, 0))
test_plotCircos <- function() {
    checkException(plotCircos(dfNameGroup, NULL, 
        initialize = TRUE, featureNames = FALSE, groupSector = FALSE, 
        groupName = FALSE, links = TRUE, highlight = FALSE))
    checkException(plotCircos(dfNameGroupRT, NULL, 
        initialize = FALSE, featureNames = TRUE, groupSector = FALSE, 
        groupName = FALSE, links = FALSE, highlight = FALSE))
    checkException(plotCircos(dfNameGroupRT, linkMat, initialize = TRUE, 
        featureNames = FALSE, groupSector = FALSE, groupName = FALSE, 
        links = TRUE, highlight = FALSE)) ## names are different
}
## END unit test for plotCircos


## START unit test for highlight
test_highlight <- function() {
    checkException(highlight(dfNameGroup, 1, NULL))
    checkException(highlight(dfNameGroup, dim(dfNameGroupRT)[1]+1,NULL))
    checkException(highlight(dfNameGroupRT, 1, linkMat))
}
## END unit test for highlight

## START unit test for circosLegend

## END unit test for circosLegend
test_circosLegend <- function() {
    checkException(circosLegend(dfNameGroup[,1], TRUE))
}

## START unit test for getLinkMatrixIndices
circos.clear()
## set circlize paramters
circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), 
           track.margin = c(0.0, 0))
plotCircos(dfNameGroup, NULL, initialize = TRUE, 
    featureNames = FALSE, groupSector = FALSE, groupName = FALSE, links = FALSE, highlight = FALSE)
test_getLinkMatrixIndices <- function() {
    checkEquals(getLinkMatrixIndices(dfNameGroup[1,], linkMat), 13)
    checkEquals(getLinkMatrixIndices(dfNameGroup[2,], linkMat), integer())
    checkEquals(getLinkMatrixIndices(dfNameGroup[3,], linkMat), 14)
    checkEquals(getLinkMatrixIndices(dfNameGroup[4,], linkMat), 15)
    checkEquals(getLinkMatrixIndices(dfNameGroup[5,], linkMat), 16)
    checkEquals(getLinkMatrixIndices(dfNameGroup[1:5,], linkMat), c(13:16))
    checkException(getLinkMatrixIndices(dfNameGroup[1,], NULL))
}
## END unit test for getLinkMatrixIndices

## START unit test for truncateName
test_truncateName <- function() {
    checkEquals(
        as.character(truncateName(dfNameGroupRT[1,], 2, nameGroup=TRUE)[2]), 
        "1398.71/1018.98")
    checkEquals(
        as.character(truncateName(dfNameGroup[1,], 2, nameGroup=FALSE)[2]), 
        "231.05/1020.69")
    checkEquals(dim(truncateName(dfNameGroup, 2, nameGroup = FALSE)), 
                dim(dfNameGroup))
    checkEquals(as.character(truncateName(dfNameGroup[1,], 2, nameGroup = TRUE)[2]),
                "NA/NA")
    checkEquals(
        as.character(truncateName(dfNameGroupRT[1,], 2, nameGroup = FALSE)[2]),
        "1/NA")
}
## END unit test truncateName

## START unit test minFragCart2Polar
degreeFeatures <- lapply(dfNameGroup$name, 
    function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
test_minFragCart2Polar <- function() {
    checkEquals(minFragCart2Polar(1,0,degreeFeatures), 74)
    checkEquals(minFragCart2Polar(0.1,0.9,degreeFeatures), 77)
    checkEquals(minFragCart2Polar(1,1, degreeFeatures), NA)
    checkEquals(minFragCart2Polar(1, 0, NULL), integer())
    checkException(minFragCart2Polar(NA, NA, degreeFeatures))
    checkException(minFragCart2Polar(1, NA, degreeFeatures))
    checkException(minFragCart2Polar(NA, 1, degreeFeatures))
}
## END unit test minFragCart2Polar


## START unit test cart2Polar
test_cart2Polar <- function() {
    checkEquals(cart2Polar(0, 0), list(r = 0, theta = 0))
    checkEquals(cart2Polar(1, 1), list(r = 1.414214, theta = 45), 
        tolerance = 0.00001)
    checkEquals(cart2Polar(0, 1), list(r = 1, theta = 90))
    checkEquals(cart2Polar(-1, 1), list(r = 1.414214, theta = 135),
        tolerance = 0.00001)
    checkEquals(cart2Polar(-1, -1), list(r = 1.414214, theta = 225), 
        tolerance = 0.00001)
    checkEquals(cart2Polar(1, -1), list(r = 1.414214, theta = 315), 
        tolerance = 0.00001)
    checkException(cart2Polar(NA, NA))
    checkException(cart2Polar(1, NA))
    checkException(cart2Polar(NA, 1))
}
## END unit test cart2Polar

## START unit test orderNames
test_orderNames <- function() {
    checkException(orderNames(dfNameGroup, NULL, "neworder"))
    checkException(orderNames(NULL, NULL, "retentionTime"))
    checkException(orderNames(dfNameGroup, NULL, "clustering"))
    checkEquals(dim(orderNames(dfNameGroup, NULL, "retentionTime")), 
                dim(dfNameGroup))
    checkTrue(is.data.frame(orderNames(dfNameGroup, NULL, "retentionTime")))
    checkTrue(is.data.frame(orderNames(dfNameGroup, NULL, "mz")))
    checkTrue(
        is.data.frame(orderNames(dfNameGroup, similarityMat, "clustering")))
}
## END unit test orderNames
