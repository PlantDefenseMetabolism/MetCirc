## create objects which will be used in unit tests
data("binnedMSP", package = "MetCirc")
## use only a selection 
binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
similarityMat <- createSimilarityMatrix(binnedMSP)
simMat <- createOrderedSimMat(similarityMat, order = "mz")
groupname <- rownames(similarityMat)
groupnameO <- rownames(simMat)
## create link mat
linkMat <- createLinkMatrix(similarityMatrix = similarityMat, 
                            threshold_low = 0.95, threshold_high = 1)

## START unit test for plotCircos
circos.clear()
circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), 
           track.margin = c(0.0, 0))
test_plotCircos <- function() {
    checkException(plotCircos(groupname, NULL, 
        initialize = TRUE, featureNames = FALSE, groupSector = FALSE, 
        groupName = FALSE, links = TRUE, highlight = FALSE, colour = NULL, 
        transparency = 0.2))
    checkException(plotCircos(groupnameO, linkMat, initialize = TRUE, 
        featureNames = FALSE, groupSector = FALSE, groupName = FALSE, 
        links = TRUE, highlight = FALSE, colour = NULL, transparency = 0.2)) ## names are different
    checkException(plotCircos(featureNames = TRUE))
    checkException(plotCircos(groupname, linkMat, cexFeatureNames = 0.3,
            initialize = TRUE, featureNames = FALSE, groupSector = FALSE, 
            groupName = FALSE, links = TRUE, highlight = FALSE, 
            transparency = 0.2, colour = c("red", "blue")))
    checkException(plotCircos(groupname, linkMat, cexFeatureNames = 0.3,
                initialize = TRUE, featureNames = FALSE, groupSector = TRUE, 
                groupName = FALSE, links = TRUE, highlight = FALSE, 
                colour = c("red", "blue"), transparency = 0.2))
}
## END unit test for plotCircos


## START unit test for highlight
test_highlight <- function() {
    checkException(highlight(groupnameO, 1, NULL, NULL, 0.4))
    checkException(highlight(groupnameO, length(groupnameO)+1,NULL, NULL, 0.4))
    checkException(highlight(groupnameO, length(groupnameO)+1,NULL, NULL, 0.4,
                             colour = "red"))
    checkException(highlight(groupname, 
        LinkMatrix = cbind("name1" = 0, "name2" = 0), ind = 1))
    checkException(highlight(groupname, 
        LinkMatrix = cbind("name1" = c(0, 2), "name2" = c(2, 0)), ind = 1))
    checkException(highlight(groupname, LinkMatrix = linkMat, links = FALSE, ind = 100))
    checkException(highlight(groupname, LinkMatrix = linkMat, links = TRUE, ind = 100))
    ## names in linkMat do not match names in groupnameO
    checkException(highlight(groupnameO, 1, linkMat, NULL, 0.4))
}
## END unit test for highlight

## START unit test for circosLegend
test_circosLegend <- function() {
    checkException(circosLegend(1, highlight = TRUE, colour = NULL))
    checkException(circosLegend(1, highlight = TRUE, colour = "red"))
    checkException(circosLegend(1, highlight = FALSE, colour = NULL))
    checkException(circosLegend(1, highlight = FALSE, colour = "red"))
}
## END unit test for circosLegend


## START unit test for getLinkMatrixIndices
circos.clear()
## set circlize paramters
circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), 
           track.margin = c(0.0, 0))
plotCircos(groupname, NULL, initialize = TRUE, 
    featureNames = FALSE, groupSector = FALSE, groupName = FALSE, links = FALSE, highlight = FALSE)
test_getLinkMatrixIndices <- function() {
    checkEquals(getLinkMatrixIndices(groupname[4], linkMat), numeric())
    checkEquals(getLinkMatrixIndices(groupname[5], linkMat), 1:3)
    checkEquals(getLinkMatrixIndices(groupname[6], linkMat), c(1, 4:5))
    checkEquals(getLinkMatrixIndices(groupname[7], linkMat), c(2, 4, 6))
    checkEquals(getLinkMatrixIndices(groupname[8], linkMat), c(3, 5, 6))
    checkEquals(getLinkMatrixIndices(groupname[4:6], linkMat), c(1:3, 1, 4, 5))
    checkException(getLinkMatrixIndices(groupname[1], NULL))
}
## END unit test for getLinkMatrixIndices

## START unit test for truncateName
test_truncateName <- function() {
    checkEquals(truncateName(groupname[1], 2), "579.35/991.14")
    checkEquals(truncateName(groupname[1], 0), "579/991")
    checkEquals(length(truncateName(groupname[1:10], 2)), 10)
}
## END unit test truncateName

## START unit test minFragCart2Polar

degreeFeatures <- lapply(groupname, 
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
