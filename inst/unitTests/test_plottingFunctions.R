## create objects which will be used in unit tests
namesPrec <- rownames(binnedMSP)
dfNameGroup <- data.frame(
    group = unlist(lapply(strsplit(namesPrec, "_"), "[[", 1)), name = namesPrec)
## order according to compartment
dfNameGroup <- dfNameGroup[order(dfNameGroup[,"group"]),] 
dfNameGroupRT <- MetabolomicTools::orderNames(dfNameGroup = dfNameGroup, 
                            similarityMatrix = NULL, order = "retentionTime")
linkMat <- MetabolomicTools::createLinkMatrix(similarityMatrix = similarityMat, 
                threshold=0.95, dfNameGroup = dfNameGroup)

## START unit test for plotCircos
circos.clear()
circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), 
           track.margin = c(0.0, 0))
test_plotCircos <- function() {
    checkException(MetabolomicTools::plotCircos(dfNameGroup, NULL, 
        initialize = TRUE, featureNames = FALSE, groupName = FALSE, 
        links = TRUE, highlight = FALSE))
    checkException(MetabolomicTools::plotCircos(dfNameGroupRT, NULL, 
        initialize = FALSE, featureNames = TRUE, groupName = FALSE,
                    links = FALSE, highlight = FALSE))
    checkException(MetabolomicTools::plotCircos(dfNameGroupRT, linkMat, 
        initialize = TRUE, featureNames = FALSE, groupName = FALSE, 
        links = TRUE, highlight = FALSE)) ## names are different
}
## END unit test for plotCircos


## START unit test for highlight
test_highlight <- function() {
    checkException(MetabolomicTools::highlight(dfNameGroup, 1, NULL))
    checkException(MetabolomicTools::highlight(dfNameGroup, dim(dfNameGroupRT)[1]+1,NULL))
    checkException(MetabolomicTools::highlight(dfNameGroupRT, 1, linkMat))
}
## END unit test for highlight

## START unit test for getLinkMatrixIndices
circos.clear()
## set circlize paramters
circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), 
           track.margin = c(0.0, 0))
MetabolomicTools::plotCircos(dfNameGroup, NULL, initialize = TRUE, 
    featureNames = FALSE, groupName = FALSE, links = FALSE, highlight = FALSE)
test_getLinkMatrixIndices <- function() {
    checkEquals(MetabolomicTools:::getLinkMatrixIndices(dfNameGroup[1,], linkMat), 1)
    checkEquals(MetabolomicTools:::getLinkMatrixIndices(dfNameGroup[2,], linkMat), integer())
    checkEquals(MetabolomicTools:::getLinkMatrixIndices(dfNameGroup[3,], linkMat), 2)
    checkEquals(MetabolomicTools:::getLinkMatrixIndices(dfNameGroup[4,], linkMat), integer())
    checkEquals(MetabolomicTools:::getLinkMatrixIndices(dfNameGroup[5,], linkMat), integer())
    checkEquals(MetabolomicTools:::getLinkMatrixIndices(dfNameGroup[1:5,], linkMat), c(1,2))
    checkException(MetabolomicTools:::getLinkMatrixIndices(dfNameGroup[1,], NULL))
}
## END unit test for getLinkMatrixIndices

## START unit test for truncateName
test_truncateName <- function() {
    checkEquals(
        as.character(MetabolomicTools:::truncateName(dfNameGroupRT[1,], 2, nameGroup=TRUE)[2]), 
        "930.49/36.47")
    checkEquals(
        as.character(MetabolomicTools:::truncateName(dfNameGroup[1,], 2, nameGroup=FALSE)[2]), 
        "132.08/74.81")
    checkEquals(dim(MetabolomicTools:::truncateName(dfNameGroup, 2, nameGroup = FALSE)), 
                dim(dfNameGroup))
    checkEquals(as.character(MetabolomicTools:::truncateName(dfNameGroup[1,], 2, nameGroup = TRUE)[2]),
                "NA/NA")
    checkEquals(
        as.character(MetabolomicTools:::truncateName(dfNameGroupRT[1,], 2, nameGroup = FALSE)[2]),
        "1/NA")
}
## END unit test truncateName

## START unit test minFragCart2Polar
degreeFeatures <- lapply(dfNameGroup$name, 
    function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
test_minFragCart2Polar <- function() {
    checkEquals(MetabolomicTools:::minFragCart2Polar(1,0,degreeFeatures), 358)
    checkEquals(
        MetabolomicTools:::minFragCart2Polar(0.1,0.9,degreeFeatures), 266)
    checkEquals(MetabolomicTools:::minFragCart2Polar(1,1, degreeFeatures), NA)
    checkEquals(MetabolomicTools:::minFragCart2Polar(1, 0, NULL), integer())
    checkException(MetabolomicTools:::minFragCart2Polar(NA, NA, degreeFeatures))
    checkException(MetabolomicTools:::minFragCart2Polar(1, NA, degreeFeatures))
    checkException(MetabolomicTools:::minFragCart2Polar(NA, 1, degreeFeatures))
}
## END unit test minFragCart2Polar


## START unit test cart2Polar
test_cart2Polar <- function() {
    checkEquals(MetabolomicTools:::cart2Polar(0, 0), list(r = 0, theta = 0))
    checkEquals(
        MetabolomicTools:::cart2Polar(1, 1), list(r = 1.414214, theta = 45), 
        tolerance = 0.00001)
    checkEquals(
        MetabolomicTools:::cart2Polar(0, 1), list(r = 1, theta = 90))
    checkEquals(
        MetabolomicTools:::cart2Polar(-1, 1), list(r = 1.414214, theta = 135),
        tolerance = 0.00001)
    checkEquals(
        MetabolomicTools:::cart2Polar(-1, -1), list(r = 1.414214, theta = 225),
        tolerance = 0.00001)
    checkEquals(
        MetabolomicTools:::cart2Polar(1, -1), list(r = 1.414214, theta = 315),
        tolerance = 0.00001)
    checkException(MetabolomicTools:::cart2Polar(NA, NA))
    checkException(MetabolomicTools:::cart2Polar(1, NA))
    checkException(MetabolomicTools:::cart2Polar(NA, 1))
}
## END unit test cart2Polar

## START unit test orderNames
test_orderNames <- function() {
    checkException(MetabolomicTools::orderNames(dfNameGroup, NULL, "neworder"))
    checkException(MetabolomicTools::orderNames(NULL, NULL, "retentionTime"))
    checkException(MetabolomicTools::orderNames(dfNameGroup, NULL, "clustering"))
    checkEquals(dim(MetabolomicTools::orderNames(dfNameGroup, NULL, "retentionTime")), 
                dim(dfNameGroup))
    checkTrue(is.data.frame(MetabolomicTools::orderNames(dfNameGroup, NULL, "retentionTime")))
    checkTrue(is.data.frame(MetabolomicTools::orderNames(dfNameGroup, NULL, "mz")))
    checkTrue(
        is.data.frame(MetabolomicTools::orderNames(dfNameGroup, similarityMat, "clustering")))
}
## END unit test orderNames
