## START unit test getPrecursorMZ
finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)

test_getPrecursorMZ <- function() {
    checkEquals(length(MetabolomicTools:::getPrecursorMZ(finalMSP)), 360)
    checkTrue(is.vector(MetabolomicTools:::getPrecursorMZ(finalMSP)))
    checkTrue(is.numeric(MetabolomicTools:::getPrecursorMZ(finalMSP)))
    checkEquals(
        MetabolomicTools:::getPrecursorMZ(sd02_deconvoluted), numeric())
}
## END unit test getPrecursorMZ

## START unit test getRT
test_getRT <- function() {
    checkEquals(length(MetabolomicTools:::getRT(finalMSP)), 360)
    checkTrue(is.vector(MetabolomicTools:::getRT(finalMSP)))
    checkTrue(is.numeric(MetabolomicTools:::getRT(finalMSP)))
    checkEquals(
        MetabolomicTools:::getRT(sd02_deconvoluted), numeric())
}
## END unit test getRT

## START unit test getBegEndIndMSP
BegEndIndMSP <- MetabolomicTools:::getBegEndIndMSP(finalMSP)
test_getBegEndIndMSP <- function() {
    checkTrue(is.list(MetabolomicTools:::getBegEndIndMSP(finalMSP)))
    checkTrue(length(BegEndIndMSP[[1]]) == length(BegEndIndMSP[[2]]))
    checkTrue(all(BegEndIndMSP[[1]] <=  BegEndIndMSP[[2]]))
    checkTrue(is.numeric(BegEndIndMSP[[1]]))
    checkTrue(is.vector(BegEndIndMSP[[1]]))
    checkTrue(is.numeric(BegEndIndMSP[[2]]))
    checkTrue(is.vector(BegEndIndMSP[[2]]))
}
## END unit test getBegIndMSP

## START unit test binning
finalMSPCut <- finalMSP[1:94,]
compartment <- c("a", "a", "b", "b", "c", "c", "d", "d")

test_binning <- function() {
    checkTrue(
        is.matrix(MetabolomicTools::binning(finalMSPCut, 0.1, compartment)))
    checkTrue(
        is.numeric(MetabolomicTools::binning(finalMSPCut, 0.1, compartment)))
    checkEquals(
        dim(MetabolomicTools::binning(finalMSPCut, 0.01, compartment)), c(8,21)) 
    checkEquals(
        dim(MetabolomicTools::binning(finalMSPCut, 0.1, compartment)), c(8,20)) 
    checkEquals(
        dim(MetabolomicTools::binning(finalMSPCut, 1, compartment)), c(8,18)) 
    checkEquals(
        dim(MetabolomicTools::binning(finalMSPCut, 10, compartment)), c(8,10))
    checkEquals(
        dim(MetabolomicTools::binning(finalMSPCut, 100, compartment)), c(8,3))
    checkEquals(
        dim(MetabolomicTools::binning(finalMSPCut, 1000, compartment)), c(8,1))
    checkEquals(
        dim(MetabolomicTools::binning(finalMSPCut, 10000, compartment)), c(8,1))
    checkException(MetabolomicTools::binning(finalMSPCut, 1, compartment[1:7]))
}
## END unit test binning
