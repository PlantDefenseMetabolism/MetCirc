

## START unit test getBegEndIndMSP
testMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", 
                       splitIndMZ = 2, splitIndRT = NULL)
testMSPmsp <- getMSP(testMSP)
BegEndIndMSP <- getBegEndIndMSP(testMSPmsp)
test_getBegEndIndMSP <- function() {
    checkTrue(is.list(getBegEndIndMSP(testMSPmsp)))
    checkTrue(length(BegEndIndMSP[[1]]) == length(BegEndIndMSP[[2]]))
    checkTrue(all(BegEndIndMSP[[1]] <=  BegEndIndMSP[[2]]))
    checkTrue(is.numeric(BegEndIndMSP[[1]]))
    checkTrue(is.vector(BegEndIndMSP[[1]]))
    checkTrue(is.numeric(BegEndIndMSP[[2]]))
    checkTrue(is.vector(BegEndIndMSP[[2]]))
}
## END unit test getBegIndMSP

## START unit test binning
compartment <- c(rep("a", 90), rep("b", 90), rep("c", 90), rep("d", 90))
## create binnedMSPs
binnedMSP001 <- binning(testMSP, 0.01, group = compartment, method = "mean")
test_binning <- function() {
    checkTrue(is.matrix(binnedMSP001))
    checkTrue(is.numeric(binnedMSP001))
    checkEquals(dim(binnedMSP001), c(360, 764)) 
    checkException(binning(finalMSP, 1, compartment[1:7]))
}
## END unit test binning
