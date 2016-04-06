## START unit test getPrecursorMZ
finalMSPmsp <- getMSP(convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2))

test_getPrecursorMZ <- function() {
    checkEquals(length(getPrecursorMZ(finalMSPmsp)), 360)
    checkTrue(is.vector(getPrecursorMZ(finalMSPmsp)))
    checkTrue(is.numeric(getPrecursorMZ(finalMSPmsp)))
    checkEquals(getPrecursorMZ(sd02_deconvoluted), numeric())
}
## END unit test getPrecursorMZ

## START unit test getRT
test_getRT <- function() {
    checkEquals(length(getRT(finalMSPmsp)), 360)
    checkTrue(is.vector(getRT(finalMSPmsp)))
    checkTrue(is.numeric(getRT(finalMSPmsp)))
    checkEquals(
        getRT(sd02_deconvoluted), numeric())
}
## END unit test getRT

## START unit test getBegEndIndMSP
BegEndIndMSP <- getBegEndIndMSP(finalMSPmsp)
test_getBegEndIndMSP <- function() {
    checkTrue(is.list(getBegEndIndMSP(finalMSPmsp)))
    checkTrue(length(BegEndIndMSP[[1]]) == length(BegEndIndMSP[[2]]))
    checkTrue(all(BegEndIndMSP[[1]] <=  BegEndIndMSP[[2]]))
    checkTrue(is.numeric(BegEndIndMSP[[1]]))
    checkTrue(is.vector(BegEndIndMSP[[1]]))
    checkTrue(is.numeric(BegEndIndMSP[[2]]))
    checkTrue(is.vector(BegEndIndMSP[[2]]))
}
## END unit test getBegIndMSP

## START unit test binning
finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
compartment <- c(rep("a", 90), rep("b", 90), rep("c", 90), rep("d", 90))
## create binnedMSPs
binnedMSP001 <- binning(finalMSP, 0.01, group = compartment)
binnedMSP1 <- binning(finalMSP, 1, group = compartment)

test_binning <- function() {
    checkTrue(is.matrix(binnedMSP001))
    checkTrue(is.numeric(binnedMSP001))
    checkEquals(dim(binnedMSP001), c(360, 691)) 
    checkEquals(dim(binnedMSP1), c(360,351))
    checkException(binning(finalMSP, 1, compartment[1:7]))
}
## END unit test binning
