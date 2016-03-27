## START unit test cutUniquePreMZ
test_cutUniquePreMZ <- function() {
    checkTrue(is.vector(cutUniquePreMZ(sd02_deconvoluted[,4],
        splitPattern = " _ ", splitInd = 2, returnCharacter = TRUE)))
    checkEquals(
        length(cutUniquePreMZ(sd02_deconvoluted[,4], 
            splitPattern = " _ ", splitInd = 2, returnCharacter = TRUE)), 360)
    checkTrue(
        is.character(cutUniquePreMZ(sd02_deconvoluted[,4], 
            splitPattern = " _ ", splitInd = 2, returnCharacter = TRUE)))
    checkTrue(
        is.numeric(cutUniquePreMZ(sd02_deconvoluted[,4], 
            splitPattern = " _ ", splitInd = 2, returnCharacter = FALSE)))
}
## END unit test cutUniquePreMZ

## START unit test convert2MSP
testMSP <- convert2MSP(sd02_deconvoluted, 
    splitPattern = " _ ", splitInd = 2)
test_convert2MSP <- function() {
    checkTrue(class(testMSP) == "MSP")
    checkEquals(length(testMSP), 360)
    checkTrue(is.data.frame(getMSP(testMSP)))
    checkEquals(dim(getMSP(testMSP)), c(7263, 2))
    checkEquals(as.numeric(
        table(getMSP(testMSP)[,1])["ADDUCTIONNAME: "]), 360)
    checkEquals(as.numeric(
        table(getMSP(testMSP)[,1])["METABOLITENAME: "]), 360)
    checkEquals(as.numeric(
        table(getMSP(testMSP)[,1])["NAME: "]), 360)
    checkEquals(as.numeric(
        table(getMSP(testMSP)[,1])["Num Peaks: "]), 360)
    checkEquals(as.numeric(
        table(getMSP(testMSP)[,1])["PRECURSORMZ: "]), 360)
    checkEquals(as.numeric(
        table(getMSP(testMSP)[,1])["RETENTIONTIME: "]), 360)
}
## END unit test convert2MSP

## START unit test msp2FunctionalLossesMSP
testMSPNL <- msp2FunctionalLossesMSP(testMSP)
test_msp2FunctionalLossesMSP <- function() {
    checkTrue(class(testMSPNL) == "MSP")
    checkEquals(length(testMSPNL), 360)
    checkTrue(is.data.frame(testMSPNL@msp))
    checkEquals(dim(testMSPNL@msp), c(7263, 2))
    checkTrue(is.data.frame(testMSPNL@msp))
    checkEquals(as.numeric(
        table(getMSP(testMSPNL)[,1])["ADDUCTIONNAME: "]), 360)
    checkEquals(as.numeric(
        table(getMSP(testMSPNL)[,1])["METABOLITENAME: "]), 360)
    checkEquals(as.numeric(
        table(getMSP(testMSPNL)[,1])["NAME: "]), 360)
    checkEquals(as.numeric(
        table(getMSP(testMSPNL)[,1])["Num Losses: "]), 360)
    checkEquals(as.numeric(
        table(getMSP(testMSPNL)[,1])["PRECURSORMZ: "]), 360)
    checkEquals(as.numeric(
        table(getMSP(testMSPNL)[,1])["RETENTIONTIME: "]), 360)
}
## END unit test msp2FunctionalLossesMSP