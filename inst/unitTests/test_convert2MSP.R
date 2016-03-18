## START unit test cutUniquePreMZ
test_cutUniquePreMZ <- function() {
    checkTrue(is.vector(MetabolomicTools:::cutUniquePreMZ(sd02_deconvoluted[,4],
        splitPattern = " _ ", splitInd = 2, returnCharacter = TRUE)))
    checkEquals(
        length(MetabolomicTools:::cutUniquePreMZ(sd02_deconvoluted[,4], 
            splitPattern = " _ ", splitInd = 2, returnCharacter = TRUE)), 360)
    checkTrue(
        is.character(MetabolomicTools:::cutUniquePreMZ(sd02_deconvoluted[,4], 
            splitPattern = " _ ", splitInd = 2, returnCharacter = TRUE)))
    checkTrue(
        is.numeric(MetabolomicTools:::cutUniquePreMZ(sd02_deconvoluted[,4], 
            splitPattern = " _ ", splitInd = 2, returnCharacter = FALSE)))
}
## END unit test cutUniquePreMZ

## START unit test convert2MSP
testMSP <- MetabolomicTools::convert2MSP(sd02_deconvoluted, 
    splitPattern = " _ ", splitInd = 2)
test_convert2MSP <- function() {
    checkTrue(class(testMSP) == "MSP")
    checkEquals(length(testMSP), 360)
    checkTrue(is.data.frame(testMSP@msp))
    checkEquals(dim(testMSP@msp), c(7263, 2))
    checkEquals(as.numeric(
        table(testMSP@msp[,1])["ADDUCTIONNAME: "]), 360)
    checkEquals(as.numeric(
        table(testMSP@msp[,1])["METABOLITENAME: "]), 360)
    checkEquals(as.numeric(
        table(testMSP@msp[,1])["NAME: "]), 360)
    checkEquals(as.numeric(
        table(testMSP@msp[,1])["Num Peaks: "]), 360)
    checkEquals(as.numeric(
        table(testMSP@msp[,1])["PRECURSORMZ: "]), 360)
    checkEquals(as.numeric(
        table(testMSP@msp[,1])["RETENTIONTIME: "]), 360)
}
## END unit test convert2MSP

## START unit test msp2FunctionalLossesMSP
testMSPNL <- MetabolomicTools::msp2FunctionalLossesMSP(testMSP)
test_msp2FunctionalLossesMSP <- function() {
    checkTrue(class(testMSPNL) == "MSP")
    checkEquals(length(testMSPNL), 360)
    checkTrue(is.data.frame(testMSPNL@msp))
    checkEquals(dim(testMSPNL@msp), c(7263, 2))
    checkTrue(is.data.frame(testMSPNL@msp))
    checkEquals(as.numeric(
        table(testMSPNL@msp[,1])["ADDUCTIONNAME: "]), 360)
    checkEquals(as.numeric(
        table(testMSPNL@msp[,1])["METABOLITENAME: "]), 360)
    checkEquals(as.numeric(
        table(testMSPNL@msp[,1])["NAME: "]), 360)
    checkEquals(as.numeric(
        table(testMSPNL@msp[,1])["Num Losses: "]), 360)
    checkEquals(as.numeric(
        table(testMSPNL@msp[,1])["PRECURSORMZ: "]), 360)
    checkEquals(as.numeric(
        table(testMSPNL@msp[,1])["RETENTIONTIME: "]), 360)
}
## END unit test msp2FunctionalLossesMSP