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
    checkTrue(is.data.frame(testMSP))
    checkEquals(dim(testMSP), c(7263, 2))
    checkTrue(is.data.frame(testMSP))
    checkEquals(as.numeric(table(testMSP[,1])["ADDUCTIONNAME: "]), 360)
    checkEquals(as.numeric(table(testMSP[,1])["METABOLITENAME: "]), 360)
    checkEquals(as.numeric(table(testMSP[,1])["NAME: "]), 360)
    checkEquals(as.numeric(table(testMSP[,1])["Num Peaks: "]), 360)
    checkEquals(as.numeric(table(testMSP[,1])["PRECURSORMZ: "]), 360)
    checkEquals(as.numeric(table(testMSP[,1])["RETENTIONTIME: "]), 360)
}
## END unit test convert2MSP

## START unit test msp2FunctionalLossesMSP
testMSPNL <- msp2FunctionalLossesMSP(testMSP)
test_msp2FunctionalLossesMSP <- function() {
    checkTrue(is.data.frame(testMSPNL))
    checkEquals(dim(testMSPNL), c(7263, 2))
    checkTrue(is.data.frame(testMSPNL))
    checkEquals(as.numeric(table(testMSPNL[,1])["ADDUCTIONNAME: "]), 360)
    checkEquals(as.numeric(table(testMSPNL[,1])["METABOLITENAME: "]), 360)
    checkEquals(as.numeric(table(testMSPNL[,1])["NAME: "]), 360)
    checkEquals(as.numeric(table(testMSPNL[,1])["Num Losses: "]), 360)
    checkEquals(as.numeric(table(testMSPNL[,1])["PRECURSORMZ: "]), 360)
    checkEquals(as.numeric(table(testMSPNL[,1])["RETENTIONTIME: "]), 360)
}
## END unit test msp2FunctionalLossesMSP