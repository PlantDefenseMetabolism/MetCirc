data("sd02_deconvoluted", package = "MetCirc")
## START unit test cutUniquePreMZ
test_cutUniquePreMZ <- function() {
    checkTrue(is.vector(cutUniquePrecursor(sd02_deconvoluted[,4],
        splitPattern = " _ ", splitInd = 2, returnCharacter = TRUE)))
    checkEquals(
        length(cutUniquePrecursor(sd02_deconvoluted[,4], 
            splitPattern = " _ ", splitInd = 2, returnCharacter = TRUE)), 360)
    checkTrue(
        is.character(cutUniquePrecursor(sd02_deconvoluted[,4], 
            splitPattern = " _ ", splitInd = 2, returnCharacter = TRUE)))
    checkTrue(
        is.numeric(cutUniquePrecursor(sd02_deconvoluted[,4], 
            splitPattern = " _ ", splitInd = 2, returnCharacter = FALSE)))
}
## END unit test cutUniquePreMZ

## START unit test convert2MSP
testMSP <- convert2MSP(sd02_deconvoluted, splitPattern = " _ ", splitIndMZ = 2)

test_convert2MSP <- function() {
    checkTrue(class(testMSP) == "MSP")
    checkEquals(length(testMSP), 360)
    checkTrue(is.data.frame(peaks(testMSP)))
    checkEquals(dim(peaks(testMSP)), c(5463, 2))
    checkEquals(length(testMSP@adduct), 360)
    checkEquals(length(testMSP@information), 360)
    checkEquals(length(testMSP@names), 360)
    checkEquals(as.numeric(
        table(peaks(testMSP)[,1])["Num Peaks: "]), 360)
    checkEquals(length(testMSP@mz), 360)
    checkEquals(length(testMSP@rt), 360)
}
## END unit test convert2MSP

## START unit test msp2FunctionalLossesMSP
testMSPNL <- msp2FunctionalLossesMSP(testMSP)
test_msp2FunctionalLossesMSP <- function() {
    checkTrue(class(testMSPNL) == "MSP")
    checkEquals(length(testMSPNL), 360)
    checkTrue(is.data.frame(testMSPNL@msp))
    checkEquals(dim(testMSPNL@msp), c(5463, 2))
    checkTrue(is.data.frame(testMSPNL@msp))
    checkEquals(length(testMSPNL@adduct), 360)
    checkEquals(length(testMSPNL@information), 360)
    checkEquals(length(testMSPNL@names), 360)
    checkEquals(as.numeric(
        table(peaks(testMSPNL)[,1])["Num Losses: "]), 360)
    checkEquals(length(testMSPNL@mz), 360)
    checkEquals(length(testMSPNL@rt), 360)
}
## END unit test msp2FunctionalLossesMSP

## START unit test MSP-class
test_MSP <- function() {
    checkEquals(is(testMSP), "MSP")
}

## END unit test MSP-class

## START unit test length-method
test_length <- function() {
    checkEquals(length(testMSP), 360)
    checkEquals(length(testMSP), length(testMSPNL))
    checkTrue(is.numeric(length(testMSP)))
}
## END unit test length-method

## START unit test show-method
test_show <- function() {
    checkTrue(is.null(show(testMSPNL)))
}
## END unit test show-method

## START unit test getPrecursorMZ
test_getPrecursorMZ <- function() {
    checkEquals(length(getPrecursorMZ(testMSP)), 360)
    checkTrue(is.vector(getPrecursorMZ(testMSP)))
    checkTrue(is.numeric(getPrecursorMZ(testMSP)))
    checkException(getPrecursorMZ(sd02_deconvoluted))
}
## END unit test getPrecursorMZ

## START unit test getRT
test_getRT <- function() {
    checkEquals(length(getRT(testMSP)), 360)
    checkTrue(is.vector(getRT(testMSP)))
    checkTrue(is.numeric(getRT(testMSP)))
    checkException(getRT(sd02_deconvoluted))
}
## END unit test getRT

## START unit test peaks-method
test_peaks <- function() {
    checkTrue(is.data.frame(peaks(testMSP)))
    checkEquals(dim(peaks(testMSP)), c(5463, 2))
    checkEquals(dim(peaks(testMSP)), dim(peaks(testMSPNL)))
}
## END unit test peaks-method

## START unit test combine-method
test_combine <- function() {
    checkEquals(length(combine(testMSP, testMSP)), 720)
    checkEquals(dim(peaks(combine(testMSP, testMSP))), c(10926, 2))
}
## END unit test combine-method

## START unit test getName-method
test_getName <- function() {
    checkTrue(all(is.character(names(testMSP))))
    checkTrue(length(names(testMSP)) == length(testMSP))
}
## END unit test getName-method


## START unit test getInformation-method
test_getMetaboliteName <- function() {
    checkTrue(all(is.character(information(testMSP))))
    checkTrue(length(information(testMSP)) == length(testMSP))
    checkException(information("x"))
}
## END unit test getMetaboliteName-method

## START unit test getMetaboliteClass-method
test_getMetaboliteClass <- function() {
    checkTrue(all(is.character(classes(testMSP))))
    checkTrue(length(classes(testMSP)) == length(testMSP))
    checkException(classes("x"))
}
## END unit test getMetaboliteClass-method

## START unit test [
test_extract <- function() {
    checkEquals(length(testMSP[1]), 1)
    checkEquals(length(testMSP[1:10]), 10)
    checkEquals(testMSP@msp[1:13,], peaks(testMSP[1:2]))
    checkTrue(is(testMSP[1:10]) == "MSP")
    checkException(testMSP[400])
}
## END unit test ]

## START unit test getBegEndIndMSP
testMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", 
                       splitIndMZ = 2, splitIndRT = NULL)
testMSPmsp <- peaks(testMSP)
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