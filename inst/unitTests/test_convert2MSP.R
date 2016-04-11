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
    checkEquals(dim(getMSP(testMSP)), c(7623, 2))
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
    checkEquals(dim(testMSPNL@msp), c(7623, 2))
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

## START unit test getMSP-method
test_getMSP <- function() {
    checkTrue(is.data.frame(getMSP(testMSP)))
    checkEquals(dim(getMSP(testMSP)), c(7623, 2))
    checkEquals(dim(getMSP(testMSP)), dim(getMSP(testMSPNL)))
}
## END unit test getMSP-method

## START unit test combine-method
test_combine <- function() {
    checkEquals(length(combine(testMSP, testMSP)), 720)
    checkEquals(dim(getMSP(combine(testMSP, testMSP))), c(15246, 2))
}
## END unit test combine-method

## START unit test getName-method
test_getName <- function() {
    checkTrue(all(is.character(getName(testMSP))))
    checkTrue(length(getName(testMSP)) == length(testMSP))
    checkException(getName("x"))
}
## END unit test getName-method

## START unit test setName-method
test_setName <- function() {
    checkException(setName(testMSP, "a"))
    checkTrue(all(
        getName(setName(testMSP, rep("a", length(testMSP)))) == rep("a", 360)))
}
## END unit test setName-method

## START unit test getMetaboliteName-method
test_getMetaboliteName <- function() {
    checkTrue(all(is.character(getMetaboliteName(testMSP))))
    checkTrue(length(getMetaboliteName(testMSP)) == length(testMSP))
    checkException(getMetaboliteName("x"))
}
## END unit test getMetaboliteName-method

## START unit test setMetaboliteName-method
test_setMetaboliteName <- function() {
    checkException(setMetaboliteName(testMSP, "a"))
    checkTrue(
        all(getMetaboliteName(
            setMetaboliteName(testMSP, rep("a", length(testMSP)))
            ) == rep("a", 360)))
}
## END unit test setMetaboliteName-method

## START unit test getMetaboliteClass-method
test_getMetaboliteClass <- function() {
    checkTrue(all(is.character(getMetaboliteClass(testMSP))))
    checkTrue(length(getMetaboliteClass(testMSP)) == length(testMSP))
    checkException(getMetaboliteClass("x"))
}
## END unit test getMetaboliteClass-method

## START unit test setMetaboliteClass-method
test_setMetaboliteClass <- function() {
    checkException(setMetaboliteClass(testMSP, "a"))
    checkTrue(
        all(getMetaboliteClass(
            setMetaboliteClass(testMSP, rep("a", length(testMSP)))
            ) == rep("a", 360)))
}
## END unit test setMetaboliteClass-method

## START unit test [
test_extract <- function() {
    checkEquals(length(testMSP[1]), 1)
    checkEquals(length(testMSP[1:10]), 10)
    checkEquals(testMSP@msp[1:25,], getMSP(testMSP[1:2]))
    checkTrue(is(testMSP[1:10]) == "MSP")
    checkException(testMSP[400])
}
## END unit test ]