## load data sets
data("sd02_deconvoluted", package = "MetCirc")
testMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", 
                       splitIndMZ = 2, splitIndRT = NULL)
testMSPmsp <- peaks(testMSP)

## START unit test binning
compartment <- c(rep("a", 90), rep("b", 90), rep("c", 90), rep("d", 90))
## create binnedMSPs
binnedMSP001 <- binning(testMSP, 0.01, group = compartment, method = "mean")
binnedMSP0 <- binning(testMSP, 0, group = compartment, method = "mean")
test_binning <- function() {
    checkTrue(is.matrix(binnedMSP001))
    checkTrue(is.numeric(binnedMSP001))
    checkEquals(dim(binnedMSP001), c(360, 764)) 
    checkException(binning(finalMSP, 1, compartment[1:7]))
    checkEquals(dim(binnedMSP0)[2], length(unique(testMSPmsp[,1])) - 2 ) 
    ## check if all unique fragments are used 
    ## (- 2 because of "Num Peaks: " and " ")
    checkException(binning(finalMSP, -0.1))
}
## END unit test binning
