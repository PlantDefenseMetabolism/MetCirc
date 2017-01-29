data("convertMSP2MSP", package = "MetCirc")
msp <- convertMSP2MSP(msp2msp)
## START unit test convertMSP2MSP
test_convertMSP2MSP <- function() {
    checkEquals(dim(msp@msp), c(6822, 2))
    checkEquals(length(msp@mz), 150)
    checkTrue(is.numeric(msp@mz))
    checkTrue(is.numeric(msp@rt))
    checkTrue(is.character(msp@names))
    checkTrue(is.character(msp@classes))
    checkTrue(is.character(msp@information))
    checkTrue(is.character(msp@adduct))
    checkException(
        convertMSP2MSP(data.frame(x = c("a", "b", "c"), y = c(1, 2, 3))))
    checkException(
        convertMSP2MSP(data.frame(x = c("Name: ", "Rt: ", "c"), y = c(1, 2, 3))))
}
## END unit test convertMSP2MSP

