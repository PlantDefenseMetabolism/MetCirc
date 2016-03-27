## START unit test allocatePrecursor2mz
testAllocatePrecursor2MZ <- allocatePrecursor2mz(sd01_outputXCMS, sd02_deconvoluted)

test_allocatePrecursor2mz <- function() {
    checkTrue(is.data.frame(testAllocatePrecursor2MZ))
    checkEquals(dim(testAllocatePrecursor2MZ), c(360, 192))
}
## END unit test allocatePrecursor2mz
