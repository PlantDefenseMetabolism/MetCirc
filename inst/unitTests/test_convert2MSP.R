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

## END 

#' @name msp2FunctionalLossesMSP
#' @title Convert MSP to MSP with functional losses
#' @description msp2FunctionalLossesMSP converts a data.frame in msp format 
#' (with fragments) into a data.frame in msp format (with neutral losses)
#' @usage msp2FunctionalLossesMSP(msp)
#' @param msp data.frame, a data.frame in msp format (with fragments)
#' @details msp2FunctionalLosses can be used when you want to calculate 
#' the similarity based on neutral losses instead of fragments
#' @return msp2FunctionalLossesMSP returns a data.frame in msp format 
#' (with neutral losses).
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{msp2FunctionalLossesMSP(msp)}
#' load(system.file("data/sd02_deconvoluted.RData", 
#'      package = "MetabolomicTools")) 
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' finalMSPNL <- msp2FunctionalLossesMSP(msp = finalMSP)
#' @export
msp2FunctionalLossesMSP <- function(msp) {
    precmz <- getPrecursorMZ(msp)
    rt <- getRT(msp)
    indices <- getBegEndIndMSP(msp)
    indBegL <- indices[[1]]
    indEndL <- indices[[2]]
    ## create data frame for MSP file
    finalMSP <- matrix(data = NA, nrow = dim(msp)[1], ncol = 2) 
    finalMSP <- as.data.frame(finalMSP)
    
    ## create MSP from 
    for (i in 1:length(precmz)) {
        
        indBeg <- indBegL[i]
        indEnd <- indEndL[i]
        
        neutralL <- (as.numeric(precmz[i]) - as.numeric(msp[indBeg:indEnd,1]))
        neutralL <- -1 * neutralL
        
        entry <- rbind(
            c("NAME: ", "Unknown"),
            c("RETENTIONTIME: ", rt[i]),
            c("PRECURSORMZ: ", precmz[i]),
            c("METABOLITENAME: ", "Unknown"),
            c("ADDUCTIONNAME: ", "Unknown"),
            c("Num Losses: ", length(indBeg:indEnd)),
            matrix(c(neutralL, msp[indBeg:indEnd,2]), ncol = 2),
            c(" ", " ")
        )
        entry <- as.matrix(entry)
        ## determine first empty line
        newstart <- which(is.na(finalMSP[,1]))[1]
        ## determine last line to write to
        newend <- newstart + dim(entry)[1] - 1
        finalMSP[newstart:newend,] <- entry
    }
    
    return(finalMSP)
}
