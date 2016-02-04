#' @name cutUniquePreMZ
#' @title Get unique precursor ions
#' @description Get unique precursor ions
#' @usage cutUniquePreMZ(precursor, splitPattern = splitPattern, 
#'      splitInd = splitInd, returnCharacter = TRUE)
#' @param precursor, character with splitPattern
#' @param splitPattern character, character vector to use for splitting, 
#'      see ?strsplit for further information
#' @param splitInd numeric, extract precursor mz at position splitInd
#' @param returnCharacter logical, if TRUE return character, if FALSE 
#'      return numeric
#' @details Internal function.
#' @return cutUniquePreMZ returns character as specified by parameters
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{cutUniquePreMZ(precursor, splitPattern = splitPattern, 
#'      splitInd = splitInd, returnCharacter = TRUE)}
cutUniquePreMZ <- function(precursor, splitPattern = splitPattern, 
                            splitInd = splitInd, returnCharacter = TRUE) {
    ## split precursors according to split pattern
    precursor <- as.character(precursor)
    splitPrecursor <- strsplit(precursor, split = splitPattern)
    ## extract precursor mz at position splitInd
    splitPrecursor <- lapply(splitPrecursor,"[", splitInd)
    PrecursorMZ <- unlist(splitPrecursor)
    lenPreMZ <- length(PrecursorMZ)
    
    ## change character to numeric
    if (!returnCharacter)
        PrecursorMZ <- as.numeric(PrecursorMZ)
    #uniquePreMZ <- unique(precursor)
    #lenUniquePreMZ <- length(uniquePreMZ)
    uniquePreMZ_cut <- unique(PrecursorMZ)
    
    return(uniquePreMZ_cut)
}

#' @name convert2MSP
#' @title Convert deconvoluted matrix into MSP format
#' @description Convert deconvoluted matrix into MSP format
#' @usage convert2MSP(mm, splitPattern = "_", splitInd = 1)
#' @param mm matrix, mm has four columns, the first column contains the m/z 
#'  value, the second column the rt, the third column the intensity, the fourth
#'  column the pcgroup_precursorMZ
#' @param splitPattern character, splitPattern is the pattern which separates 
#'      elements and precursor m/z
#' @param splitInd numeric, the position of the precursor m/z concerning 
#'      separation by splitPattern
#' @details Creates a data entry for each precursor ion. Each entry in the 
#' return object has the following information: NAME, RETENTIONTIME, 
#'      PRECURSORMZ, METABOLITENAME, ADDUCTIONNAME, Num Peaks and a list of 
#'      fragments together with their intensities.
#' @return convert2MSP returns a data.frame in .msp file format
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' load(system.file("data/sd02_deconvoluted.RData", 
#'      package = "MetabolomicTools")) 
#' convert2MSP(mm = sd01_outputXCMS, splitPattern = "_", splitInd = 1)
#' @export
convert2MSP <- function (mm, splitPattern = "_", splitInd = 1) {
    
    colNames <- colnames(mm)
    if (colNames[1] != "mz") stop("name of first colomn is not mz")
    if (colNames[2] != "rt") stop("name of second column is not rt")
    if (colNames[3] != "intensity") 
                            stop("name of third column is not intensity")
    
    
    ## if (colNames[4] != "pcgroup_precursorMZ") break
     
    precursor <- mm[,4]
    precursor <- as.character(precursor)
    
    uniquePreMZ <- unique(precursor)
    uniquePreMZ_cut <- cutUniquePreMZ(precursor = precursor, 
            splitPattern = splitPattern, splitInd = splitInd)
    lenUniquePreMZ <- length(uniquePreMZ_cut)
    
    ## add PrecursorMZ to deconvoluted idMSMS
    ## mm <- cbind(mm, PrecursorMZ)
    
    ## create data frame for MSP file
    finalMSP <- matrix(data = NA, nrow = 7 * lenUniquePreMZ + dim(mm)[1], 
            ncol = 2) ## 7 new entries + all fragment ion entries
    finalMSP <- as.data.frame(finalMSP)
    
    ## write to data frame
    for (i in 1:lenUniquePreMZ) {
        ind <- which(uniquePreMZ[i] == precursor)    
        entry <- rbind(
            c("NAME: ", "Unknown"),
            c("RETENTIONTIME: ", mean(mm[ind,"rt"])),
            c("PRECURSORMZ: ", uniquePreMZ_cut[i]),
            c("METABOLITENAME: ", "Unknown"),
            c("ADDUCTIONNAME: ", "Unknown"),
            c("Num Peaks: ", length(ind)),
            mm[ind,c(1,3)],
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
