#' @name convertMSP2MSP
#' @title Convert msp data frame into MSP format
#' @description Convert msp data frame into MSP format
#' @usage convertMSP2MSP(msp)
#' @param msp data.frame, see \code{Details} for further information. 
#' @details msp is a data frame of a .MSP file, a typical data file for 
#' MS/MS libraries. The data frame has two columns and contains in the first
#' column the entries "NAME:", 
#' "PRECURSORMZ:" (or "EXACTMASS:"), "Num Peaks:"  and information on fragments and 
#' peak areas/intensities. It may additionally contain row entries:
#' \code{convertMSP2MSP} will try to find the row entries "RETENTIONTIME:",
#' "ADDUCTIONNAME:" (or "PRECURSORTYPE:"), "CLASS:" and "INFORMATION:" and 
#' extract the respective information in the second column.
#' @return convertMSP2MSP returns an object of class MSP. 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("convertMSP2MSP", package = "MetCirc")
#' convertMSP2MSP(msp = msp2msp)
#' @export
convertMSP2MSP <- function(msp) {
    
    nameInd <- grep("NAME:", msp[,1]) ## indices of Name
    
    ## get indices of precursor mz or exactmass
    mzInd_1 <- grep("PRECURSORMZ:", msp[,1])
    mzInd_2 <- grep("EXACTMASS:", msp[,1])
    
    mzInd <- if (length(mzInd_1) >= length(mzInd_2)) {
        mzInd_1 
    } else {
        mzInd_2
    }
    
    ## get indices of "Num Peaks:"  
    numpeaksInd <- grep("Num Peaks:", msp[,1])
    
    ## get indices of "RETENTIONTIME:"
    rtInd <- grep("RETENTIONTIME:", msp[,1])
    
    ## get indices of "ADDUCTIONNAME:" or "PRECURSORTYPE:"
    adductInd_1 <- grep("ADDUCTIONNAME:", msp[,1])
    adductInd_2 <- grep("PRECURSORTYPE:", msp[,1])
    
    adductInd <- if (length(adductInd_1) >= length(adductInd_2)) {
        adductInd_1
    } else {
        adductInd_2
    }
    
    ## get indices of "CLASS:"
    classInd <- grep("CLASS:", msp[,1])
    
    ## get indices of "INFORMATION:"
    informationInd <- grep("INFORMATION:", msp[,1])
    
    
    NAMES <- as.character(msp[nameInd, 2])
    MZ <- as.numeric(as.character(msp[mzInd, 2]))
    NUMPEAKS <- as.numeric(as.character(msp[numpeaksInd, 2]))
    
    ## how many entries are in msp
    numEntries <- length(NAMES)
    
    if (numEntries != length(MZ)) 
        stop("length of precursor mz != length of names")
    
    if (numEntries != length(NUMPEAKS)) 
        stop("length of precursor mz != length of MS/MS fragments entries")
    
    
    ## create annotation vectors
    RT <- if(length(rtInd) == numEntries) {
        as.numeric(as.character(msp[rtInd, 2]))
    } else {
        rep(NaN, numEntries)   
    }
       
    ADDUCT <- if(length(adductInd) == numEntries) {
        as.character(msp[adductInd, 2])
    } else {
        rep("Unknown", numEntries)   
    }
       
    CLASS <- if(length(classInd) == numEntries) {
        as.character(msp[classInd, 2])
    } else {
        rep("Unknown", numEntries)
    }
    
    INFORMATION <- if(length(informationInd) == numEntries) {
        as.character(msp[informationInd, 2])
    } else {
        rep("Unknown", numEntries)
    }
   
    
    
    MSP <- NULL
    for (i in 1:numEntries) {
        beg <- numpeaksInd[i] + 1
        end <- numpeaksInd[i] + NUMPEAKS[i]
        fragment <- as.numeric(as.character(msp[beg:end, 1]))
        intensity <- as.numeric(as.character(msp[beg:end, 2]))
        
        ## sort fragments according to increasing fragment values
        fragment <- sort(fragment)
        intensity <- intensity[order(fragment)]
        
        ## delete double entries
        intensity <- intensity[!duplicated(fragment)]
        fragment <- fragment[!duplicated(fragment)]
        
        ## calculate percentages
        intensity <- intensity / max(intensity) * 100
        
        mspI <- rbind(
                    c("Num Peaks: ", length(intensity)),
                    cbind(fragment, intensity), 
                    c("", "")
        )
        
        MSP <- rbind(MSP, mspI)
        
    }
    
    MSP <- as.data.frame(MSP)

    msp <- new("MSP", msp = MSP, mz = MZ, rt = RT, names = NAMES, 
               classes = CLASS, information = INFORMATION, adduct = ADDUCT)
    return(msp)
}
