#' @name convertMSP2MSP
#' @title Convert msp data frame into MSP format
#' @description Convert msp data frame into MSP format
#' @usage convertMSP2MSP(msp)
#' @param msp data frame, in msp format, has the row entries "Name:", "Rt=",
#'  "Num Peaks:" and information on fragments and peak areas 
#' @details msp data frame from other sources. 
#' @return convertMSP2MSP returns an object of class MSP
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' \dontrun{convertMSP2MSP(msp = msp)}
#' @export
convertMSP2MSP <- function(msp) {
    MSP <- as.matrix(as.character(msp[,1]))
    nameInd <- grep("Name: ", MSP[,1]) ## indices of Name
    names <- unlist(lapply(strsplit(as.character(MSP[nameInd, 1]), " "), "[", 2))
    rettimeInd <- grep("Rt=", MSP[,1]) ## indices of Rt=
    ## truncate the row with Rt= such, that retention time is retrieve
    rettime <- unlist(lapply(strsplit(
        unlist(lapply(strsplit(as.character(MSP[rettimeInd,]), " "), "[", 2)), split = "="), "[", 2))

    numpeaks <- grep("Num Peaks:", MSP[,1]) ## indices of Num of peaks
    peaks <- numpeaks + 1 ## indices of fragments and intensities
    peakentry <- as.character(MSP[peaks,])
    peakentry <- strsplit(peakentry, " ")
    
    ## retrieve fragment and intensity
    peakentry_fragment_l <- lapply(peakentry, function(x) as.numeric(x[seq(1, length(x) - 1, 2)]))
    peakentry_intensity_l <- lapply(peakentry, function(x) as.numeric(x[seq(2, length(x), 2)]))
    
    ## sort peakentries according to increasing fragment values
    peakentry_fragment_l_s <- lapply(peakentry_fragment_l, sort)
    peakentry_intensity_l_s <- lapply(1:length(peakentry_intensity_l), function(x) peakentry_intensity_l[[x]][order(peakentry_fragment_l[[x]])])
    peakentry_intensity_l_s <- lapply(peakentry_intensity_l_s, function(x) x / max(x) * 100)

    ## create MSP entry
    msp_l <- lapply(1:length(names), function(x) {
        rbind(
        c("NAME: ", names[x]),
        c("RETENTIONTIME: ", rettime[x]),
        c("PRECURSORMZ: ", max(peakentry_fragment_l_s[[x]])),
        c("METABOLITENAME: ", names[x]),
        c("METABOLITECLASS: ", "Unknown"),
        c("ADDUCTIONNAME: ", "Unknown"),
        c("Num Peaks: ", strsplit(as.character(MSP[numpeaks[x],]), " ")[[1]][3]),
        cbind(peakentry_fragment_l_s[[x]], peakentry_intensity_l_s[[x]]),
        c("", ""))})

    finalMSP <- do.call(rbind.data.frame, msp_l)

    return(new("MSP", msp = finalMSP))
}
