#' @name allocatePrecursor2mz
#' @title allocatePrecursor2mz: Join two data sources
#' @description Allocates precursor ions to candidate m / z values based on 
#' minimal distance of m / z and deviance of rt based on an objective function
#' @usage allocatePrecursor2mz(sd01, sd02, kNN = 10, mzCheck = 1, rtCheck = 30, 
#'      mzVsRTbalance = 10000, splitPattern = "_", splitInd = 2)
#' @param sd01 is the output of the \code{XCMS} and \code{CAMERA} 
#' processing and statistical analysis and \code{XCMS} and \code{CAMERA} 
#' scripts (see Li et al. 2015 and vignette for further information)
#' @param sd02 is a data.frame with idMS/MS deconvoluted spectra with fragment 
#' ions (m/z, retention time, relative intensity in \%) and the corresponding 
#' principal component group with the precursor ion. sd02 
#' has four columns, the first column contains the m/z 
#' value, the second column the rt, the third column the intensity, the fourth
#' column the pcgroup_precursorMZ
#' @param kNN numerical, number of k-nearest neighbours based on deviation
#' from m/z (i.e. the k entries with the smallest deviation)
#' @param mzCheck numerical, maximum tolerated distance for m/z (strong 
#'      criterion here)
#' @param rtCheck numerical, maximum tolerated distance for retention time
#' @param mzVsRTbalance numerical, multiplicator for mz value before calculating 
#' the (euclidean) distance between two peaks, high value means that there is 
#' a strong weight on the deviation m/z value 
#' @param splitPattern character, character vector to use for splitting, see
#'  ?strsplit for further information
#' @param splitInd numeric, extract precursor mz at position splitInd
#' @details This function combines different data sources. 
#' \code{convertExampleDF} is a \code{data.frame} which comprises information 
#' on a specific metabolite per 
#' row stating the average retention time, average m/z, the name of the 
#' metabolite, the adduct ion name, the spectrum 
#' reference file name and additional information (TRIO/LVS). 
#' \code{allocatePrecursor2mz} uses \code{data.frame}s of the kind of 
#' \code{sd01\_outputXCMS} and \code{sd02\_deconvoluted} to create a 
#' \code{data.frame} of the kind of \code{convertExampleDF}. Allocation of 
#' precursor ions to candidate m/z values is based on minimal distance of m/z 
#' and deviance of retention time based on an objective function. We can specify 
#' threshold values for m/z and retention time to be used in 
#' \code{allocatePrecursor2mz}, as well as the number of neighbours based on 
#' deviation from m/z values. Also, we can specify the weight to base the 
#' selection on the m/z compared to the retention time (\code{mzVsRTbalance}). 
#' This might be useful because m/z values might differ less than the retention 
#' time in \code{sd01\_outputXCMS} and \code{sd02\_deconvoluted}. Please note, 
#' that it might be problematic to compare \code{sd01\_outputXCMS} and 
#' \code{sd02\_deconvoluted} and allocate precursor ions therefrom, 
#' especially when data were acquired under different conditions.
#' @return allocatePrecursor2mz returns a data.frame containing average 
#'      retention time, average mz, metabolite name, adduct ion name, 
#'      spectrum reference
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @references Li et al. (2015): Navigating natural variation in 
#' herbivory-induced secondary metabolism in coyote tobacco populations using 
#' MS/MS structural analysis. PNAS, 112, E4147--E4155, 10.1073/pnas.1503106112.
#' @examples
#' data("sd01_outputXCMS", package = "MetCirc")
#' data("sd02_deconvoluted", package = "MetCirc") 
#' data("convertExampleDF", package = "MetCirc")
#' allocatePrecursor2mz(sd01 = sd01_outputXCMS, sd02 = sd02_deconvoluted, 
#'      kNN = 10, mzCheck = 1, rtCheck = 30, mzVsRTbalance = 10000, splitPattern = "_", splitInd = 2)
#' @export
allocatePrecursor2mz <- function(sd01, sd02, kNN = 10, mzCheck = 1, 
    rtCheck = 30, mzVsRTbalance = 10000, splitPattern = "_", splitInd = 2) {
    ## Multiplicator for mz value before calculating the (euclidean) distance 
    ## between two peaks high value means that there is a strong weight on the 
    ## dev m/z value mzVsRTbalance
    if (kNN <= 0) stop("kNN has to be > 0")
    if (mzCheck <= 0) stop("mzCheck has to be > 0")
    if (rtCheck <= 0) stop("rtCheck has to be > 0")
    if (mzVsRTbalance <= 0) stop("mzVsRTbalance has to be > 0")
    if (!all(sort(colnames(sd02)) == c("intensity", "mz", "pcgroup_precursorMZ", "rt"))) 
        stop("colnames(sd02) have to be 'mz', 'rt', 'intensity', 'pcgroup_precursorMZ' ")
    colNamesSd01 <- colnames(sd01)
    if (!all(c("mz", "rt", "npeaks", "isotopes", "adduct", "pcgroup") %in% colNamesSd01)) 
        stop("colnames(sd01) have to have 'mz', 'rt', 'npeaks', 'isotopes', 'adduct', 'pcgroup'")
    precursor <- sd02[, "pcgroup_precursorMZ"]
    ## isolated mz values from e.g. pcgroup_precursorMZ column in 
    ## sd02_deconvoluted
    uniquePrecursor <- cutUniquePreMZ(precursor, 
                                      splitPattern = splitPattern, splitInd= splitInd)
    
    ## create finalCluster, which is the data.frame to store data
    finalCluster <- matrix(nrow = length(uniquePrecursor), ncol = dim(sd01)[2] + 7)
    colnames(finalCluster) <- c(colnames(sd01), "Metabolite Name", 
                "Adduct ion name", "Spectrum reference file name", 
                "check RT", "dev RT", "check mz", "deviation m/z")
    colnames(finalCluster)[which(colnames(finalCluster)=="mz")] <- "average mz"
    colnames(finalCluster)[which(colnames(finalCluster)=="rt")] <- "average rt"
    
    finalCluster <- as.data.frame(finalCluster)
    
    uniquePreMZ <- unique(precursor)
   
    
    ## LOOP WHICH WRITES TO finalCluster
    for (i in 1:length(uniquePrecursor)) {
        
        sd02mz <- uniquePrecursor[i]
        sd02mz <- as.numeric(sd02mz)
        sd01mz <- sd01[, "mz"]
        sd01mz <- as.character(sd01mz)
        sd01mz <- as.numeric(sd01mz)
        devmzOld <- devmz <- abs(sd02mz - sd01mz)
        
        ## use only kNN m/z dev
        sortDevMz <- sort(devmz)[1:kNN]
        
        ## get indices in sd01 of the smallest deviances
        indSortDevMZOld <- indSortDevMZ <- match(sortDevMz, devmz) 
        ## truncate devmz such that it only includes kNN m/z deviances
        devmz <- devmz[indSortDevMZ]
        
        ## check if devmz is in tolerated distance
        if (any(devmz <= mzCheck)) {
            ## truncate such that indSortDevMZ includes only indices and 
            ## devmz only m/z within the tolerance value
            indSortDevMZ <- indSortDevMZ[devmz <= mzCheck]
            devmz <- devmz[devmz <= mzCheck]
            ToleranceCheckMZ <- TRUE
        } else {
            print (c(i,"Deviation of m/z is greater than tolerance value. 
                     I won't truncate the kNN."))
            devmz <- devmz
            ToleranceCheckMZ <- FALSE
        }
        
        ## calculate fake rt from sd02 (from fragment rt values)
        ind <- which(uniquePreMZ[i] == precursor)    
        sd02rt <- sd02[ind, "rt"]
        sd02rt <- mean(sd02rt) 
        
        ## determine rt values from sd01
        sd01rt <- sd01[indSortDevMZ, "rt"]
        sd01rt <- as.character(sd01rt)
        sd01rt <- as.numeric(sd01rt)
        
        devrt <- abs(sd02rt - sd01rt)
        
        if (any(devrt <= rtCheck)) {
            ## truncate devmz and devrt that it is included in the tolerance 
            ## value
            devmz <- devmz[devrt <= rtCheck] 
            devrt <- devrt[devrt <= rtCheck]
            ToleranceCheckRT <- TRUE
            objective <- mzVsRTbalance * devmz + devrt
        } else {
            print(c(i, "Deviation of rt is greater than tolerance value. 
                    I won't use rt as a criterion."))
            ToleranceCheckRT <- FALSE
            objective <- devmz ## use only devmz
        }
        
        ## find smallest value for objective function 
        minInd <- which.min(objective) 
        ## get index in sd01 
        minInd <- which(devmz[minInd] == devmzOld) 
        
        
        ## write to finalCluster:
        ## get the entry of sd01 with the smallest value
        XCMS <- sd01[minInd,]
        replaceInd <- !is.na(match(colnames(finalCluster), colNamesSd01))
        finalCluster[i, replaceInd] <- as.vector(as.matrix(XCMS[replaceInd]))

        finalCluster[i, "average rt"] <- sd02rt
        finalCluster[i, "average mz"] <- uniquePrecursor[i]
        finalCluster[i, "Metabolite Name"] <- "Unknown"
        finalCluster[i, "Adduct ion name"] <- 
            if (nchar(as.character(XCMS[, "adduct"]) == 0)) "Unknown" else XCMS[,"adduct"]
        
        finalCluster[i, "Spectrum reference file name"] <- "Unknown"
        
        ##x <- XCMS[,which(colnames(XCMS) == "1"):which(colnames(XCMS) == "183")]
        ##x <- as.matrix(x)
        ##x <- as.vector(x)
        ##entry[, which(colnames(entry)=="1"):which(colnames(entry)=="183")] <- x
        
        finalCluster[i, "check RT"] <- ToleranceCheckRT
        finalCluster[i, "dev RT"] <- sd02rt - as.numeric(as.character(XCMS[,"rt"]))
        finalCluster[i, "check mz"] <- ToleranceCheckMZ
        finalCluster[i, "deviation m/z"] <- sd02mz - as.numeric(as.character(XCMS[,"mz"]))

    }
    
    return(finalCluster)
}

