#' @name getBegEndIndMSP
#' @title Get beginning and end indices of each entry in a data.frame in 
#' msp format
#' @description Get beginning and end indices of each entry in a data.frame in 
#' msp format
#' @usage getBegEndIndMSP(msp)
#' @param msp data.frame in msp format, see ?convert2MSP for further information
#' @details Internal use to retrieve indices when fragments start and end. 
#' @return getBegEndIndMSP returns a list of length 2 where the first entry
#' contains the start indices and the second the end indices
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", 
#'                          splitIndMZ = 2, splitIndRT = 3)
#' finalMSPdf <- getMSP(finalMSP)
#' getBegEndIndMSP(finalMSPdf)
#' @export
getBegEndIndMSP <- function(msp) {
    
    ## beginning 
    indPeaks <- which(msp[,1] == "Num Peaks: ")
    indLosses <- which(msp[,1] == "Num Losses: ")
    
    if (length(indPeaks) > length(indLosses)) indNumPeaks <- indPeaks
    ## <=: to get all possibilities
    if (length(indPeaks) <= length(indLosses)) indNumPeaks <- indLosses 
    
    indEnd <- as.numeric(msp[indNumPeaks, 2])
    indBeg <- indNumPeaks + 1 
    indEnd <- indEnd + indBeg - 1 
    
    return(list(indBeg, indEnd))
}

#' @name binning
#' @title Bin m/z values
#' @description Bin m/z values
#' @usage binning(msp, tol = 0.01, group = NULL, method = c("median", "mean"))
#' @param msp data.frame in msp format, see ?convert2MSP for further information
#' @param tol numerical, boundary value until which neighboured peaks will be 
#'      joined together
#' @param group character vector, to which group does the entry belong to
#' @param method character vector, method has to be "median" or "mean"
#' @details The functions bins fragments together by obtaining bins via 
#' calculating either mean or medians of fragments which were put in intervals 
#' according to the \code{tol} parameter. 
#' @return binning returns a matrix where rownames are precursor ions
#' (m/z / retention time) and colnames are newly calculated m/z values which 
#' were binned. 
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples data("idMSMStoMSP", package = "MetCirc")
#' ##group <- sample(c("yl", "ol", "s","r"), size = length(finalMSP), replace=TRUE) 
#' binning(msp = finalMSP, tol = 0.01, group = NULL, method = "median")
#' @export
#' @importFrom stats median
binning <- function(msp, tol = 0.01, group = NULL, method = c("median", "mean")) { 
    
    method <- match.arg(method)
    ## msp is .msp file
    ## tol is tolerance value for binning
    if (!is(msp) == "MSP") stop("msp is not of class MSP.")
    
    precmz <- getPrecursorMZ(msp)
    rt <- getRT(msp)
    
    msp <- getMSP(msp)
    
    if (is.null(group)) {
        print("argument group is not specified, will create dummy group")
        group <- rep("a", length(precmz))
    }
    ## if group is not NULL then: 
    if (length(precmz) != length(group)) 
        stop("length of precursor ions != length(group)")
    
    indices <- getBegEndIndMSP(msp)
    indBeg <- indices[[1]]
    indEnd <- indices[[2]]
    
    ## create list with indices in msp for each precursor
    l <- lapply(1:length(indBeg), function(x) c(indBeg[x]:indEnd[x]))
    
    IndFrag <- unlist(l)
    ## fragments 
    frag <- msp[IndFrag, 1] 
    frag <- as.numeric(frag)
    frag <- unique(frag)
    ## sort frag and get order
    frag_s <- sort(frag)
    ##frag_order <- order(frag)
    
    steps <- (max(frag_s) - min(frag_s)) / tol

    if (method == "median") {
        ## calculate median of values in bins
        bins <- tapply(frag_s, cut(frag_s, steps), median)
    }
    
    if (method == "mean") {
        ## calculate mean of values in bins
        bins <- tapply(frag_s, cut(frag_s, steps), mean)
    }
    ## remove bins which no not show up
    bins <- bins[!is.na(bins)]
    ## vectorise bins (do not use named vector)
    bins <- as.vector(bins)
    
    mm <- matrix(data = 0, nrow = length(precmz), ncol = length(bins))
    ## convoluted MZ is column names
    colnames(mm) <- bins
    rownamesMM <- paste(precmz, rt, sep="/")
    rownames(mm) <- paste(group, rownamesMM, sep="_")
    
    ## write each entry of msp to mm
    for (i in 1:length(l)) {
        ## get ith entry of msp
        mspI <- msp[l[[i]],]
        mspI <- data.matrix(mspI)
        ## inds hosts the column indices for the respective fragment
        inds <- lapply(mspI[,1], FUN = function(x) which.min(abs(x - bins)))
        inds <- as.numeric(unlist(inds))
        ## assign to entry which is closest to the bin value
        ## check if there are duplicated entries
        if (any(duplicated(inds))) {
            mm[i, inds[!duplicated(inds)]] <- mspI[!duplicated(inds),2]
            mm[i, inds[which(duplicated(inds))]] <- mm[i, inds[which(duplicated(inds))]] + mspI[duplicated(inds),2]
        } else {
            mm[i, inds] <- mspI[,2]
        }
    }
    
    class(mm) <- "numeric"
    return(mm)
}