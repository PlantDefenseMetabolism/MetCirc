#' @name binning
#' @title Bin m/z values
#' @description Bin m/z values
#' @usage binning(msp, tol = 0.01, group = NULL, method = c("median", "mean"), verbose = FALSE)
#' @param msp \code{MSP}-object, see ?convert2MSP for further information
#' @param tol \code{numerical}, boundary value until which neighboured peaks will be 
#'      joined together
#' @param group \code{character} vector, to which group does the entry belong to
#' @param method \code{character} vector, method has to be \code{median} or 
#'  \code{mean}
#' @param verbose \code{logical} vector, if set to TRUE information will be printed
#' if groups were not detected
#' @details The functions \code{binning} bins fragments together by obtaining bins via 
#' calculating either mean or medians of fragments which were put in intervals 
#' according to the \code{tol} parameter. 
#' @return \code{binning} returns a \code{matrix} where rownames are precursor ions
#' (m/z / retention time) and colnames are newly calculated m/z values which 
#' were binned. Entires are intensity values in %. 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples data("idMSMStoMSP", package = "MetCirc")
#' binning(msp = finalMSP, tol = 0.01, group = NULL, method = "median", verbose = FALSE)
#' @export
#' @importFrom stats median
binning <- function(msp, tol = 0.01, group = NULL, method = c("median", "mean"), verbose = FALSE) { 
    
    method <- match.arg(method)
    ## msp is .msp file
    ## tol is tolerance value for binning
    if (!is(msp) == "MSP") stop("msp is not of class MSP.")
    
    precmz <- getPrecursorMZ(msp)
    rt <- getRT(msp)
    
    msp <- peaks(msp)
    
    if (is.null(group)) {
        if (verbose) print("argument group is not specified, will create dummy group")
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
    frag <- as.character(frag)
    frag <- as.numeric(frag)
    frag <- unique(frag)
    ## sort frag and get order
    frag_s <- sort(frag)
    ##frag_order <- order(frag)
    
    
    ## three cases for tol: smaller 0, equal to 0, greater to 0
    if (tol < 0) stop("tol has to be positive ")
    
    if (tol == 0) {
        bins <- frag_s
    }
    
    if (tol > 0) {
        steps <- (max(frag_s) - min(frag_s)) / tol

        if (method == "median") {
            ## calculate median of values in bins
            if (length(frag_s) > 1) {
                bins <- tapply(frag_s, cut(frag_s, steps), median)
            } else {
                bins <- frag_s
            }
        }
    
        if (method == "mean") {
            ## calculate mean of values in bins
            if (length(frag_s) > 1) {
                bins <- tapply(frag_s, cut(frag_s, steps), mean)    
            } else {
                bins <- frag_s
            }
             
        }
        ## remove bins which no not show up
        bins <- bins[!is.na(bins)]
        ## vectorise bins (do not use named vector)
        bins <- as.vector(bins)
    }
    
    mm <- matrix(data = 0, nrow = length(precmz), ncol = length(bins))
    ## convoluted MZ is column names
    colnames(mm) <- bins
    rownamesMM <- paste(precmz, rt, sep="/")
    rownames(mm) <- paste(group, rownamesMM, sep="_")
    
    ## write each entry of msp to mm
    for (i in 1:length(l)) {
        ## get ith entry of msp
        mspI <- msp[l[[i]],]
        mspI <- as.matrix(mspI)
        class(mspI) <- "numeric"
        #mspI <- data.matrix(mspI)
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
    
    ## scale back to percent
    mm <- apply(mm, 1, function(x) {x / max(x) * 100})
    if (length(frag) > 1)  {
        mm <- t(mm)
        
    } 
    
    ## if there is only one fragment do the following
    if (length(frag) == 1) {
        mm <- as.matrix(mm)
    }
    
    colnames(mm) <- bins
    
    return(mm)
}
