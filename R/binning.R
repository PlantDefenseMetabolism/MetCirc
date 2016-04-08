#' @name getPrecursorMZ
#' @title Get precursor m/z values
#' @description Get precursor m/z values of a data.frame in msp format
#' @usage getPrecursorMZ(msp)
#' @param msp data.frame in msp format, see ?convert2MSP for further information
#' @details Internal use to retrieve precursor m/z values.
#' @return getPrecursorMZ returns a character vector with all precursor 
#' values
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' finalMSPdf <- getMSP(finalMSP)
#' getPrecursorMZ(finalMSPdf)
#' @export
getPrecursorMZ <- function (msp) {
    ## get indices with precursor mz
    IndPrecMZ <- which(msp[,1] == "PRECURSORMZ: ")
    ## get precursor mz
    precmz <- msp[IndPrecMZ,2]
    
    ## change to numeric
    precmz <- as.numeric(precmz)
    
    return(precmz)
}
#' @name getRT
#' @title Get precursor RT values
#' @description Get precursor RT values of a data.frame in msp format
#' @usage getRT(msp)
#' @param msp data.frame in msp format, see ?convert2MSP for further information
#' @details Internal use to retrieve retention time values.
#' @return getRT returns a character vector with all retention time values
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' finalMSPdf <- getMSP(finalMSP)
#' getRT(finalMSPdf)
#' @export
getRT <- function (msp) {
    ## get indices with rt
    IndRT <- which(msp[,1] == "RETENTIONTIME: ")
    ## get rt 
    rt <- msp[IndRT,2]
    
    ## change to numeric
    rt <- as.numeric(rt)
    
    return(rt)
}


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
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
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
#' @usage binning(msp, tol = 0.01, group)
#' @param msp data.frame in msp format, see ?convert2MSP for further information
#' @param tol numerical, boundary value until which neighboured peaks will be 
#'      joined together
#' @param group character vector, to which group does the entry belong to
#' @details The functions bins fragments together by calculating 
#' @return binning returns a matrix where rownames are precursor ions
#' (m/z / retention time) and colnames are newly calculated m/z values which 
#' were binned. The algorithm which is currently implemented joins the 
#' two nearest m / z values (here the m / z values can also be joined by 
#' several fragment ions) and recalculates the new m / z by weighing for the 
#' number of m / z fragments.
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples data("idMSMStoMSP", package = "MetCirc")
#' group <- sample(c("yl", "ol", "s","r"), size = length(finalMSP), replace=TRUE) 
#' binning(msp = finalMSP, tol = 0.01, group = group)
#' @export
binning <- function(msp, tol = 0.01, group = NULL) { 
    
    if (!is(msp) == "MSP") stop("msp is not of class MSP.")
    
    msp <- getMSP(msp)
    
    if (is.null(group)) {
        print("argument group is not specified, will create dummy group")
        group <- rep("a", length(getPrecursorMZ(msp)))
    }
    ## msp is .msp file
    ## tol is tolerance value for binning
    precmz <- getPrecursorMZ(msp)
    rt <- getRT(msp)
    
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
    ## sort frag and get order
    frag_s <- sort(frag)
    frag_order <- order(frag)
    
    
    ## calculate distance to neighbours  
    dist <- list(NA)
    for (i in 1:length(frag_s)) {
        if (i != length(frag_s)) {
            dist[[i]] <- c(frag_s[i], frag_s[i+1] - frag_s[i])
        } else dist[[i]] <- c(frag_s[i], Inf)
    }
    distAdapt <- dist
    
    ## get distance values
    dist2Adapt <- dist2 <- lapply(dist, "[[", 2)
    indConv <- which(unlist(dist2) == 0) ## get indices which have distance of 0
    
    mapping <- lapply(1:length(dist), function(x) x)

    ## map to the next element if it has same mz
    for (i in 1:length(indConv)) mapping[[indConv[i]]] <- indConv[i] + 1
    
    ## conv is the vector which shows convoluted mz
    conv <- numeric(length(mapping)) 
    x <- 1
    
    ## for mz values which have distance of 0 create a identifier Mx, where x 
    ## is an increasing number to be able to trace back same mz
    for (i in 1:length(mapping)) {
        if (mapping[i] != i & conv[i] == "0") {
            conv[i] <- paste("M", x, sep="")
            conv[i+1] <- paste("M", x, sep="")
            j <- i
            ## check when there is a sequence of mz which have distance 0 and 
            ## allocate then the identical mx
            while (unlist(mapping)[j+1] != (j +1) ) { 
                conv[j + 1] <- paste("M", x, sep="")
                conv[j + 2] <- paste("M", x, sep="")
                j <- j +1 
            }
            i <- j
            x <- x +1 
        }
    }
    
    ## actual binning script starts here 
    indGreater0 <- which(unlist(dist2Adapt) > 0) ## get distances greater 0
    ## find smallest distance which is greater than zero
    minDist2AdaptGreater0 <- which.min(dist2Adapt[indGreater0])
    indGreater0Min <- indGreater0[minDist2AdaptGreater0]
    
    while (distAdapt[indGreater0Min][[1]][2] < tol) {
        ## write all which have Mx value to Mx+1
        if (conv[indGreater0Min + 1 ] == 0) { ## then create new Mx
            str <- unlist(strsplit(unique(conv),split="M"))
            str <- str[which(str != "")]
            if (0 %in% str) str <- str[which(str != "0")]
            str <- max(as.numeric(str)) + 1
            str <- paste("M", str, sep="")
            if (conv[indGreater0Min] == 0) {
                conv[indGreater0Min] <- str
            } else {
                conv[which(conv[indGreater0Min] == conv)] <- str
            }
            conv[indGreater0Min + 1] <- str   
        } else { ## if conv[indGreater0min + 1] != 0, i.e. if it is Mx, 
            ## then use "old" Mx
            if (conv[indGreater0Min] == 0) {
                conv[indGreater0Min] <- conv[indGreater0Min + 1 ]
            } else {
                conv[which(conv[indGreater0Min] == conv)] <- conv[indGreater0Min + 1 ]}
        }
        ## calculate new mean for all instances with Mx+1
        indAdapt <- which(conv[indGreater0Min + 1 ] == conv)
        newMean <- mean(unlist(lapply(distAdapt[indAdapt], "[", 1)))
        ## write new mean to all instances with Mx+1
        for (i in indAdapt) distAdapt[[i]][1] <- newMean
        ## calculate new distances and write new distances
        distAdaptOld <- distAdapt
        for (i in 1:length(distAdapt)) { ## calculate for all elements in the 
            ## list (this can be changed, so that we only calculate distance for 
            ## elements before and after Mx+1)
            if (i != length(distAdapt)) {
                distAdapt[[i]] <- c(distAdaptOld[[i]][1], distAdaptOld[[i+1]][1] - distAdaptOld[[i]][1])
            } else distAdapt[[i]] <- c(distAdaptOld[[i]][1], Inf)
        }
        dist2Adapt <- lapply(distAdapt, "[[", 2)
        unlist(dist2Adapt)
        
        indGreater0 <- which(unlist(dist2Adapt) > 0) 
        minDist2AdaptGreater0 <- which.min(dist2Adapt[indGreater0])
        indGreater0Min <- indGreater0[minDist2AdaptGreater0]
    }
    ## actual binning script ends here
    
    ## write for every conv which has "0" a new Mx 
    for (i in which(conv == "0")) {
        str <- unlist(strsplit(unique(conv),split="M"))
        str <- str[which(str != "")]
        if (0 %in% str) str <- str[which(str != "0")] 
        ## find highest x and create new one (+1)
        str <- max(as.numeric(str)) + 1 
        str <- paste("M", str, sep="")
        conv[i] <- str ## allocate Mx+1 to conv[i]
    }
    
    ## find all unique bins, these will be the colnames of mm
    uniqueMZ <- unlist(lapply(distAdapt, "[", 1))
    uniqueMZ <- unique(uniqueMZ)
    
    mm <- matrix(data = 0, nrow = length(precmz), ncol = length(uniqueMZ))
    ## convoluted MZ is column names
    colnames(mm) <- uniqueMZ
    rownames(mm) <- paste(precmz, rt, sep="/")
    
    fragMM <- unlist(l)
    ## write to mm 
    
    
    ## new
    for (i in 1:length(l)) {
        entry <- l[[i]]    
        for (j in 1:length(entry)) {
            entryJ <- as.numeric(msp[entry[j],])
            colIND <- which.min(abs(as.numeric(entryJ[1]) - uniqueMZ))
            intensity <- entryJ[2]
            mm[i, colIND] <- intensity
        }
    }
    
    ## new ende
    
    ##IndPrecMZ <- match(precmz, msp[,2])
    
#     for (i in 1:length(fragMM)) {
#         ## which frag is put first? second? ...
#         indFragMM <- fragMM[frag_order][i]
#         ## get corresponding Precursor MZ
#         correspPrecMZ <- precmz[max(which(IndPrecMZ < indFragMM))]
#         ## get corresponding rt
#         correspPrecRT <- rt[max(which(IndPrecMZ < indFragMM))]
#         ## get unique row identifier
#         uniqueIden <- paste(correspPrecMZ, correspPrecRT, sep="/")
#         rowInd <- which(uniqueIden == rownames(mm))
#         ## get col index 
#         colInd <- which(distAdapt[[i]][1] == colnames(mm))
#         ## write 
#         mm[rowInd,colInd] <- msp[indFragMM, 2]
#     }

    rownames(mm) <- paste(group, rownames(mm), sep="_")
    
    class(mm) <- "numeric"
    
    #rownames(mm) <- paste(group, rNames, sep="_")
    ## convoluted MZ is column names
    
    ## was rownames(mm) <- paste(group, sprintf("%04d", 1:length(rNames)), rNames, sep="_")
    return(mm)
}




