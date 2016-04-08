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
#' @export
cutUniquePreMZ <- function(precursor, splitPattern = splitPattern, 
                            splitInd = splitInd, returnCharacter = TRUE) {
    ## split precursors according to split pattern
    precursor <- as.character(precursor)
    precursor <- unique(precursor)
    splitPrecursor <- strsplit(precursor, split = splitPattern)
    ## extract precursor mz at position splitInd
    splitPrecursor <- lapply(splitPrecursor,"[", splitInd)
    PrecursorMZ <- unlist(splitPrecursor)
    lenPreMZ <- length(PrecursorMZ)
    
    ## change character to numeric
    if (!returnCharacter)
        PrecursorMZ <- as.numeric(PrecursorMZ)
    ##uniquePreMZ <- unique(precursor)
    ##lenUniquePreMZ <- length(uniquePreMZ)
    ##uniquePreMZ_cut <- unique(PrecursorMZ)
    
    
    return(PrecursorMZ) ## return(uniquePreMZ_cut)
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
#' @return convert2MSP returns an object of class MSP
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' convert2MSP(mm = sd02_deconvoluted, splitPattern = "_", splitInd = 1)
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
    finalMSP <- matrix(data = NA, nrow = 8 * lenUniquePreMZ + dim(mm)[1], 
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
            c("METABOLITECLASS: ", "Unknown"),
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

    return(new("MSP", msp = finalMSP))
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
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' finalMSPNL <- msp2FunctionalLossesMSP(msp = finalMSP)
#' @export
msp2FunctionalLossesMSP <- function(msp) {
    
    if (!is(msp) == "MSP") stop("msp is not of class MSP")
    
    msp <- getMSP(msp)
    
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
            c("METABOLITECLASS: ", "Unknown"),
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

    
    return(MSP(msp = finalMSP))
}

#' @import methods
NULL

#' @name MSP
#' @title MSP-class
#' @aliases MSP-class
#' @description MSP class for msp data.frame. Allows easy computation of 
#' length of entries by entering length(msp), where msp is of class MSP.
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @param msp a data.frame in msp format
#' @export
MSP <- setClass("MSP", slots = c(msp = "data.frame"))

#' @name length
#' @rdname length-method
#' @aliases length,MSP-method
####usage length(x)
#' @title length method for MSP class
#' @return numerical
#' @description Gives the number of entries in the MSP object.
#' @param x object of class MSP
#' @docType methods
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' length(finalMSP)
#' @export
setMethod("length", signature = "MSP", 
          definition = function(x) {
              length(getPrecursorMZ(x@msp))
})

#' @name show
#' @rdname show-method
#' @aliases show,MSP-method
#' @title show method for MSP class
#' @return character
#' @description Prints information on the MSP class (number of entries).
#' @param object object of class MSP
#' @docType methods
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' show(finalMSP)
#' @export
setMethod("show", signature = "MSP", 
          definition = function(object) {
              cat("An object of class", class(object), "with", 
                  length(getPrecursorMZ(object@msp)), "entries.", sep = " ")
})

#' @name getMSP
#' @aliases getMSP,MSP-method
#' @title getMSP method for MSP class
#' @return data.frame
#' @description Returns the data.frame entry of an MSP object.
#' @param object object of class MSP
#' @docType methods
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' getMSP(finalMSP)
#' @export
setGeneric("getMSP", function(object) standardGeneric("getMSP"))


#' @export
#' @describeIn getMSP returns the data.frame of an MSP object
setMethod("getMSP", signature = "MSP", definition = function(object) {object@msp})

#' @name combine
#' @aliases combine,MSP-method
#' @title combine method for MSP class
#' @return MSP object
#' @description Combines two objects of class MSP.
#' @param object1 object of class MSP
#' @param object2 object of class MSP
#' @docType methods
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP1 <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' finalMSP2 <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' combine(finalMSP1, finalMSP2)
#' @export
setGeneric("combine", function(object1, object2) standardGeneric("combine"))


#' @export
#' @describeIn combine combines two MSP objects
setMethod("combine", signature = c("MSP", "MSP"), definition = function(object1, object2) {
    new("MSP", msp = rbind(object1@msp, object2@msp))})

#' @name getNames
#' @aliases getNames,MSP-method
#' @title getNames returns names of compounds in MSP object
#' @return character
#' @description getNames returns names of compounds in MSP object.
#' @param object object of class MSP
#' @docType methods
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' getNames(finalMSP)
#' @export
setGeneric("getNames", function(object) standardGeneric("getNames"))


#' @describeIn getNames returns names of compounds in MSP objects
#' @export
getNames <- function(object) {
    ## get classes of compounds in an msp object
    df <- object@msp
    ind <- which(df[,1] == "METABOLITENAME: ")
    return(df[ind,2])
}

#' @name getMetaboliteClass
#' @aliases getMetaboliteClass,MSP-method
#' @title getMetaboliteClass returns names of compounds in MSP object
#' @return character
#' @description getMetaboliteClass returns names of compounds in MSP object.
#' @param object object of class MSP
#' @docType methods
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' getMetaboliteClass(finalMSP)
#' @export
setGeneric("getMetaboliteClass", 
           function(object) standardGeneric("getMetaboliteClass"))

#' @describeIn getMetaboliteClass returns class names of compounds in MSP 
#' objects
#' @export
getMetaboliteClass <- function(object) {
    ## get classes of compounds in an msp object
    df <- object@msp
    ind <- which(df[,1] == "METABOLITECLASS: ")
    return(df[ind,2])
}

#' @name [
#' @aliases [,MSP,numeric,missing,missing-method
#' @title Extract parts of a MSP object
#' @return MSP object
#' @description [ operator acting on an MSP object to extract parts. 
#' @param x object of class MSP
#' @param i numeric
#' @param j missing
#' @param drop missing
#' @docType methods
#' @rdname extract-methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' finalMSP[1]
#' @export
setMethod("[", 
    signature(x = "MSP", i = "numeric", j = "missing", drop = "missing"), 
    definition = function(x, i, j = "missing", drop = "missing") {
        if (max(i) > length(x)) stop("max(i) greater than length(x)")
        start <- getBegEndIndMSP(x@msp)[[1]] - 7 ## which(testMSP@msp[,1] == "NAME: "), indices of 'NAME: '
        end <- getBegEndIndMSP(x@msp)[[2]] + 1
        start <- start[i]
        end <- end[i]
        inds <- lapply(1:length(start), function(j) start[j]:end[j])
        inds <- unlist(inds)
        
        return(new("MSP", msp = x@msp[inds, ]))
    }
)
