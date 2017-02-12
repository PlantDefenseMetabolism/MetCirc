#' @import methods
NULL

#' @name MSP
#' @title MSP-class
#' @aliases MSP-class
#' @description Definiton of \code{MSP}-class in \code{MetCirc}. Entries are 
#' MS/MS features including their spectra. 
#' Allows easy computation of number of entries by entering length(msp), 
#' where msp is of class MSP. The \code{MSP}-class incorporates accessors 
#' for auxiliary information of MS/MS features (names, classes, information, 
#' adduct ion name). 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @export
MSP <- setClass("MSP", 
                slots = c(msp = "data.frame", mz = "numeric", rt = "numeric", 
                          names = "character", classes = "character", 
                          information = "character", adduct = "character"),
                validity=function(object) 
                {
                    #length_msp <- which(object@msp[,1] == "Num Peaks: " | object@msp[,1] == "Num Losses: ")
                    #length_msp <- length(length_msp)
                    if(length(which(object@msp[,1] == "Num Peaks: " | object@msp[,1] == "Num Losses: ")) != length(object@mz)) {
                        return("error with length (mz)")
                    }
                    if(length(object@mz) != length(object@rt)) {
                        return("error with length (rt)")
                    }
                    if(length(object@rt) != length(object@names)) {
                        return("error with length (names)")
                    }
                    if(length(object@names) != length(object@classes)) {
                        return("error with length (classes)")
                    }
                    if(length(object@classes) != length(object@information)) {
                        return("error with length (information)")
                    }
                    if(length(object@information) != length(object@adduct)) {
                        return("error with length (adduct)")
                    }
                    
                }
)

#' @name cutUniquePrecursor
#' @title Get unique precursor ions
#' @description Get unique precursor ions
#' @usage cutUniquePrecursor(precursor, splitPattern = splitPattern, 
#'      splitInd = splitInd, returnCharacter = TRUE)
#' @param precursor  \code{character} where features are separated by 
#'  \code{splitPattern}
#' @param splitPattern \code{character}, character vector to use for splitting, 
#'      see \code{?strsplit} for further information
#' @param splitInd \code{numeric}, extract precursor mz at position 
#'  \code{splitInd}
#' @param returnCharacter \code{logical}, if \code{TRUE} return 
#'  \code{character}, if \code{FALSE} return \code{numeric}
#' @details Function for internal usage.
#' @return The function \code{cutUniquePrecursor} returns \code{character} or \code{numeric} as 
#'  specified by parameters.
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' precursor <- "A_269.0455469_-1"
#' splitPattern <- "_"
#' splitInd <- 2
#' cutUniquePrecursor(precursor, splitPattern = splitPattern, 
#'      splitInd = splitInd, returnCharacter = TRUE)
#' @export
cutUniquePrecursor <- function(precursor, splitPattern = splitPattern, 
                            splitInd = splitInd, returnCharacter = TRUE) {
    ## split precursors according to split pattern
    precursor <- as.character(precursor)
    precursor <- unique(precursor)
    splitPrecursor <- strsplit(precursor, split = splitPattern)
    ## extract precursor mz at position splitInd
    splitPrecursor <- lapply(splitPrecursor,"[", splitInd)
    Precursor <- unlist(splitPrecursor)
    lenPre <- length(Precursor)
    
    ## change character to numeric
    if (!returnCharacter)
        Precursor <- as.numeric(Precursor)
    
    return(Precursor) 
}

#' @name convert2MSP
#' @title Convert deconvoluted matrix into \code{MSP}-object
#' @description Convert deconvoluted matrix into \code{MSP}-object
#' @usage convert2MSP(mm, splitPattern = "_", splitIndMZ = 1, splitIndRT = NULL, 
#'  rt = FALSE, names = FALSE, information = FALSE, classes = FALSE, adduct = FALSE)
#' @param mm \code{matrix}, \code{mm} has to have three columns with colnames 
#'  "mz", "intensity" and "id" (order is not important). The column comprises 
#'  information about the precursor ion which will be assessed by 
#'  \code{splitPattern} and \code{splitInd}. Optionally, \code{mm} can have 
#'  colnames "rt", "names", "information", "classes" and "adduct". 
#' @param splitPattern \code{character}, \code{splitPattern} is the pattern 
#'      which separates elements and precursor m/z
#' @param splitIndMZ \code{numeric}, the position of the precursor m/z in the 
#'      character string concerning separation by \code{splitPattern}
#' @param splitIndRT \code{numeric} or \code{NULL}, the position of the 
#'      retention time in the \code{character} string concerning separation by 
#'      \code{splitPattern}, if \code{NULL} the retention time will be the mean 
#'      of all retention time values of the MS/MS feature fragments
#' @param rt \code{logical}, should retention times be retrieved? If set to 
#'  \code{TRUE}, \code{convert2MSP} will access the column "rt" in \code{mm} 
#'  which contains the retention time values for each fragment when 
#'  \code{splitIndRT} is \code{NULL}, if \code{rt} is set to \code{TRUE} and 
#'  \code{splitIndRT} is \code{numeric}, \code{convert2MSP} will access
#'  the column "id" to get the retention time at position \code{splitIndRT}
#'  when splitting with \code{splitPattern}
#' @param names \code{logical}, should names be retrieved? If set to 
#'  \code{TRUE}, \code{convert2MSP} will access the column "names" in \code{mm}
#'  which contains the names of the metabolites
#' @param information \code{logical}, should further information of metabolites 
#'  be retrieved? If set to \code{TRUE}, \code{convert2MSP} will access the 
#'  column "information" in \code{mm} which contains information about the 
#'  metabolites
#' @param classes \code{logical}, should classes of metabolites be retrieved? If 
#'  set to \code{TRUE}, \code{convert2MSP} will access the column "classes" in 
#'  \code{mm} which contains the names of the metabolites
#' @param adduct \code{logical}, should adduct ion names of metabolites be 
#'  retrieved? If set to \code{TRUE}, \code{convert2MSP} will access the column 
#'  "adduct" in \code{mm} which contains the adduct ion names of the metabolites
#' @details The function \code{convert2MSP} creates a data entry for each precursor ion. 
#' Each entry in the return object has the following information: Num Peaks and 
#' a list of fragments together with their intensities; it will further contain
#' information on m/z values of the precursor ion, the retention time, 
#' metabolite names, classes, adduct ion name and further information. 
#' \code{convert2MSP} will access the columns "rt", "names", "information", 
#' "classes" and "adduct", respectively, if arguments are set to \code{TRUE}. The 
#' column "id" has to contain a unique identifier for each MS/MS feature. It is 
#' obligatory that each element in the column "id" contains the precursor 
#' m/z value, but may contain furhter elements (e.g. peak correlation value or 
#' retention time of the precursor ion). Information about the m/z value will be 
#' assessed by \code{splitPattern} and \code{splitInd}. E.g. items in the column 
#' "id" can be in the form of "1_163.23", which has to be accessed by setting 
#' \code{splitPattern = "_"} and \code{splitInd = 2} to access the m/z value of 
#' the precursor ion (here: 162.23). If \code{rt} is set to \code{TRUE} and 
#' \code{splitIndRT} is \code{NULL}, \code{convert2MSP} will access the column 
#' "rt" to get the retention time values corresponding to each fragment and 
#' calculate the mean value, if \code{rt} is set to \code{TRUE} and 
#' \code{splitIndRT} \code{numeric}, \code{convert2MSP} will retrieve the 
#' retention time value from column "id". 
#' @return \code{convert2MSP} returns an object of class \code{MSP}
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' convert2MSP(mm = sd02_deconvoluted, splitPattern = " _ ", splitIndMZ = 2, 
#'  splitIndRT = NULL, rt = FALSE, names = FALSE, information = FALSE, 
#'  classes = FALSE, adduct = FALSE)
#' @export
convert2MSP <- function (mm, splitPattern = "_", splitIndMZ = 1, 
    splitIndRT = NULL, rt = FALSE, names = FALSE, information = FALSE, 
    classes = FALSE, adduct = FALSE) {
    
    colNames <- colnames(mm)
    if (!("mz" %in% colNames)) stop("no column 'mz' found")
    if (rt & !is.numeric(splitIndRT) & !("rt" %in% colNames)) 
        stop("no column 'rt' found")
    if (!("intensity" %in% colNames)) stop("no column 'intensity' found")
    if (names & !("names" %in% colNames)) stop("no column 'names' found")
    if (information & !("information" %in% colNames)) 
                                        stop("no column 'information' found")
    if (adduct & !("adduct" %in% colNames)) stop("no column 'adduct' found")
    if (classes & !("classes" %in% colNames)) stop ("no column 'classes' found")
    
    if (rt & (!is.numeric(splitIndRT) & !is.null(splitIndRT))) 
                                stop("splitIndRT is not numeric and not NULL")
    ## if (colNames[4] != "pcgroup_precursorMZ") break
     
    precursor <- mm[,"id"]
    precursor <- as.character(precursor)
    
    uniquePre <- unique(precursor)
    uniquePreMZ_cut <- cutUniquePrecursor(precursor = precursor, 
        splitPattern = splitPattern, splitInd = splitIndMZ)
    ## get Retention time, mz_rt_pcgroup
    #if (!is.null(splitIndRT)) {
    #    uniquePreRT_cut <- cutUniquePrecursor(precursor = precursor, 
    #                    splitPattern = splitPattern, splitInd = splitIndRT)        
    #}

    
    ## check if pcgroup_grecursorMZ is in fourth column and unique precursor
    ## mz were assessed correctly
    if (any(is.na(as.numeric(uniquePreMZ_cut))))
            stop("unique precursors are not in the fourth column or 
                splitPattern/splitInd was chosen incorrectly")
    
    #lenUniquePreMZ <- length(uniquePreMZ_cut)
    
    ## access columns names, information, classes, adduct to retrieve information 
    ## about names, information, metabolite classes or adduct
    if (names) namesMM <- mm[, "names"]
    if (classes) classesMM <- mm[, "classes"]
    if (information) informationMM <- mm[, "information"]
    if (adduct) adductMM <- mm[, "adduct"]
    if (rt) rtMM <- mm[, "rt"]

    ## create data frame for MSP file
    lenUniquePre <- length(uniquePre)
    finalMSP <- matrix(data = NA, nrow = 2 * lenUniquePre + dim(mm)[1], 
            ncol = 2) ## 7 new entries + all fragment ion entries
    finalMSP <- as.data.frame(finalMSP)
    
    
    NAMES <- INFORMATION <- CLASS <- ADDUCT <- character(length(uniquePre))
    RETTIME <- PRECMZ <- numeric(length(uniquePre))
    
    ## write to data frame
    for (i in 1:lenUniquePre) {
        ind <- which(uniquePre[i] == precursor)   
        NAMES[i] <-  if (names) {
            unique(as.character(namesMM[ind])[1])} else "Unknown"
        RETTIME[i] <- if (rt) {
            if (!is.null(splitIndRT)) {as.numeric(
                cutUniquePrecursor(precursor[ind], splitPattern, splitIndRT))[1]
            } else mean(rtMM[ind])
        } else NaN
        PRECMZ[i] <- as.numeric(
            cutUniquePrecursor(precursor[ind], splitPattern, splitInd = splitIndMZ))
        INFORMATION[i] <- if (information) {
            unique(as.character(informationMM[ind])[1])} else "Unknown"
        CLASS[i] <- if (classes) {
            unique(as.character(classesMM[ind])[1])} else "Unknown"
        ADDUCT[i] <- if (adduct) {
            unique(as.character(adductMM[ind])[1])} else "Unknown"
        
        entry <- rbind(
            c("Num Peaks: ", length(ind)),
            mm[ind, c("mz", "intensity")],
            c(" ", " ")
        )
        entry <- as.matrix(entry)
        ## determine first empty line
        newstart <- which(is.na(finalMSP[,1]))[1]
        ## determine last line to write to
        newend <- newstart + dim(entry)[1] - 1
        finalMSP[newstart:newend,] <- entry
    }
    msp <- MSP(msp=finalMSP, mz = PRECMZ, rt = RETTIME, names = NAMES, 
        classes = CLASS, information = INFORMATION, adduct = ADDUCT)

    return(msp)
}

#' @name msp2FunctionalLossesMSP
#' @title Convert MSP to MSP with functional losses
#' @description \code{msp2FunctionalLossesMSP} converts a \code{MSP}-object 
#' (with fragments) into a \code{MSP}-object with neutral losses
#' @usage msp2FunctionalLossesMSP(msp)
#' @param msp \code{MSP}-object
#' @details The function \code{msp2FunctionalLosses} can be used when calculating
#' the similarity based on neutral losses instead of fragments.
#' @return \code{msp2FunctionalLossesMSP} returns a \code{MSP-object}  
#' (with neutral losses)
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", 
#'                      splitIndMZ = 2, splitIndRT = NULL)
#' finalMSPNL <- msp2FunctionalLossesMSP(msp = finalMSP)
#' @export
msp2FunctionalLossesMSP <- function(msp) {
    
    if (!is(msp) == "MSP") stop("msp is not of class MSP")
    
    ## do the following before getMSP
    PRECMZ <- msp@mz
    RETTIME <- msp@rt
    NAMES <- msp@names
    CLASSES <- msp@classes
    INFORMATION <- msp@information
    ADDUCT <- msp@adduct
    
    msp <- msp@msp
    
    indices <- getBegEndIndMSP(msp)
    indBegL <- indices[[1]]
    indEndL <- indices[[2]]
    ## create data frame for MSP file
    finalMSP <- matrix(data = NA, nrow = dim(msp)[1], ncol = 2) 
    finalMSP <- as.data.frame(finalMSP)
    
    ## create MSP from 
    for (i in 1:length(PRECMZ)) {
        
        indBeg <- indBegL[i]
        indEnd <- indEndL[i]
        
        neutralL <- PRECMZ[i] - as.numeric(as.character(msp[indBeg:indEnd,1]))
        neutralL <- abs(neutralL)
        
        entry <- rbind(
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

    msp <- MSP(msp=finalMSP, mz = PRECMZ, rt = RETTIME, names = NAMES, 
               classes = CLASSES, information = INFORMATION, adduct = ADDUCT)
}


#' @name length
#' @rdname length-method
#' @aliases length,MSP-method
#' @title \code{length} method for \code{MSP}-class
#' @return \code{numeric}
#' @description Gives the number of entries in the \code{MSP} object.
#' @param x object of class \code{MSP}
#' @docType methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = NULL)
#' length(finalMSP)
#' @export
setMethod("length", signature = "MSP",
          definition = function(x) {length(getPrecursorMZ(x))
})

#' @name show
#' @rdname show-method
#' @aliases show,MSP-method
#' @title \code{show} method for \code{MSP}-class
#' @return \code{character}
#' @description \code{show} prints information on the \code{MSP}-object (number of entries).
#' @param object object of class \code{MSP}
#' @docType methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                      splitIndMZ = 2, splitIndRT = NULL)
#' show(finalMSP)
#' @export
setMethod("show", signature = "MSP",
          definition = function(object) {
              cat("An object of class", class(object), "with",
                  length(getPrecursorMZ(object)), "entries.", sep = " ")
})

#' @name peaks
#' @aliases peaks,MSP-method
#' @title \code{peaks} method for \code{MSP}-class
#' @return \code{data.frame}
#' @description \code{peaks} returns the \code{data.frame} entry with peak
#'  information of an \code{MSP} object.
#' @param object object of class \code{MSP}
#' @docType methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                      splitIndMZ = 2, splitIndRT = NULL)
#' peaks(finalMSP)
#' @export
setGeneric("peaks", function(object) standardGeneric("peaks"))

#' @describeIn peaks returns the \code{data.frame} of an \code{MSP}-object
#' @export
setMethod("peaks", signature = "MSP", definition = function(object) object@msp)

#' @name combine
#' @aliases combine,MSP-method
#' @title \code{combine} method for \code{MSP}-class
#' @return \code{MSP}-object
#' @description \code{combine} combines two objects of class \code{MSP}.
#' @param object1 object of class \code{MSP}
#' @param object2 object of class \code{MSP}
#' @docType methods
#' @examples data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP1 <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = NULL)
#' finalMSP2 <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = NULL)
#' combine(finalMSP1, finalMSP2)
#' @export
setGeneric("combine", function(object1, object2) standardGeneric("combine"))


#' @export
#' @describeIn combine combines two object of class \code{MSP}
setMethod("combine", signature = c("MSP", "MSP"),
        definition = function(object1, object2) {
            new("MSP", msp=rbind(object1@msp, object2@msp),
                mz = c(object1@mz, object2@mz),
                rt = c(object1@rt, object2@rt),
                names = c(object1@names, object2@names),
                classes = c(object1@classes, object2@classes),
                information = c(object1@information, object2@information),
                adduct = c(object1@adduct, object2@adduct))
            }
)

#' @name names
#' @aliases names,MSP-method
#' @title \code{names} returns names in \code{MSP}-object
#' @return \code{character}
#' @description \code{names} returns names in \code{MSP}-object.
#' @param x object of class \code{MSP}, see \code{?convert2MSP} for further
#' information
#' @docType methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = NULL)
#' names(finalMSP)
#' @export
setMethod(f = "names", signature = "MSP", definition = function(x) x@names)

#' @name names<-
#' @aliases names<-,MSP,character-method
#' @title \code{names<-} sets names in \code{MSP}-object
#' @return \code{MSP}-object
#' @description \code{names<-} sets names in \code{MSP}-object
#' @param x object of class \code{MSP}, see \code{?convert2MSP} for further
#' information
#' @param value \code{character} vector with new names
#' @docType methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = NULL)
#' names(finalMSP) <- rep("Unknown")
#' @export
setMethod(f = "names<-", signature = c("MSP", "character"),
        definition = function(x, value) {
              x@names <- value
              ##validObject(x)
              x
        }
)

#' @name information
#' @aliases information,MSP-method
#' @title \code{information} returns information of metabolites in
#'  \code{MSP}-object
#' @return \code{character}
#' @description \code{information} returns information in \code{MSP}-object.
#' @param x object of class \code{MSP}, see \code{?convert2MSP} for further
#'  information
#' @docType methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = NULL)
#' information(finalMSP)
#' @export
setGeneric("information", function(x) standardGeneric("information"))

#' @describeIn information returns information of metabolites in
#' \code{MSP}-object
#' @export
information <- function(x) x@information

#' @name information<-
#' @aliases information<-,MSP,character-method
#' @title \code{information<-} sets information in \code{MSP}-object
#' @return \code{MSP}-object
#' @description \code{information<-} sets information in \code{MSP}-object
#' @param x object of class \code{MSP}, see \code{?convert2MSP} for further
#' information
#' @param value character vector with new information
#' @docType methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = NULL)
#' information(finalMSP) <- rep("Unknown")
#' @export
setGeneric("information<-", function(x, value) standardGeneric("information<-"))

#' @describeIn information<- sets information of metabolites in
#' \code{MSP}-object
#' @export
"information<-" <- function(x, value) {
    x@information <- value
    ##validObject(x)
    x
}

#' @name classes
#' @aliases classes,MSP-method
#' @title \code{classes} returns class names of compounds in \code{MSP}-object
#' @return \code{character}
#' @description \code{classes} returns class names of compounds in \code{MSP}-object.
#' @param x object of class \code{MSP}
#' @docType methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = NULL)
#' classes(finalMSP)
#' @export
setGeneric("classes", function(x) standardGeneric("classes"))

#' @describeIn classes returns class names of metabolites in \code{MSP}-object
#' @export
classes <- function(x) x@classes

#' @name classes<-
#' @aliases classes<-,MSP,character-method
#' @title \code{classes<-} sets information in \code{MSP}-object
#' @return \code{MSP}-object
#' @description \code{classes<-} sets information in \code{MSP}-object.
#' @param x object of class \code{MSP}, see \code{?convert2MSP} for further
#' information
#' @param value character vector with new classes
#' @docType methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = NULL)
#' classes(finalMSP) <- rep("Unknown")
#' @export
setGeneric("classes<-", function(x, value) standardGeneric("classes<-"))

#' @describeIn classes<- sets information of metabolites in \code{MSP}-object
#' @export
"classes<-" <- function(x, value) {
    x@classes <- value
    ##validObject(x)
    x
}


#' @name adduct
#' @aliases adduct,MSP-method
#' @title \code{adduct} returns adduct ion names of compounds in \code{MSP}-object
#' @return \code{character}
#' @description \code{adduct} returns adduct ion names of compounds in
#' \code{MSP}-object.
#' @param x object of class \code{MSP}
#' @docType methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = NULL)
#' adduct(finalMSP)
#' @export
setGeneric("adduct", function(x) standardGeneric("adduct"))

#' @describeIn adduct returns adduct ion names of compounds in
#'  \code{MSP}-objects
#' @export
adduct <- function(x) return(x@adduct)


#' @name adduct<-
#' @aliases adduct<-,MSP,character-method
#' @title \code{adduct<-} sets adduct ion names in \code{MSP}-object
#' @return \code{MSP}-object
#' @description \code{adduct<-} sets adduct ion names in \code{MSP}-object
#' @param x object of class \code{MSP}, see \code{?convert2MSP} for further
#'  information
#' @param value \code{character} vector with new adduct ion names
#' @docType methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = NULL)
#' adduct(finalMSP) <- rep("Unknown")
#' @export
setGeneric("adduct<-", function(x, value) standardGeneric("adduct<-"))

#' @describeIn adduct<- sets adduct ion names of metabolites in
#' \code{MSP}-object
#' @export
"adduct<-" <- function(x, value) {
    x@adduct <- value
    ##validObject(x)
    x
}

#' @name getRT
#' @aliases getRT,MSP-method
#' @title \code{getRT} returns precursor RT values of an \code{MSP}-object
#' @return \code{numeric}
#' @description \code{getRT} returns a \code{numeric} vector with all retention
#'  time values
#' @param x object of class \code{MSP}
#' @docType methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                      splitIndMZ = 2, splitIndRT = NULL)
#' getRT(finalMSP)
#' @export
setGeneric("getRT", function(x) standardGeneric("getRT"))

#' @describeIn getRT returns precursor RT values of an \code{MSP}-object
#' @export
getRT <- function(x) x@rt

#' @name getPrecursorMZ
#' @aliases getPrecursorMZ,MSP-method
#' @title \code{getPrecursorMZ} returns precursor m/z values of an
#'  \code{MSP}-object
#' @return \code{numeric}
#' @description \code{getPrecursorMZ} returns a \code{numeric} vector with
#'  precursor m/z values
#' @param x object of class \code{MSP}
#' @docType methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = NULL)
#' getPrecursorMZ(finalMSP)
#' @export
setGeneric("getPrecursorMZ", function(x) standardGeneric("getPrecursorMZ"))

#' @describeIn getPrecursorMZ returns precursor m/z values of an MSP object
#' @export
getPrecursorMZ <- function(x) x@mz

#' @name [
#' @aliases [,MSP,numeric-method
#' @title Extract parts of a \code{MSP}-object
#' @return \code{MSP}-object
#' @description [ operator acting on an \code{MSP}-object to extract parts.
#' @param x object of class \code{MSP}
#' @param i \code{numeric}
#' @docType methods
#' @rdname extract-methods
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = NULL)
#' finalMSP[1]
#' @export
setMethod("[",
    signature(x = "MSP", i = "numeric"), ## j = "missing", drop = "missing"),
    definition = function(x, i) {#}, j = "missing", drop = "missing") {
        if (max(i) > length(x)) stop("max(i) greater than length(x)")
        start <- getBegEndIndMSP(x@msp)[[1]] - 1 ## which(testMSP@msp[,1] == "NAME: "), indices of 'NAME: '
        end <- getBegEndIndMSP(x@msp)[[2]] + 1
        start <- start[i]
        end <- end[i]
        inds <- lapply(1:length(start), function(j) start[j]:end[j])
        inds <- unlist(inds)

        return(MSP(msp=x@msp[inds,], mz = x@mz[i], rt = x@rt[i], names = x@names[i],
            classes = x@classes[i], information = x@information[i],
            adduct = x@adduct[i]))

        #return(new("MSP", msp = x@msp[inds, ]))
    }
)

#' @name getBegEndIndMSP
#' @title Get beginning and end indices of each entry in a \code{data.frame} in
#' \code{peaks(MSP)}-objects
#' @description Get beginning and end indices of each entry in a
#' \code{data.frame} in a \code{peaks(MSP)}-object
#' @usage getBegEndIndMSP(msp)
#' @param msp \code{data.frame} in \code{peaks(MSP)}-object, see
#'  \code{?convert2MSP} for further information
#' @details Internal use to retrieve start and end row indices for fragments of
#' MS/MS features.
#' @return \code{getBegEndIndMSP} returns a list of length 2 where the first
#' entry contains the start indices and the second the end indices
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ",
#'                          splitIndMZ = 2, splitIndRT = 3)
#' finalMSPdf <- peaks(finalMSP)
#' getBegEndIndMSP(finalMSPdf)
#' @export
getBegEndIndMSP <- function(msp) {

    ## beginning
    indPeaks <- which(msp[,1] == "Num Peaks: ")
    indLosses <- which(msp[,1] == "Num Losses: ")

    if (length(indPeaks) > length(indLosses)) indNumPeaks <- indPeaks
    ## <=: to get all possibilities
    if (length(indPeaks) <= length(indLosses)) indNumPeaks <- indLosses

    indEnd <- as.numeric(as.character(msp[indNumPeaks, 2]))
    indBeg <- indNumPeaks + 1
    indEnd <- indEnd + indBeg - 1

    return(list(indBeg, indEnd))
}
