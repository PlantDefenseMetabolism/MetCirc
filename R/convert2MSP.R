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
#' @usage convert2MSP(mm, splitPattern = "_", splitInd = 1, names = FALSE, 
#'  metNames = FALSE, class = FALSE)
#' @param mm matrix, mm has to have four columns with colnames 
#'  mz, rt, intensity (order is not important). In the fourth column there has 
#'  to information about the precursor ion which will be assessed by 
#'  splitPattern and splitInd. Optionally, mm can have colnames names, 
#'  metNames, class. 
#' @param splitPattern character, splitPattern is the pattern which separates 
#'      elements and precursor m/z
#' @param splitInd numeric, the position of the precursor m/z concerning 
#'      separation by splitPattern
#' @param names logical, should names be retrieved? If set to TRUE, convert2MSP
#'      will access the column "names" in mm which contains the names of the 
#'      metabolites
#' @param metNames logical, should names of metabolites be retrieved? 
#'      If set to TRUE, convert2MSP will access the column "metNames" in mm 
#'      which contains the names of the metabolites
#' @param class logical, should classes of metabolites be retrieved? If set to 
#'      TRUE, convert2MSP will access the column "class" in mm which contains
#'      the names of the metabolites
#' @details Creates a data entry for each precursor ion. Each entry in the 
#' return object has the following information: NAME, RETENTIONTIME, 
#'      PRECURSORMZ, METABOLITENAME, ADDUCTIONNAME, Num Peaks and a list of 
#'      fragments together with their intensities. convert2MSP will access
#'      the column name 'name', 'metNames' and 'class', respectively, 
#'      if arguments are set to TRUE. In the fourth column there has 
#'      to information about the precursor ion which will be assessed by 
#'      splitPattern and splitInd. E.g. items in the fourth column can be in 
#'      the form of '1_163.23', which has to be accessed by setting 
#'      \code{splitPattern = "_"} and \code{splitInd = 2} to access the m/z 
#'      value of the precursor ion (here: 162.23). 
#' @return convert2MSP returns an object of class MSP
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' convert2MSP(mm = sd02_deconvoluted, splitPattern = "_", splitInd = 1, 
#'  names = FALSE, metNames = FALSE, class = FALSE)
#' @export
convert2MSP <- function (mm, splitPattern = "_", splitInd = 1, names = FALSE, 
                    metNames = FALSE, class = FALSE) {
    
    colNames <- colnames(mm)
    if (!("mz" %in% colNames)) stop("no column 'mz' found")
    if (!("rt" %in% colNames)) stop("no column 'rt' found")
    if (!("intensity" %in% colNames)) stop("no column 'intensity' found")
    if (names & !("names" %in% colNames)) stop("no column 'names' found")
    if (metNames & !("metNames" %in% colNames)) 
                                        stop("no column 'metNames' found")
    if (class & !("class" %in% colNames)) stop ("no column 'class' found")
    
    ## if (colNames[4] != "pcgroup_precursorMZ") break
     
    precursor <- mm[,4]
    precursor <- as.character(precursor)
    
    uniquePreMZ <- unique(precursor)
    uniquePreMZ_cut <- cutUniquePreMZ(precursor = precursor, 
            splitPattern = splitPattern, splitInd = splitInd)
    
    ## check if pcgroup_grecursorMZ is in fourth column and unique precursor
    ## mz were assessed correctly
    if (any(is.na(as.numeric(uniquePreMZ_cut))))
            stop("unique precursors are not in the fourth column or 
                splitPattern/splitInd was chosen incorrectly")
    
    lenUniquePreMZ <- length(uniquePreMZ_cut)
    
    ## access columns names, metNames, class to retrieve information about 
    ## names, metabolite names or metabolite class
    if (names) namesMM <- mm[, "names"]
    if (metNames) metNamesMM <- mm[, "metNames"]
    if (class) classesMM <- mm[, "classes"]

    ## create data frame for MSP file
    finalMSP <- matrix(data = NA, nrow = 8 * lenUniquePreMZ + dim(mm)[1], 
            ncol = 2) ## 7 new entries + all fragment ion entries
    finalMSP <- as.data.frame(finalMSP)
    
    ## write to data frame
    for (i in 1:lenUniquePreMZ) {
        ind <- which(uniquePreMZ[i] == precursor)    
        entry <- rbind(
            c("NAME: ", if (names) {
                unique(as.character(namesMM[ind])[1])} else "Unknown"),
            c("RETENTIONTIME: ", mean(mm[ind,"rt"])),
            c("PRECURSORMZ: ", uniquePreMZ_cut[i]),
            c("METABOLITENAME: ", if (metNames) {
                unique(as.character(metNamesMM[ind])[1])} else "Unknown"),
            c("METABOLITECLASS: ", if (class) {
                unique(as.character(classesMM[ind])[1])} else "Unknown"),
            c("ADDUCTIONNAME: ", "Unknown"),
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
    
    ## do the following before getMSP
    nameMet <- getMetaboliteName(msp)
    name <- getName(msp)
    classes <- getMetaboliteClass(msp)
    precmz <- getPrecursorMZ(msp)
    rt <- getRT(msp)
    ## 
    
    msp <- getMSP(msp)
    
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
            c("NAME: ", name[i]),
            c("RETENTIONTIME: ", rt[i]),
            c("PRECURSORMZ: ", precmz[i]),
            c("METABOLITENAME: ", nameMet[i]),
            c("METABOLITECLASS: ", classes[i]),
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

    
    return(new("MSP", msp = finalMSP))
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
          definition = function(x) {length(getPrecursorMZ(x))
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
                  length(getPrecursorMZ(object)), "entries.", sep = " ")
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

#' @name getName
#' @aliases getName,MSP-method
#' @title getName returns names in MSP object
#' @return character
#' @description getName returns names in MSP object.
#' @param object object of class MSP, see ?convert2MSP for further information
#' @docType methods
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' getName(finalMSP)
#' @export
setGeneric("getName", function(object) standardGeneric("getName"))

#' @describeIn getName returns names in MSP objects
#' @export
getName <- function(object) {
    ## get names in an msp object
    df <- object@msp
    ind <- which(df[,1] == "NAME: ")
    return(df[ind,2])
}

#' @name setName
#' @aliases setName,character,MSP-method
#' @title setName sets names in MSP objects
#' @return MSP
#' @description setName sets names in MSP objects. To set names pass a vector 
#' with names to the argument \code{class}.
#' @param object object of class MSP
#' @param name character, a vector with new names
#' @docType methods
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc") 
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = "_ ", splitInd = 2)
#' setMetaboliteName(finalMSP, c(rep("unknown", 358), "name1", "name2"))
#' @export
setGeneric("setName", function(object) standardGeneric("setName"))

#' @describeIn setName sets names in MSP objects
#' @export
setName <- function(object, name) {
    df <- object@msp
    ind <- which(df[,1] == "NAME: ") 
    if (length(name) != length(ind)) stop("number of items to replace does not 
                                           match with replacement length")
    df[ind, 2] <- name
    return(new("MSP", msp = df))
}

#' @name getMetaboliteName
#' @aliases getMetaboliteName,MSP-method
#' @title getMetaboliteName returns names of metabolites in MSP object
#' @return character
#' @description getMetaboliteName returns names in MSP object.
#' @param object object of class MSP, see ?convert2MSP for further information
#' @docType methods
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' getMetaboliteName(finalMSP)
#' @export
setGeneric("getMetaboliteName", 
           function(object) standardGeneric("getMetaboliteName"))

#' @describeIn getMetaboliteName returns names of metabolites in MSP objects
#' @export
getMetaboliteName <- function(object) {
    ## get names of metabolites in an msp object
    df <- object@msp
    ind <- which(df[,1] == "METABOLITENAME: ")
    return(df[ind,2])
}

#' @name setMetaboliteName
#' @aliases setMetaboliteName,character,MSP-method
#' @title setMetaboliteName sets metabolite names in MSP objects
#' @return MSP
#' @description setMetaboliteName sets metabolite names in MSP objects. To set 
#' metabolite names pass a vector with names to the argument \code{class}.
#' @param object object of class MSP
#' @param metName character, a vector with new metabolite names
#' @docType methods
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc") 
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = "_ ", splitInd = 2)
#' setMetaboliteName(finalMSP, c(rep("unknown", 358), "met1", "met2"))
#' @export
setGeneric("setMetaboliteName", 
           function(object) standardGeneric("setMetaboliteName"))

#' @describeIn setMetaboliteName sets metabolite names in MSP objects
#' @export
setMetaboliteName <- function(object, metName) {
    df <- object@msp
    ind <- which(df[,1] == "METABOLITENAME: ") 
    if (length(metName) != length(ind)) stop("number of items to replace does 
                                    not match with replacement length")
    df[ind, 2] <- metName
    return(new("MSP", msp = df))
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

#' @name setMetaboliteClass
#' @aliases setMetaboliteClass,character,MSP-method
#' @title setMetaboliteClass sets class names of compounds in MSP objects
#' @return MSP
#' @description setMetaboliteClass sets names of class names of compounds in 
#' MSP objects. To set names pass a vector with class names to the argument 
#' \code{class}.
#' @param object object of class MSP
#' @param class character, a vector with new class names
#' @docType methods
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc") 
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = "_ ", splitInd = 2)
#' setMetaboliteClass(finalMSP, c(rep("unknown", 359), "class1"))
#' @export
setGeneric("setMetaboliteClass", 
           function(object) standardGeneric("setMetaboliteClass"))

#' @describeIn setMetaboliteClass sets class names of compounds in MSP objects
#' @export
setMetaboliteClass <- function(object, class) {
    df <- object@msp
    ind <- which(df[,1] == "METABOLITECLASS: ") 
    if (length(class) != length(ind)) stop("number of items to replace does not 
                                           match with replacement length")
    df[ind, 2] <- class
    return(new("MSP", msp = df))
}

#' @name getRT
#' @aliases getRT,MSP-method
#' @title getRT returns precursor RT values of an MSP object
#' @return numeric
#' @description getRT returns a numeric vector with all retention time values
#' @param object object of class MSP
#' @docType methods
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' getRT(finalMSP)
#' @export
setGeneric("getRT", function(object) standardGeneric("getRT"))

#' @describeIn getRT returns precursor RT values of an MSP object
#' @export
getRT <- function(object) {
    ## get classes of compounds in an msp object
    df <- object@msp
    ind <- which(df[,1] == "RETENTIONTIME: ")
    ## get rt 
    rt <- df[ind,2]
    ## change to numeric
    rt <- as.numeric(rt)
    return(rt)
}

#' @name getPrecursorMZ
#' @aliases getPrecursorMZ,MSP-method
#' @title getPrecursorMZ returns precursor m/z values of an MSP object
#' @return numeric
#' @description getPrecursorMZ returns a numeric vector with 
#'  precursor m/z values
#' @param object object of class MSP
#' @docType methods
#' @examples 
#' data("sd02_deconvoluted", package = "MetCirc")
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' getPrecursorMZ(finalMSP)
#' @export
setGeneric("getPrecursorMZ", function(object) standardGeneric("getPrecursorMZ"))

#' @describeIn getPrecursorMZ returns precursor m/z values of an MSP object
#' @export
getPrecursorMZ <- function(object) {
    ## get classes of compounds in an msp object
    df <- object@msp
    ind <- which(df[,1] == "PRECURSORMZ: ")
    ## get m/z
    mz <- df[ind,2]
    ## change to numeric
    mz <- as.numeric(mz)
    return(mz)
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
