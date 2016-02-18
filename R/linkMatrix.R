
#' @name createLink0Matrix
#' @title Create a link matrix 
#' @description Create a link matrix which links every feature in similarity
#' matrix with another. 
#' @usage createLink0Matrix(similarityMatrix, dfNameGroup, removeDuplicates)
#' @param similarityMatrix matrix, a similarity matrix that contains the 
#' NDP similarity measure between all precursors in the data set
#' @param dfNameGroup data.frame, data.frame contains column "group" and "name"
#' @param removeDuplicates logical, whether to remove duplicate entries or not
#' @details createLink0Matrix creates a matrix from a similarityMatrix which 
#' includes all connections between features in the similarityMatrix, but 
#' exclude links which have a similarity of exactly 0.
#' @return createLink0Matrix returns a matrix that gives per each row 
#' information on linked features
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' load(system.file("data/sd02_deconvoluted.RData", 
#'      package = "MetabolomicTools")) 
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' binnedMSP <- binning(msp = finalMSP, tol = 0.01)
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' namesPrec <- rownames(binnedMSP)
#' compartment <- sample(c("yl", "ol", "s","r"), size = length(namesPrec), 
#'      replace=TRUE) 
#' namesPrec <- paste(compartment, namesPrec, sep="_")
#' dfNameGroup <- data.frame(group = compartment, name = namesPrec) 
#' dfNameGroup <- dfNameGroup[order(dfNameGroup[,"name"]),] 
#' createLink0Matrix(similarityMatrix = similarityMat,  
#'      dfNameGroup = dfNameGroup)
#' @export
createLink0Matrix <- function(similarityMatrix, dfNameGroup) { ## }, removeDuplicates = TRUE) {
    
    name <- rownames(similarityMatrix)
    if (!all(colnames(similarityMatrix) == name)) {
        stop("colnames(similarityMatrix) != rownames(similarityMatrix)")
    } 
    
    if (!all(colnames(dfNameGroup) == c("group", "name"))) {
        stop("colnames of argument dfNameGroup are not group and name")
    }
    

    dfName <- dfNameGroup[,2]
    
    if (!all(sort(colnames(similarityMatrix)) == sort(dfName))) {
        stop("colnames/rownames of similarityMatrix are not identical to names in dfNameGroup")
    }
    
#     dfNameTruncate <- strsplit(dfNameGroup[,2], split = "_")
#     if (all(unlist(lapply(dfNameTruncate, function(x){ length(x) == 2}))))
#         dfNameTruncate <- unlist(lapply(dfNameTruncate, "[", 2))
#     else
#         dfNameTruncate <- unlist(dfNameTruncate)
#     
#     if (!all(sort(dfNameTruncate) == sort(name))) {
#         stop("dfNameGroup[,'name'] != colnames(similarityMatrix) ")
#     }
#     
#     dfName <- dfNameTruncate
    mat <- matrix(data = NA, ncol = 5, 
                    nrow = (length(name)^ 2 - length(name)) / 2)
    colnames(mat) <- c("group1", "name1", "group2", "name2", "NDP") 
    
    for (i in 1:length(name)) { ## columns 
        for (j in 1:length(name)) { ## rows
            if (i < j) { ## do not include which link to the same
                if (similarityMatrix[j,i] > 0) { ## do not include which have NDP of 0
                    rowIndex <- min(which(is.na(mat[,1])))
                    ## group1
                    group1 <- dfNameGroup[which(dfNameGroup[,2] == name[i]),1]
                    group1 <- as.character(group1)
                    mat[rowIndex, "group1"] <- group1
                    ## name1
                    mat[rowIndex, "name1"] <- name[i] 
                    ## group2
                    group2 <- dfNameGroup[which(dfNameGroup[,2] == name[j]),1]
                    group2 <- as.character(group2)
                    mat[rowIndex, "group2"] <- group2
                    ## name2
                    mat[rowIndex, "name2"] <- name[j]
                    ## NDP
                    mat[rowIndex, "NDP"] <- similarityMatrix[j,i]
                }
            }
        }
    }
    lastIndex <- min(which(is.na(mat[,1]))) - 1
    mat <- mat[1:lastIndex, ]
    
#     for (i in 1:length(name)) {
#         nameI <- name[i]
#         compartmentI <- dfNameGroup[which(nameI == dfName),1]
#         compartmentI <- as.character(compartmentI)
#         entryI <- matrix(data = c(compartmentI,nameI, NA, NA, NA), ncol = 5, nrow = length(name) - 1, byrow=TRUE)
#         
#         newName <- name[-which(as.character(nameI) == name)]
#         ## group of second feature
#         indComp <- match(newName, dfName)
#         entryI[,3] <- as.character(dfNameGroup[indComp,1])
#         ## m/z / retention time of second feature
#         entryI[,4] <- newName
#         ## similarity
#         entryI[,5] <- as.numeric(similarityMat[nameI,newName])
#         ## write to mat
#         mat[((length(name)-1) * (i - 1) + 1):((length(name)-1)*i),] <- entryI
#     }
#     ###########################
#     if (removeDuplicates) {
#         ## create vector which connects mz and rt for connected features
#         nameCompound <- paste(mat[,2], mat[,4], sep = "_")
#         ## create vector which connects mz and rt for connected features (reverse)
#         nameCompoundrev <- paste(mat[,4], mat[,2], sep = "_")
#         ## create list where indices refer to indices of nameCompound and 
#         ## entries link to indices of nameCompoundrev
#         double <- lapply(nameCompound, function(nameI) match(nameI, nameCompoundrev))
#         ## set double entries NA
#         
#         for (i in 1:length(double)) {
#             if (!is.na(double[[i]])) {
#                 if (!is.na(double[[double[[i]]]]) & double[[double[[i]]]] == i) 
#                     double[[i]] <- NA                
#             }
#         }
#         ## unlist double and create mat with no double values
#         double <- unlist(double)
#         mat <- mat[double[!is.na(double)],]
#     }
    return(mat)
}

#' @name thresholdLinkMatrix
#' @title Threshold a link matrix
#' @description Threshold a link matrix 
#' @usage thresholdLinkMatrix(linkMatrix, threshold)
#' @param linkMatrix matrix, a link matrix that gives per each row 
#' information on linked features
#' @param threshold numerical, threshold value for NDP values, below this value 
#' linked features will not be returned
#' @details threshold is a numerical value and filters linked precursor ions; 
#' filtering is currently based on the normalised dot product.
#' @return thresholdLinkMatrix returns a matrix that gives per each row 
#' information on linked features which are linked above a certain threshold
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' load(system.file("data/sd02_deconvoluted.RData", 
#'      package = "MetabolomicTools")) 
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' binnedMSP <- binning(msp = finalMSP, tol = 0.01)
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' namesPrec <- rownames(binnedMSP)
#' compartment <- sample(c("yl", "ol", "s","r"), size = length(namesPrec), 
#'      replace=TRUE) 
#' namesPrec <- paste(compartment, namesPrec, sep="_")
#' dfNameGroup <- data.frame(group = compartment, name = namesPrec) 
#' dfNameGroup <- dfNameGroup[order(dfNameGroup[,"name"]),] 
#' linkMatrix <- createLink0Matrix(similarityMatrix = similarityMat,  
#'      dfNameGroup = dfNameGroup)
#' thresholdLinkMatrix(linkMatrix = linkMatrix, threshold = 0.5)
#' @export
thresholdLinkMatrix <- function(linkMatrix, threshold) {
    
    if (!all(colnames(linkMatrix) == c("group1", "name1",  "group2", "name2",  "NDP")))
        stop("linkMatrix does not have right colnames")
    
    ndp <- as.numeric(linkMatrix[, "NDP"])
    
    if (threshold > max(ndp)) 
        stop("threshold greater than max NDP value in linkMatrix")
    
    ## which rows have a coefficient >= threshold?
    indThreshold <- which(ndp >= threshold)
    
    ## cut linkMatrix
    if (length(indThreshold) <= 1) {
        thresholdLinkMatrix <- matrix(NA, ncol = ncol(linkMatrix), nrow = length(indThreshold))
        thresholdLinkMatrix[1:nrow(thresholdLinkMatrix),1:ncol(thresholdLinkMatrix)] <- linkMatrix[indThreshold,]
        colnames(thresholdLinkMatrix) <- colnames(linkMatrix)    
    } else {
        thresholdLinkMatrix <- linkMatrix[indThreshold,]
        colnames(thresholdLinkMatrix) <- colnames(linkMatrix) 
    }
    
    
    return(thresholdLinkMatrix)
    
}

#' @name createLinkMatrix
#' @title Create a matrix which contains features to link (indices)
#' @description Create a matrix which contains features to link (indices)
#' @usage createLinkMatrix(similarityMatrix, threshold, dfNameGroup)
#' @param similarityMatrix matrix, a similarity matrix that contains the 
#' NDP similarity measure between all precursors in the data set
#' @param dfNameGroup data.frame, data.frame contains column "group" and "name"
#' @param threshold numerical, threshold value for NDP values, below this value 
#' linked features will not be included
#' @details threshold is a numerical value and filters linked precursor ions; 
#' filtering is currently based on the normalised dot product.
#' @return createLinkMatrix returns a matrix that gives per each row 
#' information on linked features
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' load(system.file("data/sd02_deconvoluted.RData", 
#'      package = "MetabolomicTools")) 
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' binnedMSP <- binning(msp = finalMSP, tol = 0.01)
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' namesPrec <- rownames(binnedMSP)
#' compartment <- sample(c("yl", "ol", "s","r"), size = length(namesPrec), 
#'      replace=TRUE) 
#' namesPrec <- paste(compartment, namesPrec, sep="_")
#' dfNameGroup <- data.frame(group = compartment, name = namesPrec) 
#' dfNameGroup <- dfNameGroup[order(dfNameGroup[,"name"]),] 
#' createLinkMatrix(similarityMatrix = similarityMat, threshold = 0.5, 
#'      dfNameGroup = dfNameGroup)
#' @export
createLinkMatrix <- function(similarityMatrix, dfNameGroup, threshold) {
    ## first create a link0Matrix
    linkMatrix <- createLink0Matrix(similarityMatrix = similarityMatrix, 
                        dfNameGroup = dfNameGroup)
    ## than threshold link0Matrix
    thresholdLinkMatrix <- thresholdLinkMatrix(linkMatrix = linkMatrix, 
                        threshold = threshold)
     return(thresholdLinkMatrix)
}
## deprecated
## createLinkMatrix <- function(similarityMatrix, threshold, dfNameGroup) {
##     
##     if (!is.numeric(threshold)) stop("threshold is not numeric")
##     
##     if (!all(colnames(similarityMatrix) == rownames(similarityMatrix))) {
##         stop("colnames(similarityMatrix) != rownames(similarityMatrix)")
##         
##     } 
##     
##     ##links <- which(similarityMatrix >= threshold, arr.ind = TRUE)
##     links <- links[which(links[,1] != links[,2]),] 
##     ## filter all feature which are in the diagonal
##     ## create matrix which contains linking information in the format
##     ## group1 position1 group2 position2 NDP
##     linkMat <- matrix(data = NA, ncol = 5, nrow = dim(links)[1] )
##     colnames(linkMat) <- c("group1", "name1", "group2", "name2", "NDP") 
##     colnames(similarityMatrix) <- as.character(dfNameGroup[,"name"]) # new 
##     rownames(similarityMatrix) <- as.character(dfNameGroup[,"name"]) # new
##     
##     for (i in 1:dim(links)[1]) {
##         feature <- links[i,]
##         rowFeature <- rownames(similarityMatrix)[feature][1]
##         colFeature <- colnames(similarityMatrix)[feature][2]
##         
##         dfRowFeature <- dfNameGroup[which(rowFeature == dfNameGroup[, "name"]),]
##         dfColFeature <- dfNameGroup[which(colFeature == dfNameGroup[, "name"]),]
##         iEntry <- c(as.character(dfRowFeature[["group"]]),
##                     as.character(dfRowFeature[["name"]]), 
##                     as.character(dfColFeature[["group"]]),
##                     as.character(dfColFeature[["name"]]), 
##                     similarityMatrix[feature["row"], feature["col"]])
##         linkMat[i, ] <- iEntry
##     }
##     return(linkMat)
## }


#' @name cutLinkMatrix
#' @title Create a cut LinkMatrix 
#' @description Create a cut LinkMatrix 
#' @usage cutLinkMatrix(LinkMatrix, type = c("all", "inter", "intra"))
#' @param LinkMatrix matrix, that gives per each row 
#' information on linked features
#' @param type character, one of "all", "inter" or "intra"
#' @details This function is used to cut features from LinkMatrix. If 
#' type = "all", LinkMatrix will not be changed; if type = "inter" the cut
#' LinkMatrix will only contain entries of links which are between groups and 
#' not inside groups; contrary to that, if type = "intra" the cut LinkMatrix 
#' will only contain entries of links which are inside groups and not between 
#' groups.
#' @return cutLinkMatrix returns a matrix that gives per each row 
#' information on linked features
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' load(system.file("data/sd02_deconvoluted.RData", 
#'      package = "MetabolomicTools")) 
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' binnedMSP <- binning(msp = finalMSP, tol = 0.01)
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' namesPrec <- rownames(binnedMSP)
#' compartment <- sample(c("yl", "ol", "s","r"), size = length(namesPrec), 
#'      replace=TRUE) 
#' namesPrec <- paste(compartment, namesPrec, sep="_")
#' dfNameGroup <- data.frame(group = compartment, name = namesPrec)
#' dfNameGroup <- dfNameGroup[order(dfNameGroup[,"name"]),] 
#' linkMat <- createLinkMatrix(similarityMatrix = similarityMat, threshold = 0.5, 
#'      dfNameGroup = dfNameGroup)
#' cutLinkMatrix(LinkMatrix = linkMat, type = "all")
#' @export
cutLinkMatrix <- function(LinkMatrix, type = c("all", "inter", "intra")) {
    lM <- LinkMatrix
    type <- match.arg(type)
    
    if (type == "all") 
        lM <- lM
    if (type == "inter")
        lM <- lM[which(lM[,"group1"] != lM[,"group2"]), ]
    if (type == "intra") 
        lM <- lM[which(lM[,"group1"] == lM[,"group2"]), ]
    return(lM)
}