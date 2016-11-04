#' @name createLink0Matrix
#' @title Create a link matrix 
#' @description Create a link matrix which links every feature in similarity
#' matrix with another. 
#' @usage createLink0Matrix(similarityMatrix)
#' @param similarityMatrix matrix, a similarity matrix that contains the 
#' NDP similarity measure between all precursors in the data set
#' @details createLink0Matrix creates a matrix from a similarityMatrix which 
#' includes all connections between features in the similarityMatrix, but 
#' exclude links which have a similarity of exactly 0.
#' @return createLink0Matrix returns a matrix that gives per each row 
#' information on linked features
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' data("binnedMSP", package = "MetCirc")
#' ## truncate binnedMSP
#' binnedMSP <- binnedMSP[1:28,]
#' namesPrec <- rownames(binnedMSP)
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' link0Mat <- createLink0Matrix(similarityMatrix = similarityMat)
#' @export
createLink0Matrix <- function(similarityMatrix) { 
    
    groupname <- rownames(similarityMatrix)
    if (!all(colnames(similarityMatrix) == groupname)) {
        stop("colnames(similarityMatrix) != rownames(similarityMatrix)")
    } 
    
    ## create vector with group names (e.g. compartments)
    group <- lapply(strsplit(groupname, split = "_"), "[", 1)
    group <- unlist(group)
    
    ## get matrix indices where similarity mat 
    inds <- which(similarityMatrix > 0, arr.ind = TRUE)
    indsrow <- as.vector(inds[,"row"])
    indscol <- as.vector(inds[,"col"])

    rowcol_s <- lapply(1:length(indsrow), function(x) sort(c(indsrow[x], indscol[x])))
    ##rowcol_s <- lapply(rowcol, sort)
    duplicatedRowCol <- duplicated(rowcol_s)
    inds <- inds[!duplicatedRowCol,]
    indsrow <- indsrow[!duplicatedRowCol]
    indscol <- indscol[!duplicatedRowCol]
    ## remove 1-1, 2-2, etc.
    pairwise <- which(indsrow == indscol)
    ##pairwiseCR <- which(indscol == indsrow)
    indsrow <- indsrow[-pairwise]
    indscol <- indscol[-pairwise]
    
    mat <- matrix(data = NA, ncol = 5, nrow = length(indsrow))
    colnames(mat) <- c("group1", "name1", "group2", "name2", "NDP") 

    mat[,"group1"] <- group[indsrow]
    mat[,"group2"] <- group[indscol]
    mat[,"name1"] <- groupname[indsrow]
    mat[,"name2"] <- groupname[indscol]
    ndps <- sapply(1:length(indsrow), function(x) similarityMatrix[indsrow[x], indscol[x]])
    mat[, "NDP"] <- ndps
    
    return(mat)
}

#' @name thresholdLinkMatrix
#' @title Threshold a link matrix
#' @description Threshold a link matrix 
#' @usage thresholdLinkMatrix(linkMatrix, threshold_low, threshold_high)
#' @param linkMatrix matrix, a link matrix that gives per each row 
#' information on linked features
#' @param threshold_low numerical, threshold value for NDP values, below this 
#' value linked features will not be returned
#' @param threshold_high numerical, threshold value for NDP values, above this 
#' value linked features will not be returned
#' @details threshold_low and threshold_high are numerical values and truncates 
#' similar/identical precursor ions; 
#' similarity is momentarily based on the normalised dot product.
#' @return thresholdLinkMatrix returns a matrix that gives per each row 
#' information on linked features which are linked above a certain threshold
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' data("binnedMSP", package = "MetCirc")
#' ## use only a selection 
#' binnedMSP <- binnedMSP[c(c(1:20, 29:48, 113:132, 240:259)),]
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' linkMatrix <- createLink0Matrix(similarityMatrix = similarityMat)
#' thresholdLinkMatrix(linkMatrix = linkMatrix, 
#'      threshold_low = 0.5, threshold_high=1)
#' @export
thresholdLinkMatrix <- function(linkMatrix, 
                                threshold_low = 0.75, threshold_high = 1) {
    
    if (!all(colnames(linkMatrix) == c("group1", "name1",  "group2", "name2",  "NDP")))
        stop("linkMatrix does not have right colnames")
    
    ndp <- as.numeric(linkMatrix[, "NDP"])
    if (threshold_low > threshold_high) stop("threshold_low greater than threshold_high")
    if (threshold_high > 1) stop("threshold_high greater than 1")
    if (threshold_low > max(ndp)) 
        warning("threshold greater than max NDP value in linkMatrix")
    
    ## which rows have a coefficient >= threshold?
    indThreshold <- which(ndp >= threshold_low & ndp <= threshold_high)
    
    ## cut linkMatrix
    if (length(indThreshold) <= 1) {
        thresholdLinkMatrix <- matrix(NA, ncol = ncol(linkMatrix), nrow = length(indThreshold))
        thresholdLinkMatrix[1:nrow(thresholdLinkMatrix),1:ncol(thresholdLinkMatrix)] <- linkMatrix[indThreshold,]
    } else {
        thresholdLinkMatrix <- linkMatrix[indThreshold,]
    }
    
    colnames(thresholdLinkMatrix) <- colnames(linkMatrix)  
    
    return(thresholdLinkMatrix)
}

#' @name createLinkMatrix
#' @title Create a matrix which contains features to link (indices)
#' @description Create a matrix which contains features to link (indices)
#' @usage createLinkMatrix(similarityMatrix, threshold_low, threshold_high)
#' @param similarityMatrix matrix, a similarity matrix that contains the 
#' NDP similarity measure between all precursors in the data set
#' @param threshold_low numerical, threshold value for NDP values, below this 
#' value linked features will not be included
#' @param threshold_high numerical, threshold value for NDP values, above this 
#' value linked features will not be included
#' @details threshold_low and threshold_high are numerical values and truncate 
#' similar/identical precursor ions; similarity is currently based on the 
#' normalised dot product.
#' @return createLinkMatrix returns a matrix that gives per each row 
#' information on linked features
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' data("binnedMSP", package = "MetCirc")
#' ## use only a selection 
#' binnedMSP <- binnedMSP[c(c(1:20, 29:48, 113:132, 240:259)),]
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' createLinkMatrix(similarityMatrix = similarityMat, 
#'      threshold_low = 0.5, threshold_high=1)
#' @export
createLinkMatrix <- function(similarityMatrix, threshold_low, threshold_high) {
    ## first create a link0Matrix
    linkMatrix <- createLink0Matrix(similarityMatrix = similarityMatrix)
    ## than threshold link0Matrix
    thresholdLinkMatrix <- thresholdLinkMatrix(linkMatrix = linkMatrix, 
                threshold_low = threshold_low, threshold_high = threshold_high)
     return(thresholdLinkMatrix)
}

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
#' data("binnedMSP", package = "MetCirc")
#' ## use only a selection 
#' binnedMSP <- binnedMSP[c(c(1:20, 29:48, 113:132, 240:259)),]
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' linkMat <- createLinkMatrix(similarityMatrix = similarityMat, threshold_low = 0.75, threshold_high = 1)
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