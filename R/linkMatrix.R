#' @name createLinkMatrix
#' @title Create a matrix which contains features to link (indices)
#' @description Create a matrix which contains features to link (indices)
#' @usage createLinkMatrix(similarityMatrix, threshold, dfNameGroup)
#' @param similarityMatrix matrix, a similarity matrix that contains the 
#' NDP similarity measure between all precursors in the data set
#' @param threshold numerical, threshold value for NDP values, below this value 
#' linked features will not be included
#' @param dfNameGroup data.frame, data.frame contains column "group" and "name"
#' @details
#' @value createLinkMatrix returns a matrix that gives per each row 
#' information on linked features
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{createLinkMatrix(similarityMatrix, threshold, dfNameGroup)}
#' @export
createLinkMatrix <- function(similarityMatrix, threshold, dfNameGroup) {
    
    if (!all(colnames(similarityMatrix) == rownames(similarityMatrix))) {
        print("colnames(similarityMatrix) != rownames(similarityMatrix)")
        break  
    } 
    
    links <- which(similarityMatrix >= threshold, arr.ind = TRUE)
    links <- links[which(links[,1] != links[,2]),] 
    ## filter all feature which are in the diagonal
    ## create matrix which contains linking information in the format
    ## group1 position1 group2 position2 NDP
    linkMat <- matrix(data = NA, ncol = 5, nrow = dim(links)[1] )
    colnames(linkMat) <- c("group1", "name1", "group2", "name2", "NDP") 
    colnames(similarityMatrix) <- as.character(dfNameGroup[,"name"]) # new 
    rownames(similarityMatrix) <- as.character(dfNameGroup[,"name"]) # new
    
    for (i in 1:dim(links)[1]) {
        feature <- links[i,]
        rowFeature <- rownames(similarityMatrix)[feature][1]
        colFeature <- colnames(similarityMatrix)[feature][2]
        
        dfRowFeature <- dfNameGroup[which(rowFeature == dfNameGroup[, "name"]),]
        dfColFeature <- dfNameGroup[which(colFeature == dfNameGroup[, "name"]),]
        iEntry <- c(as.character(dfRowFeature[["group"]]),
                         as.character(dfRowFeature[["name"]]), 
                         as.character(dfColFeature[["group"]]),
                         as.character(dfColFeature[["name"]]), 
                         similarityMatrix[feature["row"], feature["col"]])
        linkMat[i, ] <- iEntry
    }
    return(linkMat)
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
#' @value cutLinkMatrix returns a matrix that gives per each row 
#' information on linked features
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{cutLinkMatrix(LinkMatrix, type = c("all", "inter", "intra"))}
#' @export
cutLinkMatrix <- function(LinkMatrix, type = c("all", "inter", "intra")) {
    
    type <- match.arg(type)
    
    if (type == "all") 
        LinkMatrix <- LinkMatrix
    if (type == "inter")
        LinkMatrix <- LinkMatrix[which(LinkMatrix[,"group1"] != LinkMatrix[,"group2"]), ]
    if (type == "intra") 
        LinkMatrix <- LinkMatrix[which(LinkMatrix[,"group1"] == LinkMatrix[,"group2"]), ]
    return(LinkMatrix)
}