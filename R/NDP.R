#' @name NDP
#' @title Calculate the normalised dot product
#' @description Calculate the normalised dot product (NDP)
#' @usage NDP(matrow1, matrow2, m = 0.5, n = 2, mass)
#' @param matrow1 \code{character} or \code{numeric} vector, the entries 
#' correspond to the mass vector and contain corresponding intensities to the 
#' masses, it is the first feature to compare
#' @param matrow2 \code{character} or \code{numeric} vector, the entries 
#' correspond to the mass vector and contain corresponding intensities to the 
#' masses, it is the second feature to compare
#' @param m \code{numeric}, exponent to calculate peak intensity-based weights
#' @param n \code{numeric}, exponent to calculate peak intensity-based weights
#' @param mass \code{character} or \code{numeric} vector, vector with all masses 
#' which occur in the data set
#' @details The NDP is calculated according to the following formula: 
#'  \deqn{NDP = \frac{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2}{ \sum(W_{S1, i} ^ 2) * \sum(W_{S2, i} ^ 2) }}{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2 \sum(W_{S1, i} ^ 2) * \sum(W_{S2, i} ^ 2)},
#'  with \eqn{W = [ peak intensity] ^{m} \cdot [m/z]^n}. For further information 
#'  see Li et al. (2015): Navigating natural variation in herbivory-induced
#'  secondary metabolism in coyote tobacco populations using MS/MS structural analysis. 
#'  PNAS, E4147--E4155. NDP returns a numeric value ranging between 0 and 1, where 0 
#' indicates no similarity between the two MS/MS features, while 1 indicates 
#' that the MS/MS features are identical.
#' @return NDP returns a numeric similarity coefficient between 0 and 1
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("binnedMSP", package = "MetCirc")
#' NDP(matrow1 = binnedMSP[1,], matrow2 = binnedMSP[2,], m = 0.5, n = 2,
#'  mass = colnames(binnedMSP))
#' @export
NDP <- function(matrow1, matrow2, m = 0.5, n = 2, mass) {

    S1 <- as.numeric(matrow1)
    S2 <- as.numeric(matrow2)
    mass <- as.numeric(mass)
    
    if (length(S1) != length(S2)) 
        stop("matrow1 and matrow2 have not identical length")
    if (length(mass) != length(S1)) 
        stop("mass has not same length as matrow1 and matrow2")

    ## calculate weights 
    WS1 <- S1 ^ m * mass ^ n
    WS2 <- S2 ^ m * mass ^ n
    ## calculate NDP
    NDP <- ( sum(WS1 * WS2) ) ^ 2 / (sum(WS1 ^ 2 ) * sum(WS2 ^ 2))
    return(NDP)
}

#' @name createSimilarityMatrix
#' @title Create similarity matrix
#' @description Creates the similarity matrix by calculating the normalised dot 
#' product (NDP) between precursors
#' @usage createSimilarityMatrix(mm, m = 0.5, n = 2)
#' @param mm \code{matrix}, colnames are all fragments which occur in the 
#'  dataset, rownames are m/z / rt values, entries of \code{mm} are intensity 
#'  values corresponding to their m/z values
#' @param m \code{numeric}, see \code{?NDP} for further details
#' @param n \code{numeric}, see \code{?NDP} for further details
#' @details createSimilarityMatrix calls a function to calculate the 
#' NDP between all precursors in the data set. For further
#' information on how the NDP is calculated see \code{?NDP} and Li et al. (2015): 
#' Navigating natural variation in herbivory-induced secondary metabolism in 
#' coyote tobacco populations using MS/MS structural analysis. PNAS, 
#' E4147--E4155. Currently \code{m = 0.5} and \code{n = 2} are set as 
#' default.
#' @return \code{createSimilarityMatrix} returns a similarity matrix that 
#' contains the NDP similarity measure between all precursors in the data set
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("binnedMSP", package = "MetCirc")
#' ## truncate binnedMSP 
#' binnedMSP <- binnedMSP[1:28,]
#' createSimilarityMatrix(binnedMSP, m = 0.5, n = 2)
#' @export
createSimilarityMatrix <- function(mm, m = 0.5, n = 2) {
    
    if (!is.numeric(m)) stop("m is not numeric")
    if (!is.numeric(n)) stop("n is not numeric")
    dimMM <- dim(mm)[1]
    colNames <- colnames(mm)
    similarity <- matrix(0, nrow = dimMM, ncol = dimMM)
    groupname <- rownames(mm)
    ##orderNew <- order(groupname)
    ##mm <- mm[orderNew,]
    ##rownames(mm) <- groupname <- groupname[orderNew]
    colnames(similarity) <- rownames(similarity) <- groupname
    
    colNames <- as.numeric(colNames)
    colNames <- abs(colNames)
    ## write to similarity matrix similarity measure
    for (i in 1:dimMM) {
        for (j in 1:dimMM) {
            if (i <= j) {
                similarity[j,i] <- similarity[i,j] <- NDP(matrow1 = mm[i,], 
                                    matrow2 = mm[j,], m = m, n = n, 
                                    mass = colNames)
            }
        }
    }
    
    return(similarity)
    
}

#' @import amap
#' @name createOrderedSimMat
#' @title Update colnames and rownames of a similarity matrix according to 
#' order m/z, retention time and clustering
#' @description Internal function for shiny application. May also be used 
#' outside of shiny to reconstruct figures.
#' @usage createOrderedSimMat(similarityMatrix, order = c("retentionTime", "mz", "clustering"))
#' @param similarityMatrix \code{matrix}, \code{similarityMatrix} contains 
#' pair-wise similarity coefficients which give information about the similarity 
#' between precursors
#' @param order \code{character}, one of "retentionTime", "mz" or "clustering"
#' @details \code{createOrderSimMat} takes  a similarity matrix and a 
#' \code{character} vector
#' as arguments. It will then reorder rows and columns of 
#' the similarityMatrix object such, that it orders rows and columns of 
#' similarityMatrix according to m/z, retention time or clustering in 
#' each group. \code{createOrderSimMat} is employed in the shinyCircos 
#' function to create \code{similarityMatrix} objects which will allow to switch
#' between different types of ordering in between groups (sectors) in the 
#' circos plot. It may be used as well externally, to reproduce plots outside
#' of the reactive environment (see vignette for a workflow).
#' @return \code{createOrderedSimMat} returns a similarity matrix with ordered
#' rownames according to the \code{character} vector given to order
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("binnedMSP", package = "MetCirc")
#' data("similarityMat", package = "MetCirc")
#' ## order according to retention time 
#' createOrderedSimMat(similarityMatrix = similarityMat, order = "retentionTime")
#' @export
createOrderedSimMat <- function(similarityMatrix, order = c("retentionTime","mz", "clustering")) {
    
    order <- match.arg(order)
    groupname <- rownames(similarityMatrix)

    ## get group and name from groupname
    ## groupname is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    group <- lapply(strsplit(groupname, split = "_"), "[", 1)
    group <- unlist(group)
    name <- lapply(strsplit(groupname, split = "_"), function (x) x[length(x)])
    name <- unlist(name)
    
    ## retentionTime
    if (order == "retentionTime") {
        nameMZRT <- strsplit(name, split = "/")
        rt <- lapply(nameMZRT, "[[", 2)
        rt <- unlist(rt)
        rt <- as.numeric(rt)
        orderNew <- order(group, rt)
    }
    
    ## mz
    if (order == "mz") {
        nameMZRT <- strsplit(name, split = "/")
        mz <- lapply(nameMZRT, "[[", 1)
        mz <- unlist(mz)
        mz <- as.numeric(mz)
        orderNew <- order(group, mz)
    }
    
    ## clustering
    if (order == "clustering") {
        orderNew <- numeric(length = length(groupname))
        groupLevels <- sort(unique(group))
        ## loop in levels
        for (i in groupLevels) {
            inds <- which(group == i)
            nameGroupLevel <- groupname[inds]
            simMatI <- similarityMatrix[inds, inds]
            hClust <- amap::hcluster(simMatI, method = "spearman") 
            ## write order within groups to orderNew
            orderNew[inds] <- inds[hClust$order]
        }
    }
    
    groupNameNew <- groupname[orderNew]
    groupNameNewSplit <- strsplit(groupNameNew, "_")
    
    ## count from 1 to length(x) of unique groups
    groupNew <- unlist(lapply(groupNameNewSplit, "[", 1))
    ## length of each group
    groupNew_l <- as.vector(table(groupNew))
    ## create counter
    counter <- lapply(1:length(groupNew_l), 
        function (x) sprintf("%04d", 1:groupNew_l[x]))
    counter <- unlist(counter)
    groupNameNew <- lapply(1:length(groupname), 
        function (x) paste(groupNameNewSplit[[x]][1], 
                        counter[x], groupNameNewSplit[[x]][2], sep="_"))
    groupNameNew <- unlist(groupNameNew)
    
   ## if (order == "clustering") {
   ##     ## reorder due to clustering
   ##     orderNew <- order(groupNameNew)
##        groupNameNew <- groupNameNew[orderNew]    
   ## }
    
    simM <- similarityMatrix[orderNew, orderNew]
    colnames(simM) <- rownames(simM) <- groupNameNew
    
    return(simM)
}

