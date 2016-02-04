#' @name NDP
#' @title Calculate the normalised dot product
#' @description Calculate the normalised dot product (NDP)
#' @usage NDP(matrow1, matrow2, m = 0.5, n = 2, mass)
#' @param matrow1 character vector or numerical vector, the entries correspond 
#' to the mass vector and contain corresponding intensities to the masses, 
#' it is the first feature to compare
#' @param matrow2 character vector or numerical vector, the entries correspond 
#' to the mass vector and contain corresponding intensities to the masses, 
#' it is the second feature to compare
#' @param m numeric, exponent to calculate peak intensity-based weights
#' @param n numeric, exponent to calculate peak intensity-based weights
#' @param mass character vector or numerical vector, vector with all masses 
#' which occur in the data set
#' @details The NDP is calculated according to the following formula: 
#'  \deqn{NDP = \frac{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2}{ \sum(W_{S1, i} ^ 2) * \sum(W_{S2, i} ^ 2) }}{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2 \sum(W_{S1, i} ^ 2) * \sum(W_{S2, i} ^ 2)},
#'  with \eqn{W = [ peak intensity] ^{m} \cdot [m/z]^n}. For further information 
#'  see Li et al. (2015): Navigating natural variation in herbivory-induced
#'  secondary metabolism in coyote tobacco populations using MS/MS structural analysis. 
#'  PNAS, E4147--E4155. NDP returns a numeric value ranging between 0 and 1, where 0 
#' indicates no similarity between the two precursors, while 1 indicates 
#' a strong similarity between the two precursors.
#' @return NDP returns a numeric similarity coefficient between 0 and 1
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' load(system.file("data/sd02_deconvoluted.RData", 
#'      package = "MetabolomicTools")) 
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' binnedMSP <- binning(msp = finalMSP, tol = 0.01)
#' NDP(matrow1 = binnedMSP[1,], matrow2 = binnedMSP[2,], m = 0.5, n = 2,
#'  mass = colnames(binnedMSP))
#' @export
NDP <- function(matrow1, matrow2, m = 0.5, n = 2, mass) {

    S1 <- as.numeric(matrow1)
    S2 <- as.numeric(matrow2)
    mass <- as.numeric(mass)

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
#' @usage createSimilarityMatrix(mm)
#' @param mm matrix, colnames are all fragments which occur in the dataset, 
#'      rownames are m/z / rt values, entries of mm are intensity values 
#'      corresponding to the mass
#' @details createSimilarityMatrix calls a function to calculate the 
#' NDP between all precursors in the data set. For further
#' information on how the NDP is calculated see ?NDP and Li et al. (2015): 
#' Navigating natural variation in herbivory-induced secondary metabolism in 
#' coyote tobacco populations using MS/MS structural analysis. PNAS, 
#' E4147--E4155.
#' @return createSimilarityMatrix returns a similarity matrix that contains the 
#' NDP similarity measure between all precursors in the data set
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' load(system.file("data/sd02_deconvoluted.RData", 
#'      package = "MetabolomicTools")) 
#' finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)
#' binnedMSP <- binning(msp = finalMSP, tol = 0.01)
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' @export
createSimilarityMatrix <- function(mm) {
    n <- dim(mm)[1]
    colNames <- colnames(mm)
    similarity <- matrix(0, nrow = n, ncol = n)
    colnames(similarity) <- rownames(similarity) <- rownames(mm)

    ## write to similarity matrix similarity measure
    for (i in 1:n) {
        for (j in 1:n) {
            if (i <= j) {
                similarity[j,i] <- similarity[i,j] <- NDP(matrow1 = mm[i,], 
                                    matrow2 = mm[j,], m = 0.5, n = 2, 
                                    mass = colNames)
            }
        }
    }
    return(similarity)
    
}
    
