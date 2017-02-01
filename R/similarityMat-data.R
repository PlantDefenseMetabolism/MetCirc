#' @name similarityMat
#' @title Example data for \code{MetCirc}: \code{similarityMat}
#' @description \code{similarityMat} is a matrix containing the pair-wise
#' similarity scores derived from the \code{idMSMStissueproject} data set. 
#' See the vignette for a workflow to reproduce the object \code{similarityMat}.
#' @docType data
#' @usage similarityMat
#' @return matrix
#' @format matrix
#' @source 
#' data("binnedMSP", package = "MetCirc")
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' save(similarityMat, file = "similarityMat.RData", compress = "xz")
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com} 
NULL  
