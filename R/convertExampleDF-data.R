#' @name convertExampleDF
#' @title Example data for \code{MetCirc}: convertExampleDF 
#' @description \code{convertExampleDF} is a \code{data.frame} which comprises 
#' information on a specific metabolite per row stating the average retention 
#' time, average m/z, the name of the metabolite, the adduct ion name and the 
#' spectrum reference file name. The function \code{allocatePrecursor2mz} uses 
#' \code{data.frame}s of the kind of \code{sd01\_outputXCMS} and 
#' \code{sd02\_deconvoluted} to create a \code{data.frame} of the kind of 
#' \code{convertExampleDF}. Allocation of precursor ions to candidate
#' m/z values is based on minimal distance of m/z and deviance of retention time 
#' based on an objective function. See \code{?allocatePrecursor2mz} for further
#' information.
#' @docType data
#' @usage convertExampleDF
#' @return data.frame
#' @format data.frame
#' @source internal 
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
NULL 