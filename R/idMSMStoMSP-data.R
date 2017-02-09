#' @name idMSMStoMSP-data
#' @aliases finalMSP
#' @title Example data for \code{MetCirc}: \code{finalMSP}
#' @description \code{finalMSP} is of instance 'MSP', a container for 
#' MS/MS data. \code{finalMSP} is derived from the object \code{tissue} 
#' and \code{compartmentTissue}.
#' @docType data
#' @usage finalMSP
#' @return object of class MSP
#' @format object of class MSP
#' @source 
#' data("idMSMStissueproject", package = "MetCirc")
#' ## create vectors with precursor names present in tissue
#' tissueSPL <- compartmentTissue[compartmentTissue[,"SPL"] == TRUE, 1]
#' tissueLIM <- compartmentTissue[compartmentTissue[,"LIM"] == TRUE, 1]
#' tissueANT <- compartmentTissue[compartmentTissue[,"ANT"] == TRUE, 1]
#' tissueSTY <- compartmentTissue[compartmentTissue[,"STY"] == TRUE, 1]
#' 
#' ## truncate tissue
#' tissueSPL <- tissue[tissue[,4] %in% tissueSPL,] 
#' tissueLIM <- tissue[tissue[,4] %in% tissueLIM,]
#' tissueANT <- tissue[tissue[,4] %in% tissueANT,]
#' tissueSTY <- tissue[tissue[,4] %in% tissueSTY,]
#' 
#' ## create msp and combine msp objects of different tissues
#' finalMSP <- convert2MSP(tissueSPL, rt = TRUE)
#' finalMSP <- combine(finalMSP, convert2MSP(tissueLIM), rt = TRUE)
#' finalMSP <- combine(finalMSP, convert2MSP(tissueANT), rt = TRUE)
#' finalMSP <- combine(finalMSP, convert2MSP(tissueSTY), rt = TRUE)
#' 
#' ## write finalMSP to idMSMStoMSP.RData
#' save(finalMSP, file = "idMSMStoMSP.RData", compress = "xz")
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
NULL  
##
##