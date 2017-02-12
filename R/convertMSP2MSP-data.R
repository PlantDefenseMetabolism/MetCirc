#' @name msp2msp
#' @title Example data for \code{MetCirc}: \code{msp2msp}
#' @description \code{convertMSP2MSP} contains the object \code{msp2msp} 
#' that is a data frame in .MSP format, a typical format for MS/MS library 
#' building. Each entry consists of the metabolite name (NAME), the precursor 
#' mz (PRECURSORMZ), the retention 
#' time (RETENTIONTIME), number of peaks (Num Peaks), together with fragments and their 
#' intensity values. In the example used in the function \code{convertMSP2MSP} 
#' the \code{matrix} \code{msp2msp} is used to construct an object of class \code{MSP}. 
#' @docType data
#' @usage msp2msp
#' @return \code{data.frame}
#' @format \code{data.frame}
#' @source http://prime.psc.riken.jp/Metabolomics_Software/MS-DIAL/, truncated
#' .MSP file of GNPS MS/MS Negative (contains 22 entries): 
#' http://prime.psc.riken.jp/Metabolomics_Software/MS-DIAL/MSMS-GNPS-Curated-Neg.msp
#' @author Thomas Naake, \email{thomasnaake@googlemail.com}
NULL  
