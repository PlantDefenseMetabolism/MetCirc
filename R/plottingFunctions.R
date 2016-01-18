#' @import shiny
#' @import circlize
#' @import scales
#' @name plotCircos
#' @title Circular plot to visualise similarity
#' @description Circular plot to visualise similarity
#' @usage plotCircos(dfNameGroup, linkMat, initialize = c(TRUE, FALSE), 
#'  featureNames = c(TRUE, FALSE), cexFeatureNames = 0.2, 
#'  groupName = c(TRUE, FALSE), links = c(TRUE, FALSE), highlight = c(TRUE, FALSE))
#' @param dfNameGroup data.frame containing column "group" and "name", which is 
#' a unique identifier of the feature
#' @param linkMat data.frame containing linked features in each row, has five columns 
#'  (group1, name1, group2, name2, NDP)
#' @param initialize logical, should plot be initialized?
#' @param featureNames logical, should feature names be displayed?
#' @param cexFeatureNames numerical, size of feature names
#' @param groupName logical, should group names (e.g. compartment names or individual names) be displayed?
#' @param links logical, should links be plotted?
#' @param highlight logical, are we in highlighting mode?
#' @details Internal use for shiny app
#' @value The function will initialize a circlize plot and/or will plot 
#' features of a circlize plot. 
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{plotCircos(dfNameGroup, linkMat, initialize = TRUE, featureNames = TRUE, 
#' cexFeatureNames = 0.2, groupName = TRUE, links = TRUE)}
#' @export
plotCircos <- function (dfNameGroup, linkMat, initialize = c(TRUE, FALSE), 
        featureNames = c(TRUE, FALSE), cexFeatureNames = 0.2, 
        groupName = c(TRUE, FALSE), links = c(TRUE, FALSE), 
        highlight = c(TRUE, FALSE)) {
    
    if (!is.numeric(cexFeatureNames)) stop("cexFeatureNames is not numeric")
    if (!is.logical(initialize)) stop("initialize is not logical")
    if (!is.logical(featureNames)) stop("featureNames is not logical")
    if (!is.logical(groupName)) stop("groupName is not logical")
    if (!is.logical(links)) stop("links is not logical")
    if (!is.logical(highlight)) stop("highlight is not logical")

    dfDim <- dim(dfNameGroup)
    dfName <- as.character(dfNameGroup$name)
    if (groupName) dfGroup <- as.character(dfNameGroup$group)
    
    if (initialize) {
        circos.initialize(as.factor(dfName),
                          xlim=matrix(rep(c(0,1), dfDim[1]), ncol = 2, byrow = TRUE))
        circos.trackPlotRegion(as.factor(dfName), ylim=c(0,1))
    }
    
    ## feature names
    if (featureNames) {
        .dfNameGroup <- truncateName(dfNameGroup, nameGroup = TRUE)
        for (i in 1:dfDim[1]) {
            circos.text(x = 0.5, y = 0.5, labels = .dfNameGroup[,"name"][i],
                        sector.index = dfName[i], 
                        facing = "clockwise", cex = as.numeric(cexFeatureNames), 
                        niceFacing = TRUE)
        }
    }
    
    ## group name
    if (groupName) {
        uniqueGroup <- unique(dfGroup)
        transparency <- if (highlight) 0.05 else 0.2
        for( i in 1:length(uniqueGroup)) {
            ind <- which(uniqueGroup[i] == dfGroup)
            minInd <- min(ind)
            maxInd <- max(ind)
            circlize::highlight.sector(dfName[minInd:maxInd], 
                             col = alpha(i + 1, transparency))
            circlize::circos.text(x = 0.5, y = 1.5, labels = uniqueGroup[i], 
                        sector.index = dfName[c(minInd:maxInd)[floor(length(minInd:maxInd) / 2)]],
                        facing = "downward")    
        }
    }
    
    ## plot links
    if (links) {
        colourLink <- if (highlight) {
            rep(alpha("black", 0.1), dim(linkMat)[1]) 
            } else {
                alpha("black", alpha = (as.numeric(linkMat[,"NDP"]))^6)}
        
        for (i in 1:dim(linkMat)[1]) {
            circos.link(linkMat[i,][["name1"]], 0.5,
                        linkMat[i,][["name2"]], 0.5,
                        lwd = if (highlight) 0.5 else as.numeric(linkMat[i,][["NDP"]]),
                        ## transparency
                        col = colourLink[i])
        }
    }

}

#' @name highlight
#' @title Add links and highlight sectors
#' @description A function to add links and highlight sectors to an initialised
#'      and plotted \code{circlize} plot with one track.
#' @usage highlight(dfNameGroup, ind, LinkMatrix)
#' @param dfNameGroup data.frame , data.frame with group and unique idenfier (name)
#' @param ind numerical, indices which will be highlighted
#' @param LinkMatrix matrix, in each row there is information about features to be connected 
#' @details Internal use for shiny app.
#' @value The function will update an existing plot by highlighting a 
#'  specified sector and connected links.
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{plotCircosName(dfNameGroup, LinkMatrix, initialize = TRUE, featureNames = TRUE, groupName = TRUE, links = FALSE, highlight = FALSE); highlight(dfNameGroup, ind, LinkMatrix)}
#' @export
highlight <- function(dfNameGroup, ind, LinkMatrix) {
    dfDim <- dim(dfNameGroup)

    dfInd <- dfNameGroup[ind,]
    lMatName1 <- LinkMatrix[,"name1"]
    lMatName2 <- LinkMatrix[,"name2"]
    
    for (h in 1:length(ind)) {
        circlize::highlight.sector(sector.index = as.character(dfInd[h,"name"]), 
            col = alpha(palette()[as.numeric(dfInd[h,"group"]) + 1], 0.4))
    }
    
    ## get indices in LinkMatrix of selected features 
    #LinkMatrixInd <- which(LinkMatrix == as.character(dfInd[,"name"]), arr.ind = TRUE)[,1]
    
    LinkMatrixInd <- getLinkMatrixIndices(dfInd, LinkMatrix)

    #######
    ## plot all links
    for (i in 1:dim(LinkMatrix)[1]) {
        circos.link(lMatName1[i], 0.5,
                    lMatName2[i], 0.5,
                    lwd = 0.5,
                    ## transparency
                    col = alpha("black", 0.1))
    }
    #######
    
    ## plot highlighted links
    for (j in LinkMatrixInd) {
        circos.link(lMatName1[j], 0.5,
                    lMatName2[j], 0.5,
                    lwd = 1,
                    ## transparency
                    col = "black")
    }
}


#' @name getLinkMatrixIndices
#' @title Get indices in LinkMatrix of feature 
#' @description Gets indices in LinkMatrix of feature 
#' @usage getLinkMatrixIndices(dfInd, LinkMatrix)
#' @param dfInd row entry of data.frame dfNameGroup (selected feature)
#' @param LinkMatrix matrix, in each row there is information about features to be connected 
#' @details Internal use for function highlight.
#' @value returns indices concerning LinkMatrix to which dfInd connects
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{getLinkMatrixIndices(dfInd, LinkMatrix)}
getLinkMatrixIndices <- function(dfInd, LinkMatrix) {
    LinkMatrixInd <- lapply(as.character(dfInd[, "name"]), function(x) which(LinkMatrix == x, arr.ind = TRUE))
    LinkMatrixInd <- lapply(LinkMatrixInd, function(x) x[,1]) ## select only first column 
    LinkMatrixInd <- unlist(LinkMatrixInd) 
    return(LinkMatrixInd)
}

#' @name truncateName
#' @title Truncate names
#' @description A function to truncate names
#' @usage truncateName(dfNameGroup, roundDigits = 2, nameGroup = FALSE)
#' @param dfNameGroup data.frame with group and unique idenfier (name)
#' @param roundDigits numeric, how many digits should be displayed?
#' @param nameGroup logical, if TRUE name contains group name (e.g. "a_123/456")
#' @details The column \code{name} of \code{dfNameGroup} is a vector of \code{character} 
#'      strings of type consisting of retention time and m/z value. It is 
#'      cumbersome to display such strings. \code{truncateName} truncates 
#'      these strings by rounding retention time and m/z values by 
#'      digits given by \code{roundDigits}. \code{truncateName} is an 
#'      internal function.
#' @value returns dfNameGroup with truncated names
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{truncateName(dfNameGroup, roundDigits = 2, nameGroup = FALSE)}
truncateName <- function (dfNameGroup, roundDigits = 2, nameGroup = FALSE) {
    names <- as.character(dfNameGroup$name)
    
    if (nameGroup) {
        truncateL <- strsplit(names, split = "_")
        truncateL <- lapply(truncateL, "[", 2)
        names <- unlist(truncateL)
    }
    
    truncateL <- strsplit(names, "/")
    truncateL <- lapply(truncateL, function(x) 
        c(round(as.numeric(x[1]), roundDigits), round(as.numeric(x[2]), roundDigits)))
    newName <- lapply(truncateL, function(x) paste(x[1], x[2], sep="/"))
    newName <- unlist(newName)
    dfNameGroup$name <- newName
    return(dfNameGroup)
}

#' @name minFragCart2Polar
#' @title Calculate the nearest feature in polar coordinates given cartesian
#' coordinates
#' @description Calculates the nearest feature in polar coordinates given 
#'  cartesian coordinates
#' @usage minFragCart2Polar(x, y, degreeOfFeatures)
#' @param x cartesian x coordinate
#' @param y cartesian y coordinate
#' @param degreeOfFeatures list of positions of features
#' @details \code{minFragCart2Polar} is employed to find the feature with 
#'  the smallest distance from given cartesian coordinates. 
#' @value \code{minFragCart2Polar} returns the index of the feature that has the
#'  smallest distance to the given coordinates. As \code{minFragCart2Polar} is 
#'  used in \code{shinyCircos} for the track 1 only polar r coordinates between
#'  0.8 and 1 will be used to find the feature with smallest distance.
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{minFragCart2Polar(x, y, degreeOfFeatures)}
minFragCart2Polar <- function(x, y, degreeOfFeatures) {
    polar <- cart2Polar(x, y)
    minInd <- NA
    if (polar$r <= 1 & polar$r >= 0.8) 
        minInd <- which.min(abs(polar$theta - unlist(degreeOfFeatures)))
    return(minInd)
}

#' @name cart2Polar
#' @title Calculate polar coordinates from cartesian coordinates
#' @description cart2Polar calculates polar coordinates from cartesian coordinates
#' @usage cart2Polar(x, y)
#' @param x cartesian x coordinate
#' @param y cartesian y coordinate
#' @details \code{cart2Polar} is employed to translate cartesian coordinates 
#'  into polar coordinates especially in interactive shiny applications when
#'  using hovering and clicking features.
#' @return list of r and theta
#' @value polar coordinates r and theta
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{x <- 1; y <- 1; cart2Polar(x, y)}
cart2Polar <- function(x, y) {
    r <- sqrt( x ^ 2 + y ^ 2)
    thetaP <- atan( y/x ) * 180 / pi
    if (x == 0 & y == 0) thetaP <- 0
    if (x >= 0 & y >= 0) theta <- thetaP ## 1st quadrant
    if (x < 0 & y >= 0) theta <- thetaP + 180 ## 2nd quadrant
    if (x < 0 & y < 0) theta <- thetaP + 180 ## 3rd quadrant
    if (x >= 0 & y < 0) theta <- thetaP + 360 ## 4th quadrant
    
    return(list(r=r, theta=theta))
}
## test
cart2Polar(0, 0)
cart2Polar(1, 1)
cart2Polar(0, 1)
cart2Polar(-1, 1)
cart2Polar(-1, -1)
cart2Polar(1, -1)
##

    

#plot(x = as.numeric(msp[7:12, 1]), y = as.numeric(msp[7:12,2]), type = "h", col = "black", ylim = c(-100, 100))
#lines(x = (as.numeric(msp[7:12, 1]) + 0.0), y = - as.numeric(msp[7:12,2]), type = "h", col = "blue")
#lines(x = (as.numeric(msp[7:12, 1]) + 0.1), y = - as.numeric(msp[7:12,2]), type = "h", col = "red")
#lines(x = (as.numeric(msp[7:12, 1]) + 0.2), y = - as.numeric(msp[7:12,2]), type = "h", col = "green")
