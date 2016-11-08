#' @import shiny
#' @import circlize
#' @import scales
#' @name plotCircos
#' @title Circular plot to visualise similarity
#' @description Circular plot to visualise similarity
#' @usage plotCircos(groupname, linkMat, initialize = c(TRUE, FALSE), 
#'      featureNames = c(TRUE, FALSE), cexFeatureNames = 0.2, 
#'      groupSector = c(TRUE, FALSE), groupName = c(TRUE, FALSE), 
#'      links = c(TRUE, FALSE), highlight = c(TRUE, FALSE), colour = NULL,
#'      transparency = 0.2)
#' @param groupname vector containing "group" and "name" to display, that is 
#' a unique identifier of the features, "group" and "name" have to be separated
#' by "_" where "group" is the first and "name" is the last element
#' @param linkMat data.frame containing linked features in each row, has 
#'      five columns (group1, name1, group2, name2, NDP)
#' @param initialize logical, should plot be initialized?
#' @param featureNames logical, should feature names be displayed?
#' @param cexFeatureNames numerical, size of feature names
#' @param groupSector logical, should groups be displayed with background 
#'      colours?
#' @param groupName logical, should group names (e.g. compartment names or 
#'      individual names) be displayed?
#' @param links logical, should links be plotted? 
#' @param highlight logical, are we in highlighting mode?
#' @param colour NULL or character, colour defines the colours which are used
#'  for plotting, if NULL default colours are used
#' @param transparency numerical, defines the transparency of the colours
#' @details Internal use for shiny app
#' @return The function will initialize a circlize plot and/or will plot 
#'  features of a circlize plot. 
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' ## load binnedMSP
#' data("binnedMSP", package = "MetCirc")
#' ## use only a selection 
#' binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#' similarityMat <- createSimilarityMatrix(binnedMSP) 
#' ## order similarityMat according to retentionTime
#' simM <- createOrderedSimMat(similarityMat, order = "retentionTime")
#' ## create link matrix
#' linkMat <- createLinkMatrix(similarityMatrix = simM, 
#'      threshold_low=0.8, threshold_high=1)
#' ## cut link matrix (here: only display links between groups)
#' linkMat_cut <- cutLinkMatrix(linkMat, type = "inter")
#' ## set circlize paramters
#' circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), 
#'          track.margin = c(0.0, 0))
#' groupname <- rownames(simM)
#' ## actual plotting
#' plotCircos(groupname, linkMat_cut, initialize = TRUE, 
#'     featureNames = TRUE, cexFeatureNames = 0.2, groupSector = TRUE, 
#'      groupName = FALSE, links = FALSE, highlight = FALSE, colour = NULL, 
#'      transparency = 0.2)
#' @export
plotCircos <- function(groupname, linkMat, initialize = c(TRUE, FALSE), 
        featureNames = c(TRUE, FALSE), cexFeatureNames = 0.2, 
        groupSector = c(TRUE, FALSE), groupName = c(TRUE, FALSE), 
        links = c(TRUE, FALSE), highlight = c(TRUE, FALSE), 
        colour = NULL, transparency = 0.2) {

    ## get group and name from groupname argument
    ## groupname is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    group <- lapply(strsplit(groupname, split = "_"), "[", 1)
    group <- unlist(group)
    name <- lapply(strsplit(groupname, split = "_"), function (x) x[length(x)])
    name <- unlist(name)
    
    
    ## get length of vector groupname
    groupname_l <- length(groupname)
    
    if (!is.numeric(cexFeatureNames)) stop("cexFeatureNames is not numeric")
    if (!is.logical(initialize)) stop("initialize is not logical")
    if (!is.logical(featureNames)) stop("featureNames is not logical")
    if (!is.logical(groupSector)) stop("groupSector is not logical")
    if (!is.logical(groupName)) stop("groupName is not logical")
    if (!is.logical(links)) stop("links is not logical")
    if (!is.logical(highlight)) stop("highlight is not logical")
    if (!is.null(transparency)) {
        if(!is.numeric(transparency)) stop("transparency is not numeric")
    }
    
    
    if (initialize) {
        circos.initialize(factor(groupname),
                xlim = matrix(rep(c(0,1), groupname_l), ncol = 2, 
                byrow = TRUE) )
        circos.trackPlotRegion(groupname, ylim=c(0,1))
    }
    
    ## display feature names
    if (featureNames) {
        ##groupnameFeatName <- paste(group, name, sep="_")
        truncatedName <- truncateName(groupname)
        for (i in 1:groupname_l) {
            circos.text(x = 0.5, y = 0.5, labels = truncatedName[i],
                        sector.index = groupname[i], 
                        facing = "clockwise", cex = as.numeric(cexFeatureNames), 
                        niceFacing = TRUE)
        }
    }
    
    ## create vector with unique groups
    uniqueGroup <- unique(group)
    
    if (!is.null(colour)) {
        if (length(colour) != length(uniqueGroup)) {
            if (length(colour) != 1) {
                stop("length of colour does not match with length of group")
            }
        }
    }
    
    
    ## group sector
    if (groupSector) {

        transparency <- if (highlight) transparency - 0.1 else transparency 
        if (is.null(colour)) {
            colour <- alpha(palette()[as.numeric(as.factor(uniqueGroup))+1], transparency)
        } else {
            colour <- alpha(colour, transparency)
        }
        
        for( i in 1:length(uniqueGroup)) {
            ind <- which(uniqueGroup[i] == group)
            minInd <- min(ind)
            maxInd <- max(ind)
            circlize::highlight.sector(groupname[minInd:maxInd], 
                                       col = colour[i])
        }
    }
    
    ## group name
    if (groupName) {
        for( i in 1:length(uniqueGroup)) {
            ind <- which(uniqueGroup[i] == group)
            minInd <- min(ind)
            maxInd <- max(ind)
            circlize::circos.text(x = 0.5, y = 1.5, labels = uniqueGroup[i], 
                     sector.index = groupname[c(minInd:maxInd)[floor(length(minInd:maxInd) / 2)]],
                     facing = "downward")    
        }
    }
    
    ## plot links
    if (links) {
        ##colourLink <- rep(alpha("black", 0.05), dim(linkMat)[1]) 
        
        if (dim(linkMat)[1] != 0) {
            for (i in 1:dim(linkMat)[1]) {
                circos.link(linkMat[i,][["name1"]], 0.5,
                    linkMat[i,][["name2"]], 0.5,
                    lwd = if (highlight) 0.3 else max(0.5, as.numeric(linkMat[i,][["NDP"]])),
                    ## transparency
                    col = rep(alpha("black", 0.05))) ##colourLink[i])
            }
        }
    }
}

#' @name highlight
#' @title Add links and highlight sectors
#' @description A function to add links and highlight sectors to an initialised
#'      and plotted \code{circlize} plot with one track.
#' @usage highlight(groupname, ind, LinkMatrix, colour = NULL, transparency = 0.4)
#' @param groupname vector containing "group" and "name" to display, that is 
#' a unique identifier of the features, "group" and "name" have to be separated
#' by "_" where "group" is the first and "name" is the last element
#' @param ind numerical, indices which will be highlighted
#' @param LinkMatrix matrix, in each row there is information about features 
#'      to be connected 
#' @param colour NULL or character, colour defines the colours which are used
#'  for plotting, if NULL default colours are used
#' @param transparency numerical, defines the transparency of the colours
#' @details Internal use for shiny app.
#' @return The function will update an existing plot by highlighting a 
#'  specified sector and connected links.
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#'  ## load binnedMSP
#'  data("binnedMSP", package = "MetCirc")
#'  ## use only a selection 
#'  binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#'  similarityMat <- createSimilarityMatrix(binnedMSP)
#'  ## order similarityMat according to retentionTime and update rownames
#'  simM <- createOrderedSimMat(similarityMat, order = "retentionTime")
#'  ## create link matrix
#'  linkMat <- createLinkMatrix(similarityMatrix = simM, 
#'      threshold_low = 0.95, threshold_high = 1)
#'  ## cut link matrix (here: only display links between groups)
#'  linkMat_cut <- cutLinkMatrix(linkMat, type = "inter")
#'  ## set circlize parameters
#'  circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), 
#'          track.margin = c(0.0, 0))
#'  groupname <- rownames(simM)
#'  ## here: set selectedFeatures arbitrarily
#'  indSelected <- c(2,23,42,62)
#'  selectedFeatures <- groupname[indSelected]
#'  ## actual plotting
#'  plotCircos(groupname, linkMat_cut, initialize = TRUE, 
#'      featureNames = TRUE, cexFeatureNames = 0.2, groupSector = TRUE, 
#'      groupName = FALSE, links = FALSE, highlight = TRUE)
#'  ## highlight
#'  highlight(groupname = groupname, ind = indSelected, LinkMatrix = 
#'          linkMat_cut, colour = NULL, transparency = 0.4)
#' @export
highlight <- function(groupname, ind, LinkMatrix, colour = NULL, transparency = 0.4) {
    
    ## get group and name from groupname argument
    ## groupname is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    group <- lapply(strsplit(groupname, split = "_"), "[", 1)
    group <- unlist(group)
    name <- lapply(strsplit(groupname, split = "_"), function (x) x[length(x)])
    name <- unlist(name)
    
    ##if (length(colour))
    ## get length of vector namegroup
    groupname_l <- length(groupname)
    
    ## create vector that contains selected (ind) groupname instances
    groupnameselected <- groupname[ind]
    nameselected <- name[ind]
    
    ## retrieve name1 and name2 from LinkMatrix
    lMatName1 <- LinkMatrix[,"name1"]
    lMatName2 <- LinkMatrix[,"name2"]
    
    if (is.null(colour)) {
        colours <- alpha(palette()[as.numeric(as.factor(group))[ind]+1], transparency)   
    } else {
        colours <- alpha(colour[as.numeric(as.factor(group))[ind]], transparency)
    }
    
    for (h in 1:length(ind)) {
        highlight.sector(sector.index = as.character(groupnameselected[h]), 
            ##col = alpha(palette()[as.numeric(as.factor(group)[ind])[h] + 1], 0.4))
            col = colours[h])
    }
    
    ## get indices in LinkMatrix of selected features 
    if (dim(LinkMatrix)[1] != 0) {
        LinkMatrixInd <- getLinkMatrixIndices(groupnameselected, LinkMatrix)
    } else {LinkMatrixInd <- NULL}

    ## plot all links
    if (dim(LinkMatrix)[1] != 0) {
         for (i in 1:dim(LinkMatrix)[1]) {
             circos.link(lMatName1[i], 0.5,
                     lMatName2[i], 0.5,
                     lwd = 0.3,
                     ## transparency
                     col = alpha("black", 0.1))
         }
    }
     
    
    ## plot highlighted links
    if (!is.null(LinkMatrixInd)) {
        for (j in LinkMatrixInd) {
            circos.link(lMatName1[j], 0.5,
                    lMatName2[j], 0.5,
                    lwd = 1,
                    ## transparency
                    col = "black")
        }
    }
}


#' @name truncateName
#' @title Truncate names
#' @description A function to truncate names
#' @usage truncateName(groupname, roundDigits = 2)
#' @param groupname vector with group and unique idenfier (name)
#' @param roundDigits numeric, how many digits should be displayed?
#' @details \code{groupname} is a vector of \code{character} strings consisting 
#'      of a group, retention time and m/z value, separated by "_". It is 
#'      cumbersome to display such long strings. \code{truncateName} 
#'      truncates these strings by rounding retention time and m/z values by 
#'      digits given by \code{roundDigits}. \code{truncateName} is an 
#'      internal function.
#' @return \code{truncateName} returns groupname with truncated names 
#' without group)
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#'      groupname <- "a_100.12345/10.12345"
#'      truncateName(groupname, roundDigits = 2)
#' @export
truncateName <- function (groupname, roundDigits = 2) {
    
    ## select last element which is mz/retention time
    names <- lapply(strsplit(groupname, split = "_"), function (x) x[length(x)])
    names <- unlist(names)
    
    truncateL <- strsplit(names, "/")
    truncateL <- lapply(truncateL, function(x) 
        c(round(as.numeric(x[1]), roundDigits), 
          round(as.numeric(x[2]), roundDigits)))
    newName <- lapply(truncateL, function(x) paste(x[1], x[2], sep="/"))
    newName <- unlist(newName)
    return(newName)
}

#' @name circosLegend
#' @title Plot a legend for circos plot
#' @description circosLegend plots a legend for circos plot using group names .
#' @usage circosLegend(groupname, highlight = c(TRUE, FALSE), colour = NULL)
#' @param groupname vector containing "group" and "name" to display, that is 
#' a unique identifier of the features, "group" and "name" have to be separated
#' by "_" where "group" is the first and "name" is the last element
#' @param highlight logical, should colours be adjusted to highlight settings?
#' @param colour NULL or character, colour defines the colours which are used
#'  for plotting, if NULL default colours are used
#' @details Internal use for shiny app or outside of shiny to reproduce 
#'      figures.
#' @return The function will open a new plot and display colours together 
#'      with labels. 
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#'  ## load binnedMSP
#'  data("binnedMSP", package = "MetCirc")
#'  ## use only a selection 
#'  binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#'  similarityMat <- createSimilarityMatrix(binnedMSP)  
#'  groupname <- rownames(similarityMat)
#'  ## plot legend
#'  circosLegend(groupname, highlight = TRUE, colour = NULL)
#' @export
circosLegend <- function(groupname, highlight = c(TRUE, FALSE), colour = NULL) {
    
    ## get group and name from groupname argument
    ## groupname is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    group <- lapply(strsplit(groupname, split = "_"), "[", 1)
    group <- unlist(group)
    group <- as.factor(group)
    
    uniqNumGroup <- unique(as.numeric(group))
    
    if (is.null(colour)) {
        colours <- palette()[uniqNumGroup + 1]
    } else {
        colours <- colour[uniqNumGroup + 1]
    }
    
    plot(x=c(0,1), y=c(0,1), type="n", xlab = "", ylab = "",
         axes = FALSE, frame.plot = FALSE)
    if (highlight) {
        legend(x = c(0,1), y = c(1,0), legend = levels(group), 
               bty = "n",
               fill =  alpha(colours, 0.3),  border = alpha(colours, 0.3))
    } else { ## if not highlight
        legend(x = c(0,1), y = c(1,0), legend = levels(group), bty = "n",
               fill =  colours, border = colours)
    }
}

#' @name getLinkMatrixIndices
#' @title Get indices in LinkMatrix of feature 
#' @description Gets indices in LinkMatrix of feature 
#' @usage getLinkMatrixIndices(groupnameselected, linkMatrix)
#' @param groupnameselected vector with groupname of selected feature,
#' vector containing "group" and "name" to display, that is 
#' a unique identifier of the features, "group" and "name" have to be separated
#' by "_" where "group" is the first and "name" is the last element
#' @param linkMatrix matrix, in each row there is information about features 
#'      to be connected 
#' @details Internal use for function highlight.
#' @return \code{getLinkMatrixIndices} returns indices concerning linkMatrix to 
#'      which groupnameselected connects
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{getLinkMatrixIndices(groupnameselected, linkMatrix)}
#' @export
getLinkMatrixIndices <- function(groupnameselected, linkMatrix) {
    
    linkMatrixInd <- lapply(as.character(groupnameselected), 
                            function(x) which(linkMatrix == x, arr.ind = TRUE))
    ## select only first column
    linkMatrixInd <- lapply(linkMatrixInd, function(x) x[,1])  
    linkMatrixInd <- unlist(linkMatrixInd) 
    
    linkMatrixInd <- as.numeric(linkMatrixInd)
    
    return(linkMatrixInd)
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
#' @return \code{minFragCart2Polar} returns the index of the feature that has the
#'  smallest distance to the given coordinates. As \code{minFragCart2Polar} is 
#'  used in \code{shinyCircos} for the track 1 only polar r coordinates between
#'  0.8 and 1 will be used to find the feature with smallest distance.
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' ## load binnedMSP
#' data("binnedMSP", package = "MetCirc")
#' ## use only a selection 
#' binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#' simM <- createSimilarityMatrix(binnedMSP)
#' groupname <- rownames(simM)
#' plotCircos(groupname, NULL, initialize = TRUE, featureNames = FALSE, 
#'      groupName = FALSE, groupSector = FALSE, links = FALSE, highlight = FALSE)
#' x <- 1
#' y <- 0
#' degreeFeatures <- lapply(groupname, 
#'  function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
#' minFragCart2Polar(x, y, degreeOfFeatures = degreeFeatures)
#' @export
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
#' @return \code{cart2Polar} returns a list of colar coordinates r and theta
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' x <- 1; y <- 1
#' cart2Polar(x, y)
#' @export
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
