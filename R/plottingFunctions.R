#' @import shiny
#' @import circlize
#' @import scales
#' @name plotCircos
#' @title Circular plot to visualise similarity
#' @description Circular plot to visualise similarity
#' @usage plotCircos(dfNameGroup, linkMat, initialize = c(TRUE, FALSE), 
#'      featureNames = c(TRUE, FALSE), cexFeatureNames = 0.2, 
#'      groupSector = c(TRUE, FALSE), groupName = c(TRUE, FALSE), 
#'      links = c(TRUE, FALSE), highlight = c(TRUE, FALSE))
#' @param dfNameGroup data.frame containing column "group" and "name", which is 
#' a unique identifier of the feature
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
#' @details Internal use for shiny app
#' @return The function will initialize a circlize plot and/or will plot 
#'  features of a circlize plot. 
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#'  ## load binnedMSP
#'  data("binnedMSP", package = "MetCirc")
#'  ## use only a selection 
#'  binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#'  similarityMat <- createSimilarityMatrix(binnedMSP)#'  
#'  namesPrec <- rownames(binnedMSP)
#'  dfNameGroup <- data.frame(group = unlist(lapply(strsplit(namesPrec, "_"), "[[", 1)), 
#'      name = namesPrec) 
#'  
#'  ## order according to compartment
#'  dfNameGroup <- dfNameGroup[order(dfNameGroup[,"group"]),] 
#'  dfNameGroupRT <- orderNames(dfNameGroup = dfNameGroup, 
#'      similarityMatrix = NULL, order = "retentionTime")
#'      
#'  ## create a new similarity matrix with updated rownames
#'  simM <- createOrderedSimMat(dfNameGroupRT, similarityMat)
#'  ## create link matrix
#'  linkMat <- createLinkMatrix(similarityMatrix = simM, threshold=0.8,
#'      dfNameGroup = dfNameGroupRT)
#'  ## cut link matrix (here: only display links between groups)
#'  linkMat_cut <- cutLinkMatrix(linkMat, type = "inter")
#'  ## set circlize paramters
#'  circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), 
#'          track.margin = c(0.0, 0))
#'  ## here set selectedFeatures arbitrarily
#'  selectedFeatures <- as.character(dfNameGroupRT[c(1,101,201,301),2])
#'  
#'  ## actual plotting
#'  plotCircos(dfNameGroupRT, linkMat_cut, initialize = TRUE, 
#'      featureNames = TRUE, cexFeatureNames = 0.2, groupSector = TRUE, 
#'      groupName = FALSE, links = FALSE, highlight = FALSE)
#' @export
plotCircos <- function (dfNameGroup, linkMat, initialize = c(TRUE, FALSE), 
        featureNames = c(TRUE, FALSE), cexFeatureNames = 0.2, 
        groupSector = c(TRUE, FALSE), groupName = c(TRUE, FALSE), 
        links = c(TRUE, FALSE), highlight = c(TRUE, FALSE)) {
    
    dfNameGroup <- dfNameGroup[order(dfNameGroup[, "name"]), ]

    if (!is.numeric(cexFeatureNames)) stop("cexFeatureNames is not numeric")
    if (!is.logical(initialize)) stop("initialize is not logical")
    if (!is.logical(featureNames)) stop("featureNames is not logical")
    if (!is.logical(groupSector)) stop("groupSector is not logical")
    if (!is.logical(groupName)) stop("groupName is not logical")
    if (!is.logical(links)) stop("links is not logical")
    if (!is.logical(highlight)) stop("highlight is not logical")

    
    dfDim <- dim(dfNameGroup)
    dfName <- as.character(dfNameGroup$name)
    dfGroup <- as.character(dfNameGroup$group)
    
    if (initialize) {
        circos.initialize(dfName,
                xlim=matrix(rep(c(0,1), dfDim[1]), ncol = 2, 
                byrow = TRUE))
        circos.trackPlotRegion(dfName, ylim=c(0,1))
    }
    
    ## display feature names
    if (featureNames) {
        .dfNameGroup <- truncateName(dfNameGroup, nameGroup = TRUE)
        for (i in 1:dfDim[1]) {
            circos.text(x = 0.5, y = 0.5, labels = .dfNameGroup[,"name"][i],
                        sector.index = dfName[i], 
                        facing = "clockwise", cex = as.numeric(cexFeatureNames), 
                        niceFacing = TRUE)
        }
    }
    
    ## group sector
    if (groupSector) {
        uniqueGroup <- unique(dfGroup)
        transparency <- if (highlight) 0.1 else 0.2
        for( i in 1:length(uniqueGroup)) {
            ind <- which(uniqueGroup[i] == dfGroup)
            minInd <- min(ind)
            maxInd <- max(ind)
            circlize::highlight.sector(dfName[minInd:maxInd], 
                                       col = alpha(i + 1, transparency))
        }
    }
    
    
    ## group name
    if (groupName) {
        uniqueGroup <- unique(dfGroup)
        for( i in 1:length(uniqueGroup)) {
            ind <- which(uniqueGroup[i] == dfGroup)
            minInd <- min(ind)
            maxInd <- max(ind)
            circlize::circos.text(x = 0.5, y = 1.5, labels = uniqueGroup[i], 
                     sector.index = dfName[c(minInd:maxInd)[floor(length(minInd:maxInd) / 2)]],
                     facing = "downward")    
        }
    }
    
#     if (dendrogram) {
#         if (is.null(similarityMatrix)) stop("no similarity matrix")
#         if (rownames(similarityMatrix) != dfName) stop("no correct feature names")
#         dfGroupLevel <- unique(dfGroup)
#         for (i in dfGroupLevel) {
#             inds <- which(dfGroup == i)
#             dfNameGroupLevel <- dfNameGroup[inds,]
#             simMatI <- similarityMatrix[inds, inds]
#             hClust <- hcluster(simMatI, method = "spearman") 
#             hClust <- as.dendrogram(hClust)
#             circos.initialize(dfNameGroupLevel[,1],xlim=c(1,0))
#             circos.dendrogram(hClust, facing="outside")
#         }
#     }
# 
#     load(paste0(system.file(package = "circlize"), "/extdata/bird.orders.RData"))
#     
#     labels = hc$labels  # name of birds
#     ct = cutree(hc, 6)  # cut tree into 6 pieces
#     n = length(labels)  # number of bird species
#     dend = as.dendrogram(hc)
#     
#     circos.par(cell.padding = c(0, 0, 0, 0))
#     circos.initialize(factors = "a", xlim = c(0, n)) # only one sector
#     max_height = attr(dend, "height")  # maximum height of the trees
#     circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.3, 
#                            panel.fun = function(x, y) {
#                                for(i in seq_len(n)) {
#                                    circos.text(i-0.5, 0, labels[i], adj = c(0, 0.5), 
#                                                facing = "clockwise", niceFacing = TRUE,
#                                                col = ct[labels[i]], cex = 0.7)
#                                }
#                            })
#     
#     require(dendextend)
#     dend = color_branches(dend, k = 6, col = 1:6)
#     
#     circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA, 
#                            track.height = 0.4, panel.fun = function(x, y) {
#                                circos.dendrogram(dend, max_height = max_height)
#                            })
    ## plot links
    if (links) {
        colourLink <- if (highlight) {
            rep(alpha("black", 0.1), dim(linkMat)[1]) 
            } else {
                alpha("black", alpha = (as.numeric(linkMat[,"NDP"]))^6)}
        
        if (dim(linkMat)[1] != 0) {
            for (i in 1:dim(linkMat)[1]) {
                circos.link(linkMat[i,][["name1"]], 0.5,
                    linkMat[i,][["name2"]], 0.5,
                    lwd = if (highlight) 0.5 else as.numeric(linkMat[i,][["NDP"]]),
                    ## transparency
                    col = colourLink[i])
            }
        }
    }

}

#' @name highlight
#' @title Add links and highlight sectors
#' @description A function to add links and highlight sectors to an initialised
#'      and plotted \code{circlize} plot with one track.
#' @usage highlight(dfNameGroup, ind, LinkMatrix)
#' @param dfNameGroup data.frame , data.frame with group and unique idenfier 
#'      (name)
#' @param ind numerical, indices which will be highlighted
#' @param LinkMatrix matrix, in each row there is information about features 
#'      to be connected 
#' @details Internal use for shiny app.
#' @return The function will update an existing plot by highlighting a 
#'  specified sector and connected links.
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#'  ## load binnedMSP
#'  data("binnedMSP", package = "MetCirc")
#'  ## use only a selection 
#'  binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#'  similarityMat <- createSimilarityMatrix(binnedMSP)#'  
#'  namesPrec <- rownames(binnedMSP)
#'  dfNameGroup <- data.frame(group = unlist(lapply(strsplit(namesPrec, "_"), "[[", 1)), 
#'      name = namesPrec) 
#'  
#'  ## order according to compartment
#'  dfNameGroup <- dfNameGroup[order(dfNameGroup[,"group"]),] 
#'  dfNameGroupRT <- orderNames(dfNameGroup = dfNameGroup, 
#'      similarityMatrix = NULL, order = "retentionTime")
#'      
#'  ## create a new similarity matrix with updated rownames
#'  simM <- createOrderedSimMat(dfNameGroupRT, similarityMat)
#'  ## create link matrix
#'  linkMat <- createLinkMatrix(similarityMatrix = simM, threshold=0.95,
#'      dfNameGroup = dfNameGroupRT)
#'  ## cut link matrix (here: only display links between groups)
#'  linkMat_cut <- cutLinkMatrix(linkMat, type = "inter")
#'  ## set circlize paramters
#'  circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), 
#'          track.margin = c(0.0, 0))
#'  ## here set selectedFeatures arbitrarily
#'  selectedFeatures <- as.character(dfNameGroupRT[c(1,21,41,61),2])
#'  
#'  ## actual plotting
#'  plotCircos(dfNameGroupRT, linkMat_cut, initialize = TRUE, 
#'      featureNames = TRUE, cexFeatureNames = 0.2, groupSector = FALSE, 
#'      groupName = TRUE, links = FALSE, highlight = TRUE)
#'  indSelected <- mapply(function(x) which(x == dfNameGroupRT$name), 
#'          selectedFeatures)
#'  ## highlight
#'  highlight(dfNameGroup = dfNameGroupRT, ind = indSelected, LinkMatrix = 
#'          linkMat_cut)
#' @export
highlight <- function(dfNameGroup, ind, LinkMatrix) {
    
    dfDim <- dim(dfNameGroup)
    
    dfInd <- dfNameGroup[ind,]
    lMatName1 <- LinkMatrix[,"name1"]
    lMatName2 <- LinkMatrix[,"name2"]
    
    for (h in 1:length(ind)) {
        highlight.sector(sector.index = as.character(dfInd[h,"name"]), 
            col = alpha(palette()[as.numeric(as.factor(dfNameGroup[,"group"])[ind])[h] + 1], 0.4))
    }
    
    ## get indices in LinkMatrix of selected features 
    if (dim(LinkMatrix)[1] != 0) 
        LinkMatrixInd <- getLinkMatrixIndices(dfInd, LinkMatrix)

    #######
    ## plot all links
    if (dim(LinkMatrix)[1] != 0) {
        for (i in 1:dim(LinkMatrix)[1]) {
            circos.link(lMatName1[i], 0.5,
                    lMatName2[i], 0.5,
                    lwd = 0.5,
                    ## transparency
                    col = alpha("black", 0.1))
        }
    }
    #######
    
    ## plot highlighted links
    if (exists("LinkMatrixInd")) {
        for (j in LinkMatrixInd) {
            circos.link(lMatName1[j], 0.5,
                    lMatName2[j], 0.5,
                    lwd = 1,
                    ## transparency
                    col = "black")
        }
    }
}

#' @name getLinkMatrixIndices
#' @title Get indices in LinkMatrix of feature 
#' @description Gets indices in LinkMatrix of feature 
#' @usage getLinkMatrixIndices(dfInd, linkMatrix)
#' @param dfInd row entry of data.frame dfNameGroup (selected feature)
#' @param linkMatrix matrix, in each row there is information about features 
#'      to be connected 
#' @details Internal use for function highlight.
#' @return \code{getLinkMatrixIndices} returns indices concerning linkMatrix to 
#'      which dfInd connects
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{getLinkMatrixIndices(dfInd, linkMatrix)}
#' @export
getLinkMatrixIndices <- function(dfInd, linkMatrix) {
    
    linkMatrixInd <- lapply(as.character(dfInd[, "name"]), 
                            function(x) which(linkMatrix == x, arr.ind = TRUE))
    ## select only first column
    linkMatrixInd <- lapply(linkMatrixInd, function(x) x[,1])  
    linkMatrixInd <- unlist(linkMatrixInd) 
    
    linkMatrixInd <- as.numeric(linkMatrixInd)
    
    return(linkMatrixInd)
}

#' @name truncateName
#' @title Truncate names
#' @description A function to truncate names
#' @usage truncateName(dfNameGroup, roundDigits = 2, nameGroup = FALSE)
#' @param dfNameGroup data.frame with group and unique idenfier (name)
#' @param roundDigits numeric, how many digits should be displayed?
#' @param nameGroup logical, if TRUE name contains group name (e.g. "a_123/456")
#' @details The column \code{name} of \code{dfNameGroup} is a vector of 
#'      \code{character} strings of type consisting of retention time and m/z 
#'      value. It is cumbersome to display such strings. \code{truncateName} 
#'      truncates these strings by rounding retention time and m/z values by 
#'      digits given by \code{roundDigits}. \code{truncateName} is an 
#'      internal function.
#' @return \code{truncateName} returns dfNameGroup with truncated names
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#'      \dontrun{truncateName(dfNameGroup, roundDigits = 2, nameGroup = FALSE)}
#' @export
truncateName <- function (dfNameGroup, roundDigits = 2, nameGroup = FALSE) {
    names <- as.character(dfNameGroup$name)
    
    if (nameGroup) {
        truncateL <- strsplit(names, split = "_")
        truncateL <- lapply(truncateL, "[", 3)
        names <- unlist(truncateL)
    }
    
    if (!nameGroup) {
        truncateL <- strsplit(names, split = "_") 
        truncateL <- lapply(truncateL, "[", 2)
        names <- unlist(truncateL)
    }
    
    truncateL <- strsplit(names, "/")
    truncateL <- lapply(truncateL, function(x) 
        c(round(as.numeric(x[1]), roundDigits), 
          round(as.numeric(x[2]), roundDigits)))
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
#' namesPrec <- rownames(binnedMSP)
#' dfNameGroup <- data.frame(group = unlist(lapply(strsplit(namesPrec, "_"), 
#'      "[[", 1)), name = namesPrec)
#' ## order according to compartment
#' dfNameGroup <- dfNameGroup[order(dfNameGroup[,"group"]),] 
#' plotCircos(dfNameGroup, NULL, initialize = TRUE, featureNames = FALSE, 
#'      groupName = FALSE, groupSector = TRUE, links = FALSE, highlight = FALSE)
#' x <- 1
#' y <- 0
#' degreeFeatures <- lapply(dfNameGroup$name, 
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



#' @import amap
#' @name orderNames
#' @title Reorder names in dfNameGroup according to retention time, m/z or 
#' cluster analysis
#' @description orderNames will order the entries of the column name of the 
#' dfNameGroup object according to a given sorting instruction.
#' @usage orderNames(dfNameGroup, similarityMatrix = NULL, 
#'      order = c("retentionTime", "mz", "clustering"))
#' @param dfNameGroup data.frame containing column "group" and "name", which is 
#' a unique identifier of the feature
#' @param similarityMatrix matrix, a similarity matrix that contains the 
#' NDP similarity measure between all precursors in the data set
#' @param order character, has to be one of "retentionTime", "mz" or 
#'  "clustering"
#' @details orderNames takes dfNameGroup, similarityMatrix (only needed
#' when order = "clustering") and order as arguments. It will rearrange the 
#' order in dfNameGroup according to the retention time, m/z or 
#' cluster analysis (based on similarity between fragmentation given by 
#' entries in the similarityMatrix) and return a reordered dfNameGroup object.
#' It will manipulate the entries in the name column of dfNameGroup by 
#' adding consecutively numbered names according to the given order argument. 
#' orderNames will extract the retention time and mz from the entries of the 
#' column of dfNameGroup as these entries are composed via the pattern
#' group_mz/rt. orderNames is used in the shinyCircos function, but can also 
#' be used externally to reconstruct circos plots outside the reactive 
#' environment (see vignette for further details).
#' @return dfNameGroup object (data.frame) with updated (ordered) column 
#' name
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' ## load binnedMSP
#' data("binnedMSP", package = "MetCirc")
#' ## use only a selection 
#' binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' namesPrec <- rownames(binnedMSP)
#' dfNameGroup <- data.frame(group = unlist(lapply(strsplit(namesPrec, "_"), "[[", 1)), 
#'      name = namesPrec) 
#' ## order according to compartment
#' dfNameGroup <- dfNameGroup[order(dfNameGroup[,"group"]),] 
#' ## order according to retention time
#' dfNameGroupRT <- orderNames(dfNameGroup, similarityMat, 
#'  order = "retentionTime")
#' @export
orderNames <- function(dfNameGroup, similarityMatrix = NULL, 
                       order = c("retentionTime", "mz","clustering")) {
    
    order <- match.arg(order)
    if (order == "clustering" & is.null(similarityMatrix)) stop("no similarity matrix")
    
    dfNameGroup[,2] <- as.character(dfNameGroup[,2]) 
    dfNameGroup[,1] <- as.character(dfNameGroup[,1])
    
    dfGroup <- dfNameGroup[,"group"]
    dfGroupLevels <- unique(dfGroup)
    
    dfName <- as.character(dfNameGroup[,2])
    dfNameSplit <- strsplit(dfName, split = "_")
    
    newDfNameGroup <- dfNameGroup
    
    if (order == "retentionTime") {
        for (i in dfGroupLevels) {
            inds <- which(dfGroup == i)
            dfNameGroupLevel <- dfNameGroup[inds,]
            dfNameSplitLevel <- dfNameSplit[inds]
            dfNameSplitLevelMZRT <- lapply(dfNameSplitLevel, 
                    function(x) strsplit(x[2], split="/"))
            mzrt <- lapply(dfNameSplitLevelMZRT, "[[", 1)
            rt <- lapply(mzrt, "[[", 2)
            rt <- unlist(rt)
            rt <- as.numeric(rt)
            
            dfNameGroupI <- dfNameGroupLevel[order(rt), ]
            dfNameNew <- as.character(dfNameGroupI[,2])
            lDfNameNew <- strsplit(dfNameNew, "_")
            ordered <- sprintf("%04d", 1:length(lDfNameNew))
            dfNameNew <- lapply(1:length(ordered), function (x) 
                paste(lDfNameNew[[x]][1], ordered[x], lDfNameNew[[x]][2], sep="_"))
            dfNameGroupI[,2] <- unlist(dfNameNew)
            newDfNameGroup[inds, ] <- dfNameGroupI
            
        }
    }
    
    if (order == "mz") {
        for (i in dfGroupLevels) {
            inds <- which(dfGroup == i)
            dfNameGroupLevel <- dfNameGroup[inds,]
            dfNameSplitLevel <- dfNameSplit[inds]
            dfNameSplitLevelMZRT <- lapply(dfNameSplitLevel, 
                function(x) strsplit(x[2], split="/"))
            mzrt <- lapply(dfNameSplitLevelMZRT, "[[", 1)
            mz <- lapply(mzrt, "[[", 1)
            mz <- unlist(mz)
            mz <- as.numeric(mz)
            
            dfNameGroupI <- dfNameGroupLevel[order(mz), ]
            dfNameNew <- as.character(dfNameGroupI[,2])
            lDfNameNew <- strsplit(dfNameNew, "_")
            ordered <- sprintf("%04d", 1:length(lDfNameNew))
            dfNameNew <- lapply(1:length(ordered), function (x) 
                paste(lDfNameNew[[x]][1], ordered[x], lDfNameNew[[x]][2], sep="_"))
            dfNameGroupI[,2] <- unlist(dfNameNew)
            newDfNameGroup[inds, ] <- dfNameGroupI

        }
    }
    
    if (order == "clustering") {
        for (i in dfGroupLevels) {
            inds <- which(dfGroup == i)
            dfNameGroupLevel <- dfNameGroup[inds,]
            simMatI <- similarityMatrix[inds, inds]
            hClust <- amap::hcluster(simMatI, method = "spearman") 
            dfNameGroupI <- dfNameGroupLevel[hClust$order,]
            dfNameNew <- as.character(dfNameGroupI[,2])
            lDfNameNew <- strsplit(dfNameNew, "_")
            ordered <- sprintf("%04d", 1:length(lDfNameNew))
            dfNameNew <- lapply(1:length(ordered), function (x) 
                paste(lDfNameNew[[x]][1], ordered[x], lDfNameNew[[x]][2], sep="_"))
            dfNameGroupI[,2] <- unlist(dfNameNew)
            newDfNameGroup[inds, ] <- dfNameGroupI
        }
    }
    
    newDfNameGroup <- newDfNameGroup[order(newDfNameGroup[,2]),]
    
    
    return(newDfNameGroup)
    
}