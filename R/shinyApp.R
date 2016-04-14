#' @import grDevices
#' @import graphics
#' @name shinyCircos
#' @title Interactive visualisation of similar precursors
#' @description Visualise similar precursors.
#' @usage shinyCircos(dfNameGroup, similarityMatrix, msp, size = 400)
#' @param dfNameGroup data.frame which contains columns "group" and "name"
#' @param similarityMatrix matrix, similarityMatrix contains pair-wise 
#' similarity coefficients which give information about the similarity between
#' precursors
#' @param msp MSP, an S4 object of class 'MSP' for information about 
#' the hovered feature
#' @param size numerical, image width/height in pixels
#' @details The function is based on the shiny and circlize package. Choose
#' interactively thresholds, type of links, hover over precursors, select 
#' precursors.
#' @return shinyCircos returns a character vector with the selected 
#' precursors
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{shinyCircos(dfNameGroup, similarityMatrix, msp, size = 400)}
#' @export
shinyCircos <- function(dfNameGroup, similarityMatrix, msp = NULL, size = 400) {
    
    if (!is.numeric(size)) stop("size is not numerical")
    if (!is.null(msp)) {
        if (class(msp) != "MSP") stop("msp is not of class MSP")
        ## test if mz/rt constructer of msp and names of dfNameGroup are 
        ## identical
        mzRTdf <- sapply(
            strsplit(as.character(dfNameGroup[,"name"]), split="_"), 
            function(x) x[2])
        mzRTMSP <- paste(getPrecursorMZ(msp), getRT(msp), sep="/")
        if(!all(mzRTdf == mzRTMSP)) 
            stop("mz/rt in msp and names of dfNameGroup are not identical")
    }
    ## circlize parameters
    circos.par(gap.degree = 0, cell.padding = c(0, 0, 0, 0), 
            track.margin = c(0.0, 0))
    
    ## create plots and assign to objects by recordPlot
    ## rt
    dfNameGroupRT <- orderNames(dfNameGroup, order = "retentionTime")
    simMatRT <- createOrderedSimMat(dfNameGroupRT, similarityMatrix)
    plotCircos(dfNameGroupRT, NULL, initialize=TRUE, featureNames = TRUE, 
            groupSector = TRUE, groupName = FALSE, links = FALSE, 
            highlight = FALSE)
    PlotFilledRT <- recordPlot()
    plot.new()
     
    plotCircos(dfNameGroupRT, NULL, initialize=TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = TRUE)
    PlotHighlightRT <- recordPlot()
    plot.new()
    
    ## mz
    dfNameGroupMZ <- orderNames(dfNameGroup, order = "mz")
    simMatMZ <- createOrderedSimMat(dfNameGroupMZ, similarityMatrix)
    plotCircos(dfNameGroupMZ, NULL, initialize=TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = FALSE)
    PlotFilledMZ <- recordPlot()
    plot.new()
    
    plotCircos(dfNameGroupMZ, NULL, initialize=TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = TRUE)
    PlotHighlightMZ <- recordPlot()
    plot.new()
    
    ## clustering
    dfNameGroupCluster <- orderNames(dfNameGroup, 
        similarityMatrix = similarityMatrix, order = "clustering")
    simMatClustering <- createOrderedSimMat(dfNameGroupCluster, similarityMatrix)
    plotCircos(dfNameGroupCluster, NULL, initialize=TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = FALSE)
    PlotFilledCluster <- recordPlot()
    plot.new()
    
    plotCircos(dfNameGroupCluster, NULL, initialize=TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = TRUE)
    PlotHighlightCluster <- recordPlot()
    plot.new()
    
    
    ui <- fluidPage(
       ## fluidRow(3, 
        fluidRow( ## was: sidebarPanel(
            column(4, 
                wellPanel(
                    radioButtons("choiceLinks", "choose type of links", 
                        choices = c("all" = "all", "inter-class links" = "inter", 
                            "intra-class links" = "intra"),
                        selected = "all"),
                    sliderInput("threshold", "Threshold for similarity to display",
                        min = 0, max = 1, value = 0.75),
                    radioButtons("order", "order within groups",
                        choices = c("retention time" = "retentionTime",
                            "m/z" = "mz", "clustering" = "clustering")),
                    actionButton("resetClickIndices", "Reset features"),
                    actionButton("stop", "Stop and export \n selected features")
                )
            ),
            column(8, 
                fluidRow(    
                    column(8,
                        plotOutput("circos",
                            click = "circosClick",
                            #dblclick = "circosDblClick",
                            hover = hoverOpts(id = "circosHover", delay = 150, 
                                clip = TRUE, nullOutside = FALSE),
                            width = size, height = size)
                            ##brush = brushOpts(id = "circosBrush",
                            ##                  resetOnNew = TRUE)),
                    ),
                    column(4,  
                        plotOutput("circosLegend"))
                ), 
                        htmlOutput("hoverConnectedFeature"),
                        verbatimTextOutput("clickFeature")
            )
        )
    )
 
    
    server <- function(input, output, session) {
        
        ## use predefined similarityMatrix
        simMat <- reactive({
            if (input$order == "mz") simMat <- simMatMZ
            if (input$order == "retentionTime") simMat <- simMatRT
            if (input$order == "clustering") simMat <- simMatClustering
            simMat
        })
            
        
        ## ordering of features, use predefined dfNameGroup object
        dfNG <- reactive({
            if (input$order == "mz") dfNG <- dfNameGroupMZ
            if (input$order == "retentionTime") dfNG <- dfNameGroupRT
            if (input$order == "clustering") dfNG <- dfNameGroupCluster
            dfNG
        })

        ## get degree of features
        features <- reactive(as.character( dfNG()[,"name"] ))
        degreeFeatures <- reactive(lapply(features(), 
            function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")])))
        
        ########################
        ## calculateLink0Matrix
        link0Matrix <- reactive(createLink0Matrix(simMat(), dfNG()))
        
        ## create reactive expression for LinkMatrix
        ## create reactive expression for LinkMatrix which is cut according to 
        ## set radioButton (input$choiceLinks)
        LinkMatrix_cut <- reactive(cutLinkMatrix(link0Matrix(), 
                                                 type = input$choiceLinks))
        
        ## threshold linkMatrix_cut
        LinkMatrix_threshold <- reactive(thresholdLinkMatrix(LinkMatrix_cut(), 
                                                             input$threshold))
        
        ## reactiveValues for hover Coordinates
        CoordinatesNewHover <- reactiveValues(X = 0, Y = 0)
        CoordinatesOldHover <- reactiveValues(X = 0, Y = 0)
        
        observe({
            if (!is.null(input$circosHover$x)) {
                CoordinatesNewHover$X <- input$circosHover$x
                CoordinatesNewHover$Y <- input$circosHover$y
                CoordinatesOldHover$X <- CoordinatesNewHover$X
                CoordinatesOldHover$Y <- CoordinatesNewHover$Y
            } else {
                CoordinatesNewHover$X <- CoordinatesOldHover$X
                CoordinatesNewHover$Y <- CoordinatesOldHover$Y
            }
        })

        ## is mouse over the track 1?
        onCircle <- reactiveValues(is = NULL)
        observe({
            if (!is.null(CoordinatesNewHover$X)) {
                .dist <- sqrt(CoordinatesOldHover$X^2 + CoordinatesOldHover$Y^2)
                if (.dist >= 0.8 & .dist <= 1) {
                    onCircle$is <- TRUE 
                } else {
                    onCircle$is <- FALSE
                }
            } else onCircle$is <- FALSE
        })
        
        ## Hover: which is the current sector?
        indHover <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosHover$x))
                indHover$ind <- minFragCart2Polar(input$circosHover$x, 
                                                  input$circosHover$y, 
                                                  degreeFeatures()) 
            
        })
        
        ## click: which is the current sector?
        CoordinatesNewClick <- reactiveValues(X = 0, Y = 0)
        CoordinatesOldClick <- reactiveValues(X = 0, Y = 0)
        
        observe({
            if (!is.null(input$circosClick$x)) {
                CoordinatesNewClick$X <- input$circosClick$x
                CoordinatesNewClick$Y <- input$circosClick$y
                CoordinatesOldClick$X <- CoordinatesNewClick$X
                CoordinatesOldClick$Y <- CoordinatesNewClick$Y
            } else {
                CoordinatesNewClick$X <- CoordinatesOldClick$X
                CoordinatesNewClick$Y <- CoordinatesOldClick$Y
            }
        })
        
        
        indClick <- reactiveValues(ind = NULL, new = NULL)
        observe({
            if (!is.null(input$circosClick$x)) {
                
                indClick$new <- minFragCart2Polar(input$circosClick$x, 
                                            input$circosClick$y, 
                                            degreeFeatures()) 
                if (isolate(indClick$new %in% indClick$ind))
                    indClick$ind <- isolate(indClick$ind[-which(indClick$new == indClick$ind)])
                else 
                    indClick$ind <- isolate(c(indClick$ind, indClick$new))
            } else indClick$new <- NULL
        })
        
        observe({
            input$resetClickIndices
            isolate(indClick$ind <- NULL)
            isolate(indClick$new <- NULL)
        })

        ## plotting
        initializePlot <- reactive(plotCircos(dfNG(), NULL, initialize = TRUE, 
                featureNames = FALSE, groupName = FALSE, groupSector = FALSE,
                links = FALSE, highlight = FALSE))
        output$circos <- renderPlot({
            initializePlot()
            ##if (!is.null(PlotFilled2)) {
            if (onCircle$is) {
                if (input$order == "mz") {
                    replayPlot(PlotHighlightMZ)
                    highlight(dfNameGroupMZ, c(indHover$ind, indClick$ind), 
                              LinkMatrix_threshold())  
                }
                        
                if (input$order == "retentionTime") {
                    replayPlot(PlotHighlightRT)
                    highlight(dfNameGroupRT, c(indHover$ind, indClick$ind), 
                              LinkMatrix_threshold())    
                }
                    
                if (input$order == "clustering") {
                    replayPlot(PlotHighlightCluster)
                    highlight(dfNameGroupCluster, c(indHover$ind, indClick$ind), 
                              LinkMatrix_threshold())  
                }
            } else { ## if not onCircle$is
                if (length(indClick$ind) > 0) {
                    if (input$order == "mz") {
                        replayPlot(PlotHighlightMZ)
                        highlight(dfNameGroupMZ, c(indClick$ind), LinkMatrix_threshold()) 
                    }
                    
                    if (input$order == "retentionTime") {
                        replayPlot(PlotHighlightRT)
                        highlight(dfNameGroupRT, c(indClick$ind), LinkMatrix_threshold()) 
                    }
                    
                    if (input$order == "clustering") {
                        replayPlot(PlotHighlightCluster)
                        highlight(dfNameGroupCluster, c(indClick$ind), LinkMatrix_threshold()) 
                    }
                    
                } else {
                    if (input$order == "mz") {
                        replayPlot(PlotFilledMZ)
                        plotCircos(dfNameGroupMZ, LinkMatrix_threshold(), 
                            initialize=FALSE, featureNames = FALSE, 
                            groupSector = FALSE, groupName = FALSE, 
                            links = TRUE, highlight = FALSE)
                    }
                        
                        
                    if (input$order == "retentionTime") {
                        replayPlot(PlotFilledRT)
                        plotCircos(dfNameGroupRT, LinkMatrix_threshold(), 
                            initialize=FALSE, featureNames = FALSE, 
                            groupSector = FALSE, groupName = FALSE, 
                            links = TRUE, highlight = FALSE)
                    }
                    
                    if (input$order == "clustering") {
                        replayPlot(PlotFilledCluster)
                        plotCircos(dfNameGroupCluster, LinkMatrix_threshold(), 
                            initialize = FALSE, featureNames = FALSE,
                            groupSector = FALSE, groupName = FALSE, 
                            links = TRUE, highlight = FALSE)
                    }
                }
            }
        })
        
        output$circosLegend <- renderPlot({
            circosLegend(dfNameGroup, highlight = onCircle$is)
        })
        
        ## show when hovering the feature which connects to it
        linkMatIndsHover <- reactive({
            getLinkMatrixIndices(dfNG()[indHover$ind,], LinkMatrix_threshold())
        })
        
        output$hoverConnectedFeature <- renderUI({
            if (!is.null(onCircle$is)) {
                HTML(printInformationHover(dfNameGroup, dfNG(), msp = msp, 
                        ind = indHover$ind, lMatIndHover = linkMatIndsHover(), 
                        linkMatrixThreshold = LinkMatrix_threshold(), 
                        highlight = onCircle$is))  
            }
            
#             if (onCircle$is) {
#                 ## hoveredFeat
#                 hoveredFeat <- as.character(dfNG()[indHover$ind,"name"])
#                 ## remove first column because it contains the precursor
#                 ## connectedFeature:
#                 connFeat <- unique(as.vector(LinkMatrix_threshold()[linkMatIndsHover(), c("name1", "name2")]))[-1]
#                 if (length(linkMatIndsHover() > 0))
#                     c(hoveredFeat, "connects to", connFeat)
#                 else 
#                     c(hoveredFeat, "does not connect to any feature")
#             }
        })
        
        output$clickFeature <- renderText({
            if (length(indClick$ind) > 0)
                c("selected features: ", 
                    as.character(dfNG()[indClick$ind, "name"]))
            else "no features selected"
        })
        
        ## on exit
        observe({
            if (input$stop == 0)
                return()
            else {
                circos.clear()
                stopApp(as.character(dfNG()[indClick$ind, "name"]))
            }
        })
        
    }
    
    app <- list(ui = ui, server = server)
    runApp(app)
}
## to do
## second simmat for neutral losses

#' @name printInformationHover
#' @title Display information on connected features of hovered features
#' @description Displays information on connected features of hovered features.
#' @usage printInformationHover(dfNameGroup, dfNameGroupOrder, msp = NULL, ind, 
#'  lMatIndHover = linkMatIndsHover, linkMatrixThreshold, highlight = TRUE)
#' @param dfNameGroup data.frame which contains columns "group" and "name"
#' @param dfNameGroupOrder data.frame, an ordered data.frame derived from 
#'  dfNameGroup which contains columns "group" and "name", the column "name" 
#'  contains entries in the form of group_number_mz/rt
#' @param msp MSP, an S4 object of class 'MSP' for information about 
#'  the hovered feature
#' @param ind numeric
#' @param lMatIndsHover numeric indices of connected features
#' @param linkMatrixThreshold matrix that contains information of linked 
#'  features of a threshold or greater
#' @param highlight logical only return character when set to TRUE
#' @param similarityMatrix matrix that is used to get information on the degree 
#'  of similarity, similarityMat is an ordered version of a similarity matrix
#' @details printInformationHover is for internal use. 
#' @return character that is in HTML format
#' @examples
#' data("idMSMStoMSP", package = "MetCirc")
#' data("binnedMSP", package = "MetCirc")
#' binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' namesPrec <- rownames(binnedMSP)
#' dfNameGroup <- data.frame(group = unlist(lapply(strsplit(namesPrec, "_"),
#'                                                 "[[", 1)), name = namesPrec)
#' ## order according to compartment
#' dfNameGroup <- dfNameGroup[order(dfNameGroup[,"group"]),]
#' dfNameGroupOrder <- orderNames(dfNameGroup, order = "mz")
#' simMatO <- createOrderedSimMat(dfNameGroupOrder, similarityMat)
#' linkMat_thr <- createLinkMatrix(simMatO, dfNameGroupOrder, 0.8) 
#' ind <- 2
#' linkMatIndsHover <- getLinkMatrixIndices(dfNameGroupOrder[ind,], linkMat_thr)
#' printInformationHover(dfNameGroup = dfNameGroup, 
#'  dfNameGroupOrder = dfNameGroupOrder, msp = NULL, ind = ind, 
#'  lMatIndHover = linkMatIndsHover, linkMatrixThreshold = linkMat_thr, 
#'  highlight = TRUE, similarityMatrix = simMatO)
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
printInformationHover <- function(dfNameGroup, dfNameGroupOrder, msp = NULL, 
                            ind, lMatIndHover = linkMatIndsHover,
                            linkMatrixThreshold, highlight = c(TRUE, FALSE),
                            similarityMatrix) {
    
    if (!is.null(highlight)) {
        dfNG <- dfNameGroup
        dfNGo <- dfNameGroupOrder
        lMatThr <- linkMatrixThreshold 
        if (is.null(msp)) {
            if (highlight) {
                ## hoveredFeat
                hoveredFeat <- as.character(dfNGo[ind,"name"])
                ## remove first column because it contains the precursor
                ## connectedFeature:
                connFeat <- unique(as.vector(lMatThr[lMatIndHover, c("name1", "name2")]))[-1]
                if (length(lMatIndHover > 0)) {
                    connFeat <- paste(connFeat, collapse = " <br/>")
                    return(
                        paste(c(hoveredFeat, "connects to", "<br/>", connFeat), 
                            collapse = " "))
                } else 
                    return(
                        paste(c(hoveredFeat, "does not connect to any feature"), 
                        collapse = " "))
            }
        
        } else {
            if (highlight) {
                ## find hovered feature
                mzRTdfOrder <- sapply(strsplit(as.character(dfNGo[,"name"]), split="_"), function(x) x[3])
                mzRTMSP <- paste(getPrecursorMZ(msp), getRT(msp), sep="/")
                matchedHovMZRT <- match(mzRTdfOrder, mzRTMSP)
                hoveredFeat <- msp[matchedHovMZRT[ind]]
                hovFeat <- dfNGo[ind, "name"]
                ## connected featuures
                connect <- unique(as.vector(lMatThr[lMatIndHover, c("name1", "name2")]))[-1]
                mzRTdfcon <- sapply(strsplit(connect, split="_"), function(x) x[3])
                
                
                if (length(connect) == 0) {
                    return(paste0(hovFeat, " (", getName(hoveredFeat), ", ", 
                        getMetaboliteName(hoveredFeat), ", ", 
                        getMetaboliteClass(hoveredFeat), ") ",
                         "does not connect to any feature"))
                } else {
                    matchedConn <- match(mzRTdfcon, mzRTMSP)
                    connFeat <- msp[matchedConn]
                    connChar <- character()
                    print(degreeSimilarity <- similarityMatrix[hovFeat, ])
                    for (i in 1:length(connFeat)) {
                        connFeatI <- connFeat[i]
                       ## degSimI <- degreeSimilarity[connFeatI]
                        newFeat <- paste0(connect[i], " (", getName(connFeatI),
                                      ", ", getMetaboliteName(connFeatI), ", ", 
                                      getMetaboliteClass(connFeatI), ")",  
                                      "<br/>")
                    
                        connChar <- c(connChar, newFeat)
                    }
                    connChar <- paste(connChar, collapse=" ")
                    return(paste0(hovFeat, " (", getName(hoveredFeat), ", ", 
                            getMetaboliteName(hoveredFeat), ", ",
                            getMetaboliteClass(hoveredFeat), 
                            ") connects to ", " <br/>", connChar))
                }
            }
        }
 
    }
}


#' @name createOrderedSimMat
#' @title Update a similarity matrix according to order of name column in 
#' dfNameGroup
#' @description Internal function for shiny application. May also be used 
#' outside of shiny to reconstruct figures.
#' @usage createOrderedSimMat(dfNameGroup, similarityMatrix)
#' @param dfNameGroup data.frame which contains columns "group" and "name"
#' @param similarityMatrix matrix, similarityMatrix contains pair-wise 
#' similarity coefficients which give information about the similarity between
#' precursors
#' @details createOrderSimMat takes a dfNameGroup data.frame and a 
#' similarity matrix as arguments. It will then reorder rows and columns of 
#' the similarityMatrix object such, that it matches the order of the 
#' column name in dfNameGroup. createOrderSimMat is used in the shinyCircos 
#' function to create similarityMatrix objects which will allow to switch
#' between different types of ordering in between groups (sectors) in the 
#' circos plot. It may be used as well externally, to reproduce plots outside
#' of the reactive environment (see vignette for a workflow).
#' @return createOrderedSimMat returns a similarity matrix with ordered
#' rownames according to the name column of dfNameGroup
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples 
#' data("binnedMSP", package = "MetCirc")
#' data("similarityMat", package = "MetCirc")
#' namesPrec <- rownames(binnedMSP)
#' dfNameGroup <- data.frame(group = unlist(lapply(strsplit(namesPrec, "_"), "[[", 1)), 
#'  name = namesPrec) 
#' ## order according to compartment
#' dfNameGroup <- dfNameGroup[order(dfNameGroup[,"group"]),] 
#' ## order according to retention time and create object dfNameGroupRT
#' dfNameGroupRT <- orderNames(dfNameGroup, NULL, order = "retentionTime")
#' createOrderedSimMat(dfNameGroupRT, similarityMatrix = similarityMat)
#' @export
createOrderedSimMat <- function(dfNameGroup, similarityMatrix) {
    
    if (!("name" %in% colnames(dfNameGroup))) 
        stop("dfNameGroup does not have column 'name'")
    
    ## order according to group
    dfNameGroup <- dfNameGroup[order(dfNameGroup[,"name"]),] 
    
    dfNameGroupName <- dfNameGroup[,"name"]
    dfNameGroupName <- as.character(dfNameGroupName)
    ## crop name
    dfNameSplit <- strsplit(dfNameGroupName, "_") 
    dfNameSplit <- lapply(dfNameSplit, function(x) x[c(1, 3)])
    dfName <- lapply(dfNameSplit, 
                     function(x) paste(x[1], x[2], sep="_"))
    dfName <- unlist(dfName)
    ## order according to cropped name
    simM <- similarityMatrix[dfName, dfName]
    rownames(simM) <- colnames(simM) <- dfNameGroupName
    
    return(simM)
}
##b <- shinyCircos(dfNameGroup, similarityMat, msp = finalMSP)
