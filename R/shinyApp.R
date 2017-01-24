#' @import grDevices
#' @import graphics
#' @name shinyCircos
#' @title Interactive visualisation of similarity and navigation of MS/MS features
#' @description Visualise the similarity of MS/MS features.
#' @usage shinyCircos(similarityMatrix, msp, size = 400)
#' @param similarityMatrix matrix, similarityMatrix contains pair-wise 
#' similarity coefficients which give information about the similarity between
#' MS/MS features
#' @param msp MSP, an S4 object of class 'MSP', the MSP object will be used
#' to display information about the selected feature
#' @param size numerical, image width/height in pixels
#' @details The function is based on the shiny and circlize package. The user 
#' can choose interactively thresholds, type of links (between or within groups), 
#' display information about MS/MS features, permanently select MS/MS features
#' and export selected precursors. When running shinyCircos with the object
#' of class MSP, annotation data of selected MS/MS features will be displayed.
#' @return shinyCircos returns a character vector with the (permanently)
#' selected precursors
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("idMSMStoMSP", package = "MetCirc")
#' ## truncate files
#' finalMSP <- finalMSP[c(1:20, 29:48, 113:132, 240:259)]
#' data("binnedMSP", package = "MetCirc")
#' binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' \dontrun{shinyCircos(similarityMat, finalMSP, size = 400)}
#' @export
shinyCircos <- function(similarityMatrix, msp = NULL, size = 400) {
    
    if (!is.numeric(size)) stop("size is not numerical")
    if (!is.null(msp)) if (class(msp) != "MSP") stop("msp is not of class MSP")
    
    ## circlize parameters
    circos.par(gap.degree = 0, cell.padding = c(0, 0, 0, 0), 
            track.margin = c(0.0, 0))
    
    ## create plots and assign to objects by recordPlot
    ## rt
    simMatRT <- createOrderedSimMat(similarityMatrix, order = "retentionTime")
    groupnameRT <- rownames(simMatRT)
    plotCircos(groupnameRT, NULL, initialize=TRUE, featureNames = TRUE, 
            groupSector = TRUE, groupName = FALSE, links = FALSE, 
            highlight = FALSE)
    PlotFilledRT <- recordPlot()
    ## get group and name from groupnameRT argument
    ## groupnameRT is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    groupRT <- lapply(strsplit(groupnameRT, split = "_"), "[", 1)
    groupRT <- unlist(groupRT)
    nameRT <- lapply(strsplit(groupnameRT, split = "_"), function (x) x[length(x)])
    nameRT <- unlist(nameRT)
    ## get degree of features
    degreeFeaturesRT <- lapply(groupnameRT, 
        function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
    plot.new()
     
    plotCircos(groupnameRT, NULL, initialize=TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = TRUE)
    PlotHighlightRT <- recordPlot()
    plot.new()
    
    ## mz
    simMatMZ <- createOrderedSimMat(similarityMatrix, order = "mz")
    groupnameMZ <- rownames(simMatMZ)
    plotCircos(groupnameMZ, NULL, initialize=TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = FALSE)
    PlotFilledMZ <- recordPlot()
    ## get group and name from groupnameMZ argument
    ## groupnameMZ is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    groupMZ <- lapply(strsplit(groupnameMZ, split = "_"), "[", 1)
    groupMZ <- unlist(groupMZ)
    nameMZ <- lapply(strsplit(groupnameMZ, split = "_"), function (x) x[length(x)])
    nameMZ <- unlist(nameMZ)
    ## get degree of features
    degreeFeaturesMZ <- lapply(groupnameMZ, 
        function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
    
    plot.new()
    plotCircos(groupnameMZ, NULL, initialize=TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = TRUE)
    PlotHighlightMZ <- recordPlot()
    plot.new()
    
    ## clustering
    simMatClustering <- createOrderedSimMat(similarityMatrix, 
                                            order = "clustering")
    groupnameClustering <- rownames(simMatClustering)
    plotCircos(groupnameClustering, NULL, initialize=TRUE, 
               featureNames = TRUE, groupSector = TRUE, groupName = FALSE, 
               links = FALSE, highlight = FALSE)
    PlotFilledCluster <- recordPlot()
    ## get group and name from groupnameMZ argument
    ## groupnameMZ is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    groupClustering <- lapply(strsplit(groupnameClustering, split = "_"), "[", 1)
    groupClustering <- unlist(groupClustering)
    nameClustering <- lapply(strsplit(groupnameClustering, split = "_"), 
                             function (x) x[length(x)])
    nameClustering <- unlist(nameClustering)
    ## get degree of features
    degreeFeaturesClust <- lapply(groupnameClustering,
        function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
    plot.new()
    
    plotCircos(groupnameClustering, NULL, initialize=TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = TRUE)
    PlotHighlightCluster <- recordPlot()
    plot.new()

    
    ui <- fluidPage(
       ## fluidRow(3, 
        fluidRow( 
            column(4, 
                wellPanel(
                    radioButtons("choiceLinks", "choose type of links", 
                        choices = c("all" = "all", "inter-class links" = "inter", 
                            "intra-class links" = "intra"),
                        selected = "all"),
                    sliderInput("threshold", "Threshold for similarity to display",
                        min = 0, max = 1, value = c(0.8, 1)),
                    radioButtons("order", "order within groups",
                        choices = c("clustering" = "clustering", "m/z" = "mz",
                                    "retention time" = "retentionTime"), 
                        selected = "mz"),
                    actionButton("resetClickIndices", "Reset features"),
                    actionButton("stop", "Stop and export \n selected features")
                )
            ),
            column(8, 
                fluidRow(    
                    column(8,
                        plotOutput("circos",
                            dblclick = "circosDblClick",
                            click = "circosClick",
                            ##hover = hoverOpts(id = "circosHover", delay = 100, 
                            ##    clip = TRUE, nullOutside = FALSE),
                            width = size, height = size)
                            ##brush = brushOpts(id = "circosBrush",
                            ##                  resetOnNew = TRUE)),
                    ),
                    column(4,  
                        plotOutput("circosLegend"))
                ), 
                        htmlOutput("clickConnectedFeature"),
                        verbatimTextOutput("dblClickFeature")
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
            
        
        ## ordering of features, use predefined groupname object
        GN <- reactive({
            if (input$order == "mz") GN <- groupnameMZ
            if (input$order == "retentionTime") GN <- groupnameRT
            if (input$order == "clustering") GN <- groupnameClustering
            GN
        })

        ## get degree of features
        degreeFeatures <- reactive({
            if (input$order == "mz") degFeatures <- degreeFeaturesMZ
            if (input$order == "retentionTime") degFeatures <- degreeFeaturesRT
            if (input$order == "clustering") degFeatures <- degreeFeaturesClust
            degFeatures
        })
        
        ## calculateLink0Matrix
        link0Matrix <- reactive(createLink0Matrix(simMat()))
        
        ## create reactive expression for LinkMatrix which is cut according to 
        ## set radioButton (input$choiceLinks)
        LinkMatrix_cut <- reactive(cutLinkMatrix(link0Matrix(), 
                                                 type = input$choiceLinks))
        
        ## threshold linkMatrix_cut
        LinkMatrix_threshold <- reactive(thresholdLinkMatrix(LinkMatrix_cut(), 
                                        input$threshold[1], input$threshold[2]))
        
        ## reactiveValues for click Coordinates
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

        ## is mouse over the track 1?
        onCircle <- reactiveValues(is = NULL)
        observe({
            if (!is.null(CoordinatesNewClick$X)) {
                .dist <- sqrt(CoordinatesOldClick$X^2 + CoordinatesOldClick$Y^2)
                if (.dist >= 0.8 & .dist <= 1) {
                    onCircle$is <- TRUE 
                } else {
                    onCircle$is <- FALSE
                }
            } else onCircle$is <- FALSE
        })
        
        ## Click: which is the current sector?
        indClick <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosClick$x))
                indClick$ind <- minFragCart2Polar(input$circosClick$x, 
                                                  input$circosClick$y, 
                                                  degreeFeaturesMZ()) 
        })
        
        ## double click: which is the current sector?
        CoordinatesNewDblClick <- reactiveValues(X = 0, Y = 0)
        CoordinatesOldDblClick <- reactiveValues(X = 0, Y = 0)
        
        observe({
            if (!is.null(input$circosDblClick$x)) {
                CoordinatesNewDblClick$X <- input$circosDblClick$x
                CoordinatesNewDblClick$Y <- input$circosDblClick$y
                CoordinatesOldDblClick$X <- CoordinatesNewDblClick$X
                CoordinatesOldDblClick$Y <- CoordinatesNewDblClick$Y
            } else {
                CoordinatesNewDblClick$X <- CoordinatesOldDblClick$X
                CoordinatesNewDblClick$Y <- CoordinatesOldDblClick$Y
            }
        })
        
        ## reactive value which stores double clicked indices (inds = storage, 
        ## new = new indices)
        indDblClick <- reactiveValues(ind = NULL, new = NULL)
        observe({
            if (!is.null(input$circosDblClick$x)) {
                
                minInd <- minFragCart2Polar(input$circosDblClick$x, 
                                            input$circosDblClick$y, 
                                            degreeFeatures()) 
                if (!is.na(minInd)) {
                    GNselect <- GN()[minInd]
                    selected <- strsplit(GNselect, split = "_")[[1]]
                    groupSelected <- selected[1]
                    nameSelected <- selected[3] 
                    
                    newNG <- paste(groupSelected, nameSelected, sep = "_")
                    indDblClick$new <- newNG ## write truncated name to indDblClick$new
                } else  indDblClick$new <- NULL
            } #else indDblClick$new <- NULL
        })
        
        observe({
            input$resetClickIndices
            isolate(indDblClickMZ$ind <- NULL)
            isolate(indDblClickRT$ind <- NULL)
            isolate(indDblClickCluster$ind <- NULL)
            isolate(indDblClick$new <- NULL)
            ##isolate(indClick$ind <- NULL)
            ##isolate(indClickRT$ind <- NULL)
            ##isolate(indClickCluster$ind <- NULL)
            ##isolate(indClick$new <- NULL)
        })
        
        ## write double-clicked (truncated) names to indDblClickMZ, indDblClickRT, 
        ## indDblClickCluster
        indDblClickMZ <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosDblClick$x)) {
            if (!is.null(indDblClick$new)) {

                newMZ <- paste(groupMZ, nameMZ, sep = "_")
                newIndMZ <- match(indDblClick$new, newMZ)
                
                ##if (onCircle$is) {
                    if (isolate(newIndMZ %in% indDblClickMZ$ind)) {
                        indDblClickMZ$ind <- isolate(indDblClickMZ$ind[-which(newIndMZ == indDblClickMZ$ind)]) 
                    } else {indDblClickMZ$ind <- isolate(c(indDblClickMZ$ind, newIndMZ))}
                ##} 
            } 
            }
         })
        
        indDblClickRT <- reactiveValues(ind = NULL)
        observe({
             if (!is.null(input$circosDblClick$x)) {
                 if (!is.null(indDblClick$new)) {
         
                    newRT <- paste(groupRT, nameRT, sep = "_")
                    newIndRT <- match(indDblClick$new, newRT)
                 
                    ##if (onCircle$is) 
                    if (isolate(newIndRT %in% indDblClickRT$ind)) {
                        indDblClickRT$ind <- isolate(indDblClickRT$ind[-which(newIndRT == indDblClickRT$ind)])
                    } else {indDblClickRT$ind <- isolate(c(indDblClickRT$ind, newIndRT))}
            }
            }
        })
        # 
        indDblClickCluster <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosDblClick$x)) {
                if(!is.null(indDblClick$new)) {
                    newCl <- paste(groupClustering, nameClustering, sep = "_")
                    newIndCl <- match(indDblClick$new, newCl)
                 
                    ##if (onCircle$is) {
                    if (isolate(newIndCl %in% indDblClickCluster$ind)) {
                        indDblClickCluster$ind <- isolate(indDblClickCluster$ind[-which(newIndCl == indDblClickCluster$ind)])
                    } else {indDblClickCluster$ind <- isolate(c(indDblClickCluster$ind, newIndCl))}   
                    ##} 
                }
             }
        })
         
        ## plotting
        initializePlot <- reactive({
            circos.initialize(factor(GN()),
                xlim = matrix(rep(c(0,1), dim(similarityMatrix)[1]), ncol = 2, 
                byrow = TRUE) )
            circos.trackPlotRegion(GN(), ylim=c(0,1))  
        })
        
        output$circos <- renderPlot({
            initializePlot()
            ##if (!is.null(PlotFilled2)) {
            if (onCircle$is) {
                if (input$order == "mz") {
                    replayPlot(PlotHighlightMZ)
                    ##if (length(indDblClickMZ$ind) > 0) {
                        highlight(groupnameMZ, c(indClick$ind, indDblClickMZ$ind), 
                                  LinkMatrix_threshold())  
                    ##} else {
                    ##    plotCircos(groupnameMZ, LinkMatrix_threshold(), 
                    ##               initialize=FALSE, featureNames = FALSE, 
                    ##               groupSector = FALSE, groupName = FALSE, 
                    ##               links = TRUE, highlight = TRUE)
                    ##}
                    
                }
                        
                if (input$order == "retentionTime") {
                    replayPlot(PlotHighlightRT)
                    highlight(groupnameRT, c(indClick$ind, indDblClickRT$ind), 
                              LinkMatrix_threshold())    
                }
                    
                if (input$order == "clustering") {
                    replayPlot(PlotHighlightCluster)
                    highlight(groupnameClustering, c(indClick$ind, indDblClickCluster$ind), 
                              LinkMatrix_threshold())  
                }
            } else { ## if not onCircle$is
                if (length(indDblClickMZ$ind) > 0) {
                    if (input$order == "mz") {
                        replayPlot(PlotHighlightMZ)
                        highlight(groupnameMZ, c(indDblClickMZ$ind), LinkMatrix_threshold()) 
                    }
                    
                    if (input$order == "retentionTime") {
                        replayPlot(PlotHighlightRT)
                        highlight(groupnameRT, c(indDblClickRT$ind), LinkMatrix_threshold()) 
                    }
                    
                    if (input$order == "clustering") {
                        replayPlot(PlotHighlightCluster)
                        highlight(groupnameClustering, c(indDblClickCluster$ind), LinkMatrix_threshold())  
                    }
                    
                } else {
                    if (input$order == "mz") {
                        replayPlot(PlotFilledMZ)
                        plotCircos(groupnameMZ, LinkMatrix_threshold(), 
                            initialize=FALSE, featureNames = FALSE, 
                            groupSector = FALSE, groupName = FALSE, 
                            links = TRUE, highlight = FALSE)
                    }
                        
                        
                    if (input$order == "retentionTime") {
                        replayPlot(PlotFilledRT)
                        plotCircos(groupnameRT, LinkMatrix_threshold(), 
                            initialize=FALSE, featureNames = FALSE, 
                            groupSector = FALSE, groupName = FALSE, 
                            links = TRUE, highlight = FALSE)
                    }
                    
                    if (input$order == "clustering") {
                        replayPlot(PlotFilledCluster)
                        plotCircos(groupnameClustering, LinkMatrix_threshold(), 
                            initialize = FALSE, featureNames = FALSE,
                            groupSector = FALSE, groupName = FALSE, 
                            links = TRUE, highlight = FALSE)
                    }
                }
            }
        })
        
        output$circosLegend <- renderPlot({
            circosLegend(groupnameRT, highlight = TRUE)
        })
        
        ## show when Clicking the feature which connects to it
        linkMatIndsClick <- reactive({
            getLinkMatrixIndices(GN()[indClick$ind], LinkMatrix_threshold())
        })
        
        output$clickConnectedFeature <- renderUI({ 
            if (!is.null(onCircle$is)) {
                if (onCircle$is)
                    HTML(printInformationSelect(GN(), msp = msp, 
                        ind = indClick$ind, lMatInd = linkMatIndsClick(), 
                        linkMatrixThreshold = LinkMatrix_threshold(), 
                        similarityMatrix = simMat()))  
            }
        })
        
        output$dblClickFeature <- renderText({
            if (length(indDblClickMZ$ind) > 0) 
                c("selected features: ", 
                    paste(groupMZ[indDblClickMZ$ind], 
                          nameMZ[indDblClickMZ$ind],
                        ##sapply(strsplit(nameMZ[indDblClickMZ$ind], split="_"), function(x) x[3]),
                        sep="_"))
            else "no features selected"
        })
        
        ## on exit
        observe({
            if (input$stop == 0)
                return()
            else {
                circos.clear()
                stopApp(as.character(paste(
                    groupMZ[indDblClickMZ$ind], 
                    nameMZ[indDblClickMZ$ind],
                    ##sapply(strsplit(nameMZ[indDblClickMZ$ind], split="_"), function(x) x[3]),
                    sep="_")))
            }
        })
        
    }
    
    app <- list(ui = ui, server = server)
    runApp(app)
}
## to do
## second simmat for neutral losses

#' @name printInformationSelected
#' @title Display information on connected features of selected features
#' @description Displays information on connected features of selected features.
#' @usage printInformationSelect(groupname, msp = NULL, ind, 
#'  lMatInd, linkMatrixThreshold, similarityMatrix)
#' @param groupname vector with groupname of selected feature,
#' vector containing "group" and "name" to display, that is 
#' a unique identifier of the features, "group" and "name" have to be separated
#' by "_" where "group" is the first and "name" is the last element
#' @param msp MSP, an S4 object of class 'MSP' for information about 
#'  the selected feature
#' @param ind numeric
#' @param lMatInd numeric indices of connected features
#' @param linkMatrixThreshold matrix that contains information of linked 
#'  features of a threshold or greater
#' @param similarityMatrix matrix that is used to get information on the degree 
#'  of similarity, similarityMat is an ordered version of a similarity matrix
#' @details printInformationSelect is for internal use. 
#' @return character that is in HTML format
#' @examples
#' data("idMSMStoMSP", package = "MetCirc")
#' data("binnedMSP", package = "MetCirc")
#' ## use only a selection
#' binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' ## order similarityMat according to mz
#' simMat <- createOrderedSimMat(similarityMat, order = "mz")
#' groupname <- rownames(simMat)
#' linkMat_thr <- createLinkMatrix(simMat, 0.9, 1) 
#' ind <- 19
#' linkMatInds <- getLinkMatrixIndices(groupname[ind], linkMat_thr)
#' MetCirc:::printInformationSelect(groupname = groupname, 
#'  msp = NULL, ind = ind, lMatInd = linkMatInds, 
#'  linkMatrixThreshold = linkMat_thr, 
#'  similarityMatrix = simMat)
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
printInformationSelect <- function(groupname, msp = NULL, 
                ind, lMatInd, linkMatrixThreshold, similarityMatrix) {
    ## get group and name from groupname argument
    ## groupname is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    group <- lapply(strsplit(groupname, split = "_"), "[", 1)
    group <- unlist(group)
    name <- lapply(strsplit(groupname, split = "_"), function (x) x[length(x)])
    name <- unlist(name)
        
    lMatThr <- linkMatrixThreshold 
    if (is.null(msp)) {
        ## selected feature
        selectedFeat <- groupname[ind]
        ## get connected features
        connFeat <- unique(as.vector(lMatThr[lMatInd, c("name1", "name2")]))
        ## remove selectedFeat from connFeat
        if (selectedFeat %in% connFeat) 
            connFeat <- connFeat[-which(connFeat == selectedFeat)]
                
        if (length(lMatInd) > 0) {
            connFeat <- paste(connFeat, collapse = " <br/>")
            return(paste(c(selectedFeat, "connects to", "<br/>", connFeat), 
                            collapse = " "))
        } else 
            return(paste(c(selectedFeat, "does not connect to any feature"), 
                            collapse = " "))
            
    
    } else { ## if !is.null(msp)
            
        ## find clicked feature
        mzRTMSP <- paste(getPrecursorMZ(msp), getRT(msp), sep="/")
        matchedHovMZRT <- match(name, mzRTMSP)
        selectedFeat <- msp[matchedHovMZRT[ind]]
        selectFeat <- groupname[ind] 
        ## connected features
        connect <- unique(as.vector(lMatThr[lMatIndSelect, c("name1", "name2")]))
        ## remove duplicated hovFeat in connect
        if (selectFeat %in% connect) connect <- connect[-which(connect == selectFeat)]
        mzRTcon <- sapply(strsplit(connect, split="_"), function(x) x[3])
        
        if (length(connect) == 0) {
            return(paste0(selectFeat, " (", getName(clickedFeat), ", ", 
                getMetaboliteName(selectedFeat), ", ", 
                getMetaboliteClass(selectedFeat), ") ",
                 "does not connect to any feature"))
        } else {
            matchedConn <- match(mzRTcon, mzRTMSP)
            connFeat <- msp[matchedConn]
            connChar <- character()
            degreeSimilarity <- similarityMatrix[selectFeat, ]
            for (i in 1:length(connect)) {
                connFeatI <- connFeat[i]
                connectI <- connect[i]
                degreeSimilarityI <- round(degreeSimilarity[connectI],3)
                ## degSimI <- degreeSimilarity[connFeatI]
                newFeat <- paste0(connectI, " (", degreeSimilarityI, ", ", 
                    getName(connFeatI), ", ", getMetaboliteName(connFeatI), ", ", 
                    getMetaboliteClass(connFeatI), ")", "<br/>")
                    
                connChar <- c(connChar, newFeat)
            }
            connChar <- paste(connChar, collapse=" ")
            return(paste0(selectFeat, " (", getName(selectedFeat), ", ", 
                getMetaboliteName(selectedFeat), ", ", 
                getMetaboliteClass(selectedFeat), ") connects to ", 
                " <br/>", connChar))
        }
    }
}