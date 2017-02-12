#' @import grDevices
#' @import graphics
#' @name shinyCircos
#' @title Interactive visualisation of similarity and navigation of MS/MS features
#' @description Visualise the similarity of MS/MS features in a reactive 
#'  context. See \code{Details} the vignette for further descriptions on how to use 
#'  \code{shinyCircos}.
#' @usage shinyCircos(similarityMatrix, msp = NULL, ...)
#' @param similarityMatrix \code{matrix}, \code{similarityMatrix} contains 
#' pair-wise similarity coefficients which give information about the similarity 
#' between MS/MS features
#' @param msp \code{MSP}, an S4 object of class \code{MSP}, the 
#'  \code{MSP}-object will be used to display information about the selected 
#'  feature
#' @param ... further arguments passed to \code{shinyCircos}, e.g. 
#' \code{cexFeatureNames} to pass to \code{plotCircos} to set font size in 
#' \code{plotCircos} of feature names
#' @details The function is based on the \code{shiny} and \code{circlize} package. 
#' The user can choose interactively thresholds, type of links (between or 
#' within groups), display information about MS/MS features, permanently select 
#' MS/MS features and export selected precursors. When running 
#' \code{shinyCircos} with the object of class \code{MSP}, annotation data of 
#' selected MS/MS features will be displayed.
#' @return \code{shinyCircos} returns a \code{character} vector with the 
#' (permanently) selected precursors or an object with the entries \code{msp}
#' and \code{selectedFeatures} if a \code{MSP}-object was passed to 
#' \code{shinyCircos}
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("idMSMStoMSP", package = "MetCirc")
#' ## truncate files
#' finalMSP <- finalMSP[c(1:20, 29:48, 113:132, 240:259)]
#' data("binnedMSP", package = "MetCirc")
#' binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' \dontrun{shinyCircos(similarityMatrix = similarityMat, msp = finalMSP)}
#' @export
shinyCircos <- function(similarityMatrix, msp = NULL, ...) {
    
    if (!is.null(msp)) if (class(msp) != "MSP") stop("msp is not of class MSP")
    
    ## circlize parameters
    circos.par(gap.degree = 0, cell.padding = c(0, 0, 0, 0), 
            track.margin = c(0.0, 0))
    
    groupname <- rownames(similarityMatrix)
    ## create plots and assign to objects by recordPlot
    ## rt
    simMatRT <- createOrderedSimMat(similarityMatrix, order = "retentionTime")
    link0MatRT <- createLink0Matrix(simMatRT)
    groupnameRT <- rownames(simMatRT)
    plotCircos(groupnameRT, NULL, initialize=TRUE, featureNames = TRUE, 
         groupSector = TRUE, groupName = FALSE, links = FALSE, 
         highlight = FALSE, ...)
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
            highlight = TRUE, ...)
    PlotHighlightRT <- recordPlot()
    plot.new()
    
    ## mz
    simMatMZ <- createOrderedSimMat(similarityMatrix, order = "mz")
    link0MatMZ <- createLink0Matrix(simMatMZ)
    groupnameMZ <- rownames(simMatMZ)
    plotCircos(groupnameMZ, NULL, initialize=TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = FALSE, ...)
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
               highlight = TRUE, ...)
    PlotHighlightMZ <- recordPlot()
    plot.new()
    
    ## clustering
    simMatClustering <- createOrderedSimMat(similarityMatrix, 
                                            order = "clustering")
    link0MatClustering <- createLink0Matrix(simMatClustering)
    groupnameClustering <- rownames(simMatClustering)
    plotCircos(groupnameClustering, NULL, initialize=TRUE, 
               featureNames = TRUE, groupSector = TRUE, groupName = FALSE, 
               links = FALSE, highlight = FALSE, ...)
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
               highlight = TRUE, ...)
    PlotHighlightCluster <- recordPlot()
    plot.new()
    
    ui <- fluidPage( 
            tags$head(tags$script('
                $(document).on("shiny:connected", function(e) {
                Shiny.onInputChange("innerWidth", window.innerWidth);
                });
                $(window).resize(function(e) {
                Shiny.onInputChange("innerWidth", window.innerWidth);
                });
                '
            )),
        column(4, 
            fluidRow(
                 tabsetPanel(id = "tabs",
                      tabPanel("Main", wellPanel(
                          radioButtons("choiceLinks", "choose type of links",
                              choices = c("all" = "all", "inter-class links" = "inter",
                                  "intra-class links" = "intra"),
                              selected = "all"),
                          sliderInput("threshold",
                              "Threshold for similarity to display",
                              min = 0, max = 1, value = c(0.8, 1)),
                          radioButtons("order", "order within groups",
                              choices = c("clustering" = "clustering",
                                  "m/z" = "mz", "retention time" = "retentionTime"),
                              selected = "mz"),
                          
                          uiOutput("annotationName"),
                          uiOutput("annotationClass"),
                          uiOutput("annotationInformation"),
                          uiOutput("annotationAdduct"),
                          uiOutput("annotationButton"),
                          actionButton("resetClickIndices", "Reset features"),
                          actionButton("stop", "Stop and export \n selected features")
                      )),
                      tabPanel("Appearance", wellPanel(
                          sliderInput("plotSize", "plot size",
                              min = 0.5, max = 1.5, value = 1),
                          sliderInput("precision", "precision of numbers", value = 2, min = 0, max = 5, step = 1),
                          checkboxInput("legend", "legend", value = FALSE)
                      ))
                     )
            ),
            plotOutput("circosLegend", height = "300")
        ),
        column(8,
            fluidRow(uiOutput("sized_plot")),
            fluidRow(
                verbatimTextOutput("dimension_display"),
                htmlOutput("clickConnectedFeature"),
                verbatimTextOutput("dblClickFeature")
            )
        )
    )
    
    server <- function(input, output, session) {
        
        
        ## annotation
        if(!is.null(msp)) {mspannotation <- reactiveValues(
            names = msp@names, information = msp@information, 
            classes = msp@classes, adduct = msp@adduct)}
        
        output$annotationName <- renderUI({
            if (!is.null(msp)) {
                if (length(indClick$ind) > 0 && onCircle$is) {
                    
                    textInput("names", label = "name", 
                              value = isolate(names(MSP())[indMSP()]))
                } else NULL  
            }
        })
        
        output$annotationClass <- renderUI({
            if (!is.null(msp)) {
                if (length(indClick$ind) > 0 && onCircle$is) {
                    textInput("classes", label = "class", 
                              value = isolate(classes(MSP())[indMSP()]))
                } else NULL  
            }
        })
        
        output$annotationInformation <- renderUI({
            if (!is.null(msp)) {
                if (length(indClick$ind) > 0 && onCircle$is) {
                #if (onCircle$is) {
                    textInput("information", label = "information", 
                              value = isolate(information(MSP())[indMSP()])) 
                } else NULL  
            }
        })
        output$annotationAdduct <- renderUI({
            if (!is.null(msp)) {
                if (length(indClick$ind) > 0 && onCircle$is) {
                    textInput("adduct", label = "adduct", 
                              value = isolate(adduct(MSP())[indMSP()])) 
                } else NULL  
            }
        })
        
        
        output$annotationButton <- renderUI({
            if (!is.null(msp)) {
                if (length(indClick$ind) > 0 && onCircle$is) 
                actionButton(inputId = "annotate", label = "update annotation")
                } 
        })
        
        indMSP <- reactive({
            if (length(indClick$ind) > 0) {
            nameClick <- GN()[indClick$ind]
            trNameClick <- truncateName(nameClick, roundDigits = NULL, group = TRUE)
            which(trNameClick == groupname)
            }
        })
    
        
        indMSPAnn <- eventReactive(input$annotate, {
            indMSP()
        })
        
        ## eventReactive for input$name
        annotateNames <- eventReactive(input$annotate, {
                    as.character(input$names)
        })
        ## eventReactive for input$classes
        annotateClasses <- eventReactive(input$annotate, {
                    as.character(input$classes)
        })
        ## eventReactive for input$information
        annotateInformation <- eventReactive(input$annotate, {
                    as.character(input$information)
        })
        ## eventReactive for input$adduct
        annotateAdduct <- eventReactive(input$annotate, {
            as.character(input$adduct)
        })
        
        if (!is.null(msp)) {
        observe({
            mspannotation$names[indMSPAnn()] <- annotateNames()
        })
        observe({
            mspannotation$classes[indMSPAnn()] <- annotateClasses()
        })
        observe({
            mspannotation$information[indMSPAnn()] <- annotateInformation()
        })
        observe({
                mspannotation$adduct[indMSPAnn()] <- annotateAdduct()
        })
        }
 
        ## reactive expression for msp
        MSP <- reactive({
            if (!is.null(msp)) {
                    msp <- new("MSP", msp = msp@msp, 
                    names = mspannotation$names, 
                    classes = mspannotation$classes,
                    information = mspannotation$information,
                    adduct = mspannotation$adduct, mz = msp@mz, rt = msp@rt)
            } else {msp <- NULL}
            msp
        })
        ## end annotation

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
        link0Matrix <- reactive({
            #if (!is.null(input$order)) {
                if (input$order == "mz") link0Mat <- link0MatMZ
                if (input$order == "retentionTime") link0Mat <- link0MatRT
                if (input$order == "clustering") link0Mat <- link0MatClustering
                link0Mat  
            #}
        })
        
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
                    input$circosClick$y, degreeFeatures()) 
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
            isolate(onCircle$is <- FALSE)
            isolate(indClick$ind <- NULL)
        })
        
        ## reset indClick when changing radio button order
        observe({
            input$order
            isolate(onCircle$is <- FALSE)
            isolate(indClick$ind <- NULL)
        })
           
        
        ## write double-clicked (truncated) names to indDblClickMZ, indDblClickRT, 
        ## indDblClickCluster
        indDblClickMZ <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosDblClick$x)) {
            if (!is.null(indDblClick$new)) {

                newMZ <- paste(groupMZ, nameMZ, sep = "_")
                newIndMZ <- match(indDblClick$new, newMZ)
                
                if (isolate(newIndMZ %in% indDblClickMZ$ind)) {
                    indDblClickMZ$ind <- isolate(indDblClickMZ$ind[-which(newIndMZ == indDblClickMZ$ind)]) 
                } else {indDblClickMZ$ind <- isolate(c(indDblClickMZ$ind, newIndMZ))}
            } 
            }
         })
        
        indDblClickRT <- reactiveValues(ind = NULL)
        observe({
             if (!is.null(input$circosDblClick$x)) {
                 if (!is.null(indDblClick$new)) {
         
                    newRT <- paste(groupRT, nameRT, sep = "_")
                    newIndRT <- match(indDblClick$new, newRT)
                 
                    if (isolate(newIndRT %in% indDblClickRT$ind)) {
                        indDblClickRT$ind <- isolate(indDblClickRT$ind[-which(newIndRT == indDblClickRT$ind)])
                    } else {indDblClickRT$ind <- isolate(c(indDblClickRT$ind, newIndRT))}
            }
            }
        })
        
        indDblClickCluster <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosDblClick$x)) {
                if(!is.null(indDblClick$new)) {
                    newCl <- paste(groupClustering, nameClustering, sep = "_")
                    newIndCl <- match(indDblClick$new, newCl)
                
                    if (isolate(newIndCl %in% indDblClickCluster$ind)) {
                        indDblClickCluster$ind <- isolate(indDblClickCluster$ind[-which(newIndCl == indDblClickCluster$ind)])
                    } else {indDblClickCluster$ind <- isolate(c(indDblClickCluster$ind, newIndCl))}  
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
            if (onCircle$is) {
                if (input$order == "mz") {
                    replayPlot(PlotHighlightMZ)
                    highlight(groupnameMZ, c(indClick$ind, indDblClickMZ$ind), 
                            LinkMatrix_threshold())  
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
                        #replayPlot(PlotFilledRT)
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
        
        output$sized_plot <- renderUI({
            plotOutput("circos",
                dblclick = "circosDblClick",
                click = "circosClick",
                width = ifelse(is.null(input$innerWidth), 0, input$innerWidth*0.5*input$plotSize), 
                height = ifelse(is.null(input$innerWidth), 0, input$innerWidth*0.5*input$plotSize))
        })
        
        output$circosLegend <- renderPlot({
            if (!is.null(input$legend)) if(input$legend)
                circosLegend(groupnameRT, highlight = TRUE)
        })
        
        ## show when Clicking the feature which connects to it
        linkMatIndsClick <- reactive({
            getLinkMatrixIndices(GN()[indClick$ind], LinkMatrix_threshold())
        })
        
        output$clickConnectedFeature <- renderUI({ 
            if (!is.null(onCircle$is)) {
                if (onCircle$is)
                    HTML(printInformationSelect(groupname = groupname, 
                        msp = MSP(), ind = indMSP(), lMatInd = linkMatIndsClick(), 
                        linkMatrixThreshold = LinkMatrix_threshold(), 
                        similarityMatrix = similarityMatrix, roundDigits = input$precision))  
            }
        })
        
        output$dblClickFeature <- renderText({
            if (length(indDblClickMZ$ind) > 0) 
                c("(permanently) selected features: ", 
                        truncateName(groupnameMZ[indDblClickMZ$ind], 
                                    roundDigits = input$precision, group = TRUE)
                  )
            else "no features (permanently) selected"
        })
        
        ## on exit
        observe({
            if (input$stop == 0)
                return()
            else {
                circos.clear()
                selectedFeatures <- as.character(paste(
                    groupMZ[indDblClickMZ$ind], nameMZ[indDblClickMZ$ind], 
                    sep="_"))
                stopApp(
                    if (!is.null(msp)) {
                        list(msp = MSP(), selectedFeatures = selectedFeatures)
                    } else {selectedFeatures}
                    )
            }
        })
        
    }
    
    app <- list(ui = ui, server = server)
    runApp(app)
}
## to do
## second simmat for neutral losses

#' @name printInformationSelect
#' @title Display information on connected features of selected features
#' @description Displays information on connected features of selected features.
#' @usage printInformationSelect(groupname, msp = NULL, ind, 
#'  lMatInd, linkMatrixThreshold, similarityMatrix, roundDigits = 2)
#' @param groupname \code{character} vector with groupname of selected feature,
#' vector containing "group" and "name" to display, that is 
#' a unique identifier of the features, "group" and "name" have to be separated
#' by \code{"_"} where "group" is the first and "name" is the last element
#' @param msp \code{MSP}, an S4 object of class \code{MSP} for information about 
#'  the selected feature
#' @param ind \code{numeric}
#' @param lMatInd \code{numeric} indices of selected features
#' @param linkMatrixThreshold \code{matrix} that contains information of linked 
#'  features for given thresholds
#' @param similarityMatrix \code{matrix} that is used to get information on the 
#' degree of similarity, \code{similarityMat} is an ordered version of a 
#' similarity matrix, see \code{?createOrderedSimMat}
#' @param roundDigits \code{numeric},  how many digits should be displayed?
#' @details \code{printInformationSelect} is for internal use. 
#' @return \code{character} that is in HTML format
#' @examples
#' data("idMSMStoMSP", package = "MetCirc")
#' data("binnedMSP", package = "MetCirc")
#' ## use only a selection
#' binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' groupname <- rownames(similarityMat)
#' ## order similarityMat according to mz
#' simMat <- createOrderedSimMat(similarityMat, order = "mz") 
#' groupnameMZ <- rownames(simMat)
#' linkMat_thr <- createLinkMatrix(simMat, 0.8, 1) 
#' ind <- 2
#' indMZ <- which(groupname[ind] == truncateName(groupnameMZ, NULL, group = TRUE))
#' linkMatInds <- getLinkMatrixIndices(groupnameMZ[indMZ], linkMat_thr)
#' MetCirc:::printInformationSelect(groupname = groupname, 
#'  msp = NULL, ind = ind, lMatInd = linkMatInds, 
#'  linkMatrixThreshold = linkMat_thr, 
#'  similarityMatrix = similarityMat, roundDigits = 2)
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
printInformationSelect <- function(groupname, msp = NULL, 
                ind, lMatInd, linkMatrixThreshold, similarityMatrix, roundDigits = 2) {
    
    ## truncate name of linkMatrixThreshold
    lMatThr <- linkMatrixThreshold 
    if (is.null(msp)) {
        ## selected feature
        selectedFeat <- groupname[ind]
        
        

        
        if (length(lMatInd) >= 1) {
            ## get connected features
            connect <- unique(as.vector(lMatThr[lMatInd, c("name1", "name2")]))
            selectedFeat <- truncateName(selectedFeat, roundDigits, TRUE)
            connect <- truncateName(connect, roundDigits, TRUE)
            ## remove selectedFeat from connect
            if (selectedFeat %in% connect) 
                connect <- connect[-which(connect == selectedFeat)]
            connect <- paste(connect, collapse = " <br/>")
            
            return(paste(c(selectedFeat, "connects to", "<br/>", connect), 
                            collapse = " "))
        } else 
            return(paste(c(selectedFeat, "does not connect to any feature"), 
                            collapse = " "))
            
    
    } else { ## if !is.null(msp)
            
        ## clicked feature: create msp and get identifier
        selectFeatMSP <- msp[ind]
        selectFeat <- groupname[ind] 

        if (length(lMatInd) == 0) {
            selectFeat <- truncateName(selectFeat, roundDigits = roundDigits, 
                                       group = TRUE)
            return(paste0(selectFeat, " (", selectFeatMSP@names, ", ", 
                selectFeatMSP@information, ", ", 
                selectFeatMSP@classes, ",", selectFeatMSP@adduct,  ") ",
                 "does not connect to any feature"))
        } else {
            ## connected features: find
            connect <- unique(as.vector(lMatThr[lMatInd, c("name1", "name2")]))
            connect <- truncateName(connect, NULL, TRUE)
            ## remove duplicated hovFeat in connect
            if (selectFeat %in% connect) connect <- connect[-which(connect == selectFeat)]
            
            matchedInd <- match(connect, groupname)
            connectMSP <- msp[matchedInd]
            connChar <- character()
            degreeSimilarity <- similarityMatrix[selectFeat, ]
            
            for (i in 1:length(connect)) {
                connectMSPI <- connectMSP[i]
                connectI <- connect[i]
                degreeSimilarityI <- round(degreeSimilarity[connectI],3)
                degreeSimilarityI <- as.numeric(degreeSimilarityI)
                connectI <- truncateName(connectI, roundDigits = roundDigits, 
                                         group = TRUE)

                newFeat <- paste0(connectI, " (", degreeSimilarityI, ", ", 
                    connectMSPI@names, ", ", connectMSPI@information, ", ", 
                    connectMSPI@classes, ", ", connectMSPI@adduct, ")", "<br/>")
                    
                connChar <- c(connChar, newFeat)
            }
            
            connChar <- paste(connChar, collapse=" ")
            selectFeat <- truncateName(selectFeat, roundDigits = roundDigits, 
                                       group = TRUE)
            
            return(paste0(selectFeat, " (", selectFeatMSP@names, ", ", 
                selectFeatMSP@information, ", ", 
                selectFeatMSP@classes, ", ", selectFeatMSP@adduct, ") connects to ", 
                " <br/>", connChar))
        }
    }
}