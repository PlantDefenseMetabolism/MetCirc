#' @name shinyCircos
#' @title Interactive visualisation of similar precursors
#' @description Visualise similar precursors.
#' @usage shinyCircos(dfNameGroup, similarityMatrix, msp)
#' @param dfNameGroup data.frame which contains columns "group" and "name"
#' @param similarityMatrix matrix, similarityMatrix contains pair-wise 
#' similarity coefficients which give information about the similarity between
#' precursors
#' @param msp data.frame, a data.frame in msp file format for information about the hovered feature
#' @details The function is based on the shiny and circlize package. Choose
#' interactively thresholds, type of links, hover over precursors, select 
#' precursors.
#' @value character
#' @return shinyCircos returns a character vector with the selected 
#' precursors
#' @author Thomas Naake, \email{naake@@stud.uni-heidelberg.de}
#' @examples \dontrun{shinyCircos(dfNameGroup, similarityMatrix, msp)}
#' @export
shinyCircos <- function(dfNameGroup, similarityMatrix, msp) {
    ## circlize parameters
    circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), track.margin = c(0.0, 0))
    
    plotCircos(dfNameGroup, NULL, initialize=TRUE, featureNames = TRUE, groupName = TRUE, links = FALSE, highlight = FALSE)
    PlotFilled <- recordPlot()
    plot.new()
    
    plotCircos(dfNameGroup, NULL, initialize=TRUE, featureNames = TRUE, groupName = TRUE, links = FALSE, highlight = TRUE)
    PlotHighlight <- recordPlot()
    plot.new()
    
    ## get degree of features
    features <- as.character( dfNameGroup[,"name"] )
    degreeFeatures <- lapply(features, function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
    
    ui <- fluidPage(
        sidebarPanel(
            radioButtons("choiceLinks", "choose type of links", 
                choices = c("all" = "all", "inter-class links" = "inter", 
                            "intra-class links" = "intra"),
                selected = "all"),
            sliderInput("threshold", "Threshold for similarity to display",
                min = 0, max = 1, value = 0.75),
            actionButton("resetClickIndices", "Reset features"),
            actionButton("stop", "Stop and export \n selected features")
        ),
        mainPanel(
            plotOutput("circos",
                click = "circosClick",
                #dblclick = "circosDblClick",
                hover = hoverOpts(id = "circosHover", delay = 400, clip = TRUE,nullOutside = FALSE)),
                #brush = brushOpts(id = "circosBrush",
                #                  resetOnNew = TRUE)),
            textOutput("hoverConnectedFeature"),
            verbatimTextOutput("clickFeature"))
    )
    server <- function(input, output, session) {
        
        ## create reactive expression for LinkMatrix
        LinkMatrix <- reactive(createLinkMatrix(similarityMatrix, input$threshold, dfNameGroup))
        ## create reactive expression for LinkMatrix which is cut according to 
        ## set radioButton (input$choiceLinks)
        LinkMatrix_cut <- reactive(cutLinkMatrix(LinkMatrix(), type = input$choiceLinks))
        
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
                if (.dist >= 0.8 & .dist <= 1) onCircle$is <- TRUE else onCircle$is <- FALSE
            } else onCircle$is <- FALSE
        })
        
        ## Hover: which is the current sector?
        indHover <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosHover$x))
                indHover$ind <- minFragCart2Polar(input$circosHover$x, 
                                                  input$circosHover$y, 
                                                  degreeFeatures) 
            
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
                                            degreeFeatures) 
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
        output$circos <- renderPlot({
                if (onCircle$is) {
                   # replayPlot(PlotFilled)
                    replayPlot(PlotHighlight)
                    highlight(dfNameGroup, c(indHover$ind, indClick$ind), LinkMatrix_cut()) }
            else { ## if not onCircle$is
                if (length(indClick$ind) > 0) {
                    replayPlot(PlotHighlight)
                    highlight(dfNameGroup, c(indClick$ind), LinkMatrix_cut())
                } else {
                    replayPlot(PlotFilled)
                    plotCircos(dfNameGroup, LinkMatrix_cut(), initialize=FALSE, featureNames = FALSE, groupName = FALSE, links = TRUE, highlight = FALSE)
                }
            }
            
        })
        
        ## show when hovering the feature which connects to it
        linkMatIndsHover <- reactive({getLinkMatrixIndices(dfNameGroup[indHover$ind,], LinkMatrix_cut())})
        
        output$hoverConnectedFeature <- renderText({
            if (onCircle$is) {
                if (length(linkMatIndsHover() > 0))
                    c(as.character(dfNameGroup[indHover$ind,"name"]), "connects to", 
                        unique(as.vector(LinkMatrix_cut()[linkMatIndsHover(), c("name1","name2")]))[-1])
                else 
                    c(as.character(dfNameGroup[indHover$ind,"name"]), "does not connect to any feature")
            }
        })
        
        output$clickFeature <- renderText({
            if (length(indClick$ind) > 0)
                c("selected features: ", as.character(dfNameGroup[indClick$ind, "name"]))
            else "no features selected"
        })
        
        ## on exit
        observe({
            if (input$stop == 0)
                return()
            else {
                circos.clear()
                stopApp(as.character(dfNameGroup[indClick$ind, "name"]))
            }
        })
        
    }
    
    app <- list(ui = ui, server = server)
    runApp(app)
    
}
## to do
##order within sectors according to rt or hierarchical clustering
## second simmat for neutral losses

##b <- shinyCircos(dfNameGroup, similarityMat, msp = finalMSP)
