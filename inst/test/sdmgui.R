# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Sep 2015
# Version 1.0
# Licence GPL v3

if (!isGeneric("gui")) {
  setGeneric("gui", function(...)
    standardGeneric("gui"))
}


setMethod('gui', signature('ANY'), 
          function(...) {
            l <- c(require(shiny),require(shinythemes),require(shinyBS),require(leaflet))
            if (!all(l)) stop(paste(paste(c('shiny','shinythemes','shinyBS','leaflet')[l],collapse=', '),'are required but not installed!'))
            
            #shiny::runApp(system.file('shinyApps/sdm_base', package='sdm'))
            shinyApp(
              ui=fluidPage(theme=shinytheme("Flatly"),
                           
                           #headerPanel(img(src="r_terminal_sdm.png", height = 300, width = 900)),
                           #br(),
                           
                           HTML ('<div style="background-color:black; color:white; margin:10px; padding:10px;">
                                 <h2>> library(sdm)</h2>
                                 <p>
                                 A framework to model and simulate species distributions
                                 </p>
                                 </div>'
                           ),
                           #     
                           
                           
                           sidebarPanel(
                             bsCollapsePanel("Data Settings",
                                             radioButtons('format','Input Train Data Format',c('CSV','SHP')),
                                             selectizeInput("train","Select Train Dataset:","NULL"),
                                             hr(),
                                             selectizeInput("predictors","Select Raster Predictors file","NULL"),
                                             br(),
                                             br()
                                             
                             ),
                             # shinyFilesButton('files', label='File select', title='Please select a file', multiple=FALSE),
                             bsCollapsePanel("Models Settings",
                                             checkboxGroupInput("methods","Methods",c("GLM","GAM","BRT","RF","SVM","MARS","NNet","Bioclim","Domin","Mahalanobis","Maxent","Maxlike","Ensemble"),inline=T, selected="GLM"),
                                             hr(),
                                             checkboxGroupInput("repMethod","Replication Methods",c("Subsampling","Bootsraping","Cross-Validation"),selected="Subsampling"),
                                             numericInput("replicates","Replicates",1,min=1),
                                             numericInput("cv.folds","CV.Folds",5,min=1),
                                             sliderInput('test.percent','Test persentage',0,100,30,1)
                             ),
                             bsToggleButton("run","Run!")
                             
                             
                           ),
                           mainPanel(
                             h3(textOutput("caption")),
                             tabsetPanel(position='above',
                                         tabPanel("Plot",plotOutput('mpgPlot')),
                                         tabPanel("Summary", verbatimTextOutput("summary")),
                                         tabPanel("Map", leafletOutput("mymap"),
                                                  p(),
                                                  actionButton("recalc", "New points"))
                             )
                           )
                           ),
              server=function(input, output,session) {
                observe({
                  if (input$format == 'CSV') {
                    updateSelectizeInput(session,"train",choices=list.files(pattern='.csv$'))
                  } else if (input$format == 'SHP') {
                    updateSelectizeInput(session,"train",choices=list.files(pattern='.shp$'))
                  }
                  updateSelectizeInput(session,"predictors",choices=list.files(pattern='grd$'),selected='')
                  #if (input$run != 0) updateButton(session,'run',value=0)
                  
                })
                train <- reactive({
                  #     if (is.null(input$train)) retrun(NULL)
                  #     else {
                  #       ext <- strsplit(input$train,'[.]')[[1]]
                  #       ext <- ext[length(ext)]
                  #       if (ext == 'csv') {
                  #         return(read.csv(input$train))
                  #       } else if (ext == 'shp') {
                  #         return(readShapeSpatial(input$train))
                  #       } else return (NULL)
                  #     }
                  if (input$format == 'CSV') {
                    return(read.csv(input$train))
                  } else if (input$format == 'SHP') {
                    return(readShapeSpatial(input$train))
                  }
                  
                })
                preds <- reactive({
                  if (input$predictors != '') {
                    return(brick(input$predictors))
                  } else return (NULL)
                })
                
                method <- reactive({
                  paste(input$methods,collapse=', ')
                })
                #shinyFileChoose(input, 'files', session=session,roots=c(wd='.'), filetypes=c('', '.txt'))
                ss <- reactive({
                  if (!input$run) return (NULL)
                  if (!is.null(train)) {
                    if (is.null(preds())) d <- sdmData(train())
                    else d <- sdmData(train=train(),predictors = preds())
                  }
                  f <- as.formula(paste(d@train@Occurrence@species.name,'~.'))
                  sdm(f,d,methods = input$methods,test.perc=input$test.percent,replicate.method=input$repMethod,replicates=as.numeric(input$replicates),cv.folds=as.numeric(input$cv.folds))
                  
                  
                })
                output$caption <- renderText({
                  if (!is.null(ss())) method()
                })
                
                output$summary <- renderPrint({
                  if (!is.null(ss())) show(ss())
                })
                
                
                output$mpgPlot <- renderPlot({
                  #s <- sdm(Occurrence~.,d,methods = input$methods,test.perc=input$test.percent,replicate.method=input$repMethod)
                  #roc(s)
                  if (!is.null(ss())) roc(ss())
                })
                
                output$mymap <- renderLeaflet({
                  leaflet() %>%
                    addTiles() %>%
                    #addMarkers(data = points())
                    addMarkers(lng=174.768, lat=-36.852, popup="The birthplace of R")
                })
                
              }
              
              )
          }
)
