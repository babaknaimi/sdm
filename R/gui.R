# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2016
# Version 2.7
# Licence GPL v3



.wellP <- function(...,style="background:#27242E;color:#fff") {
  shiny::tags$div(class='well',style=style,...)
}

.gui.sdm <- function(x) {
  
  .css <- readRDS(system.file("shinyApps/css.rds", package="sdm"))
  
  wt <- sapply(x@run.info[,7:9],function(x) any(x))
  wt <- names(wt)[wt]
  
  wtSel <- c()
  
  thStats <- names(getEvaluation(x,w=x@run.info[which(x@run.info[,7])[1],1],stat=2,wtest='training'))
  thStats <- thStats[3:length(thStats)]
  
  ui <- shiny::fluidPage(
    shiny::tags$head(.css),
    
    shiny::tabsetPanel(position='above',
                       shiny::tabPanel("Summary", .wellP(shiny::verbatimTextOutput("summary"))),
                       shiny::tabPanel("Model Runs Details", .wellP(shiny::tableOutput('run_info'))),
                       shiny::tabPanel("Evaluation",
                                       .wellP(
                                         shiny::navbarPage('Evaluation results',
                                                           shiny::tabPanel('Plot',
                                                                           shiny::fluidRow(
                                                                             .wellP(
                                                                               shinyBS::bsCollapse(
                                                                                 shinyBS::bsCollapsePanel("ROC-AUC",
                                                                                                          .wellP(
                                                                                                            shiny::fluidRow(
                                                                                                              shiny::column(3,shiny::selectInput("species", label = "Species Names",choices = as.character(unique(x@run.info$species)), selected = 1)),
                                                                                                              shiny::column(3,shiny::selectizeInput("models","Modeling Methods",choices=unique(as.character(x@run.info$method)),multiple=TRUE)),
                                                                                                              shiny::column(3,shiny::selectizeInput("replications","Replications",choices=as.character(unique(x@run.info$replication)),multiple=TRUE)),
                                                                                                              shiny::column(3,shiny::selectizeInput("modelID", label = "Model IDs", choices = as.character(x@run.info$modelID), multiple = TRUE))
                                                                                ),
                                                                                shiny::fluidRow(
                                                                                  shiny::column(6,shiny::checkboxGroupInput('wtest','Which Test Dataset?',choices=wt,selected = if(length(wt) > 2) wt[1:2] else wt)),
                                                                                  shiny::column(6,shiny::radioButtons('smooth','Should Be Smoothed?',choices=c('TRUE','FALSE'),selected ='TRUE'))
                                                                                ),
                                                                                shiny::fluidRow(.wellP(shiny::plotOutput('evalPlot1')))
                                                                              ),style='primary'
                                                     ),
                                                     shinyBS::bsCollapsePanel("Calibration",
                                                                              .wellP(
                                                                                shiny::fluidRow(
                                                                                  shiny::column(3,shiny::selectInput("modelID2", label = "Model IDs", choices = as.character(x@run.info$modelID), selectize = TRUE)),
                                                                                  shiny::column(3,shiny::textInput("species_name", label = "Species Names")),
                                                                                  shiny::column(3,shiny::textInput("model_name", label = "Method")),
                                                                                  shiny::column(3,shiny::textInput("replication_name", label = "Replication"))
                                                                                ),
                                                                                shiny::fluidRow(
                                                                                  shiny::column(3,shiny::selectInput('wtest2','Which Test Dataset?',wt,selected=wt[length(wt)]))
                                                                                ),
                                                                                shiny::fluidRow(.wellP(shiny::plotOutput('evalPlot2')))
                                                                              ),style='primary'
                                                                              
                                                     ),
                                                     
                                                     shinyBS::bsCollapsePanel("Threshold Optimisations",
                                                                              .wellP(
                                                                                shiny::fluidRow(
                                                                                  shiny::column(3,shiny::selectInput("modelID3", label = "Model IDs", choices = as.character(x@run.info$modelID), selectize = TRUE)),
                                                                                  shiny::column(3,shiny::textInput("species_name2", label = "Species Names")),
                                                                                  shiny::column(3,shiny::textInput("model_name2", label = "Method")),
                                                                                  shiny::column(3,shiny::textInput("replication_name2", label = "Replication"))
                                                                                ),
                                                                                shiny::fluidRow(
                                                                                  shiny::column(3,shiny::selectInput('wtest3','Which Test Dataset?',wt,selected=wt[length(wt)])),
                                                                                  shiny::column(6,shiny::selectInput('thStat','Select Statistics',thStats,selected='TSS'))
                                                                                ),
                                                                                shiny::fluidRow(.wellP(shiny::plotOutput('evalPlot3')))
                                                                              ),style='primary'
                                                                              
                                                                              
                                                     ),
                                                     shinyBS::bsCollapsePanel("Density",
                                                                              .wellP(
                                                                                shiny::fluidRow(
                                                                                  shiny::column(3,shiny::selectInput("modelID4", label = "Model IDs", choices = as.character(x@run.info$modelID), selectize = TRUE)),
                                                                                  shiny::column(3,shiny::textInput("species_name3", label = "Species Names")),
                                                                                  shiny::column(3,shiny::textInput("model_name3", label = "Method")),
                                                                                  shiny::column(3,shiny::textInput("replication_name3", label = "Replication"))
                                                                                ),
                                                                                shiny::fluidRow(
                                                                                  shiny::column(3,shiny::selectInput('wtest4','Which Test Dataset?',wt,selected=wt[length(wt)]))
                                                                                ),
                                                                                shiny::fluidRow(.wellP(shiny::plotOutput('evalPlot4')))
                                                                              ),style='primary'
                                                     ),
                                                     shinyBS::bsCollapsePanel("Boxplot",
                                                                              .wellP(
                                                                                shiny::fluidRow(
                                                                                  shiny::column(3,shiny::selectInput("modelID5", label = "Model IDs", choices = as.character(x@run.info$modelID), selectize = TRUE)),
                                                                                  shiny::column(3,shiny::textInput("species_name4", label = "Species Names")),
                                                                                  shiny::column(3,shiny::textInput("model_name4", label = "Method")),
                                                                                  shiny::column(3,shiny::textInput("replication_name4", label = "Replication"))
                                                                                ),
                                                                                shiny::fluidRow(
                                                                                  shiny::column(3,shiny::selectInput('wtest5','Which Test Dataset?',wt,selected=wt[length(wt)]))
                                                                                ),
                                                                                shiny::fluidRow(.wellP(shiny::plotOutput('evalPlot5')))
                                                                                
                                                                              ),style='primary'
                                                     ), id='plotType'
                                                   )
                                                   
                                                 )
                                               ),icon=shiny::icon('line-chart',"fa-2x")
                                      ),
                                      shiny::tabPanel('threshold-based',
                                               .wellP(
                                                 shiny::fluidRow(
                                                   shiny::column(3,shiny::selectInput("modelID6", label = "Model IDs", choices = as.character(x@run.info$modelID), selectize = TRUE)),
                                                   shiny::column(3,shiny::textInput("species_name5", label = "Species Names")),
                                                   shiny::column(3,shiny::textInput("model_name5", label = "Method")),
                                                   shiny::column(3,shiny::textInput("replication_name5", label = "Replication"))
                                                 ),
                                                 shiny::fluidRow(
                                                   shiny::column(3,shiny::selectInput('wtest6','Which Test Dataset?',wt,selected=wt[length(wt)]))
                                                 )
                                               ),
                                               .wellP(shiny::tableOutput('thTable')),
                                               icon=shiny::icon('arrow-circle-down',"fa-2x")
                                      ),
                                      shiny::tabPanel('threshold-independent',
                                               .wellP(
                                                 shiny::fluidRow(
                                                   shiny::column(3,shiny::selectInput("modelID7", label = "Model IDs", choices = as.character(x@run.info$modelID), selectize = TRUE)),
                                                   shiny::column(3,shiny::textInput("species_name6", label = "Species Names")),
                                                   shiny::column(3,shiny::textInput("model_name6", label = "Method")),
                                                   shiny::column(3,shiny::textInput("replication_name6", label = "Replication"))
                                                 ),
                                                 shiny::fluidRow(
                                                   shiny::column(3,shiny::selectInput('wtest7','Which Test Dataset?',wt,selected=wt[length(wt)]))
                                                 )
                                               ),
                                               .wellP(shiny::tableOutput('thTable2')),icon=shiny::icon('arrows-h',"fa-2x")
                                      )
                                      
                                      ,id='Statistic'
                                      
                           )
                         )),
                       shiny::tabPanel("Variable Importance",
                                       .wellP(
                                         shiny::fluidRow(
                                           shiny::column(3,shiny::selectInput("modelID8", label = "Model IDs", choices = as.character(x@run.info$modelID), selectize = TRUE)),
                                           shiny::column(3,shiny::textInput("species_name7", label = "Species Names")),
                                           shiny::column(3,shiny::textInput("model_name7", label = "Method")),
                                           shiny::column(3,shiny::textInput("replication_name7", label = "Replication"))
                                         ),
                                         shiny::fluidRow(
                                           shiny::column(3,shiny::selectInput('wtest8','Which Test Dataset?',wt,selected=wt[length(wt)])),
                                           shiny::column(3,shiny::selectInput('varStat','VarImportance Metric',c('Correlation test','AUC test'),selected=1))
                                         ),
                                         shiny::fluidRow(.wellP(shiny::plotOutput('evalPlot6')))
                                       )
                       ),id='main_tabs'
    )
  )
  
  #------------------------
  server <- function(input,output,session) {
    
    single_model_info <- shiny::reactive({
      w <- NULL
      if (!is.null(input$plotType)) {
        if (input$plotType == 'Calibration') {
          if (!is.null(input$modelID2)) w <- which(x@run.info$modelID == input$modelID2)
        } else if (input$plotType == 'Threshold Optimisations') {
          if (!is.null(input$modelID3)) w <- which(x@run.info$modelID == input$modelID3)
        } else if (input$plotType == 'Density') {
          if (!is.null(input$modelID4)) w <- which(x@run.info$modelID == input$modelID4)
        } else if (input$plotType == 'Boxplot') {
          if (!is.null(input$modelID5)) w <- which(x@run.info$modelID == input$modelID5)
        } 
      } 
      if (!is.null(w)) c(species=as.character(x@run.info[w,"species"]), method=as.character(x@run.info[w,"method"]),replication=as.character(x@run.info[w,"replication"]))
      
    })
    
    
    single_model_info2 <- shiny::reactive({
      w <- NULL
      if (input$main_tabs == "Variable Importance") {
        if (!is.null(input$modelID8)) w <- which(x@run.info$modelID == input$modelID8)
      }
      if (!is.null(w)) c(species=as.character(x@run.info[w,"species"]), method=as.character(x@run.info[w,"method"]),replication=as.character(x@run.info[w,"replication"]))
      
    })
    
    modelID <- shiny::reactive({
      getModelInfo(x,species=input$species,method=input$models,replication=input$replications)$modelID
    })
    
    shiny::observe({
      
      output$summary <- shiny::renderPrint({
        show(x)
      })
      
      
      
      if (input$Statistic == 'Plot') {
        
        if (!is.null(input$plotType)) {
          
          if (input$plotType == 'ROC-AUC') {
            
            shiny::updateSelectizeInput(session,"modelID",choices=modelID())
            
            output$evalPlot1 <- shiny::renderPlot({
              if (!is.null(input$modelID)) {
                roc(x,p=input$modelID,wtest=input$wtest,smooth=as.logical(input$smooth))
              } else {
                roc(x,p=modelID(),wtest=input$wtest,smooth=as.logical(input$smooth))
              } 
            })
            
          } else if (input$plotType == 'Calibration') {
            shiny::updateTextInput(session,'species_name',value=as.character(single_model_info()[1]))
            shiny::updateTextInput(session,'model_name',value=as.character(single_model_info()[2]))
            shiny::updateTextInput(session,'replication_name',value=as.character(single_model_info()[3]))
            shiny::updateSelectInput(session,'modelID2',selected = input$modelID2)
            
            output$evalPlot2 <- shiny::renderPlot({
              sinfo <- as.character(single_model_info())
              if (length(sinfo) != 0) {
                ev <- x@models[[sinfo[1]]][[sinfo[2]]][[as.character(input$modelID2)]]@evaluation[[as.character(input$wtest2)]]
                plot(calibration(ev))
              }
              
            })
            
          } else if (input$plotType == 'Threshold Optimisations') {
            
            
            shiny::updateTextInput(session,'species_name2',value=as.character(single_model_info()[1]))
            shiny::updateTextInput(session,'model_name2',value=as.character(single_model_info()[2]))
            shiny::updateTextInput(session,'replication_name2',value=as.character(single_model_info()[3]))
            shiny::updateSelectInput(session,'modelID3',selected = input$modelID3)
            
            
            output$evalPlot3 <- shiny::renderPlot({
              sinfo <- as.character(single_model_info())
              if (length(sinfo) != 0) plot(x@models[[sinfo[1]]][[sinfo[2]]][[as.character(input$modelID3)]]@evaluation[[as.character(input$wtest3)]],as.character(input$thStat))
            })
            
          } else if (input$plotType == 'Density') {
            shiny::updateTextInput(session,'species_name3',value=as.character(single_model_info()[1]))
            shiny::updateTextInput(session,'model_name3',value=as.character(single_model_info()[2]))
            shiny::updateTextInput(session,'replication_name3',value=as.character(single_model_info()[3]))
            shiny::updateSelectInput(session,'modelID4',selected = input$modelID4)
            
            output$evalPlot4 <- shiny::renderPlot({
              sinfo <- as.character(single_model_info())
              if (length(sinfo) != 0) {
                ev <- x@models[[sinfo[1]]][[sinfo[2]]][[as.character(input$modelID4)]]@evaluation[[as.character(input$wtest4)]]
                density(ev)
              }
              
            })
            
          } else if (input$plotType == 'Boxplot') {
            shiny::updateTextInput(session,'species_name4',value=as.character(single_model_info()[1]))
            shiny::updateTextInput(session,'model_name4',value=as.character(single_model_info()[2]))
            shiny::updateTextInput(session,'replication_name4',value=as.character(single_model_info()[3]))
            shiny::updateSelectInput(session,'modelID5',selected = input$modelID5)
            
            
            output$evalPlot5 <- shiny::renderPlot({
              sinfo <- as.character(single_model_info())
              if (length(sinfo) != 0) {
                ev <- x@models[[sinfo[1]]][[sinfo[2]]][[as.character(input$modelID5)]]@evaluation[[as.character(input$wtest5)]]
                boxplot(ev)
              }
              
            })
            
          }
        }
        
      } else if (input$Statistic == 'threshold-based') {
        if (!is.null(input$modelID6)) {
          w <- which(x@run.info$modelID == input$modelID6)
          sinfo <- as.character(c(species=as.character(x@run.info[w,"species"]), method=as.character(x@run.info[w,"method"]),replication=as.character(x@run.info[w,"replication"])))
          shiny::updateTextInput(session,'species_name5',value=as.character(sinfo[1]))
          shiny::updateTextInput(session,'model_name5',value=as.character(sinfo[2]))
          shiny::updateTextInput(session,'replication_name5',value=as.character(sinfo[3]))
          shiny::updateSelectInput(session,'modelID6',selected = input$modelID6)
          
          output$thTable <- shiny::renderTable({
            getEvaluation(x,w=as.numeric(input$modelID6),stat=2,wtest=as.character(input$wtest6))
          })
        }
        
      } else if (input$Statistic == 'threshold-independent') {
        if (!is.null(input$modelID7)) {
          w <- which(x@run.info$modelID == input$modelID7)
          sinfo <- as.character(c(species=as.character(x@run.info[w,"species"]), method=as.character(x@run.info[w,"method"]),replication=as.character(x@run.info[w,"replication"])))
          shiny::updateTextInput(session,'species_name6',value=as.character(sinfo[1]))
          shiny::updateTextInput(session,'model_name6',value=as.character(sinfo[2]))
          shiny::updateTextInput(session,'replication_name6',value=as.character(sinfo[3]))
          shiny::updateSelectInput(session,'modelID7',selected = input$modelID7)
          
          s1 <- getEvaluation(x,w=as.numeric(input$modelID7),stat=1,wtest=as.character(input$wtest7))
          
          s2 <- NULL
          
          if (length(sinfo) != 0) {
            ev <- x@models[[sinfo[1]]][[sinfo[2]]][[as.character(input$modelID7)]]@evaluation[[as.character(input$wtest7)]]
            s2 <- calibration(ev)@statistic
          }
          df <- data.frame(matrix(nrow=5,ncol=2))
          colnames(df) <- c('Statistic','Value')
          df[,1] <- c(names(s1),'Calibrartion')
          df[,2] <- c(s1[[1]],s1[[2]],s1[[3]][1],s1[[4]],s2)
          output$thTable2 <- shiny::renderTable({
            df
          })
        }
        
      }
      
      
      if (input$main_tabs == "Variable Importance") {
        shiny::updateTextInput(session,'species_name7',value=as.character(single_model_info2()[1]))
        shiny::updateTextInput(session,'model_name7',value=as.character(single_model_info2()[2]))
        shiny::updateTextInput(session,'replication_name7',value=as.character(single_model_info2()[3]))
        shiny::updateSelectInput(session,'modelID8',selected = input$modelID8)
        
        
        output$evalPlot6 <- shiny::renderPlot({
          sinfo <- as.character(single_model_info2())
          if (length(sinfo) != 0) {
            if (input$varStat == 'Correlation test') vst <- 'corTest'
            else vst <- 'AUCtest'
            plot(getVarImp(x,id=as.numeric(input$modelID8),wtest = input$wtest8),vst)
          }
          
        })
        
      }
      
    })
    
    
    
    
    
    shiny::observeEvent(input$wtest,{
      if (length(input$wtest) > 2) {
        ad <- input$wtest[!input$wtest %in% wtSel]
        wtSel <<- wtSel[-1]
        wtSel <<- c(wtSel,ad)
        shiny::updateCheckboxGroupInput(session,'wtest',selected=wtSel)
        
      } else {
        if (length(input$wtest) == 1) wtSel <<- c(wtSel,as.character(input$wtest))
        else wtSel <<- c(wtSel,input$wtest[!input$wtest %in% wtSel])
      }
      
    })
    
    output$run_info <- shiny::renderTable({
      x@run.info
    })
  }
  
  shiny::runApp(shiny::shinyApp(ui=ui,server=server))
}

if (!isGeneric("gui")) {
  setGeneric("gui", function(x,...)
    standardGeneric("gui"))
}


setMethod('gui', signature('sdmModels'), 
          function(x,...) {
            l <- c(.loadLib(c('shiny','shinyBS')))
            if (!all(l)) stop(paste(paste(c('shiny','shinyBS')[l],collapse=', '),'are not installed; install them or use installAll() function to install all the functions that may be required by some functions in the package...'))
            
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            .gui.sdm(x)
          }
)
