.trainCheckVars <- function(train) {
  nnFact <- nnf <- nnxy <- nnt <- NULL
  w <- sdm:::.speciesDetect(train)
  if (!is.null(w$nFact)) {
    nFact <- w$nFact
    nnFact <- paste(paste('f(',w$nFact,')',sep=''),collapse='+')
  }
  if (!is.null(w$nf)) {
    nf <- w$nf
    nnf <- paste(w$nf,collapse='+')
  }
  if (!is.null(w$nxy)) {
    nxy <- w$nxy
    nnxy <- paste(paste('coords(',paste(w$nxy,collapse='+'),')',sep=''),collapse='+')
  }
  if (!is.null(w$nt)) {
    nt <- w$nt
    nnt <- paste(paste('time(',w$nt,')',sep=''),collapse='+')
  }
  formula <- as.formula(paste(paste(w$nsp,collapse="+"),'~',paste(c(nnf,nnFact,nnxy,nnt),collapse='+')),env = parent.frame())
  c(formula=formula,w)
}


#.trainCheckVars(data.frame(sp=c(1,1,1,0,0),a1=c(1,2,3,4,5),b=c('a','b','c','d','e')))



server <- function(input, output,session) {
  #----- ShiyFiles--------
  volumes <- c(wd='.',getVolumes()())
  
  shinyFileChoose(input, 'trainFile',session=session,roots=volumes, filetypes=c('', 'csv','txt','shp','sdd','rds'))
  shinyFileChoose(input, 'predFiles',session=session,roots=volumes, filetypes=c('','grd', 'tif','img','bil','asc','rst','nc','envi'))
  shinyFileChoose(input, 'testFile',session=session,roots=volumes, filetypes=c('', 'csv','txt','shp'))
  
  observeEvent(input$trainFile$files, {
    trainName <- as.character(parseFilePaths(volumes,input$trainFile)[1,4])
    updateTextInput(session,"trainFileName",value=trainName)
    ext <- strsplit(trainName,'[.]')[[1]]
    ext <- ext[length(ext)]
    if (ext == 'csv') {
      train <- read.csv(trainName)
    } else if (ext == 'shp') {
      train <- shapefile(trainName)
    } else if (ext == 'txt') {
      train <- read.table(trainName,header = TRUE,sep='\t')
    } else if (ext == 'sdd') {
      train <- read.sdm(trainName)
    } else if (ext == 'rds') {
      train <- readRDS(trainName)
    }
    if (class(train) == 'data.frame') {
      # output$dynSpecies <- renderUI({
      #   checkboxGroupInput("chk","Select species",colnames(train))
      # })
      # output$dynPreds <- renderUI({
      #   checkboxGroupInput("chk","Select predictors",colnames(train))
      # })
      output$dynSpecies <- renderUI({
        # If missing input, return to avoid error later in function
        # Create the checkboxes and select them all by default
        checkboxGroupInput("columns", "Choose columns", 
                           choices  = names(train),
                           selected = names(train))
      })
      
      # output$dynPreds <- renderTable({
      #   # If missing input, return to avoid error later in function
      #   # Keep the selected columns
      #   train <- train[, input$columns, drop = FALSE]
      #   
      #   # Return first 20 rows
      #   head(train, 20)
      # })
      # 
      output$tbl <- DT::renderDataTable({
        df <- data.frame(
          varNames=colnames(train),
          Linear=rep(TRUE,length(colnames(train))),
          Quadratic=rep(FALSE,length(colnames(train))),
          Cubic=rep(FALSE,length(colnames(train)))
        )
        addCheckboxButtons <- paste0('<input type="checkbox" name="row', df$Linear, '" value="', df$Linear, '" checked="1">',"")
        
        DT::datatable(cbind(select=addCheckboxButtons, df), options = list(searching=FALSE, paging=FALSE, ordering=FALSE,
                                                                           columnDefs=list(list(className='dt-center', targets=c(0,2)))),
                      callback = "function(table) {
                      table.on('change.dt', 'tr td input:checkbox', function() {
                      setTimeout(function () {
                      Shiny.onInputChange('varRows', $(this).add('tr td input:checkbox:checked').parent().siblings(':nth-child(2)').map(function() {
                      return $(this).text();
                      }).get())
                      }, 10);
                      }); 
      }", escape=FALSE) 
    })
      
    } else if (class(train) == 'SpatialPointsDataFrame') {
      output$dynSpecies <- renderUI({
        checkboxGroupInput("chk","Select species",colnames(train))
      })
    }
  }
  )
  
  observeEvent(input$predFiles$files, {
    updateTextInput(session,"predFileNames",value=as.character(parseFilePaths(volumes,input$predFiles)[,4]))
  }
  )
  
  observeEvent(input$trainFile$files, {
    updateTextInput(session,"testFileName",value=as.character(parseFilePaths(volumes,input$testFile)[1,4]))
  }
  )
  #---------------
  #################
  
  
  
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
    #if (!is.null(ss())) roc(ss())
  })
  
  output$mymap <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      #addMarkers(data = points())
      addMarkers(lng=174.768, lat=-36.852, popup="The birthplace of R")
  })
  
}