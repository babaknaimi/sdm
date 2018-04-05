library(shiny)
library(shinyFiles)
library(DT)
library(leaflet)
ui <- fluidPage(
  tags$head(
   tags$style(HTML('.container { background-color: #fff; border:3px solid #ccc; width:300px; height: 200px; overflow-y: scroll; }
                   .checkboxes {
    columns: 4 8em;
                   }'))
  ),
  
  img(src="sdm_design.png",height=300,width="100%"),
  fluidRow(
    .collapsePanel("Data input",
                   column(4,
                          wellPanel(
                            tags$p('Select the training file:'),
                            textInput('trainFileName','Training'),
                            shinyFilesButton("trainFile","File select","training file",FALSE)
                          )
                          ),
                   column(4,
                          wellPanel(
                            tags$p('Select the predictor file(s) (optional):'),
                            textInput('predFileNames','Predictors'),
                            shinyFilesButton("predFiles","File(s) select","predictor files",TRUE)
                          )
                   ),
                   column(4,
                          wellPanel(
                            tags$p('Select the test file (optional):'),
                            textInput('testFileName','Test'),
                            shinyFilesButton("testFile","File select","test file",FALSE)
                          )
                   )
                  )
   
    
  ),
  fluidRow(
    .collapsePanel("Data Settings",
                   column(4,wellPanel(style = "background-color: #ffffff;",tags$div(uiOutput("dynSpecies"))))
    ),
    #column(5,wellPanel(style = "background-color: #ffffff;",tags$div(uiOutput("dynCheck"),class='container')))
    column(5,wellPanel(style = "background-color: #003333; color: #ffffff",tags$div(uiOutput("dynPreds")))),
    column(5,wellPanel(style = "background-color: #003333; color: #ffffff",tags$div(DT::dataTableOutput('tbl'))))
    
    # column(5,wellPanel(HTML('<div class="container">
    # <input type="checkbox" /> This is checkbox <br />
    #                         <input type="checkbox" /> This is checkbox <br />
    #                         <input type="checkbox" /> This is checkbox <br />
    #                         <input type="checkbox" /> This is checkbox <br />
    #                         <input type="checkbox" /> This is checkbox <br />
    #                         <input type="checkbox" /> This is checkbox <br />
    #                         <input type="checkbox" /> This is checkbox <br />
    #                         <input type="checkbox" /> This is checkbox <br />
    #                         <input type="checkbox" /> This is checkbox <br />
    #                         <input type="checkbox" /> This is checkbox <br />
    #                         </div>')))
  ),
  fluidRow(HTML('<table>
<tbody>
                <!-- What should I use here to add a question with the instructions? -->
                <tr>
                <th>watches</th>
                <th>jewerly</th>
                <th>jewerly</th>
                </tr>
                <tr>
                <td>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>    
                <input type="checkbox" name="utType1[]" value="e" />Something </br> 
                </td>
                <td>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>    
                <input type="checkbox" name="utType1[]" value="e" />Something </br> 
                </td>
                <td>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>    
                <input type="checkbox" name="utType1[]" value="e" />Something </br> 
                </td>
                </tr>
                <!-- What should I use here to add a question with the instructions? -->
                <tr>
                <th>watches</th>
                <th>jewerly</th>
                <th>jewerly</th>
                </tr>
                <tr>
                <td>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>    
                <input type="checkbox" name="utType1[]" value="e" />Something </br> 
                </td>
                <td>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>    
                <input type="checkbox" name="utType1[]" value="e" />Something </br> 
                </td>
                <td>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>
                <input type="checkbox" name="utType1[]" value="e" />Something </br>    
                <input type="checkbox" name="utType1[]" value="e" />Something </br> 
                </td>
                </tr>
                </tbody>
                </table>')),
  sidebarPanel(
    radioButtons('format','Input Train Data Format',c('CSV','SHP')),
    hr(),
    selectizeInput("train","Select Train Dataset:","NULL"),
    selectizeInput("predictors","Select Raster Predictors file","NULL"),
    checkboxGroupInput("methods","Methods",c("GLM","GAM","BRT","RF","SVM","MARS","NNet","Bioclim","Domin","Mahalanobis","Maxent","Maxlike","Ensemble"),inline=T, selected="GLM"),
    hr(),
    checkboxGroupInput("repMethod","Replication Methods",c("Subsampling","Bootsraping","Cross-Validation"),selected="Subsampling"),
    numericInput("replicates","Replicates",1,min=1),
    numericInput("cv.folds","CV.Folds",5,min=1),
    sliderInput('test.percent','Test persentage',0,100,30,1),
    #bsCollapsePanel("test")
    #bsCollapsePanel("Data Settings",
    #                 radioButtons('format','Input Train Data Format',c('CSV','SHP')),
    #                 hr(),
    #                 selectizeInput("train","Select Train Dataset:","NULL"),
    #                 selectizeInput("predictors","Select Raster Predictors file","NULL"),
    #                 br(),
    #                 br()
    #                 
    # ),
    # bsCollapsePanel('test',
    #                 sliderInput('a','aa',min = 0,max = 10,1)
    #                 )
    #HTML(paste0(p1,p2,p3))
    #bsCollapsePanel('ttt')
    #shinyFilesButton('files', label='File select', title='Please select a file', multiple=FALSE)
    .collapsePanel("Models Settings",
                   checkboxGroupInput("methods1","Methods",c("GLM","GAM","BRT","RF","SVM","MARS","NNet","Bioclim","Domin","Mahalanobis","Maxent","Maxlike","Ensemble"),inline=T, selected="GLM"),
                   hr(),
                   checkboxGroupInput("repMethod1","Replication Methods",c("Subsampling","Bootsraping","Cross-Validation"),selected="Subsampling"),
                   numericInput("replicates1","Replicates",1,min=1),
                   numericInput("cv.folds1","CV.Folds1",5,min=1),
                   sliderInput('test.percent1','Test persentage',0,100,30,1)
    )
    # bsButton("run","Run!")
    # 
    
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
)
#library(shinyFiles)
