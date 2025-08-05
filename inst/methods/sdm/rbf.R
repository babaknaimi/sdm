# Author: Babak Naimi, naimi.b@gmail.com
# Date :  April 2018
# last update: April 2018
# Version 1.0
# Licence GPL v3

#-------------
methodInfo <- list(name=c('rbf','RBF','nnet.rbf','nnetRbf'),
                   packages='RSNNS',
                   modelTypes = c('pa','pb','ab','n'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame',v='sdmVariables'),
                   fitSettings = list(size=10,maxit=500),
                   fitFunction = function(formula,data,v,...) {
                     x <- .getData.sdmMatrix(formula,data,normalize=TRUE,frame=v@varInfo$numeric,scale=TRUE)
                     y <- .getData.sdmY(formula,data)
                     rbf(x=x,y=y,...)
                   },
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object='model',formula='standard.formula',newx='sdmDataFrame',v='sdmVariables'),
                   predictSettings=NULL,
                   predictFunction=function(object,formula,newx,v) {
                     newx <- .getData.sdmMatrix(formula,newx,normalize=TRUE,frame=v@varInfo$numeric,scale=TRUE)
                     .p <- predict(object,newx)[,1]
                     #-----
                     # If the range of predicted values is not in [0, 1], they are modified:
                     .r <- range(.p,na.rm = TRUE)
                     if (.r[1] < 0) .p <- .p + abs(.r[1])
                     
                     if (max(.p,na.rm = TRUE) > 1) {
                       .p <- .p / max(.p,na.rm = TRUE)
                     }
                     #-------
                     .p
                     
                   },
                   #------ metadata (optional):
                   title='radial basis function (RBF) network',
                   creator='Babak Naimi',
                   authors=c("Christoph Bergmeir et al. (for the package RSNNS)"), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Manual',title = "SNNS Stuttgart Neural Network Simulator User Manual,",
                                          author = as.person("A. Zell, [aut]"),
                                          year='1998',
                                          publisher="University of Stuttgart and WSI, University of Tübingen"
                   )
                   ),
                   description="The RBF is a type of Neural Network, performs a linear combination of n basis functions that are radially symmetric around a center/prototype."
)