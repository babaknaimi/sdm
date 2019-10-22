# Author: Babak Naimi, naimi.b@gmail.com
# Date :  April 2018
# last update: March 2019
# Version 1.1
# Licence GPL v3

#-------------
methodInfo <- list(name=c('mlp','MLP','nnet.mlp','nnetMLP'),
                   packages='RSNNS',
                   modelTypes = c('pa','pb','ab','n'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(size=10,maxit=500),
                   fitFunction = function(formula,data,...) {
                     x <- .getData.sdmMatrix(formula,data,normalize=TRUE)
                     y <- .getData.sdmY(formula,data)
                     mlp(x=x,y=y,...)
                   },
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object='model',formula='standard.formula',newx='sdmDataFrame'),
                   predictSettings=NULL,
                   predictFunction=function(object,formula,newx) {
                     newx <- .getData.sdmMatrix(formula,newx,normalize=TRUE)
                     predict(object,newx)[,1]
                   },
                   #------ metadata (optional):
                   title='multilayer perceptron (MLP) network',
                   creator='Babak Naimi',
                   authors=c("Christoph Bergmeir et al. (for the package RSNNS)"), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry(bibtype='Manual',title = "SNNS Stuttgart Neural Network Simulator User Manual",
                                          author = person("A. Zell, [aut]"),
                                          year='1998',
                                          organization="University of Stuttgart and WSI, University of Tübingen"
                   )
                   ),
                   description="MLP is a type of Neural Network, a fully connected feedforward networks, and probably the most common network architecture in use."
)