# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  Jan 2021
# Version 1.0
# Licence GPL v3

#-------------
methodInfo <- list(name=c('ranger','rangerRF','rangerForest'),
                   packages='ranger',
                   modelTypes = c('pa','pb','ab','n'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(num.trees=1000,
                                      mtry=NULL,
                                      importance='none',
                                      probability=TRUE,
                                      quantreg=FALSE,
                                      keep.inbag = FALSE,
                                      num.threads=1,
                                      verbose = FALSE
                   ),
                   fitFunction = 'ranger',
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object='model',data='sdmDataFrame'),
                   predictSettings=list(type='response',num.threads=1),
                   predictFunction=function(object,data,type,num.threads=num.threads) {
                     predict(object,data=data,type=type)$predictions
                   },
                   #------ metadata (optional):
                   title='Random Forest (Ranger)',
                   creator='Babak Naimi',
                   authors=c('Marvin N. Wright'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = " ranger: A fast implementation of random forests for high dimensional data in C++ and R",
                                          author = as.person("M.N. Wright [aut]"),
                                          year = "2017",
                                          journal = "J Stat Softw",
                                          number="77",
                                          pages="1-17"
                   )
                   ),
                   description="a fast implementation of random forests (Breiman 2001) or recursive partitioning, particularly suited for high dimensional data."
)