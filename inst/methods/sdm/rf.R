# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  July 2017
# Version 1.1
# Licence GPL v3

#-------------
methodInfo <- list(name=c('rf','RF','randomForest','rforest'),
                   packages='randomForest',
                   modelTypes = c('pa','pb','ab','n'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(ntree=1000,
                                      replace=TRUE,
                                      importance=TRUE
                   ),
                   fitFunction = 'randomForest',
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object='model',newdata='sdmDataFrame'),
                   predictSettings=list(type='response'),
                   predictFunction='predict',
                   #------ metadata (optional):
                   title='Random Forest',
                   creator='Babak Naimi',
                   authors=c('Andy Liaw','Matthew Wiener'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = "Random Forests",
                                          author = as.person("L. Breiman [aut]"),
                                          year = "2001",
                                          journal = "Machine Learning",
                                          number="45(1)",
                                          pages="5-32"
                   )
                   ),
                   description="implements Breiman's random forest algorithm (based on Breiman and Cutler's original Fortran code) for classification and regression."
)