# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  July 2017
# Version 1.1
# Licence GPL v3

#-------------
methodInfo <- list(name=c('bioclim','Bioclim'),
                   packages=NULL,
                   modelTypes = c('po'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(c=2,weights=NULL),
                   fitFunction = '.bioclimFit',
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object='model',newdata='sdmDataFrame'),
                   predictSettings=NULL,
                   predictFunction='predict',
                   #------ metadata (optional):
                   title='Model occurrence probability using presence-only data',
                   creator='Babak Naimi',
                   authors=c('Babak Naimi'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = "sdm: a reproducible and extensible R platform for species distribution modelling",
                                          author = as.person("B. Naimi [aut], M.B. Araujo [aut]"),
                                          year='2016',
                                          journal = "Ecography"
                                          
                   )
                   ),
                   description="Estimates habitat suitability using presence-only data"
)