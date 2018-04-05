# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update) :  July 2017
# Version 1.2
# Licence GPL v3

#-------------
methodInfo <- list(name=c('maxent','Maxent','MAXENT','entropy','maxentropy'),
                   packages=c('dismo','rJava'),
                   modelTypes = c('pb'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(beta=1,prevalence=0.5,feat='autofeature',args=NULL),
                   fitFunction = '.maxent',
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object='model',x='sdmDataFrame'),
                   predictSettings=NULL,
                   predictFunction='predict',
                   #------ metadata (optional):
                   title='Maxent',
                   creator='Babak Naimi',
                   authors=c('Steven Phillips and Robert J. Hijmans'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = "Maximum entropy modeling of species geographic distributions",
                                          author = as.person("S. J. Phillips [aut], R.P. Anderson [aut], R.E. Schapire [aut]"),
                                          year='2006',
                                          number='190',
                                          pages='231-259',
                                          journal = "Ecological Modelling"
                                          
                   )
                   ),
                   description="Build a 'MaxEnt' (Maximum Entropy) species distribution model."
)