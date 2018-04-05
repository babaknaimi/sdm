# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update) :  Nov. 2016
# Version 1.1
# Licence GPL v3

#-------------
methodInfo <- list(name=c('bioclim.dismo','BioclimDismo','bioclimDismo','bioclim.dis','bioclimD'),
                   packages=c('dismo'),
                   modelTypes = c('po'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = NULL,
                   fitFunction = '.bioclimDismo',
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object='model',x='sdmDataFrame'),
                   predictSettings=NULL,
                   predictFunction='predict',
                   #------ metadata (optional):
                   title='Bioclim',
                   creator='Babak Naimi',
                   authors=c('Robert J. Hijmans'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = "Testing the ability of climate envelope models to predict the effect of climate change on species distributions",
                                          author = as.person("R. J. Hijmans [aut], C.H. Graham [aut]"),
                                          year='2006',
                                          number='12',
                                          pages='2272-2281',
                                          journal = "Global change biology"
                                          
                   )
                   ),
                   description="Build a 'Bioclim' (climate-envelope-model) species distribution model."
)