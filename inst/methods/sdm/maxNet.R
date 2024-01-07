# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Jan 2024
# last update: Jan 2024
# Version 1.0
# Licence GPL v3
#-------------

methodInfo <- list(name=c('mexNet','maxnet','mnet','mNet'),
                   packages='maxnet',
                   modelTypes = c('pb'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(f=NULL,regfun = NULL ,regmult = 1,addsamplestobackground=TRUE),
                   fitFunction = '.maxNet',
                   settingRules = NULL,
                   tuneParams = NULL,
                   predictParams=list(object='model',newdata='sdmDataFrame'),
                   predictSettings=list(type='cloglog'),
                   predictFunction='predict',
                   #------ metadata (optional):
                   title='Model Maxent over glmnet',
                   creator='Babak Naimi',
                   authors=c('Steven Phillips'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = "Opening the black box: An open‐source release of Maxent",
                                          author = as.person("S. Phillips [aut] et al."),
                                          year='2017',
                                          journal = "Ecography"
                                          
                   )
                   ),
                   description="Procedures to fit species distributions models from occurrence records and environmental variables, using 'glmnet' for model fitting."
)