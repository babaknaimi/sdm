# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  July 2017
# Version 1.1
# Licence GPL v3

#-------------
methodInfo <- list(name=c('svm','SVM','ksvm'),
                   packages='kernlab',
                   modelTypes = c('pa','pb','ab','n'),
                   fitParams = list(x='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(type='eps-svr',kernel='rbfdot',epsilon=0.1,prob.model=FALSE,tol=0.001,shrinking=TRUE),
                   fitFunction = 'ksvm',
                   settingRules = function(x='sdmVariables',f='fitSettings') {
                     if (x@distribution == 'multinomial') f[['type']] <- 'C-svc'
                     list(fitSettings=f)
                   },
                   tuneParams = NULL,
                   predictParams=list(object='model',newdata='sdmDataFrame'),
                   predictSettings=list(type='response'),
                   predictFunction='predict',
                   #------ metadata (optional):
                   title='Support Vector Machines',
                   creator='Babak Naimi',
                   authors=c('Alexandros Karatzoglou'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = "LIBSVM: a library for Support Vector Machines",
                                          author = as.person("C. Chih-Chung [aut], L. Chih-Jen [aut]"),
                                          year='2015',
                                          journal = "http://www.csie.ntu.edu.tw/~cjlin/libsvm"
                   )
                   ),
                   description='Support Vector Machines are an excellent tool for classification, novelty detection, and regression'
)