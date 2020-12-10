# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  Dec. 2020
# Version 1.0
# Licence GPL v3

#-------------
methodInfo <-list(name=c('glmpoly','glmpolynomial','glmp'),
                  packages='stats',
                  modelTypes = c('pa','pb','ab','n'),
                  fitParams = list(v='sdmVariables',data='sdmDataFrame'),
                  fitSettings = list(family=binomial(link='logit'),degree=3,weights=NULL,model=FALSE),
                  fitFunction = function(v,data,degree=3,...) {
                    .f <- .getFormula.glmPoly(n=c(v@response,v@variables$numeric),nFact = v@variables$nFact,degree=degree)
                    glm(formula = .f, data = data,...)
                  },
                  settingRules = function(x='sdmVariables',f='fitSettings') {
                    if (x@distribution %in% c('poisson','multinomial')) {
                      f[['family']] <- x@distribution
                    }
                    list(fitSettings=f)
                  },
                  tuneParams = NULL,
                  predictParams=list(object='model',newdata='sdmDataFrame'),
                  predictSettings=list(type='response'),
                  predictFunction='predict.glm',
                  #------ metadata (optional):
                  title='Generalized Linear Model (Polynomial)',
                  creator='Babak Naimi',
                  authors=c('R Core team'), # authors of the main method
                  email='naimi.b@gmail.com',
                  url='http://r-gis.net',
                  citation=list(bibentry('book',title = "Generalized linear models",
                                         author = as.person("P. McCullagh [aut], J. A. Nelder [aut]"),
                                         year = "1989",
                                         publisher = "Chapman and Hall",
                                         address = "London")
                  ),
                  description='glm is used to fit generalized linear models, specified by giving a symbolic description of the linear predictor and a description of the error distribution [see the help for glm function in stats package]'
)
#------------