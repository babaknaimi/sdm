# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  August 2025
# Version 7.2
# Licence GPL v3


setClassUnion("listORcharacter", c("list","character")) 
setClassUnion("characterORnull", c("character", "NULL"))
#setClassUnion("CRSorNULL", c("CRS", "NULL"))
setClassUnion("formulaORnull", c("formula", "NULL"))
setClassUnion("numericORnull", c("numeric", "NULL"))
setClassUnion("characterORmissing", c("character", "missing"))
setClassUnion("listORnull", c("list", "NULL"))
setClassUnion("integerORnull", c("integer", "NULL"))
setClassUnion("functionORnull", c("function", "NULL"))
setClassUnion("matrixORnull", c("matrix", "NULL"))
setClassUnion("data.frameORnull", c("data.frame", "NULL"))
setClassUnion("data.frameORmatrix", c("data.frame", "matrix"))
setClassUnion("functionORcharacter", c("function", "character"))
setClassUnion("environmentORnull", c("environment", "NULL"))
setClassUnion("expressionORnull", c("expression", "NULL"))




setClass(".Metadata",
         representation(
           title='characterORnull',
           creators='listORnull',
           authors='listORnull',
           email='characterORnull',
           description='characterORnull',
           date='characterORnull',
           Help='characterORnull',
           url='characterORnull',
           citations='listORcharacter',
           licence='characterORnull'
         )
)

setClassUnion(".MetadataORnull", c(".Metadata", "NULL"))

setClass(".methodTemplate",
         representation(
           name='character',
           aliases='characterORnull',
           arguments='character',
           user.arguments='characterORnull',
           user.argument.values='listORnull',
           Help='characterORnull',
           Function='function',
           metadata='.MetadataORnull'
         )
)
#----------
#--- variable parameters 
#------- numeric keeps params like min, max, mean, sd for continuous variable
#------- categorical keeps keeps params for categorical variables (levels, count)
#------- others is a data.frame with list of other variables (time, group, info, xy)
setClass('.variables',
         representation(
           names='character',
           species='characterORnull',
           numeric='data.frameORnull',
           categorical='listORnull',
           others='data.frameORnull'
         )
)


################################
setClass('sdmFormula',
         representation(
           formula='formula',
           vars='.variables',
           data.terms='listORnull',
           model.terms='listORnull'
         )
)

#-------

#########----- terms in a formula are converted to appropriate classes:

# setClass('.pcaSetting',
#          representation(
#            vars='character',
#            scale='logical',
#            n='numericORnull'
#          )
# )
#----

setClass('.pcaSetting',
         representation(
           vars='character',
           n='ANY',
           term="call"
         )
)
#----

setClass('.scaleSetting',
         representation(
           vars='character',
           method='character',
           term='call'
         )
)
#----

# if a model in the formula is used as a term:
# ----- m and r function in formula returns model with .prediction and .residual as output
setClass('.nestedModel',
         representation(response='character',
                        predictors='character',
                        method='characterORnull',
                        distribution='character',
                        setting='listORnull',
                        output='character',
                        term='call'
         )
)
#-------
# a simple variable term:
setClass('.var',
         representation(
           name='character'
         )
)

# select function
# vars: variables to which the selection is applied
# keep: variables that are excluded from the selection procedure (should be kept in model)
# n: number to keep
# stat: selection method (either based on collinearity test, or cross-entropy, or ... )

setClass('.selectFrame',
         representation(vars='characterORnull',
                        method='character',
                        th='numericORnull',
                        keep='characterORnull',
                        term='call'
         )
)



# name of the variable (column) based on which the records are grouped
setClass('.grouping',
         representation(
           group.var='character',
           term='call'
         )
)


# name of all variables:
setClass('.all.vars',
         representation(
           names='character'
         )
)


# coordinate columns:
setClass('.coord.vars',
         representation(
           xy='character'
         )
)

# .time keeps data/time info
setClass('.time',
         representation(
           name='character',
           terms='list',
           term='call'
         )
)

# .Info keeps names of variables contain information (not for using in model)
setClass('.Info',
         representation(
           names='character'
         )
)
# --------------

setClass('.formulaFunction',
         representation(cls='call',
                        name='character',
                        args='character',
                        setFrame='functionORnull',
                        getFeature='functionORnull'
         )
)

#------- container of functions in formula:
setRefClass('.formulaFunctions',
            fields=list(
              funcNames='character',
              funcs='list'
            ),
            methods=list(
              initialize=function() {
                .self$funcs <- list()
              },
              add=function(x) {
                if (!inherits(x,'.formulaFunction')) stop('the definition of the formula function is not appropriate!')
                .self$funcs[[x@name[1]]] <- x
                .self$funcNames=unique(c(.self$funcNames,x@name[1]))
              },
              getNames=function(alt=FALSE) {
                if (alt) {
                  unique(unlist(lapply(names(.self$funcs),function(x) .self$funcs[[x]]@name)))
                } else .self$funcNames
              },
              getFuncs=function(n) {
                if (missing(n)) .self$funcs[getNames()]
                else {
                  mn <- getNames()
                  if (!all(n %in% mn)) {
                    w <- which(!n %in% mn)
                    names(mn) <- mn
                    mnlist <- lapply(mn,function(x) .self$funcs[[x]]@name)
                    for (i in w) {
                      u <- unlist(lapply(mnlist,function(x) n[i] %in% x))
                      if (any(u)) n[i] <- names(u)[which(u)]
                      else {
                        u <- unlist(lapply(mn,function(x) any(!is.na(pmatch(n[i],x)))))
                        if (any(u)) n[i] <- names(u)[which(u)]
                        else warning(paste(n[i],'is not a registered formula function!'))
                      }
                    }
                  }
                  if (any(n %in% mn)) {
                    w <- which(n %in% mn)
                    nw <- n[w]
                    names(nw) <- nw
                    .self$funcs[nw]
                  } 
                }
              },
              setClasses=function(n=getNames()) {
                try(a <- lapply(n,function(x) eval(.self$funcs[[x]]@cls)),silent=TRUE)
                rm(a)
              },
              show=function(...) {
                cat('container class                 :' , class(.self), '\n')
                cat('=====================================================','\n')
                cat('number of methods               : ' ,length(.self$funcs), '\n')
                cat('name of methods                 : ' , paste(getNames(),collapse=', '),'\n')
                cat('-----------------------------------------------------\n')
                
              }
            )
)
#-----------
setOldClass("princomp")

setClass(".pcaObject",
         representation(
           data='ANY',
           pcaObject='princomp',
           scaled='logical'
         )
)
#------
#------------


setClass('featuresFrame',
         representation(
           vars='character',
           feature.types='list',
           response.specific='listORnull'
         )
)
# setClass('featuresFrame',
#          representation(
#            vars='character',
#            feature.types='list',
#            resonse.specific='listORnull',
#            model.specific='listORnull'
#          )
# )

# setClass('.featureFrame',
#          representation(
#            var='character',
#            feature.name='character',
#            type='character',
#            params='listORnull',
#            response='characterORnull'
#          )
# )
setClass(".featureFrame",
         representation(
           responses='character',
           predictors='character',
           numeric='data.frameORnull',
           categorical='listORnull',
           transformers = 'functionORnull', # list of functions for data transform (scale, pca, etc.)
           featureGenerator='functionORnull'
         )
)
###########################


#---- classes corresponding to sdmdata

setClass('.group',
         representation(
           name='character',
           values='data.frame',
           indices='list'
         )
)

setClass(".info",
         representation(
           info='data.frameORnull',
           time='data.frameORnull',
           coords='matrixORnull',
           crs='ANY',
           metadata='.MetadataORnull'
         )
)
#------
setClassUnion(".infoORnull", c(".info", "NULL"))
#-------
setClass('.species.data',
         representation(
           name='character',
           type='character',
           presence='numericORnull',
           absence='numericORnull',
           background='numericORnull',
           abundance='data.frameORnull',
           numerical='data.frameORnull',
           Multinomial='data.frameORnull'
         )
)
#----------


setClass('sdmdata',
         representation(
           species.names='character',
           species='list',
           features='data.frameORnull',
           features.name='characterORnull',
           factors='characterORnull',
           info='.infoORnull',
           groups='list',
           sdmFormula='sdmFormula',
           errorLog='list'
         )
)
#-----

.methods <- setRefClass('.methods',
                        fields=list(Methods="list",
                                    arguments="vector",
                                    outputs="list",
                                    test.values="list",
                                    template="function",
                                    help='character'),
                        methods=list(
                          initialize=function() {
                            
                          },
                          #----
                          getMethodNames=function(alt=FALSE) {
                            if (alt) {
                              mn <- names(.self$Methods)
                              names(mn) <- mn
                              lapply(mn,function(x) unique(c(x,.self$Methods[[x]]@aliases)))
                            } else names(.self$Methods)
                          },
                          #---
                          getHelp=function() {cat(.self$help)},
                          #---
                          whichMethod=function(n) {
                            if (length(n) > 1) n <- n[1]
                            mn <- getMethodNames(alt=TRUE)
                            if (!n %in% names(mn)) {
                              w <- which(unlist(lapply(mn,function(x) n %in% x)))
                              if (length(w) == 0) w <- which(unlist(lapply(tolower(names(mn)),function(x) any(!is.na(pmatch(tolower(n),x))))))
                              if (length(w) > 0) names(mn)[w[1]]
                            } else n
                          },
                          #---
                          getFunctions=function(n=getMethodNames()) {
                            mn <- getMethodNames()
                            if (!all(n %in% mn)) {
                              w <- which(!n %in% mn)
                              names(mn) <- mn
                              mnlist <- lapply(mn,function(x) .self$Methods[[x]]@aliases)
                              for (i in w) {
                                u <- unlist(lapply(mnlist,function(x) n[i] %in% x))
                                if (any(u)) n[i] <- names(u)[which(u)]
                                else {
                                  u <- unlist(lapply(mn,function(x) any(!is.na(pmatch(n[i],x)))))
                                  if (any(u)) n[i] <- names(u)[which(u)]
                                  else warning(paste(n[i],'is not a registered method!'))
                                }
                              }
                            }
                            if (any(n %in% mn)) {
                              w <- which(n %in% mn)
                              nw <- n[w]
                              names(nw) <- nw
                              lapply(nw, function(x) .self$Methods[[x]]@Function)
                            } 
                          },
                          #---
                          addMethod=function(x,echo=TRUE) {
                            if (x@name %in% unlist(getMethodNames(alt=TRUE))) stop('a method with the same name (or alternative name) does exist,\n Use different name; or use updateMethod to change the existing method')
                            else {
                              w <- unlist(lapply(names(x@arguments),function(n) n %in% names(.self$arguments)))
                              if (!all(w)) {
                                w <- which(!w)
                                x@user.arguments <- x@arguments[w]
                                x@arguments <- x@arguments[-w]
                                x@user.argument.values <- x@user.argument.values[names(x@user.arguments)]
                              }
                              x@Function <- .templateMatch(x@Function,.self$template)
                              if (.testMethod(x,template=.self$template,arguments=.self$arguments,outputs=.self$outputs,test.args=.self$test.values)) {
                                .self$Methods[[x@name]] <- x
                                if (echo) cat('method',x@name,'is successfully added to the',class(.self),' object.\n')
                              } else cat('Error: Method is not added...!')
                              
                            }
                          },
                          
                          #---
                          updateMethod=function(name,alt=NULL,args=NULL,Help=NULL,f=NULL) {
                            if (!name %in% unlist(getMethodNames(alt=TRUE))) stop('the specified method does not exist!')
                            else {
                              name <- whichMethod(name)
                              x <- .self$Methods[[name]]
                              if (!is.null(alt)) x@aliases <- alt
                              if (!is.null(args)) x@arguments <- args
                              if (!is.null(Help)) x@help <- Help
                              if (!is.null(f)) {
                                x@Function <- .templateMatch(f,.self$template)
                                if (.testMethod(x,template=.self$template,arguments=.self$arguments,outputs=.self$outputs,test.args=.self$test.values)) {
                                  .self$Methods[[x@name]] <- x
                                  cat('method',x@name,'is successfully updated.\n')
                                } else stop('Error: Method is not updated!')
                              } else {
                                .self$Methods[[x@name]] <- x
                                cat('method',x@name,'is successfully updated.\n')
                              }
                            }
                          },
                          #-----
                          deleteMethod=function(name) {
                            if (!name %in% unlist(getMethodNames(alt=TRUE))) stop('a method with the specified name does exist!')
                            name <- whichMethod(name)
                            .self$Methods <- .self$Methods[-which(getMethodNames() == name)]
                            cat('Method',name,'is successfully deleted.\n')
                          },
                          show=function(...) {
                            cat('container class                 :' , class(.self), '\n')
                            cat('=====================================================','\n')
                            cat('number of methods               : ' ,length(.self$Methods), '\n')
                            cat('name of methods                 : ' , paste(getMethodNames(),collapse=', '),'\n')
                            cat('reserved argument names         : ' ,paste(names(.self$arguments),collapse=', ') , '\n')
                            cat('-----------------------------------------------------','\n')
                            
                          },
                          example=function(name) {
                            name <- whichMethod(name)
                            if (!is.null(name)) {
                              x <- .self$Methods[[name]]
                              test.args <- c(.self$test.values,x@user.argument.values)
                              o <- try(do.call(x@Function,test.args),TRUE)
                              if(!inherits(o, "try-error")) o
                            }
                          }
                        )
)



# sdm methods:
setClass("sdmCorrelativeMethod",
         representation(
           name='character',
           aliases='characterORnull',
           packages='characterORnull',
           modelTypes='characterORnull',
           dataArgument.names='listORnull',
           fitParams='list',
           fitSettings='listORnull',
           settingRules='functionORnull',
           fitFunction='function',
           tuneParams='listORnull',
           predictParams='listORnull',
           predictSettings='listORnull',
           predictFunction='functionORnull',
           metadata='.MetadataORnull',
           .temp.env='environmentORnull'
         )
)
#--------

#---------------------------------------------------------
setRefClass('.sdmMethodsContainer',
            fields=list(
              MethodDefinitions='data.frame',
              Methods="list",
              test.data="list",
              userFunctions='environment',
              help='character'),
            methods=list(
              initialize=function() {
                d <- data.frame(matrix(nrow=0,ncol=6))
                colnames(d) <- c('name','apprach','type','dataType','formulaType','inTempEnv')
                for (i in 1:5) d[,i] <- as.character(d[,i])
                .self$MethodDefinitions <- d
              },
              #----
              getMethodNames=function(alt=FALSE) {
                if (alt) {
                  mn <- names(.self$Methods)
                  names(mn) <- mn
                  lapply(mn,function(x) unique(c(x,.self$Methods[[x]]@aliases)))
                } else names(.self$Methods)
              },
              #---
              getHelp=function() {cat(.self$help)},
              #---
              whichMethod=function(n) {
                if (length(n) > 1) n <- n[1]
                mn <- getMethodNames(alt=TRUE)
                if (!n %in% names(mn)) {
                  w <- which(unlist(lapply(mn,function(x) n %in% x)))
                  if (length(w) == 0) w <- which(unlist(lapply(tolower(names(mn)),function(x) any(!is.na(pmatch(tolower(n),x))))))
                  if (length(w) > 0) names(mn)[w[1]]
                } else n
              },
              #---
              fixNames=function(n) {
                mn <- getMethodNames()
                if (!all(n %in% mn)) {
                  w <- which(!n %in% mn)
                  names(mn) <- mn
                  mnlist <- lapply(mn,function(x) .self$Methods[[x]]@aliases)
                  for (i in w) {
                    u <- unlist(lapply(mnlist,function(x) n[i] %in% x))
                    if (any(u)) n[i] <- names(u)[which(u)]
                    else {
                      u <- unlist(lapply(mn,function(x) any(!is.na(pmatch(n[i],x)))))
                      if (any(u)) n[i] <- names(u)[which(u)]
                      else warning(paste(n[i],'is not a registered sdm method!'))
                    }
                  }
                }
                n
              },
              #---
              getFitFunctions=function(n=getMethodNames()) {
                mn <- getMethodNames()
                n <- fixNames(n)
                if (any(n %in% mn)) {
                  w <- which(n %in% mn)
                  nw <- n[w]
                  names(nw) <- nw
                  lapply(nw, function(x) .self$Methods[[x]]@fitFunction)
                } 
              },
              #---
              getFitArguments=function(n=getMethodNames()) {
                mn <- getMethodNames()
                n <- fixNames(n)
                if (any(n %in% mn)) {
                  w <- which(n %in% mn)
                  nw <- n[w]
                  names(nw) <- nw
                  lapply(nw, function(x) list(params=.self$Methods[[x]]@fitParams,settings=.self$Methods[[x]]@fitSettings))
                } else stop('none of the specified methods are registered sdm Methods!')
              },
              #---
              getPredictFunctions=function(n=getMethodNames()) {
                mn <- getMethodNames()
                n <- fixNames(n)
                
                if (any(n %in% mn)) {
                  w <- which(n %in% mn)
                  nw <- n[w]
                  names(nw) <- nw
                  lapply(nw, function(x) .self$Methods[[x]]@predictFunction)
                } 
              },
              #---
              getPredictArguments=function(n=getMethodNames()) {
                mn <- getMethodNames()
                n <- fixNames(n)
                if (any(n %in% mn)) {
                  w <- which(n %in% mn)
                  nw <- n[w]
                  names(nw) <- nw
                  lapply(nw, function(x) list(params=.self$Methods[[x]]@predictParams,settings=.self$Methods[[x]]@predictSettings))
                } else stop('none of the specified methods are registered sdm Methods!')
              },
              #---
              getPackageNames=function(m=getMethodNames()) {
                for (i in seq_along(m)) m[i] <- .self$whichMethod(m[i])
                names(m) <- m
                lapply(m,function(x) {.self$Methods[[x]]@packages})
              },
              #---
              addMethod=function(x,echo=TRUE) {
                if (x@name %in% unlist(getMethodNames(alt=TRUE))) stop('a method with the same name (or alternative name) does exist,\n Use different name; or use updateMethod to change the existing method')
                else {
                  .self$Methods[[x@name]] <- x
                  i <- nrow(.self$MethodDefinitions)+1
                  .self$MethodDefinitions[i,1] <- x@name
                  if (inherits(x,'sdmCorrelativeMethod'))  {
                    .self$MethodDefinitions[i,2] <- 'correlative'
                    
                    if (!is.null(x@modelTypes)) .self$MethodDefinitions[i,3] <- paste(x@modelTypes,collapse='.')
                    else .self$MethodDefinitions[i,3] <- 'all'
                    
                    if (!is.null(x@dataArgument.names)) .self$MethodDefinitions[i,4] <- paste(x@dataArgument.names,collapse=';')
                    else {
                      nfit <- npred <- nt <- NULL
                      w <- which(x@fitParams %in% c('data.frame','sdmDataFrame','sdmX','sdmDataFrame.norm','sdmX.norm','sdmY','sdmRaster','sdmMatrix','sdmMatrix.norm'))
                      if (length(w) > 0) {
                        nt <- x@fitParams[w]
                        nfit <- names(x@fitParams)[w]
                      }
                      
                      w <- which(x@predictParams %in% c('data.frame','sdmDataFrame','sdmX','sdmDataFrame.norm','sdmX.norm','sdmY','sdmRaster','sdmMatrix','sdmMatrix.norm'))
                      if (length(w) > 0) {
                        npred <- names(x@predictParams)[w]
                        nt <- c(nt,x@predictParams[w])
                      }
                      n <- c(nfit,npred)
                      if (!is.null(n)) {
                        .self$MethodDefinitions[i,4] <- paste(unique(nt),collapse=';')
                        .self$Methods[[x@name]]@dataArgument.names <- list(fit=nfit,predict=npred)
                      }
                      
                    }
                    if (!is.null(x@fitParams$formula)) .self$MethodDefinitions[i,5] <- x@fitParams$formula
                    if ('.temp' %in% x@packages) {
                      .self$userFunctions <- .movEnv(x@.temp.env,.self$userFunctions)
                      #e <- .self$userFunctions
                      .self$MethodDefinitions[i,6] <- TRUE
                      #.movEnv2sdm(e)
                      #rm(e)
                      #x@packages <- x@packages[-which(x@packages == '.temp')]
                    } else .self$MethodDefinitions[i,6] <- FALSE
                  }
                  
                  ## -other types of models needs to be checked and included here
                  
                  if (echo) cat('method',x@name,'is successfully added to the',class(.self),' object.\n')
                }
              },
              
              #---
              updateMethod=function(x=NULL,...) {
                if (!is.null(x)) {
                  if (!x@name %in% unlist(getMethodNames(alt=FALSE))) stop('the specified method does not exist!')
                  .self$Methods[[x@name]] <- x
                } else {
                  name <- list(...)[['name']]
                  if (is.null(name)) stop('the name of method to update is not specified!')
                  if (!name %in% unlist(getMethodNames(alt=TRUE))) stop('the specified method does not exist!')
                  else {
                    name <- whichMethod(name)
                    x <- .self$Methods[[name]]
                    xx <- .update.sdmCorrelativeMethod(x,...)
                    .self$Methods[[name]] <- xx
                    cat('method',x@name,'is successfully updated.\n')
                  }
                }
                
              },
              #-----
              deleteMethod=function(name,echo=TRUE) {
                if (!name %in% unlist(getMethodNames(alt=TRUE))) stop('a method with the specified name does exist!')
                name <- whichMethod(name)
                .self$Methods <- .self$Methods[-which(getMethodNames() == name)]
                if (echo) cat('Method',name,'is successfully deleted.\n')
              },
              getDataArgumentNames=function(name) {
                if (!name %in% unlist(getMethodNames(alt=TRUE))) stop('a method with the specified name does exist!')
                name <- whichMethod(name)
                .self$Methods[[name]]@dataArgument.names
              },
              show=function(...) {
                cat('container class                 :' , class(.self), '\n')
                cat('=====================================================','\n')
                cat('number of methods               : ' ,length(.self$Methods), '\n')
                cat('name of methods                 : ' , paste(getMethodNames(),collapse=', '),'\n')
                cat('-----------------------------------------------------\n')
                
              },
              example=function() {
                cat('example is not implemented...\n')
                cat('template:\n')
                print(.self$template)
              }
            )
)

#----------------------

# container class of replicate methods
setRefClass('.ReplicateMethods',
            contain='.methods',
            methods=list(
              initialize=function() {
                .self$arguments=c(x='numeric',replicates='numeric',nfolds='numeric',test.percent='numeric',family='character',stratify='logical')
                .self$test.values=list(x=c(1,1,1,1,0,0,0,0),replicates=2,nfolds=2,test.percent=20,family='binomial',stratify=TRUE)
                .self$outputs= list(c('numeric','list'),'matrix')
                .self$template=function(x,...) {
                  list()
                }
                #----
                .self$help='The ReplicateMethods object is a container of the resampling methods that partition the main dataset into training and test.
                A user can add a new method by supplying a function in which the arguments are selected from the reserved list (new arguments can also be included) and the output of the function should be the same as the defined output type.
                
                Following is the reserved arguments as well as the output type.
                # inputs:
                -  x: numeric vector e.g. species occurrence: c(1,1,0,0,1,0,1)
                - replicates= number of replicates
                - nfolds : number of folds in cross-validation procedure
                - family : distribution family of values in x
                - test.percent: a proportion of data that should be used as a test dataset
                - stratify: for binomial data, specifies whether the resampling should be stratified based on presence/absence
                
                # output: a list with two items:
                #-----: [[1]] a numeric vector with the same length as x 
                #             including values of 1 or 2, specifies whether the 
                #             corresponding item in x should be used for train or test
                #-----: [[2]] a matrix (nrows=length(x), ncol= number of replicates)
                #             each column include values ranging between 1:length(x),
                #-----------------
                
                run example function for this object to see the output of the existing methods for a simple example.
                '
                
                
                #----
                
              },
              example=function(name) {
                name <- whichMethod(name)
                if (!is.null(name)) {
                  x <- .self$Methods[[name]]
                  test.args <- c(.self$test.values,x@user.argument.values)
                  o <- try(do.call(x@Function,test.args),TRUE)
                  if(!inherits(o, "try-error")) o
                }
              }
            )
)
#----------------

setClass('.sdmCorModel',
         representation(
           mID='numeric',
           method='character',
           response='character',
           object='ANY',
           evaluation='list',
           varImportance='list',
           errorLog='list'
         )
)
#-------

setClass('.parallelSetting',
         representation(
           cl='ANY',
           hosts='characterORnull',
           ncore='numeric',
           fork='logical',
           doParallel='expressionORnull',
           type='characterORnull',
           method='character',
           strategy='characterORnull'
         )
)
#----
setClassUnion(".parallelSettingORnull", c(".parallelSetting", "NULL"))
#------------
# 
# setClass('.sdmCorSetting',
#          representation(
#            methods='character',
#            sdmFormula='sdmFormula',
#            featuresFrame='featuresFrame',
#            distribution='character',
#            interaction.depth='numericORnull',
#            test.percentage='numericORnull',
#            replicate='characterORnull',
#            n.replicates='numericORnull',
#            cv.folds='numericORnull',
#            pseudo.absence.methods='characterORnull',
#            n.pseudo.absence='numericORnull',
#            varImportance.methods='characterORnull',
#            var.selection='logical',
#            response.curve='logical',
#            modelSettings='listORnull',
#            seed='numericORnull',
#            parallelSettings='.parallelSettings',
#            errorLog='list'
#          )
# )
setClass('.sdmCorSetting',
         representation(
           methods='character',
           sdmFormula='sdmFormula',
           featureFrame='.featureFrame',
           distribution='character',
           interaction.depth='numericORnull',
           test.percentage='numericORnull',
           replicate='characterORnull',
           n.replicates='numericORnull',
           cv.folds='numericORnull',
           pseudo.absence.methods='characterORnull',
           n.pseudo.absence='numericORnull',
           varImportance.methods='characterORnull',
           var.selection='logical',
           response.curve='logical',
           modelSettings='listORnull',
           seed='numericORnull',
           parallelSetting='.parallelSetting',
           errorLog='list'
         )
)
setClass('.sdmVariables',
         representation(
           response='character',
           variables='list',
           distribution='character',
           features.numeric='characterORnull',
           features.factor='characterORnull',
           number.of.records='numeric',
           n.presence='numericORnull',
           varInfo='listORnull'
         )
)
#-------------------
setClass("sdmModels",
         representation(
           data='sdmdata',
           recordIDs='list',
           setting='.sdmCorSetting',
           run.info='data.frame',
           replicates='list',
           models='list'
         )
)

#-------

setRefClass(".workload",
            fields=list(
              data='sdmdata',
              setting='.sdmCorSetting',
              frame=".featureFrame",
              train='list',
              test='listORnull',
              sdmVariables='list',
              params='list',
              arguments='list',
              dataObject.names='list',
              funs='list',
              replicates='list',
              settingRules='list',
              tuneParams='list',
              recordIDs='list',
              tasks='ANY',
              ncore='numericORnull',
              filename='characterORnull'
            ),
            methods=list(
              fit=function(woL=.self,species,models,runs,hasTest,.parMethod=.self$setting@parallelSetting@method,.hostnames=.self$setting@parallelSetting@hosts,.fork=.self$setting@parallelSetting@fork,filename=.self$filename) {
                .fit(woL=woL,species=species,runs=runs,hasTest = hasTest,.parMethod=.parMethod,.hostnames = .hostnames,.fork = .fork,filename = filename)
              },
              getSdmVariables=function(sp,nf,nFact) {
                if (length(.self$sdmVariables) > 0 && !is.null(.self$sdmVariables[[sp]])) .self$sdmVariables[[sp]]
                else {
                  if (missing(nFact) || is.null(nFact)) {
                    nFact <- .self$setting@sdmFormula@vars@categorical
                  }
                  
                  if (missing(nf) || is.null(nFact)) {
                    nf <- .self$setting@sdmFormula@vars@numeric$names
                  }
                  
                  .self$sdmVariables[[sp]] <- new('.sdmVariables',response=sp,variables=list(numeric=nf,factors=nFact),distribution=.self$setting@distribution[[sp]],features.numeric=.excludeVector(colnames(.self$train[[sp]]$sdmDataFrame),c(sp,nFact)),features.factor=nFact,
                                                  number.of.records=if (is.null(.self$test)) nrow(.self$train[[sp]]$sdmDataFrame) else c(train=nrow(.self$train[[sp]]$sdmDataFrame),test=nrow(.self$test[[sp]]$sdmDataFrame)),
                                                  n.presence=if (.self$setting@distribution[[sp]] == 'binomial') length(which(.self$train[[sp]]$sdmDataFrame[,1] == 0)) else NULL,varInfo=list(numeric=.self$frame@numeric,categorical=.self$frame@categorical))
                  .self$sdmVariables[[sp]]
                }
              },
              getDataFrame=function(id=NULL,train=TRUE,sp=NULL,grp=NULL) {
                if (is.null(id)) {
                  if (train) as.data.frame(.self$data,grp='train',sp=sp)
                  else as.data.frame(.self$data,grp=grp,sp=sp)
                } else {
                  as.data.frame(.self$data,ind=id,sp=sp,grp=grp)
                }
              },
              normalize=function(x,except=NULL,scale=FALSE) {
                
                .fr <- .self$frame@numeric
                
                w <- !.where(is.factor,x)
                if (!is.null(except)) {
                  w[except] <- FALSE
                }
                #---
                if (any(w)) {
                  if (all(names(w[w]) %in% .fr$names)) {
                    .fr <- .fr[.fr$names %in% names(w[w]),]
                  } else .fr <- NULL
                  
                  
                  if (!is.null(.fr)) {
                    if (scale) .sc <- .getScaleFunction(.fr,scl='minmax')
                    else .sc <- .getScaleFunction(.fr,scl='center')
                    #-------
                    x <- .sc(x)
                  } else {
                    for (i in seq_along(w)) {
                      if (w[i]) {
                        if (scale) {
                          .min <- min(x[,i],na.rm=TRUE)
                          .max <- max(x[,i],na.rm=TRUE)
                          x[,i] <- (x[,i] - .min) / (.max - .min)
                        } else {
                          .mean <- mean(x[,i],na.rm=TRUE)
                          .sd <- sd(x[,i],na.rm=TRUE)
                          x[,i] <- x[,i] - .mean
                          if (.sd != 0) x[,i] <- x[,i] / .sd
                        }
                      }
                    }
                  }
                }
                x
              },
              generateParams=function(n,sp,train=TRUE,data=TRUE) {
                # if data=FALSE, the type of data is returned rather than data object (i.e., sdmDataFrame)
                for (i in seq_along(n)) {
                  if (n[[i]] == 'sdmDataFrame') {
                    if (data) {
                      if (train) n[[i]] <- .self$train[[sp]]$sdmDataFrame
                      else n[[i]] <- .self$test[[sp]]$sdmDataFrame
                    } #else n[[i]] <- 'sdmDataFrame'
                  } else if (n[[i]] == 'data.frame') {
                    if (data) {
                      if (train) n[[i]] <- .self$getDataFrame(train=TRUE,sp=sp)
                      else n[[i]] <- .self$getDataFrame(train=FALSE,grp='test',sp=sp)
                    } #else n[[i]] <- 'sdmDataFrame'
                  } else if (n[[i]] == 'sdmVariables') {
                    n[[i]] <- getSdmVariables(sp)
                  } else if (n[[i]] == 'standard.formula') {
                    n[[i]] <- .getFormula(colnames(.self$train[[sp]]$sdmDataFrame),env=parent.frame(2))
                  } else if (n[[i]] == 'gam.mgcv.formula') {
                    sv <- .self$getSdmVariables(sp)
                    n[[i]] <- .getFormula.gammgcv(c(sp,sv@features.numeric),sv@features.factor,env=parent.frame(2))
                  } else if (n[[i]] == 'sdmX') {
                    if (data) {
                      if (train) n[[i]] <- .self$train[[sp]]$sdmDataFrame[,colnames(.self$train[[sp]]$sdmDataFrame) != sp,drop=FALSE]
                      else n[[i]] <- .self$test[[sp]]$sdmDataFrame[,colnames(.self$test[[sp]]$sdmDataFrame) != sp,drop=FALSE]
                    } #else 'sdmX'
                  } else if (n[[i]] == 'sdmY') {
                    if (data) {
                      if (train) n[[i]] <- .self$train[[sp]]$sdmDataFrame[,sp]
                      else n[[i]] <- .self$test[[sp]]$sdmDataFrame[,sp]
                    } #else 'sdmY'
                  } else if (n[[i]] == 'sdmX.norm') {
                    if (data) {
                      if (train) n[[i]] <- .self$normalize(.self$train[[sp]]$sdmDataFrame[,colnames(.self$train[[sp]]$sdmDataFrame) != sp,drop=FALSE])
                      else n[[i]] <- .self$normalize(.self$test[[sp]]$sdmDataFrame[,colnames(.self$test[[sp]]$sdmDataFrame) != sp,drop=FALSE])
                    } #else 'sdmX.norm'
                  } else if (n[[i]] == 'sdmDataFrame.norm') {
                    if (data) {
                      if (train) n[[i]] <- .self$normalize(.self$train[[sp]]$sdmDataFrame,except=sp)
                      else n[[i]] <- .self$normalize(.self$test[[sp]]$sdmDataFrame,except=sp)
                    } #else 'sdmDataFrame.norm'
                  }  else if (n[[i]] == 'sdmDataFrame.scaled') {
                    if (data) {
                      if (train) n[[i]] <- .self$normalize(.self$train[[sp]]$sdmDataFrame,except=sp,scale = TRUE)
                      else n[[i]] <- .self$normalize(.self$test[[sp]]$sdmDataFrame,except=sp,scale = TRUE)
                    } #else 'sdmDataFrame.norm'
                  } else if (n[[i]] == 'sdmMatrix') {
                    if (data) {
                      if (train) n[[i]] <- model.matrix(.getFormula(colnames(.self$train[[sp]]$sdmDataFrame),env=parent.frame(2)),.self$train[[sp]]$sdmDataFrame)[,-1]
                      else n[[i]] <- model.matrix(.getFormula(colnames(.self$train[[sp]]$sdmDataFrame),env=parent.frame(2)),.self$test[[sp]]$sdmDataFrame)[,-1]
                    } #else 'sdmMatrix'
                  } else if (n[[i]] == 'sdmMatrix.norm') {
                    if (data) {
                      if (train) {
                        #w <- .where(is.factor,.self$train[[sp]]$sdmDataFrame)
                        n[[i]] <- model.matrix(.getFormula(colnames(.self$train[[sp]]$sdmDataFrame),env=parent.frame(2)),.self$normalize(.self$train[[sp]]$sdmDataFrame,except=sp))[,-1]
                      } else {
                        n[[i]] <- model.matrix(.getFormula(colnames(.self$train[[sp]]$sdmDataFrame),env=parent.frame(2)),.self$normalize(.self$test[[sp]]$sdmDataFrame,except=sp))[,-1]
                      }
                    }
                  } else if (n[[i]] == 'sdmMatrix.scaled') {
                    if (data) {
                      if (train) {
                        n[[i]] <- model.matrix(.getFormula(colnames(.self$train[[sp]]$sdmDataFrame),env=parent.frame(2)),.self$normalize(.self$train[[sp]]$sdmDataFrame,except=sp,scale = TRUE))[,-1]
                      } else {
                        n[[i]] <- model.matrix(.getFormula(colnames(.self$train[[sp]]$sdmDataFrame),env=parent.frame(2)),.self$normalize(.self$test[[sp]]$sdmDataFrame,except=sp,scale = TRUE))[,-1]
                      }
                    }
                  } else if (n[[i]] == 'sdmdata') {
                    if (data) n[[i]] <- .self$data
                  } else if (n[[i]] == 'sdmSetting') {
                    if (data) n[[i]] <- .self$setting
                  } else if (n[[i]] == 'sdmNrRecords') {
                    if (data) {
                      if (is.null(names(.self$replicates))) {
                        n[[i]] <- list(train=nrow(.self$train[[sp]]$sdmDataFrame),test=if (is.null(.self$test)) 0 else nrow(.self$test[[sp]]$sdmDataFrame),replicates=NULL)
                      } else {
                        .re <- list()
                        for (.r in .self$setting@replicate) {
                          .w <- which(unlist(lapply(.self$replicates[[sp]],function(x) x$method== .r)))[1]
                          .re[[.r]] <- c(train=length(.self$replicates[[sp]][[.w]]$train),test=length(.self$replicates[[sp]][[.w]]$test))
                        }
                        n[[i]] <- list(train=nrow(.self$train[[sp]]$sdmDataFrame),test=if (is.null(.self$test)) 0 else nrow(.self$test[[sp]]$sdmDataFrame),replicates=.re)
                      }
                    }
                  } else if (n[[i]] == 'sdmRaster') {
                    #####
                    if (data) {
                      ###
                    } #else 'sdmRaster'
                  } else if (n[[i]] %in% names(.self$params)) {
                    n[[i]] <- do.call(.self$params[[n[[i]]]],generateParams(.CharVector2List(names(formals(.self$params[[n[[i]]]]))),sp)) 
                  }
                }
                n
              },
              # w is either 2 or 3 (train, test)
              getData=function(sp,run=NULL,w=2,d='sdmDataFrame',train=TRUE) {
                if (train) {
                  if (!is.null(run)) {
                    .self$train[[sp]][[d]][.self$replicates[[sp]][[run]][[w]],]
                  } else .self$train[[sp]][[d]]
                } else .self$test[[sp]][[d]]
              },
              getReseved.names=function() {
                c('sdmdata','data.frame','sdmDataFrame','sdmX','sdmY','sdmRaster','sdmVariables','standard.formula','gam.mgcv.furmula','sdmMatrix','sdmMatrix.norm','sdmMatrix.scaled','sdmDataFrame.norm','sdmDataFrame.scaled','sdmX.norm','sdmNrRecords','sdmSetting')
              },
              getFitArgs=function(sp,mo) {
                o <- list()
                pa <- .self$arguments$fit[[mo]]$params
                n <- names(pa)
                ww <- which(names(pa) %in% .self$dataObject.names[[mo]][['fit']])
                if (length(ww) == 0) stop('data object required by the fit function is not recognised!')
                
                #for (nn in n) o[[nn]] <- .self$generateParams(pa[[nn]],sp,data=FALSE)
                o <- .self$generateParams(pa[n],sp,data=FALSE)
                #o[[n[ww]]] <- pa[[ww]]
                o <- c(o,.self$arguments$overriden_settings$fit[[mo]][[sp]])
                o
              },
              getPredictArgs=function(sp,mo) {
                # return a list in which the first element is reserved for 'model'
                # and the second element is reserved for data e.g., 'newdata'
                # these two elements will be updated before putting in the predict function
                o <- list()
                pa <- .self$arguments$predict[[mo]]$params
                n <- names(pa)
                ww <- which(pa == 'model')
                o[[n[ww]]] <- 'model'
                n <- n[-ww]
                pa<- pa[-ww]
                ww <- which(names(pa) %in% .self$dataObject.names[[mo]][['predict']])
                if (length(ww) == 0) stop('data object required by the predict function is not recognised!')
                o[[n[ww]]] <- pa[[ww]]
                n <- n[-ww]
                pa<- pa[-ww]
                
                if (length(n) > 0) {
                  #for (nn in n) o[[nn]] <- .self$generateParams(pa[[nn]],sp)
                  o <- c(o,.self$generateParams(pa[n],sp))
                }
                
                o <- c(o,.self$arguments$overriden_settings$predict[[mo]][[sp]])
                o
              },
              setRules=function(mo,sp) {
                se <- .self$settingRules[[mo]]
                if (inherits(se,'function')) {
                  fo <- as.list(formals(se))
                  n <- names(fo)
                  if (any(n == '...')) {
                    .w <- which(n == '...')
                    fo <- fo[-.w]
                    n <- n[-.w]
                  }
                  .w <- .pmatch(unlist(fo),c(.self$getReseved.names(),'fitSettings','predictSettings'))
                  if (!all(is.na(.w))) fo[which(!is.na(.w))] <- .w[which(!is.na(.w))]
                  
                  .ww <- which(fo %in% c('fitSettings','predictSettings'))
                  if (length(.ww) > 0) {
                    .w <- which(fo %in% c('fitSettings'))
                    if (length(.w) > 0) {
                      fo[[.w]] <- .self$arguments$fit[[mo]]$settings
                    }
                    .w <- which(fo %in% c('predictSettings'))
                    if (length(.w) > 0) {
                      fo[[.w]] <- .self$arguments$predict[[mo]]$settings
                    }
                    fo[-.ww] <- .self$generateParams(fo[-.ww],sp=sp)
                    do.call(se,fo)
                  } else NULL
                } else NULL
              },
              show=function(...) {
                cat('container class                 :' , class(.self), '\n')
                cat('=====================================================','\n')
                cat('species name                    : ' ,paste(names(.self$data),collapse=', '),'\n')
                cat('number of methods               : ' ,length(.self$setting@methods), '\n')
                cat('name of methods                 : ' , paste(.self$setting@methods,collapse=', '),'\n')
                cat('-----------------------------------------------------\n')
                
              }
              
            )
)
#----------




# 
# setRefClass(".workloadP",
#             fields=list(
#               obj='list', # list of models (@models) from sdmModels
#               #data='data.frame', # data used to fit the models
#               newdata='list',
#               modelFrame='data.frameORnull',
#               params='list',
#               arguments='list',
#               dataObject.names='character',
#               funs='list',
#               settingRules='list',
#               runTasks='data.frame',
#               ncore='numericORnull'
#             ),
#             methods=list(
#               predictMID=function(IDs) {
#                 options(warn=-1)
#                 IDs <- which(.self$runTasks$modelID %in% IDs)
#                 m <- lapply(IDs,function(i) {
#                   p <- .self$getPredictArgs(.self$runTasks$species[i],.self$runTasks$method[i])
#                   p[[1]] <- .self$obj[[.self$runTasks$speciesID[i]]][[.self$runTasks$methodID[i]]][[.self$runTasks$mIDChar[i]]]@object
#                   
#                   if (is.null(.self$modelFrame$specis_specific)) p[[2]] <- .self$modelFrame$features
#                   else p[[2]] <- cbind(.self$modelFrame$features,.self$modelFrame$specis_specific[[.self$runTasks$species[i]]])
#                   
#                   m <- try(.self$funs[[.self$runTasks$methodID[i]]](p),silent=TRUE)
#                   m
#                   
#                 })
#                 options(warn=0)
#                 m
#               },
#               predictID=function(i) {
#                 options(warn=-1)
#                 i <- which(.self$runTasks$modelID == i)
#                 p <- .self$getPredictArgs(.self$runTasks$species[i],.self$runTasks$method[i])
#                 p[[1]] <- .self$obj[[.self$runTasks$species[i]]][[.self$runTasks$method[i]]][[.self$runTasks$mIDChar[i]]]@object
#                 
#                 if (is.null(.self$modelFrame$specis_specific)) p[[2]] <- .self$modelFrame$features
#                 else p[[2]] <- cbind(.self$modelFrame$features,.self$modelFrame$specis_specific[[.self$runTasks$species[i]]])
#                 
#                 m <- try(.self$funs[[.self$runTasks$method[i]]](p),silent=TRUE)
#                 options(warn=0)
#                 m
#               },
#               getFeatures=function(sp) {
#                 if (!is.null(.self$modelFrame$species_specific)) .self$modelFrame$features
#                 else cbind(.self$modelFrame$features,.self$modelFrame$species_specific[[sp]])
#               },
#               generateParams=function(n,sp=NULL) {
#                 if (n == 'sdmDataFrame') {
#                   if (!is.null(.self$modelFrame$specis_specific)) {
#                     if (is.null(sp)) stop('species should be specified!')
#                     cbind(.self$modelFrame$features,.self$modelFrame$specis_specific[[sp]])
#                   } else .self$modelFrame$features
#                 } else if (n %in% c('newdata','data.frame')) {
#                   .self$newdata$data.frame
#                 } else if (n == 'standard.formula') {
#                   .getFormula(c(sp,colnames(.self$generateParams('sdmDataFrame',sp))),env=parent.frame(2))
#                 } else if (n == 'sdmX') {
#                   #####
#                   
#                 } else if (n == 'sdmY') {
#                   #####
#                   
#                 } else if (n == 'sdmRaster') {
#                   #####
#                   
#                 } else if (n %in% names(.self$params)) {
#                   do.call(.self$params[[n]],generateParams(.CharVector2List(names(formals(.self$params[[n]]))),sp)) 
#                 }
#               },
#               getReseved.names=function() {
#                 c('sdmDataFrame','sdmX','sdmY','sdmRaster','sdmVariables','standard.formula','gam.mgcv.furmula')
#               },
#               getPredictArgs=function(sp,mo) {
#                 # return a list in which the first element is reserved for 'model'
#                 # and the second element is reserved for data e.g., 'newdata'
#                 # these two elements will be updated before putting in the predict function
#                 o <- list()
#                 pa <- .self$arguments[[mo]]$params
#                 n <- names(pa)
#                 ww <- which(pa == 'model')
#                 o[[n[ww]]] <- 'model'
#                 n <- n[-ww]
#                 pa<- pa[-ww]
#                 ww <- which(names(pa) %in% .self$dataObject.names)
#                 if (length(ww) == 0) stop('data object required by the predict function is not recognised!')
#                 o[[n[ww]]] <- pa[[ww]]
#                 n <- n[-ww]
#                 pa<- pa[-ww]
#                 
#                 if (length(n) > 0) {
#                   for (nn in n) o[[nn]] <- .self$generateParams(pa[[nn]],sp)
#                 }
#                 
#                 o <- c(o,.self$arguments[[mo]]$settings)
#                 o
#               },
#               setRules=function() {
#                 # check if any rule is defined as a function for each method,
#                 # run the function to change the setting
#                 
#               }
#             )
# )


setRefClass(".workloadPredict",
            fields=list(
              obj='list', # list of models (@models) from sdmModels
              featFrame='.featureFrame',
              params='list',
              arguments='list',
              dataObject.names='character',
              funs='list',
              settingRules='list',
              runTasks='data.frame',
              parallelSetting='.parallelSettingORnull'
            ),
            methods=list(
              predictMID=function(IDs,.frame) {
                options(warn=-1)
                IDs <- which(.self$runTasks$modelID %in% IDs)
                m <- lapply(IDs,function(i) {
                  p <- .self$getPredictArgs(.self$runTasks$species[i],.self$runTasks$method[i],.frame)
                  p[[1]] <- .self$obj[[.self$runTasks$species[i]]][[.self$runTasks$method[i]]][[.self$runTasks$mIDChar[i]]]@object
                  
                  p[[2]] <- .frame
                  
                  m <- try(.self$funs[[.self$runTasks$method[i]]](p),silent=TRUE)
                  if (!inherits(m,'try-error')) as.numeric(m)
                  else m
                })
                options(warn=0)
                m
              },
              getSdmVariables=function(sp=NULL,nFact=NULL,nf=NULL) {
                # need to be revised (in predict, we need featureFrame, so it needs to be adjusted!)
                if (missing(nFact) || is.null(nFact)) {
                  nFact <- names(.self$featFrame@categorical)
                }
                
                if (missing(nf) || is.null(nFact)) {
                  nf <- .self$featFrame@numeric$names
                }
                
                new('.sdmVariables',response=sp,variables=list(numeric=nf,factors=nFact),distribution="binomial",features.numeric=nf,features.factor=nFact,
                    number.of.records=0,
                    n.presence=0,varInfo=list(numeric=.self$featFrame@numeric,categorical=.self$featFrame@categorical))
                
              },
              predictID=function(i,.frame) {
                options(warn=-1)
                i <- which(.self$runTasks$modelID == i)
                p <- .self$getPredictArgs(.self$runTasks$species[i],.self$runTasks$method[i],.frame)
                p[[1]] <- .self$obj[[.self$runTasks$species[i]]][[.self$runTasks$method[i]]][[.self$runTasks$mIDChar[i]]]@object
                
                p[[2]] <- .frame
                
                m <- try(.self$funs[[.self$runTasks$method[i]]](p),silent=TRUE)
                options(warn=0)
                if (!inherits(m,'try-error')) as.numeric(m)
                else m
              },
              generateParams=function(n,sp=NULL,.frame=NULL) {
                if (n == 'sdmDataFrame') {
                  .frame
                } else if (n %in% c('newdata','data.frame')) {
                  .frame
                } else if (n == 'standard.formula') {
                  .getFormula(c(sp,colnames(.self$generateParams('sdmDataFrame',sp,.frame))),env=parent.frame(2))
                } else if (n == 'sdmVariables') {
                  getSdmVariables(sp)
                } else if (n[[i]] == 'sdmX') {
                  .frame[,colnames(.frame) != sp,drop=FALSE]
                } else if (n[[i]] == 'sdmY') {
                  .frame[,colnames(.frame) == sp]
                } else if (n[[i]] == 'sdmX.norm') {
                  
                  .normalize(.frame[,colnames(.frame) != sp,drop=FALSE],frame=.self$featFrame@numeric)
                  
                } else if (n[[i]] == 'sdmDataFrame.norm') {
                  .normalize(.frame,except=sp,frame=.self$featFrame@numeric)
                  
                }  else if (n[[i]] == 'sdmDataFrame.scaled') {
                  .normalize(.frame,except=sp,frame=.self$featFrame@numeric,scale=TRUE)
                } else if (n[[i]] == 'sdmMatrix') {
                  .f <- .self$generateParams('standard.formula',sp,.frame)
                  model.matrix(.f,.frame)[,-1]
                  
                } else if (n[[i]] == 'sdmMatrix.norm') {
                  .f <- .self$generateParams('standard.formula',sp,.frame)
                  model.matrix(.f,.self$generateParams('sdmDataFrame.norm',sp,.frame))[,-1]
                  
                } else if (n[[i]] == 'sdmMatrix.scaled') {
                  .f <- .self$generateParams('standard.formula',sp,.frame)
                  model.matrix(.f,.self$generateParams('sdmDataFrame.scaled',sp,.frame))[,-1]
                  
                } else if (n[[i]] == 'sdmdata') {
                  .self$obj@data
                } else if (n %in% names(.self$params)) {
                  do.call(.self$params[[n]],generateParams(.CharVector2List(names(formals(.self$params[[n]]))),sp)) 
                }
              },
              getReseved.names=function() {
                c('sdmDataFrame','sdmX','sdmY','sdmRaster','sdmVariables','standard.formula','gam.mgcv.furmula','sdmMatrix.scaled',
                  'sdmMatrix.norm','sdmDataFrame.scaled','sdmDataFrame.norm')
              },
              getPredictArgs=function(sp,mo,.frame=NULL) {
                # return a list in which the first element is reserved for 'model'
                # and the second element is reserved for data e.g., 'newdata'
                # these two elements will be updated before putting in the predict function
                o <- list()
                pa <- .self$arguments[[mo]]$params
                n <- names(pa)
                ww <- which(pa == 'model')
                o[[n[ww]]] <- 'model'
                n <- n[-ww]
                pa<- pa[-ww]
                ww <- which(names(pa) %in% .self$dataObject.names)
                if (length(ww) == 0) stop('data object required by the predict function is not recognised!')
                o[[n[ww]]] <- pa[[ww]]
                n <- n[-ww]
                pa<- pa[-ww]
                
                if (length(n) > 0) {
                  for (nn in n) o[[nn]] <- .self$generateParams(pa[[nn]],sp,.frame)
                }
                
                o <- c(o,.self$arguments[[mo]]$settings)
                o
              },
              setRules=function() {
                # check if any rule is defined as a function for each method,
                # run the function to change the setting
                
              }
            )
)


setRefClass(".workloadP",
            fields=list(
              obj='list', # list of models (@models) from sdmModels
              #data='data.frame', # data used to fit the models
              newdata='list',
              modelFrame='data.frameORnull',
              params='list',
              arguments='list',
              dataObject.names='character',
              funs='list',
              settingRules='list',
              runTasks='data.frame',
              ncore='numericORnull'
            ),
            methods=list(
              predictMID=function(IDs) {
                options(warn=-1)
                IDs <- which(.self$runTasks$modelID %in% IDs)
                m <- lapply(IDs,function(i) {
                  p <- .self$getPredictArgs(.self$runTasks$species[i],.self$runTasks$method[i])
                  p[[1]] <- .self$obj[[.self$runTasks$species[i]]][[.self$runTasks$method[i]]][[.self$runTasks$mIDChar[i]]]@object
                  
                  p[[2]] <- .self$modelFrame
                  
                  m <- try(.self$funs[[.self$runTasks$method[i]]](p),silent=TRUE)
                  m
                })
                options(warn=0)
                m
              },
              predictID=function(i) {
                options(warn=-1)
                i <- which(.self$runTasks$modelID == i)
                p <- .self$getPredictArgs(.self$runTasks$species[i],.self$runTasks$method[i])
                p[[1]] <- .self$obj[[.self$runTasks$species[i]]][[.self$runTasks$method[i]]][[.self$runTasks$mIDChar[i]]]@object
                
                p[[2]] <- .self$modelFrame
                
                m <- try(.self$funs[[.self$runTasks$method[i]]](p),silent=TRUE)
                options(warn=0)
                m
              },
              getFeatures=function(sp) {
                .self$modelFrame
                
              },
              generateParams=function(n,sp=NULL) {
                if (n == 'sdmDataFrame') {
                  .self$modelFrame
                } else if (n %in% c('newdata','data.frame')) {
                  .self$newdata$data.frame
                } else if (n == 'standard.formula') {
                  .getFormula(c(sp,colnames(.self$generateParams('sdmDataFrame',sp))),env=parent.frame(2))
                } else if (n == 'sdmX') {
                  #####
                  
                } else if (n == 'sdmY') {
                  #####
                  
                } else if (n == 'sdmRaster') {
                  #####
                  
                } else if (n %in% names(.self$params)) {
                  do.call(.self$params[[n]],generateParams(.CharVector2List(names(formals(.self$params[[n]]))),sp)) 
                }
              },
              getReseved.names=function() {
                c('sdmDataFrame','sdmX','sdmY','sdmRaster','sdmVariables','standard.formula','gam.mgcv.furmula')
              },
              getPredictArgs=function(sp,mo) {
                # return a list in which the first element is reserved for 'model'
                # and the second element is reserved for data e.g., 'newdata'
                # these two elements will be updated before putting in the predict function
                o <- list()
                pa <- .self$arguments[[mo]]$params
                n <- names(pa)
                ww <- which(pa == 'model')
                o[[n[ww]]] <- 'model'
                n <- n[-ww]
                pa<- pa[-ww]
                ww <- which(names(pa) %in% .self$dataObject.names)
                if (length(ww) == 0) stop('data object required by the predict function is not recognised!')
                o[[n[ww]]] <- pa[[ww]]
                n <- n[-ww]
                pa<- pa[-ww]
                
                if (length(n) > 0) {
                  for (nn in n) o[[nn]] <- .self$generateParams(pa[[nn]],sp)
                }
                
                o <- c(o,.self$arguments[[mo]]$settings)
                o
              },
              setRules=function() {
                # check if any rule is defined as a function for each method,
                # run the function to change the setting
                
              }
            )
)
#-------------
#-------------
setClass("sdmEvaluate",
         representation(
           observed='numeric',
           predicted='numeric',
           statistics='list',
           threshold_based='data.frameORnull'
         )
)
#----------


setClass(".maxlikeModel",
         representation(
           fit='list'
         )
)

#----------

setRefClass(".sdmOptions",
            fields=list(
              options='list'
            ),
            methods=list(
              addOption=function(n,v) {
                .self$options[[n]] <- v
              },
              getOption=function(n) {
                .self$options[[n]]
              },
              getOptions=function() {
                .self$options
              },
              deleteOption=function(n) {
                if (n %in% names(.self$options)) {
                  .self$options <- .self$options[names(.self$options) != n]
                }
              }
            )
)
.sdmOptions <- new('.sdmOptions')

#-------

setClass(".sdmCalibration",
         representation(
           statistic='numeric',
           calibration='data.frame'
         )
)

#--------
setClass(".varImportance",
         representation(
           variables='character',
           varImportance='data.frame'
         )
)
#--------
setClass(".varImportanceList",
         representation(
           variables='character',
           varImportanceList='list',
           varImportanceMean='list'
         )
)
#-----------

setClass('.bioclimModel',
         representation(
           features='character',
           min='numeric',
           max='numeric',
           median='numeric',
           q25='numeric',
           q75='numeric',
           c='numeric',
           weights='numericORnull'
         )
)
#--------
setClass(".responseCurve",
         representation(
           variables='character',
           categorical='characterORnull',
           multi='logical',
           response='list'
         )
)
#------

setClass(".nicheRaster",
         representation(
           names='characterORnull',
           nicheRaster='RasterLayer',
           scaleParams='data.frame'
         )
)
#----------
setClass(".nicheSpatRaster",
         representation(
           names='characterORnull',
           nicheRaster='SpatRaster',
           scaleParams='data.frame'
         )
)
#-----
setClass(".envSpace",
         representation(
           names='character',
           coords='data.frameORnull',
           scaledVariables='list',
           scaleParams='data.frame'
         )
)
#----------