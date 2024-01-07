# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  Jan 2024
# Version 5.1
# Licence GPL v3
#--------


.methodFix <- function(n) {
  .addMethods()
  for (i in seq_along(n)) {
    nx <- .sdmMethods$whichMethod(n[i])
    if (!is.null(nx)) n[i] <- nx
    else n[i] <- NA
  }
  n
}
#----------
.replicate.methodFix <- function(n) {
  .addMethods()
  for (i in seq_along(n)) {
    nx <- .replicateMethods$whichMethod(n[i])
    if (!is.null(nx)) n[i] <- nx
    else n[i] <- NA
  }
  n
}
#----------
.getSpeciesDistribution <- function(data,sp=NULL) {
  if (!is.null(sp)) sp <- data@species[sp]
  else sp <- data@species
  
  o <- lapply(sp,function(x) {
    if (!is.null(x@presence)) return('binomial')
    else if (!is.null(x@abundance)) return('poisson')
    else if (!is.null(x@Multinomial)) return('multinomial')
    else return(NA)
  })
  n <- names(o)
  o <- as.character(o)
  names(o) <- n
  o
}
#-------------

.getFormula.rhs <- function(n,env=parent.frame()) {
  as.formula(paste('~',paste(n,collapse='+'),sep=''),env = env)
}

.getFormula <- function(n,env=parent.frame()) {
  as.formula(paste(n[1],'~',paste(n[-1],collapse='+'),sep=''),env = env)
}

.getFormula.glmPoly <- function(n,nFact=NULL,degree=3,env=parent.frame()) {
  if (!is.null(nFact)) as.formula(paste(n[1],'~',paste(c(paste(paste0('poly(',n[-1],sep=''),', degree = ',degree,')',sep=''),nFact),collapse='+'),sep=''),env = env)
  else as.formula(paste(n[1],'~',paste(paste(paste('poly(',n[-1],sep=''),', degree = ',degree,')',sep=''),collapse='+'),sep=''),env = env)
}


.getFormula.gammgcv.rhs <- function(n,nFact=NULL,k=-1,bs='tp',env=parent.frame()) {
  # if (!is.null(nFact)) as.formula(paste('~',paste(c(paste(paste('s(',n,sep=''),')',sep=''),nFact),collapse='+'),sep=''),env = env)
  # else as.formula(paste('~',paste(paste(paste('s(',n,sep=''),')',sep=''),collapse='+'),sep=''),env = env)
  
  if (length(k) > 1) {
    if (length(k) != length(n)) stop('the number of items in k (smothing parameter for the GAM formula) should be either one or equal to the number of continuous variables!')
  } 
  
  if (!is.null(nFact)) as.formula(paste('~',paste(c(paste(paste0('s(',n,sep=''),', k = ',k,', bs = "',bs,'")',sep=''),nFact),collapse='+'),sep=''),env = env)
  else as.formula(paste('~',paste(paste(paste('s(',n,sep=''),', k = ',k,', bs = "',bs,'")',sep=''),collapse='+'),sep=''),env = env)
  
}

.getFormula.gammgcv <- function(n,nFact=NULL,k=-1,bs='tp',env=parent.frame()) {
  if (length(k) > 1) {
    if (length(k) != length(n[-1])) stop('the number of items in k (smothing parameter for the GAM formula) should be either one or equal to the number of continuous variables!')
  } 
  
  if (!is.null(nFact)) as.formula(paste(n[1],'~',paste(c(paste(paste0('s(',n[-1],sep=''),', k = ',k,', bs = "',bs,'")',sep=''),nFact),collapse='+'),sep=''),env = env)
  else as.formula(paste(n[1],'~',paste(paste(paste('s(',n[-1],sep=''),', k = ',k,', bs = "',bs,'")',sep=''),collapse='+'),sep=''),env = env)
  
}
#----------

.addLHSformula <- function(f,lhs,env=parent.frame()) {
  # ~ x1 + x2; add lhs to such formula
  as.formula(paste(lhs,'~',as.character(f)[2]),env = env )
} 


#----------------

.mahal <- function(d1,d2) {
  co <- solve(cov(d1))
  mahalanobis(d2,colMeans(d1,na.rm=TRUE),co,inverted=TRUE)
}
#----------
.checkFactor <- function(f1,f2) {
  f2 <- factor(f2)
  f1 <- factor(f1)
  l <- levels(f2)[!levels(f2) %in% levels(f1)]
  ww <- NULL
  if (length(l) > 0) {
    w1 <- w2 <- ww <- list()
    for (i in seq_along(l)) w1[[l[i]]] <- which(f2 == l[i])
    l2 <- levels(f2)[levels(f2) %in% levels(f1)]
    for (i in seq_along(l2)) w2[[l2[i]]] <- which(f2 == l2[i])
    ww[['p']] <- w1
    ww[['np']] <- w2
  }
  ww
}


#----------------
.eqFactLevel <- function(data1,data2) {
  nFact <- colnames(data2)[.where(is.factor,data2)]
  for (i in seq_along(nFact)) data2[,nFact[i]] <- factor(data2[,nFact[i]],levels = levels(data1[,nFact[i]]))
  data2
}
#-------------


.is.windows <- function() {
  s <- Sys.info()
  if (!is.null(s)) s[['sysname']] == 'Windows'
  else FALSE
}
#---------
.getRunID <- function(ri,sp,m) {
  # from sdmModels, extract the modelID and runID (replicates) to be assigned to the model objects generated through fitting
  w1 <- ri[,2] == sp
  w2 <- ri[,3] == m
  w1 <- w1 & w2
  list(mID=ri[w1,1],rID=ri[w1,5])
}

#-----------

.getRecordID <- function(x,sp,id,train) {
  # x is recordID list
  # it finds the record ID of observation by the rowID used to generate train or dependent test, or in independent test
  # train = FALSE means id belongs to independent test
  if (train) {
    x[[sp]]$train$rID[x[[sp]]$train$rowID %in% id]
  } else {
    x[[sp]]$test$rID[x[[sp]]$test$rowID %in% id]
  }
}
#--------------
# convert a character vector to a list and use the values as the names of the list item unless the names are specified through n
.charVector2List <- function(x,n=NULL) {
      xx <- as.list(x)
      if (!is.null(n) && length(n) == length(x)) names(xx) <- n
      else names(xx) <- x
      xx
}
#------------
.pkgLoad <- function(m) {
  # m is method names (this loads the packages that are required for the methods)
  # written to be used for evaluating in clusters (parallelisation...)
  pkgs <- .sdmMethods$getPackageNames(m)
  .sdm...temp <- NULL; rm(.sdm...temp)
  pos <- 1
  
  w <- sapply(pkgs, function(x) '.temp' %in% x)
  
  if (any(w)) {
    if (".sdm...temp" %in% ls(pattern='^.sdm..',pos=1,all.names = TRUE)) {
      ww <- ls(.sdmMethods$userFunctions)
      rm(list=ww,pos=1)
      rm(.sdm...temp,pos=1)
    }
    #----
    for (i in which(w)) pkgs[[i]] <- pkgs[[i]][which(pkgs[[i]] != '.temp')]
    #----
    w <- ls(.sdmMethods$userFunctions)
    
    if (length(w) > 0) {
      assign('.sdm...temp',c(),envir = as.environment(pos))
      for (ww in w) {
        if (!exists(ww,where=1)) {
          assign(ww,.sdmMethods$userFunctions[[ww]],envir = as.environment(pos))
          .sdm...temp <<- c(.sdm...temp,ww)
        }
      }
    }
  }
  
  ww <- .loadLib(pkgs)
  all(ww)
}

#----------
.generateWL <- function(d,s,filename=NULL) {
  pkgs <- .sdmMethods$getPackageNames(s@methods)
  .sdm...temp <- NULL; rm(.sdm...temp)
  pos <- 1
  
  w <- sapply(pkgs, function(x) '.temp' %in% x)
  
  if (any(w)) {
    if (".sdm...temp" %in% ls(pattern='^.sdm..',pos=1,all.names = TRUE)) {
      ww <- ls(.sdmMethods$userFunctions)
      rm(list=ww,pos=1)
      rm(.sdm...temp,pos=1)
    }
    #----
    for (i in which(w)) pkgs[[i]] <- pkgs[[i]][which(pkgs[[i]] != '.temp')]
    #----
    w <- ls(.sdmMethods$userFunctions)
    
    if (length(w) > 0) {
      assign('.sdm...temp',c(),envir = as.environment(pos))
      for (ww in w) {
        if (!exists(ww,where=1)) {
          assign(ww,.sdmMethods$userFunctions[[ww]],envir = as.environment(pos))
          .sdm...temp <<- c(.sdm...temp,ww)
        }
      }
    }
  }
  
  ww <- .loadLib(pkgs)
  if (!all(ww)) {
    if (!any(ww)) {
      cat('Some methods are ignored as they depend on package(s) which is/are not installed on your machine!\n')
      cat('you can use the installAll() function to  install all packages required by methods in sdm\n')
      stop(paste('The packages rquired by the selected methods are not installed -> package names:',paste(unlist(pkgs),collapse=', ')))
    } else {
      cat('Some methods are ignored as they depend on package(s) which is/are not installed on your machine!\n')
      cat('you can use the installAll() function to  install all packages required by methods in sdm\n')
      warning(paste('The packages rquired by the following methods are not installed:',paste(s@methods[!ww],collapse=', '),'; These methods are excluded! The packages need to be installed for these methods:',paste(unlist(pkgs[!ww]),collapse=', ')))
      s@methods <- s@methods[ww]
    }
  }
  #------------
  if (!is.null(s@seed)) {
    if (is.numeric(s@seed)) set.seed(s@seed[1])
    else if (is.logical(s@seed) && s@seed) {
      s@seed <- sample(1000000,1)
      set.seed(s@seed)
    }
  }
  #-----------
  if (length(s@parallelSetting@ncore) == 0) s@parallelSetting@ncore <- 1
  
  w <- new('.workload',ncore=s@parallelSetting@ncore,data=d,setting=s,frame=s@featureFrame,filename=filename)
  #-----
  hasTest <- 'test' %in% d@groups$training@values[,2]
  nFact <- NULL
  if (!is.null(d@factors) > 0 && any(d@factors %in% s@sdmFormula@vars@names)) nFact <- d@factors[d@factors %in% s@sdmFormula@vars@names]
  nf <- .excludeVector(d@features.name,nFact)
  nf <- nf[nf %in% s@sdmFormula@vars@names]
  nFact <- nFact[nFact %in% s@sdmFormula@vars@names]
  #---------
  for (sp in s@sdmFormula@vars@species) {
    dt <- as.data.frame(d,sp=sp,grp='train')
    w$recordIDs[[sp]]$train <- data.frame(rID=dt[,1],rowID=1:nrow(dt))
    f <- w$frame@featureGenerator(dt)
    w$train[[sp]]$sdmDataFrame  <- cbind(dt[,sp,drop=FALSE],f)
    
    if (hasTest) {
      dt <- as.data.frame(d,sp=sp,grp='test')
      w$recordIDs[[sp]]$test <- data.frame(rID=dt[,1],rowID=1:nrow(dt))
      f <- w$frame@featureGenerator(dt)
      w$test[[sp]]$sdmDataFrame <- cbind(dt[,sp,drop=FALSE],f)
    }
    #--------
    w$getSdmVariables(sp=sp,nf=nf,nFact = nFact)
  }
  #------------
  w$funs[['fit']] <- .sdmMethods$getFitFunctions(s@methods)
  w$arguments[['fit']] <- .sdmMethods$getFitArguments(s@methods)
  w$funs[['predict']] <- .sdmMethods$getPredictFunctions(s@methods)
  w$arguments[['predict']] <- .sdmMethods$getPredictArguments(s@methods)
  #---
  mo <- s@methods
  names(mo) <- mo
  w$dataObject.names <- lapply(mo, .sdmMethods$getDataArgumentNames)
  w$settingRules <- lapply(mo,function(x) .sdmMethods$Methods[[x]]@settingRules)
  #-----------
  
  # .w <- sapply(w$settingRules,is.null)
  # if (any(!.w)) {
  #   .w <- names(.w[!.w])
  #   for (sp in s@sdmFormula@vars@species) {
  #     for (mo in .w) {
  #       f <- w$settingRules[[mo]](w$sdmVariables[[sp]],w$arguments$fit[[mo]]$settings)
  #       if ('fitSettings' %in% names(f)) w$arguments$fit[[mo]]$settings <- f$fitSettings
  #       
  #       if ('predictSettings' %in% names(f)) w$arguments$predict[[mo]]$settings <- f$predictSettings
  #       
  #     }
  #   }
  # }
  #------ applying the setting rules & overriding the user settings:
  for (mo in w$setting@methods) {
    .u <- NULL
    if (!is.null(w$setting@modelSettings) && mo %in% names(w$setting@modelSettings)) {
      .u <- w$setting@modelSettings[[mo]]
    }
    for (sp in names(w$train)) {
      w$arguments[['overriden_settings']][['fit']][[mo]][[sp]] <- w$arguments$fit[[mo]]$settings
      w$arguments[['overriden_settings']][['predict']][[mo]][[sp]] <- w$arguments$predict[[mo]]$settings
      #-------
      .set <- w$setRules(mo,sp)
      if (!is.null(.set)) {
        if ('fitSettings' %in% names(.set)) {
          w$arguments[['overriden_settings']][['fit']][[mo]][[sp]][names(.set[['fitSettings']])] <- .set[['fitSettings']]
        }
        
        if ('predictSettings' %in% names(.set)) {
          w$arguments[['overriden_settings']][['predict']][[mo]][[sp]][names(.set[['predictSettings']])] <- .set[['predictSettings']]
        }
      }
      #-------
      if (!is.null(.u)) {
        if (any(names(.u) %in% c('fitSettings','fitSetting','fitsetting','fitset','FitSetting','FitSettings','predictSettings','predictSetting','predictsetting','predictset','PredictSetting','PredictSettings'))) {
          .w <- which(names(.u) %in% c('fitSettings','fitSetting','fitsetting','fitset','FitSetting','FitSettings'))
          if (length(.w) == 1) {
            if (is.list(.u[[.w]])) {
              .ww <- names(.u[[.w]]) %in% names(w$arguments[['overriden_settings']][['fit']][[mo]][[sp]])
              if (any(.ww)) {
                w$arguments[['overriden_settings']][['fit']][[mo]][[sp]][names(.u[[.w]])[.ww]] <- .u[[.w]][.ww]
                if (any(!.ww)) w$arguments[['overriden_settings']][['fit']][[mo]][[sp]] <- c(w$arguments[['overriden_settings']][['fit']][[mo]][[sp]],.u[[.w]][!.ww])
              } else w$arguments[['overriden_settings']][['fit']][[mo]][[sp]] <- c(w$arguments[['overriden_settings']][['fit']][[mo]][[sp]],.u[[.w]])
            }
            .u <- .u[-.w]
          }
          #------
          .w <- which(names(.u) %in% c('predictSettings','predictSetting','predictsetting','predictset','PredictSetting','PredictSettings'))
          if (length(.w) == 1) {
            if (is.list(.u[[.w]])) {
              .ww <- names(.u[[.w]]) %in% names(w$arguments[['overriden_settings']][['predict']][[mo]][[sp]])
              if (any(.ww)) {
                w$arguments[['overriden_settings']][['predict']][[mo]][[sp]][names(.u[[.w]])[.ww]] <- .u[[.w]][.ww]
                if (any(!.ww)) w$arguments[['overriden_settings']][['predict']][[mo]][[sp]] <- c(w$arguments[['overriden_settings']][['predict']][[mo]][[sp]],.u[[.w]][!.ww])
              } else w$arguments[['overriden_settings']][['predict']][[mo]][[sp]] <- c(w$arguments[['overriden_settings']][['predict']][[mo]][[sp]],.u[[.w]])
            }
            .u <- .u[-.w]
          }
        }
        #------
        if (length(.u) > 0) {
          
          .ww <- names(.u) %in% names(w$arguments[['overriden_settings']][['fit']][[mo]][[sp]])
          if (any(.ww)) {
            w$arguments[['overriden_settings']][['fit']][[mo]][[sp]][names(.u)[.ww]] <- .u[.ww]
            if (any(!.ww)) w$arguments[['overriden_settings']][['fit']][[mo]][[sp]] <- c(w$arguments[['overriden_settings']][['fit']][[mo]][[sp]],.u[!.ww])
          } else w$arguments[['overriden_settings']][['fit']][[mo]][[sp]] <- c(w$arguments[['overriden_settings']][['fit']][[mo]][[sp]],.u)
        }
      }
    }
  }
  
  #-------------
  
 
  for (mo in s@methods) {
    wc <- unlist(lapply(w$arguments$fit[[mo]]$params,function(x) is.character(x)))
    if (any(!wc)) {
      if (!all(unlist(lapply(w$arguments$fit[[mo]]$params[!wc],function(x) is.function(x))))) {
        warning(paste('parameter definition for the model',mo,'in the model container is not correctly defined!'))
        cat('\nmethod',mo,'is removed because the definition for the fit parameters',paste(names(w$arguments$fit[[mo]]$params)[!wc],collapse=', '),' is not correct!\n')
        w$setting@methods <- .excludeVector(w$setting@methods,mo)
        w$arguments$fit <- w$arguments$fit[-which(names(w$arguments$fit) == mo)]
        w$arguments$predict <- w$arguments$predict[-which(names(w$arguments$predict) == mo)]
        w$funs$fit <- w$funs$fit[-which(names(w$funs$fit) == mo)]
        w$funs$predict <- w$funs$predict[-which(names(w$funs$predict) == mo)]
        w$dataObject.names <- w$dataObject.names[-which(names(w$dataObject.names) == mo)]
      } else {
        for (n in names(w$arguments$fit[[mo]]$params)[!wc]) {
          .nn <- paste0(n,'_',mo,'_fit') # a specific name is generated and added to w$params to which the function is assigned
          # the assigned function will be called when the w$generateParams is executed to generate the customised parameters defined through the function in the model definition
          w$params[[.nn]] <- w$arguments$fit[[mo]]$params[[n]]
          w$arguments$fit[[mo]]$params[[n]] <- .nn
        }
      }
    }
    #-----
    if (mo %in% w$setting@methods) {
      wc <- unlist(lapply(w$arguments$predict[[mo]]$params,function(x) is.character(x)))
      if (any(!wc)) {
        if (!all(unlist(lapply(w$arguments$predict[[mo]]$params[!wc],function(x) is.function(x))))) {
          warning(paste('parameter definition for the model',mo,'in the model container is not correctly defined!'))
          cat('\nmethod',mo,'is removed because the definition for the predict parameters',paste(names(w$arguments$predict[[mo]]$params)[!wc],collapse=', '),' is not correct!\n')
          w$setting@methods <- .excludeVector(w$setting@methods,mo)
          w$arguments$fit <- w$arguments$fit[-which(names(w$arguments$fit) == mo)]
          w$arguments$predict <- w$arguments$predict[-which(names(w$arguments$predict) == mo)]
          w$funs$fit <- w$funs$fit[-which(names(w$funs$fit) == mo)]
          w$funs$predict <- w$funs$predict[-which(names(w$funs$predict) == mo)]
          w$dataObject.names <- w$dataObject.names[-which(names(w$dataObject.names) == mo)]
        } else {
          for (n in names(w$arguments$predict[[mo]]$params)[!wc]) {
            .nn <- paste0(n,'_',mo,'_predict') 
            w$params[[.nn]] <- w$arguments$predict[[mo]]$params[[n]]
            w$arguments$predict[[mo]]$params[[n]] <- .nn
          }
        }
      }
    }
  }
  #-----------------
  
  if (!is.null(s@replicate)) {
    f <- .replicateMethods$getFunctions(s@replicate)
    for (sp in names(w$train)) {
      if (d@species[[sp]]@type %in% c('Presence-Absence','Presence-Background')) family <- 'binomial'
      else if (d@species[[sp]]@type %in% c('Abundance')) family <- 'poisson'
      else if (d@species[[sp]]@type %in% c('Multinomial')) family <- 'Multinomial'
      else family <- 'Numerical'
      
      # sdmDataFrame! Leter should be checked for other types of data!
      for (ff in f) {
        w$replicates[[sp]] <- c(w$replicates[[sp]],ff(x=w$train[[sp]][['sdmDataFrame']][,sp],family=family,stratify=TRUE,test.percent=s@test.percentage,nfolds=s@cv.folds,replicates=s@n.replicates))
      }
    }
  }
  #-----
  w
}
#----------------
if (!isGeneric("sdm")) {
  setGeneric("sdm", function(formula,data,methods,...)
    standardGeneric("sdm"))
}

setMethod('sdm', signature(formula='ANY',data='sdmdata',methods='character'), 
          function(formula,data,methods,...) {
            a <- c('interaction.depth','n','replication','cv.folds','test.percent','bg','bg.n','var.importance','response.curve','var.selection','setting','ncore','modelSettings','seed','parallelSetting','filename')
            .sdm...temp <- NULL; rm(.sdm...temp)
            dot <- list(...)
            ndot <- names(dot)
            if (length(ndot) > 0) {
              ndot <- .pmatch(ndot,a)
              w <- !is.na(ndot)
              ndot <- ndot[w]
              dot <- dot[w]
              names(dot) <- ndot
            }
            
            if ('filename' %in% ndot) {
              w <- which(ndot == 'filename')
              filename <- dot[[w]]
              ndot <- ndot[-w]
              dot <- dot[-w]
            } else filename <- NULL
            
            if (missing(formula)) formula <- NULL
            
            dot$data <- data
            dot$formula <- formula
            dot$methods <- methods
            
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            
            s <- do.call('sdmSetting',dot)
            w <- .generateWL(data,s,filename = filename)
            w <- w$fit()
            #if (".sdmMethods$userFunctions" %in% search()) detach('.sdmMethods$userFunctions')
            if (".sdm...temp" %in% ls(pattern='^.sdm..',pos=1,all.names = TRUE)) {
              ww <- ls(.sdmMethods$userFunctions)
              rm(list=ww,pos=1)
              rm(.sdm...temp,pos=1)
            }
            w
          }
)

setMethod('sdm', signature(formula='ANY',data='sdmdata',methods='.sdmCorSetting'), 
          function(formula,data,methods,...) {
            a <- c('interaction.depth','n','replication','cv.folds','test.percent','bg','bg.n','var.importance','response.curve','var.selection','setting','ncore','modelSettings','seed','parallelSetting','filename')
            .sdm...temp <- NULL; rm(.sdm...temp)
            dot <- list(...)
            ndot <- names(dot)
            if (length(ndot) > 0) {
              ndot <- .pmatch(ndot,a)
              w <- !is.na(ndot)
              ndot <- ndot[w]
              dot <- dot[w]
              names(dot) <- ndot
            }
            
            if ('filename' %in% ndot) {
              w <- which(ndot == 'filename')
              filename <- dot[[w]]
              ndot <- ndot[-w]
              dot <- dot[-w]
            } else filename <- NULL
            
            s <- methods
            
            if (missing(formula) && !all(s@sdmFormula@vars@names %in% data@features.name)) s@sdmFormula <- data@sdmFormula
            else if (inherits(formula,'sdmFormula')) s@sdmFormula <- formula
            else if (inherits(formula,'formula')) {  
              s@sdmFormula <- .exFormula(formula,as.data.frame(data)[,-1])
            }
            
            if (length(dot) > 0) {
                if (length(ndot) > 0) {
                  for (nd in ndot) {
                    if (nd == 'interaction.depth') interaction.depth <- dot[[nd]]
                    else if (nd == 'ncore') s@parallelSetting@ncore <- dot[[nd]]
                    else if (nd == 'replication') s@replicate <- dot[[nd]]
                    else if (nd == 'cv.folds') s@cv.folds <- dot[[nd]]
                    else if (nd == 'test.percent') s@test.percentage <- dot[[nd]]
                    else if (nd == 'bg') s@pseudo.absence.methods <- dot[[nd]]
                    else if (nd == 'bg.n') s@n.pseudo.absence <- dot[[nd]]
                    else if (nd == 'var.importance') s@varImportance.methods <- dot[[nd]]
                    else if (nd == 'response.curve'&& is.logical(dot[[nd]])) s@response.curve <- dot[[nd]]
                    else if (nd == 'var.selection' && is.logical(dot[[nd]])) s@var.selection <- dot[[nd]]
                    else if (nd == 'modelSettings') s@modelSettings <- dot[[nd]]
                    else if (nd == 'seed') s@seed <- dot[[nd]]
                    else if (nd == 'parallelSetting' && is.list(dot[[nd]])) {
                      parallelSetting <- dot[[nd]]
                      nparallel <- names(parallelSetting)
                      a <- c('ncore','doParallel','method','cluster','hosts','fork','type')
                      nparallel <- .pmatch(nparallel,a)
                      w <- which(!is.na(nparallel))
                      if (length(w) > 0) {
                        parallelSetting <- parallelSetting[w]
                        nparallel <- nparallel[w]
                        names(parallelSetting) <- nparallel
                      }
                      #--
                      if ('ncore' %in% nparallel) s@parallelSetting@ncore <- parallelSetting$ncore
                      else s@parallelSetting@ncore <- max(c(floor(parallel::detectCores() * 0.5),1))
                      #--
                      if ('method' %in% nparallel) {
                        if (parallelSetting$method %in% c('parallel','foreach','future')) s@parallelSetting@method <- parallelSetting$method
                        else {
                          warning('parallelisation method is not recognised; the default value ("parallel") is used!')
                          s@parallelSetting@method <- 'parallel'
                        }
                      } else s@parallelSetting@method <- 'parallel'
                      #--
                      if ('strategy' %in% nparallel) {
                        parallelSetting$strategy <- tolower(parallelSetting$strategy)[1]
                        if (!parallelSetting$strategy %in% c('species','method','replicate','simple','auto')) {
                          warning('The parallel strategy is not recognised (should be one of c("auto","species","method","replicate","simple")); the default, "auto", is used!')
                          s@parallelSetting@strategy <- 'auto'
                        } else s@parallelSetting@strategy <- parallelSetting$strategy
                      } else s@parallelSetting@strategy <- 'auto'
                      #---
                      if ('fork' %in% nparallel) {
                        if (is.logical(parallelSetting$fork)) {
                          if (parallelSetting$fork && .is.windows()) {
                            warning('"fork" in parallelisation setting cannot be TRUE on Windows Operating Systems; It is changed to FALSE!')
                            s@parallelSetting@fork <- FALSE
                          } else s@parallelSetting@fork <- parallelSetting$fork
                        } else {
                          warning('"fork" in parallelisation setting should be logical; the default value is used!')
                          s@parallelSetting@fork <- !.is.windows()
                        }
                      } else s@parallelSetting@fork <- !.is.windows()
                      #--
                      if ('type' %in% nparallel) s@parallelSetting@type <- parallelSetting$type
                      #--
                      if ('doParallel' %in% nparallel && is.expression(parallelSetting$doParallel)) s@parallelSetting@doParallel <- parallelSetting$doParallel
                      #--
                      if ('cluster' %in% nparallel && inherits(parallelSetting$cluster,'cluster')) s@parallelSetting@cl <- parallelSetting$cluster
                      #--
                      if ('hosts' %in% nparallel && is.character(parallelSetting$hosts)) s@parallelSetting@hosts <- parallelSetting$hosts
                    }
                  }
                }
              }
            
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            
            w <- .generateWL(data,s, filename = filename)
            w <- w$fit(woL=w)
            #if (".sdmMethods$userFunctions" %in% search()) detach('.sdmMethods$userFunctions')
            if (".sdm...temp" %in% ls(pattern='^.sdm..',pos=1,all.names = TRUE)) {
              ww <- ls(.sdmMethods$userFunctions)
              rm(list=ww,pos=1)
              rm(.sdm...temp,pos=1)
            }
            w
          }
)


setMethod('sdm', signature(formula='sdmdata',data='.sdmCorSetting',methods='ANY'), 
          function(formula,data,methods,...) {
            a <- c('interaction.depth','n','replication','cv.folds','test.percent','bg','bg.n','var.importance','response.curve','var.selection','setting','ncore','modelSettings','seed','parallelSetting','filename')
            .sdm...temp <- NULL; rm(.sdm...temp)
            dot <- list(...)
            ndot <- names(dot)
            if (length(ndot) > 0) {
              ndot <- .pmatch(ndot,a)
              w <- !is.na(ndot)
              ndot <- ndot[w]
              dot <- dot[w]
              names(dot) <- ndot
            }
            
            if ('filename' %in% ndot) {
              w <- which(ndot == 'filename')
              filename <- dot[[w]]
              ndot <- ndot[-w]
              dot <- dot[-w]
            } else filename <- NULL
            
            s <- data
            d <- formula
            
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            
            if (!missing(methods) && is.character(methods)) {
              m <- .methodFix(methods)
              if (any(is.na(m))) warning(paste('methods',paste(methods[is.na(m)],collapse=', '),'do not exist!'))
              if (!all(is.na(m))) {
                m <- unique(m[!is.na(m)])
                s@methods <- m
              }
            }
            #-----------
            if (length(dot) > 0) {
                if (length(ndot) > 0) {
                  for (nd in ndot) {
                    if (nd == 'interaction.depth') interaction.depth <- dot[[nd]]
                    else if (nd == 'ncore') s@parallelSetting@ncore <- dot[[nd]]
                    else if (nd == 'replication') s@replicate <- dot[[nd]]
                    else if (nd == 'cv.folds') s@cv.folds <- dot[[nd]]
                    else if (nd == 'test.percent') s@test.percentage <- dot[[nd]]
                    else if (nd == 'bg') s@pseudo.absence.methods <- dot[[nd]]
                    else if (nd == 'bg.n') s@n.pseudo.absence <- dot[[nd]]
                    else if (nd == 'var.importance') s@varImportance.methods <- dot[[nd]]
                    else if (nd == 'response.curve'&& is.logical(dot[[nd]])) s@response.curve <- dot[[nd]]
                    else if (nd == 'var.selection' && is.logical(dot[[nd]])) s@var.selection <- dot[[nd]]
                    else if (nd == 'modelSettings') s@modelSettings <- dot[[nd]]
                    else if (nd == 'seed') s@seed <- dot[[nd]]
                    else if (nd == 'parallelSetting' && is.list(dot[[nd]])) {
                      parallelSetting <- dot[[nd]]
                      nparallel <- names(parallelSetting)
                      a <- c('ncore','doParallel','method','cluster','hosts','fork','type')
                      nparallel <- .pmatch(nparallel,a)
                      w <- which(!is.na(nparallel))
                      if (length(w) > 0) {
                        parallelSetting <- parallelSetting[w]
                        nparallel <- nparallel[w]
                        names(parallelSetting) <- nparallel
                      }
                      #--
                      if ('ncore' %in% nparallel) s@parallelSetting@ncore <- parallelSetting$ncore
                      else s@parallelSetting@ncore <- max(c(floor(parallel::detectCores() * 0.5),1))
                      
                      #--
                      if ('method' %in% nparallel) {
                        if (parallelSetting$method %in% c('parallel','foreach')) s@parallelSetting@method <- parallelSetting$method
                        else {
                          warning('parallelisation method is not recognised; the default value ("parallel") is used!')
                          s@parallelSetting@method <- 'parallel'
                        }
                      } else s@parallelSetting@method <- 'parallel'
                      #--
                      if ('strategy' %in% nparallel) {
                        parallelSetting$strategy <- tolower(parallelSetting$strategy)[1]
                        if (!parallelSetting$strategy %in% c('species','method','replicate','simple','auto')) {
                          warning('The parallel strategy is not recognised (should be one of c("auto","species","method","replicate","simple")); the default, "auto", is used!')
                          s@parallelSetting@strategy <- 'auto'
                        } else s@parallelSetting@strategy <- parallelSetting$strategy
                      } else s@parallelSetting@strategy <- 'auto'
                      #---
                      if ('fork' %in% nparallel) {
                        if (is.logical(parallelSetting$fork)) {
                          if (parallelSetting$fork && .is.windows()) {
                            warning('"fork" in parallelisation setting cannot be TRUE on Windows Operating Systems; It is changed to FALSE!')
                            s@parallelSetting@fork <- FALSE
                          } else s@parallelSetting@fork <- parallelSetting$fork
                        } else {
                          warning('"fork" in parallelisation setting should be logical; the default value is used!')
                          s@parallelSetting@fork <- !.is.windows()
                        }
                      } else s@parallelSetting@fork <- !.is.windows()
                      #--
                      if ('type' %in% nparallel) s@parallelSetting@type <- parallelSetting$type
                      #--
                      if ('doParallel' %in% nparallel && is.expression(parallelSetting$doParallel)) s@parallelSetting@doParallel <- parallelSetting$doParallel
                      #--
                      if ('cluster' %in% nparallel && inherits(parallelSetting$cluster,'cluster')) s@parallelSetting@cl <- parallelSetting$cluster
                      #--
                      if ('hosts' %in% nparallel && is.character(parallelSetting$hosts)) s@parallelSetting@hosts <- parallelSetting$hosts
                    }
                  }
                }
              }
            
            w <- .generateWL(d,s,filename = filename)
            w <- w$fit(woL=w)
            #if (".sdmMethods$userFunctions" %in% search()) detach('.sdmMethods$userFunctions')
            if (".sdm...temp" %in% ls(pattern='^.sdm..',pos=1,all.names = TRUE)) {
              ww <- ls(.sdmMethods$userFunctions)
              rm(list=ww,pos=1)
              rm(.sdm...temp,pos=1)
            }
            w
          }
)

#-------------
###########################
#####################################################################






.getRunTaskTable <- function(r,species,models,runs=NULL,nc=1) {
  ######- revised!
  # r is run.info table
  # returns a list of tasks for parallelising
  #---------------
  # for type=species: the parallelisation is on the first level of the list which is species name
  # for type=runs: the parallelisation is on the third level of the list which is replicateID
  # for type=methods: the parallelisation is on the second level of the list which is method
  # for type=species.noRun: the parallelisation is on the first level of the list which is species name
  # for type=methods.noRun: the parallelisation is on the second level of the list which is method
  # simple & simple.noRun is used when parallelisation should be running over the main task list including all tasks
  #################################
  o <- list()
  for (sp in species) {
    o[[sp]] <- list()
    for (i in models) o[[sp]][[i]] <- r[r$species == sp & r$method == i,]
  }
  if (!is.null(runs)) {
    if (length(species) > length(runs)) {
      if (length(species) >= nc) .type <- 'species'
      else .type <- 'simple'
    } else {
      if (length(runs) >= length(models)) {
        if (length(runs) >= nc) .type <- 'runs'
        else .type <- 'simple'
      } else {
        .type <- 'methods'
        if (length(models) >= nc) .type <- 'methods'
        else .type <- 'simple'
      }
    }
  } else {
    if (length(species) >= length(models)) {
      if (length(species) >= nc) .type <- 'species.noRun'
      else .type <- 'simple.noRun'
    } else {
      if (length(models) >= nc) .type <- 'methods.noRun'
      else .type <- 'simple.noRun'
    }
  }
  list(type=.type,tasks=o)
}

#################

.fit_and_3evals <- function(mID,method,sp,fit,fit.par,pred,pred.par,rID,n,dt,w) {
  # fit +  3 evaluations: training +  both dependent and independent tests!
  # n is the name of data argument
  options(warn=-1)
  mo <- new('.sdmCorModel',method=method,mID=mID,response=sp)
  
  fit.par[[n]] <- dt[w$replicates[[sp]][[rID]][[2]],]
  
  mo@object <- try(fit(fit.par),silent=TRUE)
  dtype <- pred.par[[2]]
  if (!inherits(mo@object, "try-error")) {
    pred.par[[1]] <-mo@object
    pred.par[[2]] <- fit.par[[n]]
    ev <- try(evaluates(fit.par[[n]][,sp],pred(pred.par),distribution=w$setting@distribution[sp]),silent=TRUE)
    if (!inherits(ev, "try-error")) mo@evaluation[['training']] <- ev
    else mo@errorLog[['evaluation']][['training']] <- ev
    
    vi <- try(._varImp(pred.par,pred=pred,sp=sp,nsim = 5),silent=TRUE)
    if (!inherits(vi, "try-error")) mo@varImportance[['training']] <- vi
    else mo@errorLog[['varImportance']][['training']] <- vi
    
    
    pred.par[[2]] <- dt[w$replicates[[sp]][[rID]][[3]],]
    ev <- try(evaluates(pred.par[[2]][,sp],pred(pred.par),distribution=w$setting@distribution[sp]),silent=TRUE)
    if (!inherits(ev, "try-error")) mo@evaluation[['test.dep']] <- ev
    else mo@errorLog[['evaluation']][['test.dep']] <- ev
    
    vi <- try(._varImp(pred.par,pred=pred,sp=sp,nsim = 5),silent=TRUE)
    if (!inherits(vi, "try-error")) mo@varImportance[['test.dep']] <- vi
    else mo@errorLog[['varImportance']][['test.dep']] <- vi
    
    
    pred.par[[2]] <- w$generateParams(list(dtype),sp,train=FALSE)[[1]]
    ev <- try(evaluates(pred.par[[2]][,sp],pred(pred.par),distribution=w$setting@distribution[sp]),silent=TRUE)
    if (!inherits(ev, "try-error")) mo@evaluation[['test.indep']] <- ev
    else mo@errorLog[['evaluation']][['test.indep']] <- ev
    
    vi <- try(._varImp(pred.par,pred=pred,sp=sp,nsim = 5),silent=TRUE)
    if (!inherits(vi, "try-error")) mo@varImportance[['test.indep']] <- vi
    else mo@errorLog[['varImportance']][['test.indep']] <- vi
  } else mo@errorLog[['fit']] <- mo@object
  options(warn=0)
  mo
}
#-----------

.fit_and_2evals_indep <- function(mID,method,sp,fit,fit.par,pred,pred.par,n,dt=dt,rID=NULL,w) {
  # fit +  2 evaluations: training +  independent tests!
  # n is the name of data argument (rID is not used, just to make it consistent with the other methods!)
  options(warn=-1)
  mo <- new('.sdmCorModel',method=method,mID=mID,response=sp)
  
  fit.par[[n]] <- dt
  
  mo@object <- try(fit(fit.par),silent=TRUE)
  dtype <- pred.par[[2]]
  if (!inherits(mo@object, "try-error")) {
    pred.par[[1]] <-mo@object
    pred.par[[2]] <- dt
    ev <- try(evaluates(dt[,sp],pred(pred.par),distribution=w$setting@distribution[sp]),silent=TRUE)
    if (!inherits(ev, "try-error")) mo@evaluation[['training']] <- ev
    else mo@errorLog[['evaluation']][['training']] <- ev
    
    vi <- try(._varImp(pred.par,pred=pred,sp=sp,nsim = 10),silent=TRUE)
    if (!inherits(vi, "try-error")) mo@varImportance[['training']] <- vi
    else mo@errorLog[['varImportance']][['training']] <- vi
    
    
    pred.par[[2]] <- w$generateParams(list(dtype),sp,train=FALSE)[[1]]
    ev <- try(evaluates(pred.par[[2]][,sp],pred(pred.par),distribution=w$setting@distribution[sp]),silent=TRUE)
    if (!inherits(ev, "try-error")) mo@evaluation[['test.indep']] <- ev
    else mo@errorLog[['evaluation']][['test.indep']] <- ev
    
    vi <- try(._varImp(pred.par,pred=pred,sp=sp,nsim = 10),silent=TRUE)
    if (!inherits(vi, "try-error")) mo@varImportance[['test.indep']] <- vi
    else mo@errorLog[['varImportance']][['test.indep']] <- vi
    
  } else mo@errorLog[['fit']] <- mo@object
  options(warn=0)
  mo
}
#--------
.fit_and_2evals_dep  <- function (mID,method,sp,fit,fit.par,pred,pred.par,rID,n,dt=dt,w) {
  # fit +  2 evaluations: training +  dependent tests!
  # n is the name of data argument
  options(warn=-1)
  mo <- new('.sdmCorModel',method=method,mID=mID,response=sp)
  
  fit.par[[n]] <- dt[w$replicates[[sp]][[rID]][[2]],]
  
  mo@object <- try(fit(fit.par),silent=TRUE)
  if (!inherits(mo@object, "try-error")) {
    pred.par[[1]] <-mo@object
    pred.par[[2]] <- fit.par[[n]]
    ev <- try(evaluates(fit.par[[n]][,sp],pred(pred.par),distribution=w$setting@distribution[sp]),silent=TRUE)
    if (!inherits(ev, "try-error")) mo@evaluation[['training']] <- ev
    else mo@errorLog[['evaluation']][['training']] <- ev
    
    vi <- try(._varImp(pred.par,pred=pred,sp=sp,nsim = 10),silent=TRUE)
    if (!inherits(vi, "try-error")) mo@varImportance[['training']] <- vi
    else mo@errorLog[['varImportance']][['training']] <- vi
    
    
    pred.par[[2]] <- dt[w$replicates[[sp]][[rID]][[3]],]
    ev <- try(evaluates(pred.par[[2]][,sp],pred(pred.par),distribution=w$setting@distribution[sp]),silent=TRUE)
    if (!inherits(ev, "try-error")) mo@evaluation[['test.dep']] <- ev
    else mo@errorLog[['evaluation']][['test.dep']] <- ev
    
    vi <- try(._varImp(pred.par,pred=pred,sp=sp,nsim = 10),silent=TRUE)
    if (!inherits(vi, "try-error")) mo@varImportance[['test.dep']] <- vi
    else mo@errorLog[['varImportance']][['test.dep']] <- vi
    
  } else mo@errorLog[['fit']] <- mo@object
  options(warn=0)
  mo
}
#----------
.fit_and_1eval <- function(mID,method,sp,fit,fit.par,pred,pred.par,n,dt=dt,rID=NULL,w) {
  # fit4: when no test is available (evaluation would be only on training data)
  options(warn=-1)
  mo <- new('.sdmCorModel',method=method,mID=mID,response=sp)
  
  fit.par[[n]] <- dt
  
  mo@object <- try(fit(fit.par),silent=TRUE)
  if (!inherits(mo@object, "try-error")) {
    pred.par[[1]] <-mo@object
    pred.par[[2]] <- dt
    ev <- try(evaluates(dt[,sp],pred(pred.par),distribution=w$setting@distribution[sp]),silent=TRUE)
    if (!inherits(ev, "try-error")) mo@evaluation[['training']] <- ev
    else mo@errorLog[['evaluation']][['training']] <- ev
    
    vi <- try(._varImp(pred.par,pred=pred,sp=sp,nsim = 10),silent=TRUE)
    if (!inherits(vi, "try-error")) mo@varImportance[['training']] <- vi
    else mo@errorLog[['varImportance']][['training']] <- vi
    
  } else mo@errorLog[['fit']] <- mo@object
  options(warn=0)
  mo
}
#------------

.fitSimpleLapply <- function(id, .Fun,w, cl,...) {
  # it works inside the .fit function given the workload (w), and some objects (e.g., .tasks) generated there!
  mID <- w$tasks[id,1]
  sp <- w$tasks[id,2]
  x <- w$tasks[id,3]
  f <- w$funs$fit[[x]]
  f.par <- w$getFitArgs(sp,x)
  n <- w$dataObject.names[[x]][['fit']]
  dt <- w$generateParams(f.par[n],sp)[[1]]
  p <- w$funs$predict[[x]]
  p.par <- w$getPredictArgs(sp,x)
  rID <- w$tasks[id,5]
  .Fun(mID=mID,sp = sp,method = x,fit = f,fit.par = f.par,pred = p,pred.par = p.par,rID = rID,n = n,dt=dt,w=w)
}



.fitLapply <- function(sp, .Fun, .l1=lapply, .l2=lapply,w,cl, ...) {
  # it works inside the .fit function given the workload (w), and some objects (e.g., .tasks) generated there!
  r <- w$tasks[[sp]]
  .m <- names(r)
  names(.m) <- .m
  .l1(X=.m, FUN=function(x,...) {
    .r <- r[[x]]
    f <- w$funs$fit[[x]]
    f.par <- w$getFitArgs(sp,x)
    n <- w$dataObject.names[[x]][['fit']]
    dt <- w$generateParams(f.par[n],sp)[[1]]
    p <- w$funs$predict[[x]]
    p.par <- w$getPredictArgs(sp,x)
    IDs <- as.list(.r[,c(1,5)])
    id <- seq_along(IDs$modelID)
    names(id) <- IDs$modelID
    
    .l2(X=id,FUN=function(i,...) {
      .Fun(mID=IDs$modelID[i],sp = sp,method = x,fit = f,fit.par = f.par,pred = p,pred.par = p.par,rID = IDs$replicationID[i],n = n,dt=dt,w=w)
    },f=f,f.par=f.par,p.par=p.par,n=n,IDs=IDs,dt=dt,w=w,...)
  },w=w,cl=cl)
}


.fitlapply <- function(sp, .Fun, w, ...) { # the same as .fitLapply, but .l1 and .l2 are both lapply, so removed,
  # it works inside the .fit function given the workload (w), and some objects (e.g., .tasks) generated there!
  r <- w$tasks[[sp]]
  .m <- names(r)
  names(.m) <- .m
  lapply(X=.m, FUN=function(x,...) {
    .r <- r[[x]]
    f <- w$funs$fit[[x]]
    f.par <- w$getFitArgs(sp,x)
    n <- w$dataObject.names[[x]][['fit']]
    dt <- w$generateParams(f.par[n],sp)[[1]]
    p <- w$funs$predict[[x]]
    p.par <- w$getPredictArgs(sp,x)
    IDs <- as.list(.r[,c(1,5)])
    id <- seq_along(IDs$modelID)
    names(id) <- IDs$modelID
    
    lapply(X=id,FUN=function(i,...) {
      .Fun(mID=IDs$modelID[i],sp = sp,method = x,fit = f,fit.par = f.par,pred = p,pred.par = p.par,rID = IDs$replicationID[i],n = n,dt=dt,w=w)
    },f=f,f.par=f.par,p.par=p.par,n=n,IDs=IDs,dt=dt,w=w,...)
  },w=w)
}
#---------

###############
.parLapply <- function(X, FUN, cl, ...) {
  parLapply(cl=cl,X=X,fun=FUN,...)
}
#--------
.feLapply <- function(X, FUN,cl=NULL,...) {
  # wrapper to foreach (cl is not used, just to make it consistent with the other functions!)
  xx <- eval(expression({foreach(i=X,.errorhandling='pass') %dopar% FUN(i,...)}))
  if (!is.null(names(X))) names(xx) <- names(X)
  xx
}
#---------
.fuLapply <- function(X, FUN,cl=NULL,...) {
  # wrapper to future_lapply (cl is not used, just to make it consistent with the other functions!)
  xx <- eval(expression({future_lapply(X=X,FUN=FUN,future.seed = TRUE,...)}))
  if (!is.null(names(X))) names(xx) <- names(X)
  xx
}
#-----
.Lapply <- function(X, FUN, ...) {
  # wrapper to lapply (cl is not used, just to make it consistent with the other functions!)
  lapply(X, FUN, ...)
}


.memory_used <- function () {
  # report the memory used by R objects in the session (based on the mem_used function in pryr)
  sum(gc()[, 1] * c(8L * .Machine$sizeof.pointer, 8))  / (1024*1024)
}
#----------


.fit=function(woL,species,models,runs,hasTest,.parMethod='parallel',.hostnames=NULL,.fork=TRUE,filename=NULL,compress=TRUE) {
  if (missing(species)) {
    if (length(woL$train) > 0) species <- names(woL$train)
    else stop('no species is available!')
  }
  
  if (missing(models)) models <- names(woL$funs$fit)
  names(models) <- models
  if (missing(runs)) {
    if (length(woL$replicates) == 0) {
      runs <- NULL
      n.total <- length(models)*length(species)
    } else {
      runs <- 1:length(woL$replicates[[1]])
      n.total <- length(models)*length(species)*length(runs)
    }
  } else {
    if (is.null(runs)) n.total <- length(models)*length(species)
    else {
      runs <- runs[runs %in% c(1:length(woL$replicates[[1]]))]
      if (length(runs) == 0) stop('all the specified runs index do not exist in the replications')
      n.total <- length(models)*length(species)*length(runs)
    }
  }
  
  if (missing(hasTest)) hasTest <- !is.null(woL$test)
  
  if (is.null(woL$ncore)) nc <- 1L
  else {
    #.require('parallel')
    #nc <- detectCores()
    #if (woL$ncore < nc) nc <- woL$ncore
    nc <- woL$ncore
    if (.is.windows() || 'maxent' %in% models) .fork <- FALSE
  }
  #----####################
  .run.info <- data.frame(matrix(NA,ncol=9,nrow=n.total))
  colnames(.run.info) <- c('modelID','species','method','replication','replicationID','success','training','test.dep','test.indep')
  .run.info[,1] <- 1:n.total
  
  
  if (!is.null(runs)) {
    .run.info[,c(5,3,2)] <- expand.grid(runs,models,species)
    .run.info[,4] <- rep(unlist(lapply(woL$replicates[[1]],function(x) x$method)),length.out=n.total)
  } else {
    .run.info[,3:2] <- expand.grid(models,species)
  }
  sm <- new('sdmModels',replicates=woL$replicates,data=woL$data,setting=woL$setting,recordIDs=woL$recordIDs,run.info=.run.info)
  #--------- RUN....:
  
  
  .memo <- NULL
  # if (.require("mraster")) {
  #   .memo <- eval(expression({memory(session=TRUE,echo=FALSE)}))
  #   if (is.numeric(.memo)) {
  #     .wf <- eval(parse(text = "mraster:::.change_unit(utils::object.size(woL$data@features)*150,'B','M')[[1]]")) # guessing the size of modelObj
  #     .mf <- floor(0.75 * ( .memo[2] /  .wf) )[[1]] # how many modelObj fits into 75% of the available memory!
  #     if (nrow(.run.info) > .mf) .ch <- ceiling(nrow(.run.info) / .mf)
  #     else .ch <- NULL
  #   } else .ch <- NULL 
  #   
  # } else {
  #   if (!is.null(filename)) warning('the package mraster is not installed (use installAll() to get it installed); this package helps when the output object size is big!')
  #   .ch <- NULL
  # }
  .ch <- NULL # temporary:=> until the mraster package is fixed and the above lines are modified based on!
  
  if (is.null(filename) && !is.null(.ch)) {
    #if (nrow(.run.info) > floor(.memo[1] / .wf)[[1]])
    if (.ch > 20) warning('It seems that the output of sdm would be a big object, and therefore, it may cause the error that it does not fit in the memory!\nYou can specify a filename to avoid the error!')
    .ch <- NULL
  }
  
  .tasks <- .getRunTaskTable(.run.info,species,models,runs,nc=nc)
  .tasks.type <- .tasks[[1]]
  if (.tasks.type %in% c('simple','simple.noRun')) {
    .tasks <- .run.info
    .tasks$species <- as.character(.tasks$species)
    .tasks$method <- as.character(.tasks$method)
  } else .tasks <- .tasks[[2]]
  
  if (!is.null(filename) && !is.null(.ch)) {
    if (grepl('.noRun',.tasks.type)) .tasks.type <- 'simple.noRun'
    else .tasks.type <- 'simple'
    .tasks <- .run.info
  }
  
  if (.tasks.type %in% c('species','runs','methods','simple')) {
    if (hasTest) .fun <- .fit_and_3evals
    else .fun <- .fit_and_2evals_dep
  } else if (.tasks.type %in% c('species.noRun','methods.noRun','simple.noRun')) {
    if (hasTest) .fun <- .fit_and_2evals_indep
    else .fun <- .fit_and_1eval
  }
  #----------
  cl <- NULL
  woL$tasks <- .tasks
  
  if (nc > 1) {
    #####################
    .require('parallel')
    
    if (.parMethod == 'future') {
      if (!.require('future.apply')) {
        warning('To use "future" as the parallel method, packages of "future" and "future.apply" are needed but they are not available, so the method is changed to "parallel"')
        .parMethod <- "parallel"
      }
    } else if (.parMethod == 'foreach') {
      if (!.require('foreach')) {
        warning('To use "foreach" as the parallel method, the "foreach" package is needed but it is not available, so the method is changed to "parallel"')
        .parMethod <- "parallel"
      }
    }
    #----
    if (.parMethod == 'future') {
      if (.fork) {
        .eval('plan(multicore,workers=nc,gc=TRUE)',environment())
      } else {
        .eval('plan(multisession,workers=nc,gc=TRUE)',environment())
      }
      .lapply <- .fuLapply
    } else {
      if (.fork) {
        cl <- makeForkCluster(nc)
        .lapply <- .parLapply
      } else {
        if (!is.null(.hostnames)) {
          cl <- try(makePSOCKcluster(.hostnames),silent = TRUE)
          if (inherits(cl,'try-error')) {
            cat('\n Error in connecting to remote servers:' ,print(cl))
            cat('\n cores on the local machine is considered!')
            if (nc == 1) {
              nc <- length(which(.hostnames == 'localhost'))
              nc <- max(nc,1)
              #if (nc > detectCores()) nc <- detectCores()
              cl <- makePSOCKcluster(nc)
            }
          }
        } else cl <- makePSOCKcluster(nc) # I should work more on providing the hostnames (on windows the connection needs to be through plink or PUTTY...)
        #--------
        clusterExport(cl,c('woL','.Lapply','.fun','.tasks','.fitLapply','.fitSimpleLapply','.fitlapply'),envir = environment())
        .w <- try(clusterEvalQ(cl,{
          library(sdm)
          sdm:::.addMethods()
          sdm:::.pkgLoad(woL$setting@methods)
        }),silent = TRUE)
        
        if (inherits(.w,'try-error')) {
          nc <- 1
          stopCluster(cl)
          warning('for some reasons, the clusters (for parallel processing) are not working!')
          cat('\n ncore is changed to 1!')
          .lapply <- .Lapply
        } else .lapply <- .parLapply
      }
      #---
      if (.parMethod == 'foreach' && .require('foreach') && (.require('doParallel') | .eval('getDoParRegistered()',environment()))) {
        .lapply <- .feLapply
        if (!.eval('getDoParRegistered()',environment())) {
          eval(expression({registerDoParallel(cl,nc)}))
        }
        
        if (eval(expression({getDoParWorkers() != nc}))) {
          eval(expression({warning(paste0('the number of workers registered for the foreach backend parallelisation is :',getDoParWorkers(),', which is different than the n.cores specified in the function! (n.cores is changed to this number)'))}))
          eval(expression({nc <- getDoParWorkers()}))
        }
      }
    }
  } else .lapply <- .Lapply
  #-------
  names(species) <- species
  
  if (is.null(.ch)) {
    if (.tasks.type %in% c('species','species.noRun')) {
      .mo <- .lapply(species, function(x,...) .fitlapply(x,...),.Fun=.fun,w=woL,cl=cl)
    } else if (.tasks.type == 'runs') {
      .mo <- lapply(species, function(x,...) .fitLapply(x,...),.Fun=.fun,.l1=lapply,.l2=.lapply,w=woL,cl=cl)
    } else if (.tasks.type %in% c('methods','methods.noRun')) {
      .mo <- lapply(species, function(x,...) .fitLapply(x,...),.Fun=.fun,.l1=.lapply,.l2=lapply,w=woL,cl=cl)
    } else if (.tasks.type %in% c('simple','simple.noRun')) {
      .w <- .lapply(.tasks$modelID, function(id, ...) .fitSimpleLapply(id,...), .Fun=.fun,w=woL, cl=cl)
      .mo <- vector('list',length=length(species))
      names(.mo) <- species
      .ww <- vector('list',length=length(models))
      names(.ww) <- models
      for (sp in species) {
        .mo[[sp]] <- .ww
        for (.wm in models) {
          .mo[[sp]][[.wm]] <- list()
          .www <- which(.tasks$species == sp & .tasks$method == .wm)
          for (i in .www) {
            .mo[[sp]][[.wm]][[as.character(.tasks$modelID[i])]] <- .w[[i]]
          }
        }
      }
      rm (.w, .ww, .www, .wm, sp, i); gc()
    }
    #------
    if (!is.null(filename)) {
      for (sp in species) {
        for (.wm in models) {
          .w <- names(.mo[[sp]][[.wm]])
          for (n in .w) {
            .n <- paste0('m_',n,'_',sp,'_',.wm,'.sdm')
            saveRDS(.mo[[sp]][[.wm]][[n]],paste0(filename,'/',.n),compress=compress)
            .mo[[sp]][[.wm]][[n]] <- .n
          }
        }
      }
    }
  } else {
    .mo <- vector('list',length=length(species))
    names(.mo) <- species
    .ww <- vector('list',length=length(models))
    names(.ww) <- models
    for (sp in species) {
      .mo[[sp]] <- .ww
      for (.wm in models) {
        .mo[[sp]][[.wm]] <- list()
      }
    }
    
    .w <- ceiling(nrow(.tasks) / .ch)
    .idList <- c()
    for (i in 1:(.ch-1)) {
      .idList[[i]] <- .tasks$modelID[(i-1)*.w + 1:.w]
    }
    i <- i+1
    .idList[[i]] <- .tasks$modelID[((i-1)*.w + 1):nrow(.run.info)]
    
    for (i in 1:length(.idList)) {
      .w <- .lapply(.idList[[i]], function(id, ...) .fitSimpleLapply(id,...), .Fun=.fun,w=woL, cl=cl)
      for (j in seq_along(.idList[[i]])) {
        .wid <-  which(.tasks$modelID == .idList[[i]][j])
        sp <- .tasks$species[.wid]
        .wm <- .tasks$method[.wid]
        n <- as.character(.tasks$modelID[.wid])
        .n <- paste0('m_',n,'_',sp,'_',.wm,'.sdm')
        saveRDS(.w[[j]],paste0(filename,'/',.n),compress=compress)
        .mo[[sp]][[.wm]][[n]] <- .n
      }
      rm (.w); gc()
    }
    rm (.wid, .wm, sp, i, j, .idList); gc()
  }
  #----------
  
  try(stopCluster(cl),silent = TRUE)
  sm@models <- .mo
  
  #---- update the run.info with whether the fitting and evaluations were successfully done!
  for (i in sm@run.info[,1]) {
    o <- sm@models[[sm@run.info[i,2]]][[sm@run.info[i,3]]][[as.character(i)]]
    sm@run.info[i,6:9] <- c(!is.null(o@object) && !inherits(o@object,"try-error"),inherits(o@evaluation[['training']],"sdmEvaluate"),inherits(o@evaluation[['test.dep']],"sdmEvaluate"),inherits(o@evaluation[['test.indep']],"sdmEvaluate"))
  }
  sm
}
#----------