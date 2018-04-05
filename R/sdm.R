# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  April 2018
# Version 3.4
# Licence GPL v3
#--------


.methodFix <- function(n) {
  for (i in seq_along(n)) {
    nx <- .sdmMethods$whichMethod(n[i])
    if (!is.null(nx)) n[i] <- nx
    else n[i] <- NA
  }
  n
}
#----------
.replicate.methodFix <- function(n) {
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

.getFormula.gammgcv.rhs <- function(n,nFact=NULL,env=parent.frame()) {
  if (!is.null(nFact)) as.formula(paste('~',paste(c(paste(paste('s(',n,sep=''),')',sep=''),nFact),collapse='+'),sep=''),env = env)
  else as.formula(paste('~',paste(paste(paste('s(',n,sep=''),')',sep=''),collapse='+'),sep=''),env = env)
}
.getFormula.gammgcv <- function(n,nFact=NULL,env=parent.frame()) {
  if (!is.null(nFact)) as.formula(paste(n[1],'~',paste(c(paste(paste('s(',n[-1],sep=''),')',sep=''),nFact),collapse='+'),sep=''),env = env)
  else as.formula(paste(n[1],'~',paste(paste(paste('s(',n[-1],sep=''),')',sep=''),collapse='+'),sep=''),env = env)
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

.factorFix <- function(data1,data2,nFact,nf) {
  # assign the problematic factors to the more similar factors according to continuous variables
  # if no continus variable does exist, then it is assigned to a dominant class.
  if (missing(nFact) || is.null(nFact)) nFact <- colnames(data2)[.where(is.factor,data2)]
  if (missing(nf) || is.null(nf)) {
    nf <- colnames(data2)[!colnames(data2) %in% nFact]
    if (length(nf) == 0) nf <- NULL
  }
  
  dd <- data2
  
  for (i in seq_along(nFact)) {
    data1[,nFact[i]] <- factor(data1[,nFact[i]])
    data2[,nFact[i]] <- factor(data2[,nFact[i]])
    dd[,nFact[i]] <- as.character(dd[,nFact[i]])
    
    fc <- .checkFactor(data1[,nFact[i]],data2[,nFact[i]])
    if (!is.null(fc)) {
      p <- names(fc$p)
      np <- names(fc$np)
      for (j in seq_along(p)) {
        w <- fc[['p']][[p[[j]]]]
        if (!is.null(nf)) {
          m <- rep(NA,length(np))
          d2 <- data2[w,which(colnames(data2) %in% nf)]
          if (length(nf) > 1) {
            options(warn=-1)
            for (k in seq_along(np)) {
              ww <- fc[['np']][[np[[k]]]]
              d1 <- data2[ww,which(colnames(data2) %in% nf)]
              m[k] <- mean(try(.mahal(d1,d2),silent=TRUE),na.rm=TRUE)
            }
            options(warn=0)
          } else {
            for (k in seq_along(np)) {
              ww <- fc[['np']][[np[[k]]]]
              d1 <- data2[ww,which(colnames(data2) %in% nf)]
              m[k] <- abs(mean(d2,na.rm=TRUE) - mean(d1,na.rm=TRUE))
            }
          }
          
          ww <- which.min(m)
          if (length(ww) > 0) dd[w,nFact[i]] <- np[ww]
          else {
            dom.class <- summary(data1[,nFact[i]])
            dom.class <- attr(sort(dom.class, decreasing = TRUE)[1], "names")
            dd[w,nFact[i]] <- dom.class
          }
        } else {
          dom.class <- summary(data1[,nFact[i]])
          dom.class <- attr(sort(dom.class, decreasing = TRUE)[1], "names")
          dd[w,nFact[i]] <- dom.class
        }
        
      }
      
    }
  }
  for (i in seq_along(nFact)) dd[,nFact[i]] <- factor(dd[,nFact[i]])
  .eqFactLevel(data1,dd)
}
#--------

.factorFixW <- function(data1,data2,nFact,nf) {
  # just generates a list to report which class should be assigned to which class!
  if (missing(nFact) || is.null(nFact)) nFact <- colnames(data2)[.where(is.factor,data2)]
  if (missing(nf) || is.null(nf)) {
    nf <- colnames(data2)[!colnames(data2) %in% nFact]
    if (length(nf) == 0) nf <- NULL
  }
  
  dd <- data2
  o <- list()
  for (i in seq_along(nFact)) {
    data1[,nFact[i]] <- factor(data1[,nFact[i]])
    data2[,nFact[i]] <- factor(data2[,nFact[i]])
    dd[,nFact[i]] <- as.character(dd[,nFact[i]])
    
    fc <- .checkFactor(data1[,nFact[i]],data2[,nFact[i]])
    if (!is.null(fc)) {
      p <- names(fc$p)
      np <- names(fc$np)
      for (j in seq_along(p)) {
        w <- fc[['p']][[p[[j]]]]
        if (!is.null(nf)) {
          m <- rep(NA,length(np))
          d2 <- data2[w,which(colnames(data2) %in% nf)]
          if (length(nf) > 1) {
            options(warn=-1)
            for (k in seq_along(np)) {
              ww <- fc[['np']][[np[[k]]]]
              d1 <- data2[ww,which(colnames(data2) %in% nf)]
              m[k] <- mean(try(.mahal(d1,d2),silent=TRUE),na.rm=TRUE)
            }
            if (any(is.na(m))) {
              m <- matrix(NA,nrow=length(nf),ncol=length(np))
              for (k in seq_along(np)) {
                ww <- fc[['np']][[np[[k]]]]
                d1 <- data2[ww,which(colnames(data2) %in% nf)]
                for (nfi in seq_along(nf)) {
                  m[nfi,k] <- abs(mean(d2[[nf[nfi]]],na.rm=TRUE) - mean(d1[[nf[nfi]]],na.rm=TRUE))
                }
              }
              m <- abs(apply(t(apply(m,1,function(x) (x - mean(x)) / sd(x))),2,mean))
            }
            options(warn=0)
          } else {
            for (k in seq_along(np)) {
              ww <- fc[['np']][[np[[k]]]]
              d1 <- data2[ww,which(colnames(data2) %in% nf)]
              m[k] <- abs(mean(d2,na.rm=TRUE) - mean(d1,na.rm=TRUE))
            }
          }
          
          ww <- which.min(m)
          if (length(ww) > 0) {
            o <- c(o,list(c(field=nFact[i],old=unique(dd[w,nFact[i]]),new=np[ww])))
          } else {
            dom.class <- summary(data1[,nFact[i]])
            dom.class <- attr(sort(dom.class, decreasing = TRUE)[1], "names")
            o <- c(o,list(c(field=nFact[i],old=unique(dd[w,nFact[i]]),new=dom.class)))
          }
        } else {
          dom.class <- summary(data1[,nFact[i]])
          dom.class <- attr(sort(dom.class, decreasing = TRUE)[1], "names")
          o <- c(o,list(c(field=nFact[i],old=unique(dd[w,nFact[i]]),new=dom.class)))
        }
      }
    }
  }
  o
}

#--------

.factorFix.bm <- function(data1,data2,nFact,nf) {
  # problematic factors are moved from data2 to data1
  if (missing(nFact) || is.null(nFact)) nFact <- colnames(data2)[.where(is.factor,data2)]
  er <- FALSE
  
  for (i in seq_along(nFact)) {
    data1[,nFact[i]] <- factor(data1[,nFact[i]])
    data2[,nFact[i]] <- factor(data2[,nFact[i]])
    fc <- .checkFactor(data1[,nFact[i]],data2[,nFact[i]])
    if (!is.null(fc)) {
      er <- TRUE
      p <- names(fc$p)
      for (j in seq_along(p)) {
        ww <- fc[['p']][[p[[j]]]]
        ww <- sample(ww,1)
        data1 <- rbind(data1,data2[ww,])
        data2 <- data2[-ww,]
        #dd <- dd[-ww,]
      }
    }
  }
  if (er) {
    for (i in seq_along(nFact)) {
      data1[,nFact[i]] <- factor(data1[,nFact[i]])
      data2[,nFact[i]] <- factor(data2[,nFact[i]])
    }
    return(list(train=data1,test=data2,IDs=fc$p))
  }
}
#--------

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
.require <- function(x) {
  # based on simplifying the code of the reqiure function in the base package
  loaded <- paste("package", x, sep = ":") %in% search()
  if (!loaded) {
    value <- tryCatch(library(x,character.only = TRUE, logical.return = TRUE, warn.conflicts = FALSE, quietly = TRUE), error = function(e) e)
    if (inherits(value, "error")) {
      return(FALSE)
    }
    if (!value) return(FALSE)
  } else value <- TRUE
  value
}
#----------
.loadLib <- function(pkgs) {
  options(warn=-1)
  return(unlist(lapply(pkgs,function(x) {
    all(unlist(lapply(x,function(p) {.require(p)})))
  })))
  options(warn=0)
}
#---------
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
  # if (!all(ww)) {
  #   if (!any(ww)) {
  #     cat('some methods are removed because they depend on some packages that are not installed on this machine!\n')
  #     cat('you can use installAll() function to simply install all the packages that may be required by some functions in the sdm package!\n')
  #     stop(paste('There is no installed packages rquired by the selected methods. Package names:',paste(unlist(pkgs),collapse=', ')))
  #   } else {
  #     cat('some methods are removed because they depend on some packages that are not installed on this machine!\n')
  #     cat('you can use installAll() function to simply install all the packages that may be required by some functions in the sdm package!\n')
  #     warning(paste('There is no installed packages rquired by the methods:',paste(s@methods[!ww],collapse=', '),'; These methods are excluded! The packages need to be installed for these methods:',paste(unlist(pkgs[!ww]),collapse=', ')))
  #     s@methods <- s@methods[ww]
  #   }
  # }
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
      cat('some methods are removed because they depend on some packages that are not installed on this machine!\n')
      cat('you can use installAll() function to simply install all the packages that may be required by some functions in the sdm package!\n')
      stop(paste('There is no installed packages rquired by the selected methods. Package names:',paste(unlist(pkgs),collapse=', ')))
    } else {
      cat('some methods are removed because they depend on some packages that are not installed on this machine!\n')
      cat('you can use installAll() function to simply install all the packages that may be required by some functions in the sdm package!\n')
      warning(paste('There is no installed packages rquired by the methods:',paste(s@methods[!ww],collapse=', '),'; These methods are excluded! The packages need to be installed for these methods:',paste(unlist(pkgs[!ww]),collapse=', ')))
      s@methods <- s@methods[ww]
    }
  }
  #------------
  if (!is.null(s@seed)) set.seed(s@seed)
  #-----------
  #fr <- .getFeaturetype(d,s@sdmFormula)
  w <- new('.workload',ncore=s@parallelSettings@ncore,data=d,setting=s,frame=s@featuresFrame,filename=filename)
  hasTest <- 'test' %in% d@groups$training@values[,2]
  nFact <- NULL
  if (!is.null(d@factors) > 0 && any(d@factors %in% s@sdmFormula@vars)) nFact <- d@factors[d@factors %in% s@sdmFormula@vars]
  nf <- .excludeVector(d@features.name,nFact)
  nf <- nf[nf %in% s@sdmFormula@vars]
  nFact <- nFact[nFact %in% s@sdmFormula@vars]
  
  for (sp in s@sdmFormula@species) {
    dt <- as.data.frame(d,sp=sp,grp='train')
    w$recordIDs[[sp]]$train <- data.frame(rID=dt[,1],rowID=1:nrow(dt))
    f <- .getModelFrame(w$frame,dt,response=sp)
    if (!is.null(f$specis_specific)) {
      w$train[[sp]]$sdmDataFrame <- cbind(dt[,sp],f$features,f$specis_specific)
    } else w$train[[sp]]$sdmDataFrame  <- cbind(dt[,sp],f$features)
    colnames(w$train[[sp]]$sdmDataFrame)[1] <- sp
    
    if (hasTest) {
      dt <- as.data.frame(d,sp=sp,grp='test')
      w$recordIDs[[sp]]$test <- data.frame(rID=dt[,1],rowID=1:nrow(dt))
      f <- .getModelFrame(w$frame,dt,response=sp)
      if (!is.null(f$specis_specific)) {
        w$test[[sp]]$sdmDataFrame <- cbind(dt[,sp],f$features,f$specis_specific)
      } else w$test[[sp]]$sdmDataFrame <- cbind(dt[,sp],f$features)
      colnames(w$test[[sp]]$sdmDataFrame)[1] <- sp
      
      if (!is.null(nFact)) {
        
        for (nF in nFact) {
          fc <- .checkFactor(w$train[[sp]]$sdmDataFrame[,nF],w$test[[sp]]$sdmDataFrame[,nF])
          if (!is.null(fc)) {
            p <- names(fc$p)
            for (j in seq_along(p)) {
              ww <- fc[['p']][[p[[j]]]]
              if (length(ww) > 1) ww <- sample(ww,1)
              w$train[[sp]]$sdmDataFrame <- rbind(w$train[[sp]]$sdmDataFrame,w$test[[sp]]$sdmDataFrame[ww,])
              w$test[[sp]]$sdmDataFrame[ww,] <- w$test[[sp]]$sdmDataFrame[-ww,]
              w$data <- .updateGroup(w$data,c(.getRecordID(w$recordIDs,sp=sp,id = ww,train=FALSE),'test','train'))
            }
          }
        }
        
        #f <- .factorFix.bm(w$train[[sp]]$sdmDataFrame,w$test[[sp]]$sdmDataFrame,nFact)
        #if (!is.null(f)) {
        #  w$train[[sp]]$sdmDataFrame <- f$train
        #  w$test[[sp]]$sdmDataFrame <- f$test
        #}
      }
    }
  }
  #------------
  w$funs[['fit']] <- .sdmMethods$getFitFunctions(s@methods)
  w$arguments[['fit']] <- .sdmMethods$getFitArguments(s@methods)
  w$funs[['predict']] <- .sdmMethods$getPredictFunctions(s@methods)
  w$arguments[['predict']] <- .sdmMethods$getPredictArguments(s@methods)
  #---
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
  
  #w$dataObject.names <- unique(unlist(lapply(s@methods, .sdmMethods$getDataArgumentNames)))
  mo <- s@methods
  names(mo) <- mo
  w$dataObject.names <- lapply(mo, .sdmMethods$getDataArgumentNames)
  w$settingRules <- lapply(mo,function(x) .sdmMethods$Methods[[x]]@settingRules)
  #-----------
  
  #reserved.names <- w$getReseved.names()
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
      else family <- 'xxx'
      
      # sdmDataFrame! Leter should be checked for other types of data!
      for (ff in f) {
        w$replicates[[sp]] <- c(w$replicates[[sp]],ff(x=w$train[[sp]][['sdmDataFrame']][,sp],family=family,stratify=TRUE,test.percent=s@test.percentage,nfolds=s@cv.folds,replicates=s@n.replicates))
      }
    }
  }
  
  # for each replicatios, checks whether the factor level is going to be problematic
  # if so, move the record from test to train
  if (!is.null(nFact)) {
    for (nF in nFact) {
      for (sp in names(w$replicates)) {
        for (i in 1:length(w$replicates[[sp]])) {
          #fc <- .checkFactor(w$train[[sp]]$sdmDataFrame[w$runtasks$runIndex[[sp]][[i]]$train,nF],w$train[[sp]]$sdmDataFrame[w$runtasks$runIndex[[sp]][[i]]$test,nF])
          fc <- .checkFactor(w$train[[sp]]$sdmDataFrame[w$replicates[[sp]][[i]]$train,nF],w$train[[sp]]$sdmDataFrame[w$replicates[[sp]][[i]]$test,nF])
          if (!is.null(fc)) {
            p <- names(fc$p)
            for (j in seq_along(p)) {
              ww <- fc[['p']][[p[[j]]]]
              if (length(ww) > 1) ww <- sample(ww,1)
              w$replicates[[sp]][[i]]$train <- c(w$replicates[[sp]][[i]]$train,w$replicates[[sp]][[i]]$test[ww])
              w$replicates[[sp]][[i]]$test <- w$replicates[[sp]][[i]]$test[-ww]
            }
          }
        }
      }
    }
  }
  #-----
  w
}

#----------------------------------------
# if (!isGeneric("sdmSetting")) {
#   setGeneric("sdmSetting", function(formula,data,methods,interaction.depth=1,n=1,replication=NULL,
#                                     cv.folds=NULL,test.percent=NULL,bg=NULL,bg.n=NULL,var.importance=NULL,response.curve=TRUE,
#                                     var.selection=FALSE,ncore=1L,...)
#     standardGeneric("sdmSetting"))
# }
# 
# setMethod('sdmSetting', signature(formula='ANY','sdmdata','character'), 
#           function(formula,data,methods,interaction.depth=1,n=1,replication=NULL,
#                    cv.folds=NULL,test.percent=NULL,bg=NULL,bg.n=NULL,var.importance=NULL,response.curve=TRUE,
#                    var.selection=FALSE,ncore=1L,...) {
#             
#             if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
#             
#             dot <- list(...)
#             sobj <- NULL
#             if (length(dot) > 0) {
#               ndot <- names(dot)
#               if ('' %in% ndot) {
#                 for (i in seq_along(which(ndot == ''))) {
#                   if (inherits(dot[[i]],'.sdmCorSetting')) {
#                     sobj <- dot[[i]]
#                     break
#                   }
#                 }
#                 dot <- dot[-which(ndot == '')]
#                 ndot <- names(dot)
#               }
#               
#               a <- c('interaction.depth','replication','cv.folds','test.percent','bg','bg.n','var.importance','response.curve','var.selection','ncore')
#               ndot <- .pmatch(ndot,a)
#               w <- !is.na(ndot)
#               if (length(w) > 0) {
#                 dot <- dot[w]
#                 ndot <- ndot[w]
#                 names(dot) <- ndot
#               }
#               
#               
#               if ('setting' %in% names(dot) && inherits(dot[['setting']],'.sdmCorSetting')) {
#                 sobj <- dot[['setting']]
#                 dot <- dot[-which(ndot == 'setting')]
#                 ndot <- names(dot)
#               }
#               
#               if (length(dot) > 0) {
#                 if (length(ndot) > 0) {
#                   for (nd in ndot) {
#                     if (nd == 'interaction.depth' && interaction.depth == 1) interaction.depth <- dot[[nd]]
#                     else if (nd == 'ncore' && ncore == 1L) ncore <- dot[[nd]]
#                     else if (nd == 'replication' && is.null(replication)) replication <- dot[[nd]]
#                     else if (nd == 'cv.folds' && is.null(cv.folds)) cv.folds <- dot[[nd]]
#                     else if (nd == 'test.percent' && is.null(test.percent)) test.percent <- dot[[nd]]
#                     else if (nd == 'bg' && is.null(bg)) bg <- dot[[nd]]
#                     else if (nd == 'bg.n' && is.null(bg.n)) bg.n <- dot[[nd]]
#                     else if (nd == 'var.importance' && is.null(var.importance)) var.importance <- dot[[nd]]
#                     else if (nd == 'response.curve' && response.curve && is.logical(dot[[nd]])) response.curve <- dot[[nd]]
#                     else if (nd == 'var.selection' && !var.selection && is.logical(dot[[nd]])) var.selection <- dot[[nd]]
#                   }
#                 }
#               }
#             }
#             #--------
#             
#             m <- .methodFix(methods)
#             if (any(is.na(m))) stop(paste('methods',paste(methods[is.na(m)],collapse=', '),'do not exist!'))
#             m <- unique(m)
#             
#             s <- new('.sdmCorSetting',methods=m)
#             s@distribution <- .getSpeciesDistribution(data)
#             
#             if (missing(formula)) {
#               if (!is.null(sobj)) {
#                 if (all(sobj@sdmFormula@vars %in% data@features.name)) s@sdmFormula <- sobj@sdmFormula
#                 else s@sdmFormula <- data@sdmFormula
#               } else s@sdmFormula <- data@sdmFormula
#               
#             } else if (inherits(formula,'sdmFormula')) s@sdmFormula <- formula
#             else if (inherits(formula,'formula')) {  
#               s@sdmFormula <- .exFormula(formula,as.data.frame(data)[,-1])
#             } else if (inherits(formula,'.sdmCorSetting')) {
#               sobj <- formula
#               if (all(sobj@sdmFormula@vars %in% data@features.name)) s@sdmFormula <- sobj@sdmFormula
#               else s@sdmFormula <- data@sdmFormula
#             } else {
#               if (!is.null(sobj)) {
#                 if (all(sobj@sdmFormula@vars %in% data@features.name)) s@sdmFormula <- sobj@sdmFormula
#                 else s@sdmFormula <- data@sdmFormula
#               } else s@sdmFormula <- data@sdmFormula
#             }
#             
#             s@featuresFrame <- .getFeaturetype(data,s@sdmFormula)  
#               
#             if (!is.null(test.percent)) s@test.percentage <- test.percent
#             else {
#               if (!is.null(sobj)) {
#                 if (!is.null(sobj@test.percent)) s@test.percentage <- sobj@test.percent
#               }
#             }
#             
#             s@interaction.depth <- interaction.depth
#             if (interaction.depth ==1 && !is.null(sobj) && !is.null(sobj@interaction.depth)) s@interaction.depth <- sobj@interaction.depth
#             
#             s@ncore <- ncore
#             if (ncore == 1L && !is.null(sobj) && length(sobj@ncore) == 1) s@ncore <- sobj@ncore
#             
#             if (!is.null(replication)) {
#               nx <- .replicate.methodFix(replication)
#               if (any(is.na(nx))) warning(paste(paste(replication[is.na(nx)],collapse=', '),'methods in replication are not found [They are ignored!]'))
#               replication <- nx[!is.na(nx)]
#               s@replicate <- replication
#             } else {
#               if (!is.null(sobj)) {
#                 if (!is.null(sobj@replicate)) s@replicate <- sobj@replicate
#               }
#               if (is.null(s@replicate) && !is.null(s@test.percentage)) {
#                 s@replicate <- "subsampling"
#               }
#             }
#             
#             s@n.replicates <- n
#             if (!is.null(sobj) && !is.null(sobj@n.replicates)) s@n.replicates <- sobj@n.replicates
#             
#             if ("subsampling" %in% s@replicate) {
#               if (is.null(s@test.percentage)) s@test.percentage <- 30
#             }
#             
#             if (!is.null(cv.folds)) s@cv.folds <- cv.folds
#             else {
#               if (!is.null(sobj) && !is.null(sobj@cv.folds)) s@cv.folds <- sobj@cv.folds
#               if (is.null(s@cv.folds) && "cross_validation" %in% s@replicate) s@cv.folds <- 5
#             }
#             
#             if (!is.null(s@cv.folds) && !"cross_validation" %in% s@replicate) {
#               s@replicate <- c("cross_validation",s@replicate)
#             }
#             
#             if (!is.null(bg)) s@pseudo.absence.methods <- bg
#             else {
#               if (!is.null(sobj)) {
#                 if (!is.null(sobj@pseudo.absence.methods)) s@pseudo.absence.methods <- sobj@pseudo.absence.methods
#               }
#             }
#             if (!is.null(bg.n)) s@n.pseudo.absence <- bg.n
#             else {
#               if (!is.null(sobj)) {
#                 if (!is.null(sobj@n.pseudo.absence)) s@n.pseudo.absence <- sobj@n.pseudo.absence
#               }
#               if (is.null(s@n.pseudo.absence) && !is.null(s@pseudo.absence.methods)) {
#                 s@n.pseudo.absence <- 1000
#               }
#             }
#             if (!is.null(var.importance)) s@varImportance.methods <- var.importance
#             else {
#               if (!is.null(sobj)) {
#                 if (!is.null(sobj@varImportance.methods)) s@varImportance.methods <- sobj@varImportance.methods
#               }
#             }
#             if (response.curve) s@response.curve <- TRUE
#             else {
#               if (!is.null(sobj)) {
#                 if (!is.null(sobj@response.curve) && sobj@response.curve) s@response.curve <- sobj@response.curve
#               } else s@response.curve <- FALSE
#             }
#             
#             if (var.selection) s@var.selection <- TRUE
#             else {
#               if (!is.null(sobj)) {
#                 if (!is.null(sobj@var.selection) && sobj@var.selection) s@var.selection <- sobj@var.selection
#               } else s@var.selection <- FALSE
#             }
#             
#             if (!is.null(interaction.depth)) s@interaction.depth <- interaction.depth
#             else {
#               if (!is.null(sobj)) {
#                 if (!is.null(sobj@interaction.depth)) s@interaction.depth <- sobj@interaction.depth
#               }
#             }
#             s
#           }
# )
# #----------------
if (!isGeneric("sdmSetting")) {
  setGeneric("sdmSetting", function(formula,data,methods,interaction.depth=1,n=1,replication=NULL,
                                    cv.folds=NULL,test.percent=NULL,bg=NULL,bg.n=NULL,var.importance=NULL,response.curve=TRUE,
                                    var.selection=FALSE,ncore=1L,modelSettings=NULL,seed=NULL,parallelSettings=NULL,...)
    standardGeneric("sdmSetting"))
}

setMethod('sdmSetting', signature(formula='ANY','sdmdata','character'), 
          function(formula,data,methods,interaction.depth=1,n=1,replication=NULL,
                   cv.folds=NULL,test.percent=NULL,bg=NULL,bg.n=NULL,var.importance=NULL,response.curve=TRUE,
                   var.selection=FALSE,ncore=1L,modelSettings=NULL,seed=NULL,parallelSettings=NULL,...) {
            
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            
            dot <- list(...)
            sobj <- NULL
            if (length(dot) > 0) {
              ndot <- names(dot)
              if ('' %in% ndot) {
                for (i in seq_along(which(ndot == ''))) {
                  if (inherits(dot[[i]],'.sdmCorSetting')) {
                    sobj <- dot[[i]]
                    break
                  }
                }
                dot <- dot[-which(ndot == '')]
                ndot <- names(dot)
              }
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
              a <- c('interaction.depth','replication','cv.folds','test.percent','bg','bg.n','var.importance','response.curve','var.selection','ncore','modelSettings','seed','setting','parallelSettings')
              ndot <- .pmatch(ndot,a)
              w <- !is.na(ndot)
              if (length(w) > 0) {
                dot <- dot[w]
                ndot <- ndot[w]
                names(dot) <- ndot
              }
              
              if ('setting' %in% names(dot) && inherits(dot[['setting']],'.sdmCorSetting')) {
                sobj <- dot[['setting']]
                dot <- dot[-which(ndot == 'setting')]
                ndot <- names(dot)
              }
              
              if (length(dot) > 0) {
                if (length(ndot) > 0) {
                  for (nd in ndot) {
                    if (nd == 'interaction.depth' && interaction.depth == 1) interaction.depth <- dot[[nd]]
                    else if (nd == 'ncore' && ncore == 1L) ncore <- dot[[nd]]
                    else if (nd == 'replication' && is.null(replication)) replication <- dot[[nd]]
                    else if (nd == 'cv.folds' && is.null(cv.folds)) cv.folds <- dot[[nd]]
                    else if (nd == 'test.percent' && is.null(test.percent)) test.percent <- dot[[nd]]
                    else if (nd == 'bg' && is.null(bg)) bg <- dot[[nd]]
                    else if (nd == 'bg.n' && is.null(bg.n)) bg.n <- dot[[nd]]
                    else if (nd == 'var.importance' && is.null(var.importance)) var.importance <- dot[[nd]]
                    else if (nd == 'response.curve' && response.curve && is.logical(dot[[nd]])) response.curve <- dot[[nd]]
                    else if (nd == 'var.selection' && !var.selection && is.logical(dot[[nd]])) var.selection <- dot[[nd]]
                    else if (nd == 'modelSettings' && is.null(modelSettings)) modelSettings <- dot[[nd]]
                    else if (nd == 'seed' && is.null(seed)) seed <- dot[[nd]]
                    else if (nd == 'parallelSettings' && is.null(parallelSettings) && is.list(dot[[nd]])) parallelSettings <- dot[[nd]]
                  }
                }
              }
            }
            #--------
            
            m <- .methodFix(methods)
            if (any(is.na(m))) stop(paste('methods',paste(methods[is.na(m)],collapse=', '),'do not exist!'))
            m <- unique(m)
            #---------
            s <- new('.sdmCorSetting',methods=m)
            #---------
            if (missing(formula)) {
              if (!is.null(sobj)) {
                if (all(sobj@sdmFormula@vars %in% data@features.name)) s@sdmFormula <- sobj@sdmFormula
                else s@sdmFormula <- data@sdmFormula
              } else s@sdmFormula <- data@sdmFormula
              
            } else if (inherits(formula,'sdmFormula')) s@sdmFormula <- formula
            else if (inherits(formula,'formula')) {  
              s@sdmFormula <- .exFormula(formula,as.data.frame(data)[,-1])
            } else if (inherits(formula,'.sdmCorSetting')) {
              sobj <- formula
              if (all(sobj@sdmFormula@vars %in% data@features.name)) s@sdmFormula <- sobj@sdmFormula
              else s@sdmFormula <- data@sdmFormula
            } else {
              if (!is.null(sobj)) {
                if (all(sobj@sdmFormula@vars %in% data@features.name)) s@sdmFormula <- sobj@sdmFormula
                else s@sdmFormula <- data@sdmFormula
              } else s@sdmFormula <- data@sdmFormula
            }
            
            s@featuresFrame <- .getFeaturetype(data,s@sdmFormula)  
            #---------
            s@distribution <- .getSpeciesDistribution(data,sp=s@sdmFormula@species)
            #---------
            if (!is.null(test.percent)) s@test.percentage <- test.percent
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@test.percentage)) s@test.percentage <- sobj@test.percentage
              }
            }
            #---------
            if (!is.null(parallelSettings) && is.list(parallelSettings)) {
              nparallel <- names(parallelSettings)
              a <- c('ncore','doParallel','method','cluster','hosts','fork','type')
              nparallel <- .pmatch(nparallel,a)
              w <- !is.na(nparallel)
              if (length(w) > 0) {
                parallelSettings <- parallelSettings[w]
                nparallel <- nparallel[w]
                names(parallelSettings) <- nparallel
              }
              #--
              if ('ncore' %in% nparallel) s@parallelSettings@ncore <- parallelSettings$ncore
              else {
                s@parallelSettings@ncore <- ncore
                if (ncore == 1L && !is.null(sobj) && length(sobj@parallelSettings@ncore) == 1) s@parallelSettings@ncore <- sobj@parallelSettings@ncore
              }
              #--
              if ('method' %in% nparallel) {
                if (parallelSettings$method %in% c('parallel','foreach')) s@parallelSettings@method <- parallelSettings$method
                else {
                  warning('parallelisation method is not recognised; the default value ("parallel") is used!')
                  s@parallelSettings@method <- 'parallel'
                }
              } else s@parallelSettings@method <- 'parallel'
              #--
              if ('fork' %in% nparallel) {
                if (is.logical(parallelSettings$fork)) {
                  if (parallelSettings$fork && .is.windows()) {
                    warning('"fork" in parallelisation setting cannot be TRUE on Windows Operating Systems; It is changed to FALSE!')
                    s@parallelSettings@fork <- FALSE
                  } else s@parallelSettings@fork <- parallelSettings$fork
                } else {
                  warning('"fork" in parallelisation setting should be logical; the default value is used!')
                  s@parallelSettings@fork <- !.is.windows()
                }
              } else s@parallelSettings@fork <- !.is.windows()
              #--
              if ('type' %in% nparallel) s@parallelSettings@type <- parallelSettings$type
              #--
              if ('doParallel' %in% nparallel && is.expression(parallelSettings$doParallel)) s@parallelSettings@doParallel <- parallelSettings$doParallel
              #--
              if ('cluster' %in% nparallel && inherits(parallelSettings$cluster,'cluster')) s@parallelSettings@cluster <- parallelSettings$cluster
              #--
              if ('hosts' %in% nparallel && is.character(parallelSettings$hosts)) s@parallelSettings@hosts <- parallelSettings$hosts
              
            } else {
              if (!is.null(sobj)) s@parallelSettings <- sobj@parallelSettings
              else {
                s@parallelSettings@ncore <- ncore
                s@parallelSettings@fork <- !.is.windows()
              }
            }
            
            #---------
            if (!is.null(replication)) {
              nx <- .replicate.methodFix(replication)
              if (any(is.na(nx))) warning(paste(paste(replication[is.na(nx)],collapse=', '),'methods in replication are not found [They are ignored!]'))
              replication <- nx[!is.na(nx)]
              s@replicate <- replication
            } else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@replicate)) s@replicate <- sobj@replicate
              }
              if (is.null(s@replicate) && !is.null(s@test.percentage)) {
                s@replicate <- "subsampling"
              }
            }
            
            s@n.replicates <- n
            if (!is.null(sobj) && !is.null(sobj@n.replicates)) s@n.replicates <- sobj@n.replicates
            
            if ("subsampling" %in% s@replicate) {
              if (is.null(s@test.percentage)) s@test.percentage <- 30
            }
            
            if (!is.null(cv.folds)) s@cv.folds <- cv.folds
            else {
              if (!is.null(sobj) && !is.null(sobj@cv.folds)) s@cv.folds <- sobj@cv.folds
              if (is.null(s@cv.folds) && "cross_validation" %in% s@replicate) s@cv.folds <- 5
            }
            
            if (!is.null(s@cv.folds) && !"cross_validation" %in% s@replicate) {
              s@replicate <- c("cross_validation",s@replicate)
            }
            #---------
            if (!is.null(bg)) s@pseudo.absence.methods <- bg
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@pseudo.absence.methods)) s@pseudo.absence.methods <- sobj@pseudo.absence.methods
              }
            }
            if (!is.null(bg.n)) s@n.pseudo.absence <- bg.n
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@n.pseudo.absence)) s@n.pseudo.absence <- sobj@n.pseudo.absence
              }
              if (is.null(s@n.pseudo.absence) && !is.null(s@pseudo.absence.methods)) {
                s@n.pseudo.absence <- 1000
              }
            }
            #---------
            if (!is.null(var.importance)) s@varImportance.methods <- var.importance
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@varImportance.methods)) s@varImportance.methods <- sobj@varImportance.methods
              }
            }
            #---------
            if (response.curve) s@response.curve <- TRUE
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@response.curve) && sobj@response.curve) s@response.curve <- sobj@response.curve
              } else s@response.curve <- FALSE
            }
            #---------
            if (var.selection) s@var.selection <- TRUE
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@var.selection) && sobj@var.selection) s@var.selection <- sobj@var.selection
              } else s@var.selection <- FALSE
            }
            #---------
            #s@interaction.depth <- interaction.depth
            #if (interaction.depth ==1 && !is.null(sobj) && !is.null(sobj@interaction.depth)) s@interaction.depth <- sobj@interaction.depth
            #---------
            
            if (!is.null(interaction.depth)) s@interaction.depth <- interaction.depth
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@interaction.depth)) s@interaction.depth <- sobj@interaction.depth
              }
            }
            #---------
            if (!is.null(modelSettings) && inherits(modelSettings,'list')) {
              .ms <- names(modelSettings)
              if (!is.null(.ms)) {
                .ms <- .methodFix(.ms)
                if (!all(.ms %in% s@methods)) warning(paste('the models in the modelSettings:',paste0(names(modelSettings)[!.ms %in% s@methods],collapse = ', '),'are not selected in the methods, or do not exitst!'))
                w <- which(.ms %in% s@methods)
                if (length(w) > 0) {
                  .ms <- .ms[w]
                  modelSettings <- modelSettings[w]
                  names(modelSettings) <- .ms
                  ww <- c()
                  for (i in seq_along(.ms)) {
                    if(!inherits(modelSettings[[.ms[i]]],'list')) ww <- c(ww,i)
                  }
                  
                  if (length(ww) > 0) {
                    if (length(ww) < length(modelSettings)) {
                      warning(paste('the modelSettings for the items:',paste(.ms[ww],collapse = ','),'are not a list, and so they are ignored!'))
                      modelSettings <- modelSettings[-ww]
                    } else {
                      warning('the arguments for each method in the modelSettings should be introduced using a list; modelSettings is ignored!')
                      modelSetting <- NULL
                    }
                  }
                } else modelSetting <- NULL
                
              } else warning('modelSettings is not in the right structure, so it is ignored!')
              
              if (!is.null(modelSettings)) {
                s@modelSettings <- modelSettings
              }
            } else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@modelSettings)) s@modelSettings <- sobj@modelSettings
              }
            }
            #---------
            if (!is.null(seed)) {
              if (is.logical(seed)) seed <- sample(100000,1)
              else if (!is.numeric(seed)) seed <- NULL
              s@seed <- seed
            } else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@seed)) s@seed <- sobj@seed
              }
            }
            #-------------
            s
          }
)
#----------------
if (!isGeneric("sdm")) {
  setGeneric("sdm", function(formula,data,methods,...)
    standardGeneric("sdm"))
}

setMethod('sdm', signature(formula='ANY',data='sdmdata',methods='character'), 
          function(formula,data,methods,...) {
            a <- c('interaction.depth','n','replication','cv.folds','test.percent','bg','bg.n','var.importance','response.curve','var.selection','setting','ncore','modelSettings','seed','parallelSettings','filename')
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
            a <- c('interaction.depth','n','replication','cv.folds','test.percent','bg','bg.n','var.importance','response.curve','var.selection','setting','ncore','modelSettings','seed','parallelSettings','filename')
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
            
            if (missing(formula) && !all(s@sdmFormula@vars %in% data@features.name)) s@sdmFormula <- data@sdmFormula
            else if (inherits(formula,'sdmFormula')) s@sdmFormula <- formula
            else if (inherits(formula,'formula')) {  
              s@sdmFormula <- .exFormula(formula,as.data.frame(data)[,-1])
            }
            
            if (length(dot) > 0) {
                if (length(ndot) > 0) {
                  for (nd in ndot) {
                    if (nd == 'interaction.depth') interaction.depth <- dot[[nd]]
                    else if (nd == 'ncore') s@parallelSettings@ncore <- dot[[nd]]
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
                    else if (nd == 'parallelSettings' && is.list(dot[[nd]])) {
                      parallelSettings <- dot[[nd]]
                      nparallel <- names(parallelSettings)
                      a <- c('ncore','doParallel','method','cluster','hosts','fork','type')
                      nparallel <- .pmatch(nparallel,a)
                      w <- !is.na(nparallel)
                      if (length(w) > 0) {
                        parallelSettings <- parallelSettings[w]
                        nparallel <- nparallel[w]
                        names(parallelSettings) <- nparallel
                      }
                      #--
                      if ('ncore' %in% nparallel) s@parallelSettings@ncore <- parallelSettings$ncore
                      #--
                      if ('method' %in% nparallel) {
                        if (parallelSettings$method %in% c('parallel','foreach')) s@parallelSettings@method <- parallelSettings$method
                        else {
                          warning('parallelisation method is not recognised; the default value ("parallel") is used!')
                          s@parallelSettings@method <- 'parallel'
                        }
                      } else s@parallelSettings@method <- 'parallel'
                      #--
                      if ('fork' %in% nparallel) {
                        if (is.logical(parallelSettings$fork)) {
                          if (parallelSettings$fork && .is.windows()) {
                            warning('"fork" in parallelisation setting cannot be TRUE on Windows Operating Systems; It is changed to FALSE!')
                            s@parallelSettings@fork <- FALSE
                          } else s@parallelSettings@fork <- parallelSettings$fork
                        } else {
                          warning('"fork" in parallelisation setting should be logical; the default value is used!')
                          s@parallelSettings@fork <- !.is.windows()
                        }
                      } else s@parallelSettings@fork <- !.is.windows()
                      #--
                      if ('type' %in% nparallel) s@parallelSettings@type <- parallelSettings$type
                      #--
                      if ('doParallel' %in% nparallel && is.expression(parallelSettings$doParallel)) s@parallelSettings@doParallel <- parallelSettings$doParallel
                      #--
                      if ('cluster' %in% nparallel && inherits(parallelSettings$cluster,'cluster')) s@parallelSettings@cluster <- parallelSettings$cluster
                      #--
                      if ('hosts' %in% nparallel && is.character(parallelSettings$hosts)) s@parallelSettings@hosts <- parallelSettings$hosts
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
            a <- c('interaction.depth','n','replication','cv.folds','test.percent','bg','bg.n','var.importance','response.curve','var.selection','setting','ncore','modelSettings','seed','parallelSettings','filename')
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
                    else if (nd == 'ncore') s@parallelSettings@ncore <- dot[[nd]]
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
                    else if (nd == 'parallelSettings' && is.list(dot[[nd]])) {
                      parallelSettings <- dot[[nd]]
                      nparallel <- names(parallelSettings)
                      a <- c('ncore','doParallel','method','cluster','hosts','fork','type')
                      nparallel <- .pmatch(nparallel,a)
                      w <- !is.na(nparallel)
                      if (length(w) > 0) {
                        parallelSettings <- parallelSettings[w]
                        nparallel <- nparallel[w]
                        names(parallelSettings) <- nparallel
                      }
                      #--
                      if ('ncore' %in% nparallel) s@parallelSettings@ncore <- parallelSettings$ncore
                      
                      #--
                      if ('method' %in% nparallel) {
                        if (parallelSettings$method %in% c('parallel','foreach')) s@parallelSettings@method <- parallelSettings$method
                        else {
                          warning('parallelisation method is not recognised; the default value ("parallel") is used!')
                          s@parallelSettings@method <- 'parallel'
                        }
                      } else s@parallelSettings@method <- 'parallel'
                      #--
                      if ('fork' %in% nparallel) {
                        if (is.logical(parallelSettings$fork)) {
                          if (parallelSettings$fork && .is.windows()) {
                            warning('"fork" in parallelisation setting cannot be TRUE on Windows Operating Systems; It is changed to FALSE!')
                            s@parallelSettings@fork <- FALSE
                          } else s@parallelSettings@fork <- parallelSettings$fork
                        } else {
                          warning('"fork" in parallelisation setting should be logical; the default value is used!')
                          s@parallelSettings@fork <- !.is.windows()
                        }
                      } else s@parallelSettings@fork <- !.is.windows()
                      #--
                      if ('type' %in% nparallel) s@parallelSettings@type <- parallelSettings$type
                      #--
                      if ('doParallel' %in% nparallel && is.expression(parallelSettings$doParallel)) s@parallelSettings@doParallel <- parallelSettings$doParallel
                      #--
                      if ('cluster' %in% nparallel && inherits(parallelSettings$cluster,'cluster')) s@parallelSettings@cluster <- parallelSettings$cluster
                      #--
                      if ('hosts' %in% nparallel && is.character(parallelSettings$hosts)) s@parallelSettings@hosts <- parallelSettings$hosts
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
# .getModel.info <- function(x,w,...) {
#   if (missing(w) || is.null(w)) {
#     a <- c('species','method','replication','run')
#     w1 <- w2 <- w3 <- w4 <- TRUE
#     dot <- list(...)
#     if (length(dot) > 0) {
#       ndot <- names(dot)
#       ndot <- .pmatch(ndot,a)
#       w <- !is.na(ndot)
#       ndot <- ndot[w]
#       dot <- dot[w]
#       names(dot) <- ndot
#       for (nd in ndot) {
#         if (nd == 'species' && !is.null(dot[[nd]])) {
#           dot[[nd]] <- .pmatch(dot[[nd]],unique(as.character(x@run.info[,2])))
#           w1 <- x@run.info[,2] %in% dot[[nd]]
#         } else if (nd == 'method' && !is.null(dot[[nd]])) {
#           dot[[nd]] <- .methodFix(dot[[nd]])
#           w2 <- x@run.info[,3] %in% dot[[nd]]
#         }
#         else if (nd == 'replication' && !is.null(dot[[nd]])) {
#           if (length(x@replicates) != 0) {
#             dot[[nd]] <- .replicate.methodFix(dot[[nd]])
#             w3 <- x@run.info[,4] %in% dot[[nd]]
#           }
#         } else if (nd == 'run' && !is.null(dot[[nd]])) {
#           if (!is.null(dot[[nd]])) {
#             if (length(x@replicates) != 0) {
#               r <- unlist(lapply(x@replicates[[1]],function(x) x$method))
#               ru <- unique(r)
#               names(ru) <- ru
#               rID <- lapply(ru,function(x) which(r == x))
#               w4 <- c()
#               for (i in 1:length(rID)) {
#                 w4 <- c(w4,rID[[i]][c(1:length(rID[[i]])) %in% dot[[nd]]])
#               }
#               w4 <- x@run.info[,5] %in% w4
#             }
#           }
#         }
#       }
#       x@run.info[w1 & w2 & w3 & w4,]
#     } else x@run.info
#   } else x@run.info[x@run.info[,1] %in% w,]
# }
#--------
.getModel.info <- function(x,w=NULL,species=NULL,method=NULL,replication=NULL,run=NULL) {
  if (missing(w) || is.null(w)) {
    
    w1 <- w2 <- w3 <- w4 <- TRUE
    if (!is.null(species)) {
      species <- .pmatch(species,unique(as.character(x@run.info[,2])))
      w1 <- x@run.info[,2] %in% species
    }
    
    if (!is.null(method)) {
      method <- .methodFix(method)
      w2 <- x@run.info[,3] %in% method
    }
    
    
    if (!is.null(replication)) {
      if (length(x@replicates) != 0) {
        replication <- .replicate.methodFix(replication)
        w3 <- x@run.info[,4] %in% replication
      }
    }
    
    if (!is.null(run)) {
      if (!is.null(run)) {
        if (length(x@replicates) != 0) {
          r <- unlist(lapply(x@replicates[[1]],function(x) x$method))
          ru <- unique(r)
          names(ru) <- ru
          rID <- lapply(ru,function(x) which(r == x))
          w4 <- c()
          for (i in 1:length(rID)) {
            w4 <- c(w4,rID[[i]][c(1:length(rID[[i]])) %in% run])
          }
          
          w4 <- x@run.info[,5] %in% w4
        }
      }
    }
    x@run.info[w1 & w2 & w3 & w4,]
  } else x@run.info[x@run.info[,1] %in% w,]
}

#--------


.getModel.info2 <- function(x,w=NULL,species=NULL,method=NULL,replication=NULL,run=NULL,wtest=NULL) {
  # comparing to .getModel.info: In this, only one species is allowed!
  # x: sdmModels
  if (!is.null(w)) {
    mi <- .getModel.info(x,w)
  } else {
    mi <- .getModel.info(x)
    u <- as.character(unique(mi[,2]))
    m <- as.character(unique(mi[,3]))
    r <- unique(mi[,4])
    if (!is.null(species)) {
      if (length(species) > 1) {
        species <- species[1]
        warning('only the first species is considered!')
      }
      if (is.numeric(species)) {
        if (length(u) <= species) species <- u[species]
        else stop('The specified species is not recognised!')
      } else {
        species <- .pmatch(species,u)
        if (is.na(species)) stop('The specified species is not recognised!')
      }
    } else {
      if (length(u) > 1) stop('This object contains models for more than one species; in species argument spcify the name of species!')
      else species <- u
    }
    
    if (!is.null(method)) {
      method <- .sdmMethods$fixNames(method)
      wm <- method %in% m
      if (any(!wm)) {
        if (all(!wm)) stop('the specified methods do not exist in the object!')
        warning(paste('Methods',paste(method[!wm],collapse=', '),'do not exsit in the object, and are excluded!'))
        method <- method[wm]
      }
    } else method <- m
    
    if (!is.null(replication)) replication <- .replicate.methodFix(replication)
    else replication <- r
    
    mi <- .getModel.info(x,species=species,method=method,replication=replication,run=run)
  }
  
  mi
}
#---------
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
    .require('parallel')
    nc <- detectCores()
    if (woL$ncore < nc) nc <- woL$ncore
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
  if (.require("mraster")) {
    .memo <- eval(expression({memory(session=TRUE,echo=FALSE)}))
    .wf <- eval(parse(text = "mraster:::.change_unit(utils::object.size(woL$data@features)*150,'B','M')[[1]]")) # guessing the size of modelObj
    .mf <- floor(0.75 * ( .memo[2] /  .wf) )[[1]] # how many modelObj fits into 75% of the available memory!
    if (nrow(.run.info) > .mf) .ch <- ceiling(nrow(.run.info) / .mf)
    else .ch <- NULL
  } else {
    if (!is.null(filename)) warning('the package mraster is not installed (use installAll() to get it installed); this package helps when the output object size is big!')
    .ch <- NULL
  }
  
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
            if (nc > detectCores()) nc <- detectCores()
            cl <- makePSOCKcluster(nc)
          }
        }
      } else cl <- makePSOCKcluster(nc) # I should work more on providing the hostnames (on windows the connection needs to be throuy plink or PUTTY...)
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
        warning('for some reason, the clusters (for parallel processing) are not working!')
        cat('\n ncore is changed to 1!')
        .lapply <- .Lapply
      } else .lapply <- .parLapply
    }
  } else .lapply <- .Lapply
  #-------
  if (eval(expression({nc > 1 && .parMethod == 'foreach' && .require('foreach') && (.require('doParallel') | getDoParRegistered())}))) {
    .lapply <- .feLapply
    if (!eval(parse(text='getDoParRegistered()'))) {
      eval(expression({registerDoParallel(cl,nc)}))
    }
    
    if (eval(expression({getDoParWorkers() != nc}))) {
      eval(expression({warning(paste0('the number of workers registered for the foreach backend parallelisation is :',getDoParWorkers(),', which is different than the n.cores specified in the function! (n.cores is changed to this number)'))}))
      eval(expression({nc <- getDoParWorkers()}))
    }
  }
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