# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Oct. 2016
# Last Update :  July 2024
# Version 3.8
# Licence GPL v3
#-------------------


.ent <- function(x,...) {
  if (all(is.na(x))) return(NA)
  n1 <- length(which(x == 1))
  n0 <- length(which(x == 0))
  x <- c(n0,n1)/(n1+n0)
  if (any(x == 0)) {
    return (0)
  } else return(- sum(x*log(x)) / log(2))
}



#--------

.getWeights <- function(x, mi, wtest, id, stat, opt) {
  mi <- mi[mi$modelID %in% id,]
  if (all(mi[,wtest[1]])) {
    weight <- getEvaluation(x,id = id,wtest = wtest[1],stat = stat,opt = opt)[,2]
  } else {
    id1 <- mi$modelID[mi[,wtest[1]]]
    if (length(wtest) > 1) {
      w <- which(!mi[,wtest[1]])
      if (any(mi[w,wtest[2]])) {
        id2 <- mi$modelID[w][mi[w,wtest[2]]]
      } else id2 <- NULL
      if (length(c(id1,id2)) < length(id)) {
        if (length(wtest) == 3) {
          w <- id[which(!id %in% c(id1,id2))]
          w <- which(mi$modelID %in% w)
          if (any(mi[w,wtest[3]])) {
            id3 <- mi$modelID[w][mi[w,wtest[3]]]
          } else id3 <- NULL
        } else id3 <- NULL
      } else id3 <- NULL
    } else id2 <- id3 <- NULL
    
    id <- id[which(id %in% c(id1,id2,id3))]
    mi <- mi[mi$modelID %in% id,]
    id <- mi$modelID
    w1 <- getEvaluation(x,w = id1,wtest = wtest[1],stat = stat,opt = opt)[,2]
    if (!is.null(id2)) w2 <- getEvaluation(x,w = id2,wtest = wtest[2],stat = stat,opt = opt)[,2]
    else w2 <- NULL
    if (!is.null(id3)) w3 <- getEvaluation(x,w = id3,wtest = wtest[2],stat = stat,opt = opt)[,2]
    else w3 <- NULL
    idw <- c(id1,id2,id3)
    w <- c(w1,w2,w3)
    weight <- w[sapply(id,function(x) which(idw == x))]
  }
  
  list(weight=weight,mi=mi)
  
}
#-------

.Ensemble_Method <- function(method) {
  method <- tolower(method)
  .m <- c()
  for (i in seq_along(method)) {
    if (method[i] %in% c('weighted','unweighted','median','mad','pa','mean-weighted','mean-unweighted','median-weighted','median-unweighted','entropy','uncertainty','cv','ci','stdev')) {
      if (method[i] =='uncertainty') .m <- c(.m,'entropy')
      else .m <- c(.m,method[i])
    } else {
      if (method[i] %in% c('weight','weighted','weited','w','wei','weig','weigh','wait','waited')) .m <- c(.m,'weighted')
      else if (method[i] %in% c('unweight','unweighted','unweited','unw','unwaited','unwait','unweigh','u','uw','mean')) .m <- c(.m,'unweighted')
      else if (method[i] %in% c('median','med','meidan','medean','medan')) .m <- c(.m,'median')
      else if (method[i] %in% c('mad')) .m <- c(.m,'mad')
      else if (method[i] %in% c('pamean','pa','pa_mean','pa-mean')) .m <- c(.m,'pa')
      else if (method[i] %in% c('mean-weighted','mean_weighted','m-w','mean-w','u-w','unweighted-weighted','unweighted_weighted')) .m <- c(.m,'mean-weighted')
      else if (method[i] %in% c('mean-unweighted','mean_unweighted','m-u','mean-u','u-u','unweighted-unweighted','unweighted_unweighted','mean-mean','m-m','m_u')) .m <- c(.m,'mean-unweighted')
      else if (method[i] %in% c('median-weighted','median_weighted','med-w','median-w','med-weighted','med_w','med_weighted','median-weited','med-weited')) .m <- c(.m,'median-weighted')
      else if (method[i] %in% c('median-unweighted','median_unweighted','med-u','median-u','median-mean')) .m <- c(.m,'median-unweighted')
      else if (method[i] %in% c('cv','cvar','covar','coef.var','c.v')) .m <- c(.m,'cv')
      else if (method[i] %in% c('ci','c.i','confidence.interval','confidence')) .m <- c(.m,'ci')
      else if (method[i] %in% c('sd','sdev','stdev','s.d','s.dev','st.dev')) .m <- c(.m,'stdev')
      else if (method[i] %in% c('uncertainty','entropy','uncert','uncertain')) .m <- c(.m,'entropy')
      else {
        method[i] <- .pmatch(method[i],c('weighted','unweighted','median','pa','mean-weighted','mean-unweighted','median-weighted','median-unweighted','entropy','uncertainty','cv','ci','stdev'))
        if (!is.na(method[i]) && method[i] =='uncertainty') .m <- c(.m,'entropy')
        else .m <- c(.m,method[i])
      }
    }
  }
  #----
  .m <- .m[!is.na(.m)]
  
  .m
}
#===
.getID_From_ModelNames <- function(n) {
  o <- rep(NA,length(n))
  for (i in seq_along(n)) {
    a <- regexpr('([[:digit:]]+)',n[i])
    
    if (a > 0) {
      o[i] <- substr(n[i],a[[1]],a[[1]]-1+attributes(a)[[1]])
    }
  }
  as.numeric(o)
  
}
#---
# remove space from the beginning of items in a character vector
.removeSpace <- function(x) {
  sapply(x,function(x) {
    if (substr(x,1,1) == " ") substr(x,2,9999) 
    else x
  })
}

#-------
# To evaluate the expr argument in setting of ensemble and get modelIDs match condition specified by user:
.evalExpr <- function(x,m,opt=2,wtest=NULL) {
  if (is.character(x)) a <- str2lang(x)
  else if (is.expression(x)) a <- as.call(x)[[1]]
  
  stopifnot("expr in setting should be a chatater, or expression, or call (created by quote)"=is.call(a))
  
  b <- deparse(a)
  b <- paste(b,collapse = '')
  b <- gsub('\"',"'",b)
  #---
  
  .distribution <- m@setting@distribution[1]
  
  if (.distribution == 'binomial') .metrics <- c('AUC','COR','Deviance','sensitivity','specificity','TSS','Kappa','NMI','phi','ppv','npv','ccr','threshold')
  else .metrics <- c('RMSE','COR','MAE','Deviance')
  
  .st <- strsplit(b,' ')[[1]]
  .stw <- .pmatch(.st,.metrics)
  .w <- which(.stw %in% .metrics)
  if (length(.w) > 0) {
    ev <- getEvaluation(m,stat = .st[.w],wtest=wtest,opt = opt)
    for (i in .w) {
      colnames(ev)[which(colnames(ev) == .stw[i])] <- .st[i]
    }
    ev[eval(str2lang(b),ev),"modelID"]
  }
  
}
#----
# weighted mean for data.frame (length of w should be equal to ncol(df) and sum(w) == 1)
.wm <- function(df,w) {
  colSums(t(df) * w,na.rm = TRUE)
}
#----

# w: weight (for weighted) and threshold (for pa, entropy); index (for mean/median-unweighted) and a list(index,weight) for mean/median-weighted
.getEnsembleR <- function(x,w=NULL,alpha=0.05,method,filename='',overwrite=FALSE,...) {
  
  n <- list(...)
  
  if (method == 'weighted') {
    weighted.mean(x,w,na.rm=TRUE,filename=filename,overwrite=overwrite,...)
  } else if (method == 'unweighted') {
    if (length(n) > 0 && 'wopt' %in% names(n)) mean(x,na.rm=TRUE,filename=filename,overwrite=overwrite,wopt=n[['wopt']])
    else mean(x,na.rm=TRUE,filename=filename,overwrite=overwrite)
  } else if (method == 'median') {
    if (length(n) > 0 && 'wopt' %in% names(n)) median(x,na.rm=TRUE,filename=filename,overwrite=overwrite,wopt=n[['wopt']])
    else median(x,na.rm=TRUE,filename=filename,overwrite=overwrite)
  } else if (method == 'pa') {
    x <- ifel(x >= w, 1, 0)
    if (length(n) > 0 && 'wopt' %in% names(n)) mean(x,na.rm=TRUE,filename=filename,overwrite=overwrite,wopt=n[['wopt']])
    else mean(x,na.rm=TRUE,filename=filename,overwrite=overwrite)
  } else if (method == 'entropy') {
    x <- ifel(x >= w, 1, 0)
    if (length(n) > 0 && 'wopt' %in% names(n)) app(x,fun=.ent,filename=filename,overwrite=overwrite,wopt=n[['wopt']])
    else app(x,fun=.ent,filename=filename,overwrite=overwrite)
  } else if (method == 'mean-weighted') {
    x <- tapp(x, index=as.numeric(w[[1]]),fun='mean',na.rm=TRUE)
    weighted.mean(x,as.numeric(w[[2]]),na.rm=TRUE,filename=filename,overwrite=overwrite,...)
  } else if (method == 'mean-unweighted') {
    x <- tapp(x, index=as.numeric(w),fun='mean',na.rm=TRUE)
    if (length(n) > 0 && 'wopt' %in% names(n)) mean(x,na.rm=TRUE,filename=filename,overwrite=overwrite,wopt=n[['wopt']])
    else mean(x,na.rm=TRUE,filename=filename,overwrite=overwrite)
  } else if (method == 'median-weighted') {
    x <- tapp(x, index=as.numeric(w[[1]]),fun='median',na.rm=TRUE)
    weighted.mean(x,as.numeric(w[[2]]),na.rm=TRUE,filename=filename,overwrite=overwrite,...)
  } else if (method == 'median-unweighted') {
    x <- tapp(x, index=as.numeric(w),fun='median',na.rm=TRUE)
    if (length(n) > 0 && 'wopt' %in% names(n)) mean(x,na.rm=TRUE,filename=filename,overwrite=overwrite,wopt=n[['wopt']])
    else mean(x,na.rm=TRUE,filename=filename,overwrite=overwrite)
  } else if (method == 'cv') {
    .m <- mean(x,na.rm=TRUE)
    .s <- stdev(x, na.rm = TRUE)
    .cv <- .s / .m
    if (filename != '') {
      if (length(n) > 0 && 'wopt' %in% names(n)) {
        writeRaster(.cv,filename=filename,overwrite=overwrite,wopt=n[['wopt']])
      } else writeRaster(.cv,filename=filename,overwrite=overwrite)
    } else .cv
  } else if (method == 'ci') {
    .n <- nlyr(x)
    #.m <- mean(x,na.rm=TRUE)
    .se <- stdev(x, na.rm = TRUE) / sqrt(.n)
    .tsc = qt(p=alpha/2, df=.n-1,lower.tail=F)
    #.lo <- .m - (.tsc * .se)
    #.up <- .m + (.tsc * .se) 
    #names(.lo) <- 'ci_lower'
    #names(.up) <- 'ci_upper'
    .ci <- (.tsc * .se) * 2 # .upper - .lower (confidence interval length)
    names(.ci) <- 'CI'
    if (filename != '') {
      if (length(n) > 0 && 'wopt' %in% names(n)) {
        writeRaster(.ci,filename=filename,overwrite=overwrite,wopt=n[['wopt']])
      } else writeRaster(.ci,filename=filename,overwrite=overwrite)
    } else .ci
  } else if (method == 'stdev') {
    if (length(n) > 0 && 'wopt' %in% names(n)) stdev(x,na.rm=TRUE,filename=filename,overwrite=overwrite,wopt=n[['wopt']])
    else stdev(x,na.rm=TRUE,filename=filename,overwrite=overwrite)
  } else if (method == 'mad') {
    if (length(n) > 0 && 'wopt' %in% names(n)) app(x,fun='mad',filename=filename,overwrite=overwrite,wopt=n[['wopt']])
    else app(x,fun='mad',filename=filename,overwrite=overwrite)
  }
  
}
#-------
# w: weight (for weighted) and threshold (for pa, entropy); index (for mean/median-unweighted) and a list(index,weight) for mean/median-weighted
.getEnsembleDF <- function(x,w=NULL,alpha=0.05,method) {
  
  if (method == 'weighted') {
    .wm(x,w)
  } else if (method == 'unweighted') {
    rowMeans(x,na.rm=TRUE)
  } else if (method == 'median') {
    apply(x,1,'median',na.rm=TRUE)
  } else if (method == 'pa') {
    colMeans(ifelse(t(x) >= w,1,0),na.rm=TRUE)
  } else if (method == 'entropy') {
    x <- ifelse(t(x) >= w, 1, 0)
    apply(x,2,.ent)
  } else if (method == 'mean-weighted') {
    x <- t(apply(x,1,function(x,v) {
      tapply(x, INDEX=v,FUN='mean',na.rm=TRUE)
    },v=as.numeric(w[[1]])))
   
    .wm(x,as.numeric(w[[2]]))
  } else if (method == 'mean-unweighted') {
    x <- t(apply(x,1,function(x,v) {
      tapply(x, INDEX=v,FUN='mean',na.rm=TRUE)
    },v=as.numeric(w)))
    rowMeans(x,na.rm = TRUE)
    
  } else if (method == 'median-weighted') {
    x <- t(apply(x,1,function(x,v) {
      tapply(x, INDEX=v,FUN='median',na.rm=TRUE)
    },v=as.numeric(w[[1]])))
    
    .wm(x,as.numeric(w[[2]]))
  } else if (method == 'median-unweighted') {
    x <- t(apply(x,1,function(x,v) {
      tapply(x, INDEX=v,FUN='median',na.rm=TRUE)
    },v=as.numeric(w)))
    
    rowMeans(x,na.rm = TRUE)
  } else if (method == 'cv') {
    .m <- rowMeans(x,na.rm = TRUE)
    .s <- apply(x,1,'median',na.rm=TRUE)
    
    .s / .m
  } else if (method == 'ci') {
    .n <- apply(x,1,function(x)  length(x[!is.na(x)]))
    if (any(.n == 0)) .n[.n == 0] <- NA
    .m <- rowMeans(x,na.rm = TRUE)
    .s <- apply(x,1,'sd',na.rm=TRUE)
    .se <- .s / sqrt(.n)
    .un <- unique(.n)
    .un <- .un[!is.na(.un)]
    .tsc <- .n
    for (u in .un) {
      w <- which(.n == u)
      .tsc[w] <- rep(qt(p=alpha / 2, df=u-1,lower.tail=FALSE),length(w))
    }
    
    .lo <- .m - (.tsc * .se)
    .up <- .m + (.tsc * .se) 
    .ci <- data.frame(.lo,.up)
    colnames(.ci) <- c('ci_lower','ci_upper')
    .ci
  } else if (method == 'stdev') {
    apply(x,1,'sd',na.rm=TRUE)
  }  else if (method == 'mad') {
    apply(x,1,'mad',na.rm=TRUE)
  }
  
}
#-------
.ensSetting <- function(x, newdata,setting) {
  .distribution <- x@setting@distribution[1]
  
  if (missing(setting)) {
    if (.distribution == 'binomial') setting <- list(method='weighted',stat='AUC',power=1,expr=NULL,wtest=NULL)
    else setting <- list(method='weighted',stat='RMSE',power=1,expr=NULL,wtest=NULL)
  }
  
  if (!is.list(setting)) stop('setting should be defined as a list object!')
  
  n <- names(setting)
  if (is.null(n) || '' %in% n) stop('one (or more) argument in "setting" has no name!')
  
  n <- .pmatch(n,c('method','stat','opt','id','wtest','weight','power','expr','alpha'))
  if (any(is.na(n))) stop('At least one argument in setting is unkown!')
  names(setting) <- n
  #------------
  
  mi <- getModelInfo(x)
  mi <- mi[mi$success,]
  #-----------id:
  if ('id' %in% n) {
    id <- setting[['id']]
    if (!is.numeric(id)) {
      warning('id (in setting) should be numeric; it is ignored')
      id <- mi$modelID
    } else {
      w <- which(mi$modelID %in% id)
      if (length(w) < length(id)) {
        if (length(w) < 2) stop('At least 2 models are required for ensemble; (is the "id" in setting correct?)!')
        warning(paste(length(id) - length(w),'of the specified IDs (id) are not available in the model object, and are excluded'))
      }
      
      mi <- mi[w,]
      id <- mi$modelID
    }
  } else id <- mi$modelID
  
  #---------- opt
  if ('opt' %in% n) {
    opt <- setting[['opt']][1]
    if (!is.null(opt)) {
      if (is.numeric(opt)) {
        if (!opt %in% 1:10) {
          opt <- 2
          warning('opt (the criteria for optimum threshold) should be a number between 1:10; opt=2 is considered (i.e., max(se+sp))')
        }
      } else {
        th.criteria <- c("sp=se","max(se+sp)","min(cost)","minROCdist","max(kappa)","max(ppv+npv)","ppv=npv","max(NMI)","max(ccr)","prevalence")
        opt <- .pmatch(opt,th.criteria)
        opt <- opt[!is.na(opt)]
        if (length(opt) == 0) {
          warning('opt (the criteria for optimum threshold) is not understood! "max(se+sp)" is considered...')
          opt <- 2
        } else {
          opt <- which(th.criteria == opt)
        }
      } 
    } else opt <- 2
  } else opt <- 2
  #--------- wtest:
  if ('wtest' %in% n) {
    wtest <- setting[[which(n == 'wtest')]]
    
    if (!is.null(wtest)) {
      wtest <- tolower(wtest)
      wtest <- .pmatch(wtest,c('training','test.dep','test.indep'))[1]
    }
    if (is.null(wtest) || is.na(wtest) || length(wtest) == 0) {
      wtest <- colnames(mi)[9:7][which(c(length(which(mi[,9])),length(which(mi[,8])),length(which(mi[,7]))) > 0)]
    }
  } else wtest <- colnames(mi)[9:7][which(c(length(which(mi[,9])),length(which(mi[,8])),length(which(mi[,7]))) > 0)]
  
  #------
  
  if ('expr' %in% n && !is.null(setting[['expr']])) {
    if (is.character(setting[['expr']]) || is.expression(setting[['expr']]) || is.call(setting[['expr']])) {
      w <- try(.evalExpr(setting[['expr']],m=x,wtest = wtest,opt=opt),silent = TRUE)
      if (!inherits(w,'try-error')) {
        if (is.null(w)) warning('expr in setting seems not right (ignored)..!')
        else if (length(w) == 0) warning('expr issue: none of models are selected by evaluating the expr condition in setting so it is ignored...!')
        else if (length(w) == 1) warning('expr issue: only 1 model is selected by evaluating the expr condition in setting which is not enough for ensemble, so it is ignored...!')
        else {
          id <- w[w %in% id]
          mi <- mi[mi$modelID %in% id,]
        }
      } else warning('expr in setting is not correct (ignored)..!')
    } else warning('expr in setting should be either a character, an expression, or a call (created by the quote function), so it is ignored..!')
  }
  #-------
  
  .species <- as.character(unique(mi[,2]))
  .mi <- NULL
  if (length(.species) > 1) {
    .mi <- list()
    for (.sp in .species) {
      w <- mi[mi$species == .sp,]
      if (nrow(w) > 1) .mi[[.sp]] <- w
      else {
        mi <- mi[-which(mi$species == .sp),]
        .species <- .excludeVector(.species,.sp)
        warning(paste0('The model corresponding to species ',.sp,' is excluded because only one model was available for the species!'))
      }
    }
    #---
    if (length(.mi) == 0) stop('Not enough models (at least 2) are available for none of species!')
    else {
      if (length(.mi) == 1) rm(.mi)
      id <- mi$modelID
    }
  } else {
    if (nrow(mi) == 0) stop('There is no successfully fitted model in the input object!')
    if (nrow(mi) < 2) stop('sdmModels object contains 1 model but ensemble requires at least 2 models...!')
  }
  #-------- method:
  
  if ('method' %in% n) {
    method <- .Ensemble_Method(setting[['method']])
    if (length(method) == 0) {
      warning('the specified method is not recognised; default ("weighted") is considered')
      method <- 'weighted'
    }
  } else method <- 'weighted'
  
  #---- alpha:
  
  if ('ci' %in% method && 'alpha' %in% n) {
    .alpha <- setting[['alpha']][1]
    if (is.numeric(.alpha)) {
      if (.alpha < 1 && .alpha > 0.8) .alpha <- 1 - .alpha
      else if (.alpha > 1 && .alpha < 15) .alpha <- .alpha / 100
      else if (.alpha > 85 && .alpha < 100) .alpha <- 1 - (.alpha / 100)
      else if (.alpha >= 100 || .alpha > 0.15) {
        .alpha <- 0.05
        warning('alpha does not seem reasonable; default = 0.05 is considered!')
      }
    } else .alpha <- 0.05
  } else .alpha <- 0.05
  #--------- stat:
  if ('stat' %in% n) {
    if (.distribution == 'binomial') s <- c('AUC','COR','Deviance','sensitivity','specificity','TSS','Kappa','NMI','phi','ppv','npv','ccr')
    else s <- c('RMSE','MAE','COR','Deviance')
    stat <- tolower(setting[[which(n == 'stat')]])
    if (length(stat) > 1) {
      warning('only one statistic should be used. The first statistic in stat is used')
      stat <- stat[1]
    }
    stat <- .pmatch(stat,s)
    if (is.na(stat)) stop('stat is unknown!')
  } else {
    if (.distribution == 'binomial') stat <- 'auc'
    else stat <- 'rmse'
  }
  
  #---------- weight
  if ('weight' %in% n) {
    weight <- setting[['weight']]
    if (!is.numeric(weight)) stop('weight should be numeric!')
    if (length(weight) != length(id)) stop('the weight has a different length than the number of models in the input object (or the length of provided ids)!')
  } else weight <- NULL
  #--------------
  if ('power' %in% n) {
    .power <- setting[['power']]
    if (!is.numeric(.power) && any(method %in% c('weighted','mean-weighted','median-weighted'))) {
      warning('power should be numeric. It is ignored, and the default (1) is considered!')
      .power <- 1
    }
  } else .power <- 1
  #---------
  list(id=id,mi=mi,.mi=.mi,method=method,stat=stat,weight=weight,.species=.species,wtest=wtest,.power=.power,opt=opt,.alpha=.alpha)
  
}

#--------
if (!isGeneric("ensemble")) {
  setGeneric("ensemble", function(x,newdata,filename,setting,overwrite=FALSE,pFilename="",...)
    standardGeneric("ensemble"))
}	




setMethod('ensemble', signature(x='sdmModels',newdata='Raster'), 
          function(x, newdata, filename,setting,overwrite=FALSE,pFilename="",...) {
            if (missing(setting)) setting <- list(method='weighted',stat='AUC',wtest=NULL,power=1)
            
            if (missing(filename)) filename <- ''
            if (missing(pFilename)) pFilename <- ''
            
            if (missing(overwrite)) overwrite <- FALSE
            
            if (filename != '' && !overwrite) {
              if (file.exists(filename)) stop('The specified filename does exist. You may use overwrite = TRUE or use a different filename...!')
            }
            
            
            if (!is.list(setting)) stop('setting should be defined as a list object!')
            
            newdata <- rast(newdata)
            ens <- ensemble(x,newdata,filename=filename,setting=setting,overwrite = overwrite,pFilename = pFilename,...)
            
            as(ens,'Raster')
          }
)

#-------
setMethod('ensemble', signature(x='sdmModels',newdata='data.frame'), 
          function(x, newdata, filename,setting,overwrite=FALSE,pFilename="",...) {
            if (missing(setting)) setting <- list(method='weighted',stat='AUC',wtest=NULL,power=1)
            if (missing(filename)) filename=''
            if (missing(pFilename)) pFilename=''
            
            if (missing(overwrite)) overwrite <- FALSE
            
            if (filename != '' && !overwrite) {
              if (file.exists(filename)) stop('The specified filename does exist. You may use overwrite = TRUE or use a different filename...!')
            } 
            wtest <- id <- stat <- opt <- .power <- .mi <- .species <- NULL
            w <- .ensSetting(x,newdata,setting)
            for (n in names(w)) assign(n,w[[n]])
            
            .n <- colnames(newdata)
            p <- NULL
            if (!all(x@setting@featureFrame@predictors %in% .n)) {
              if (length(.n) == length(id)) {
                if (!all(grepl('id_([[:digit:]]+)__sp',.n))) {
                  stop('the newdata does not contain some or all the predictor variables required by the models...')
                  
                } else {
                  if (!all(.getID_From_ModelNames(.n) %in% id)) {
                    warning('It seems that the predictions (in newdata) for some modelIDs (id) are not available..!')
                  }
                }
                
                p <- newdata
              } else {
                if (all(grepl('id_([[:digit:]]+)__sp',.n))) {
                  .id <- .getID_From_ModelNames(.n)
                  if (!all(id %in% .id)) {
                    w <- which(id %in% .id)
                    if (length(w) > 1) {
                      warning(paste0('predictions of ',length(id) - length(w),' models are not available in newdata and so they don\'t contribute in the ensemble procedure!'))
                      id <- id[w]
                    } else stop('The predictions in newdata do not correspond to the selected modelIDs (id)!')
                  }
                  p <- newdata
                  
                  if (!all(.id %in% id)) {
                    .id <- .id[.id %in% id]
                    if (!all(.id == id)) {
                      if (all(sort(.id) == id)) .id <- sort(.id)
                      else warning('It seems that the id of predicted rasters (newdata) do not match with the modelIDs required for ensemble!')
                    }
                    if (length(.id) < 2) stop('At least 2 predictions are required for ensemble that are not available in newdata (the predictions)')
                    p <- p[,.id]
                  }
                  
                } else if (any(grepl('id_([[:digit:]]+)__sp',.n))) {
                  w <- which(grepl('id_([[:digit:]]+)__sp',.n))
                  
                  if (length(w) > 2) {
                    p <- newdata[,w]
                    .n <- .n[w]
                    
                    .id <- .getID_From_ModelNames(.n)
                    if (!all(id %in% .id)) {
                      w <- which(id %in% .id)
                      if (length(w) > 1) {
                        warning(paste0('predictions of ',length(id) - length(w),' models are not available in newdata and so they don\'t contribute in the ensemble procedure!'))
                        id <- id[w]
                      } else stop('The predictions in newdata do not correspond to the selected modelIDs (id)!')
                    }
                    
                    
                    if (!all(.id %in% id)) {
                      .id <- .id[.id %in% id]
                      if (!all(.id == id)) {
                        if (all(sort(.id) == id)) .id <- sort(.id)
                        else warning('It seems that the id of predicted rasters (newdata) do not match with the modelIDs required for ensemble!')
                      }
                      if (length(.id) < 2) stop('At least 2 predictions are required for ensemble that are not available in newdata (the predictions)')
                      p <- p[,.id]
                    }
                  } else stop('At least 2 predictions are required for ensemble that are not available in newdata (the predictions)')
                  
                } else if (all(grepl('sp_(.*)__m_',.n))) {
                  # here I should write the codes for the methods of mean-weighted or mean-unweighted
                  stop('the newdata does not contain the predictor variables required by the model...')
                } else stop('the newdata does not contain some or all the predictor variables required by the model...')
                
              }
            }
            #-----------------------------------
            
            if (is.null(p)) {
              
              p <- predict(x, newdata,id=id,filename=pFilename,overwrite=overwrite)
              #----
              if (ncol(p) == 1) stop('prediction for only one model is generated (by the predict function) but ensemble requires at least two...!')
              if (ncol(p) < length(id)) {
                warning(paste0('predictions for ',length(id)-ncol(p),' models are not generated correctly, so they are excluded from the ensemble procedure'))
                id <- .getID_From_ModelNames(colnames(p))
              }
            }
            
            ######################
            
            if ('ci' %in% method) {
              w <- which(method == 'ci')
              .method <- method[-w]
              .method <- c(.method,c('ci_lower','ci_upper'))
              .ens <- data.frame(matrix(NA,nrow=nrow(p),ncol=length(.species) * length(.method)))
              w <-c()
              for (.sp in .species) {
                for (j in .method) {
                  w <- c(w,paste0('ensemble_',.sp,'_',j))
                }
              }
            } else {
              .ens <- data.frame(matrix(NA,nrow=nrow(p),ncol=length(.species) * length(method)))
              if (ncol(.ens) > 1) {
                w <-c()
                for (.sp in .species) {
                  for (j in method) {
                    w <- c(w,paste0('ensemble_',.sp,'_',j))
                  }
                }
              } else {
                w <- paste0('ensemble_',method)
              }
            }
            
            colnames(.ens) <- w
            
            
            ####################
            if (is.null(weight)) {
              weight <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = stat, opt = opt)
              mi <- weight$mi
              id <- mi$modelID
              weight <- weight$weight
              if (stat %in% c('Deviance','RMSE','MAE')) {
                w <- max(weight,na.rm=TRUE)
                if (w == 0) weight <- rep(1,length(weight))
                else weight <- (w + 0.01*w) - weight
              }
            }
            #-----
            
            .id <- .getID_From_ModelNames(colnames(p))
            
            weight <- weight ^ .power
            weight <- weight / sum(weight)
            
            
            if (any(method %in% c('mean-weighted','mean-unweighted','median-weighted','median-unweighted'))) {
              method2 <- method[method %in% c('mean-weighted','mean-unweighted','median-weighted','median-unweighted')]
              method <- method[!method %in% c('mean-weighted','mean-unweighted','median-weighted','median-unweighted')]
              if (length(method) == 0) method <- NULL
            } else method2 <- NULL
            #------------
            if (!is.null(method) && 'ci' %in% method) {
              method <- method[method != 'ci']
              if (length(method) == 0) method <- NULL
              .ci <- TRUE
            } else .ci <- FALSE
            #------------
            if (!is.null(method)) {
              if (ncol(.ens) > 1) {
                for (.m in method) {
                  if (.m %in% c('pa','entropy')) {
                    .we <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = 'threshold', opt = opt)
                    .we <- .we$weight
                  } else if (.m == 'weighted') {
                    .we <- weight
                  } else .we <- NULL
                  
                  
                  if (length(.species) > 1) {
                    for (.sp in .species) {
                      w <- which(id %in% .mi[[.sp]]$modelID)
                      if (length(w) > 1) {
                        .idw <- which(.id %in% id[w])
                        .ens[paste0('ensemble_',.sp,'_',.m)] <- .getEnsembleDF(p[,.idw],w=.we[w],method = .m)
                      } else {
                        .ens <- .ens[,-which(colnames(.ens) == paste0('ensemble_',.sp,'_',.m)),drop=FALSE]
                      }
                    }
                  } else {
                    if (length(.id) > length(id)) {
                      .idw <- which(.id %in% id)
                      .ens[paste0('ensemble_',.species,'_',.m)] <- .getEnsembleDF(p[,.idw],w=.we,method = .m)
                    } else {
                      .ens[paste0('ensemble_',.species,'_',.m)] <- .getEnsembleDF(p,w=.we,method = .m)
                    }
                  }
                }
              } else {
                if (length(.id) > length(id)) {
                  .idw <- which(.id %in% id)
                  .ens[[1]] <- .getEnsembleDF(p[[.idw]],w=weight,method = method)
                } else {
                  .ens[[1]] <- .getEnsembleDF(p,w=weight,method = method)
                }
                
              }
            }
            #-----------
            #-----------
            if (!is.null(method2)) {
              if (ncol(.ens) > 1) {
                for (.m in method2) {
                  if (length(.species) > 1) {
                    for (.sp in .species) {
                      w <- which(id %in% .mi[[.sp]]$modelID)
                      if (length(w) > 1) {
                        .idw <- which(.id %in% id[w])
                        .ind <- as.numeric(as.factor(mi$method[w]))
                        if (.m %in% c('mean-weighted','median-weighted')) {
                          .we <- tapply(weight[w],.ind,'mean')
                          .we <- .we ^ .power
                          .we <- .we / sum(.we)
                          .ind <- list(.ind,.we)
                        }
                        .ens[[paste0('ensemble_',.sp,'_',.m)]] <- .getEnsembleDF(p[,.idw],w=.ind,method = .m)
                      } else {
                        .ens <- .ens[,-which(names(.ens) == paste0('ensemble_',.sp,'_',.m)),drop=FALSE]
                      }
                    }
                  } else {
                    .ind <- as.numeric(as.factor(mi$method))
                    if (.m %in% c('mean-weighted','median-weighted')) {
                      .we <- tapply(weight,.ind,'mean')
                      .we <- .we ^ .power
                      .we <- .we / sum(.we)
                      .ind <- list(.ind,.we)
                    }
                    
                    if (length(.id) > length(id)) {
                      .idw <- which(.id %in% id)
                      .ens[[paste0('ensemble_',.species,'_',.m)]] <- .getEnsembleDF(p[,.idw],w=.ind,method = .m)
                    } else {
                      .ens[[paste0('ensemble_',.species,'_',.m)]] <- .getEnsembleDF(p,w=.ind,method = .m)
                    }
                  }
                }
              } else {
                .ind <- as.numeric(as.factor(mi$method))
                if (method2 %in% c('mean-weighted','median-weighted')) {
                  .we <- tapply(weight,.ind,'mean')
                  .we <- .we ^ .power
                  .we <- .we / sum(.we)
                  .ind <- list(.ind,.we)
                }
                if (length(.id) > length(id)) {
                  .idw <- which(.id %in% id)
                  .ens[[1]] <- .getEnsembleDF(p[[.idw]],w=.ind,method = method2)
                } else {
                  .ens[[1]] <- .getEnsembleDF(p,w=.ind,method = method2)
                }
                
              }
            }
            #---
            if (.ci) {
              if (length(.species) > 1) {
                for (.sp in .species) {
                  w <- which(id %in% .mi[[.sp]]$modelID)
                  if (length(w) > 1) {
                    .idw <- which(.id %in% id[w])
                    .ens[paste0('ensemble_',.sp,'_',c('ci_lower','ci_upper'))] <- .getEnsembleDF(p[,.idw],w=.we[w],method = 'ci')
                    
                  } else {
                    .ens <- .ens[,-which(colnames(.ens) %in% paste0('ensemble_',.sp,'_',c('ci_lower','ci_upper'))),drop=FALSE]
                  }
                }
              } else {
                if (length(.id) > length(id)) {
                  .idw <- which(.id %in% id)
                  .ens[paste0('ensemble_',.species,'_',c('ci_lower','ci_upper'))] <- .getEnsembleDF(p[,.idw],method = 'ci')
                } else {
                  .ens[paste0('ensemble_',.species,'_',c('ci_lower','ci_upper'))] <- .getEnsembleDF(p,method = 'ci')
                }
              }
            }
            #-----------
            
            if (filename != '') write.csv(.ens,filename,row.names = FALSE)
            
            .ens
          }
)

#########

setMethod('ensemble', signature(x='sdmModels',newdata='SpatRaster'), 
          function(x, newdata, filename,setting,overwrite=FALSE,pFilename="",...) {
            if (missing(setting)) setting <- list(method='weighted',stat='AUC',power=1,expr=NULL,wtest=NULL)
            
            if (missing(filename)) filename <- ''
            if (missing(pFilename)) pFilename=''
            
            if (missing(overwrite)) overwrite <- FALSE
            
            if (filename != '' && !overwrite) {
              if (file.exists(filename)) stop('The specified filename does exist. You may use overwrite = TRUE or use a different filename...!')
            }
            
            wtest <- id <- stat <- opt <- .power <- .mi <- .species <- NULL
            
            w <- .ensSetting(x,newdata,setting)
            for (n in names(w)) assign(n,w[[n]])
            
            .n <- names(newdata)
            p <- NULL
            if (!all(x@setting@featureFrame@predictors %in% .n)) {
              if (length(.n) == length(id)) {
                if (!all(grepl('id_([[:digit:]]+)__sp',.n))) {
                  cat('\n ......... the input Raster object is considered as the predicted probabilities...\n')
                  names(newdata) <- paste0('id_',id,'__sp_',mi$species,'__m_',mi$method)
                } else {
                  if (!all(.getID_From_ModelNames(.n) %in% id)) {
                    warning('It seems that the predicted rasters (in newdata) for some modelIDs (id) are not available..!')
                  }
                }
                
                p <- newdata
              } else {
                if (all(grepl('id_([[:digit:]]+)__sp',.n))) {
                  .id <- .getID_From_ModelNames(.n)
                  if (!all(id %in% .id)) {
                    w <- which(id %in% .id)
                    if (length(w) > 1) {
                      warning(paste0('predictions of ',length(id) - length(w),' models are not available in newdata and so they don\'t contribute in the ensemble procedure!'))
                      id <- id[w]
                    } else stop('The predicted layers in newdata do not correspond to the selected modelIDs (id)!')
                  }
                  p <- newdata
                  
                  if (!all(.id %in% id)) {
                    .id <- .id[.id %in% id]
                    if (!all(.id == id)) {
                      if (all(sort(.id) == id)) .id <- sort(.id)
                      else warning('It seems that the id of predicted rasters (newdata) do not match with the modelIDs required for ensemble!')
                    }
                    p <- p[[.id]]
                  }
                  
                } else if (all(grepl('sp_(.*)__m_',.n))) {
                  # here I should write the codes for the methods of mean-weighted or mean-unweighted
                  stop('the newdata does not contain some or all the predictor variables required by the model...')
                } else stop('the newdata does not contain some or all the predictor variables required by the model...')
                
              }
            }
            #-----------------------------------
            
            if (is.null(p)) {
              
              p <- predict(x, newdata,id=id,filename=pFilename,overwrite=overwrite)
              #----
              if (nlyr(p) == 1) stop('Only 1 raster is generated by the predict function which is not enough for the ensemble procedure!')
              if (nlyr(p) < length(id)) {
                warning(paste0('predictions for ',length(id)-nlyr(p),' models are not generated correctly, so they are excluded from the ensemble procedure'))
                id <- .getID_From_ModelNames(names(p))
              }
            }
            
            ######################
            .ens <- rast(newdata,nlyrs=length(.species) * length(method))
            if (nlyr(.ens) > 1) {
              w <-c()
              for (.sp in .species) {
                for (j in method) {
                  w <- c(w,paste0('ensemble_',.sp,'_',j))
                }
              }
              set.values(.ens,1,NA,layer=1)
            } else {
              w <- paste0('ensemble_',method)
            }
            names(.ens) <- w
            
            ####################
            if (is.null(weight)) {
              weight <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = stat, opt = opt)
              mi <- weight$mi
              id <- mi$modelID
              weight <- weight$weight
              if (stat %in% c('Deviance','RMSE','MAE')) {
                w <- max(weight,na.rm=TRUE)
                if (w == 0) weight <- rep(1,length(weight))
                else weight <- (w + 0.01*w) - weight
              }
            }
            #-----
            
            .id <- .getID_From_ModelNames(names(p))
            
            weight <- weight ^ .power
            weight <- weight / sum(weight)
            
            if (any(method %in% c('mean-weighted','mean-unweighted','median-weighted','median-unweighted'))) {
              method2 <- method[method %in% c('mean-weighted','mean-unweighted','median-weighted','median-unweighted')]
              method <- method[!method %in% c('mean-weighted','mean-unweighted','median-weighted','median-unweighted')]
              if (length(method) == 0) method <- NULL
            } else method2 <- NULL
            #------------
            
            if (!is.null(method)) {
              if (nlyr(.ens) > 1) {
                for (.m in method) {
                  if (.m %in% c('pa','entropy')) {
                    .we <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = 'threshold', opt = opt)
                    .we <- .we$weight
                  } else if (.m == 'weighted') {
                    .we <- weight
                  } else .we <- NULL
                  
                  
                  if (length(.species) > 1) {
                    for (.sp in .species) {
                      w <- which(id %in% .mi[[.sp]]$modelID)
                      if (length(w) > 1) {
                        .idw <- which(.id %in% id[w])
                        .ens[[paste0('ensemble_',.sp,'_',.m)]] <- .getEnsembleR(p[[.idw]],w=.we[w],method = .m)
                      } else {
                        .ens <- .ens[[-which(names(.ens) == paste0('ensemble_',.sp,'_',.m))]]
                      }
                    }
                  } else {
                    if (length(.id) > length(id)) {
                      .idw <- which(.id %in% id)
                      .ens[[paste0('ensemble_',.species,'_',.m)]] <- .getEnsembleR(p[[.idw]],w=.we,method = .m)
                    } else {
                      .ens[[paste0('ensemble_',.species,'_',.m)]] <- .getEnsembleR(p,w=.we,method = .m)
                    }
                  }
                }
              } else {
                if (length(.id) > length(id)) {
                  .idw <- which(.id %in% id)
                  if (filename == '') .ens <- .getEnsembleR(p[[.idw]],w=weight,method = method)
                  else .ens <- .getEnsembleR(p[[.idw]],w=weight,method = method,filename = filename,overwrite = overwrite,...)
                } else {
                  if (filename == '') .ens <- .getEnsembleR(p,w=weight,method = method)
                  else .ens <- .getEnsembleR(p,w=weight,method = method,filename = filename,overwrite = overwrite,...)
                }
                names(.ens) <- paste0('ensemble_',method)
              }
            }
            #-----------
            #-----------
            if (!is.null(method2)) {
              if (nlyr(.ens) > 1) {
                for (.m in method2) {
                  if (length(.species) > 1) {
                    for (.sp in .species) {
                      w <- which(id %in% .mi[[.sp]]$modelID)
                      if (length(w) > 1) {
                        .idw <- which(.id %in% id[w])
                        .ind <- as.numeric(as.factor(mi$method[w]))
                        if (.m %in% c('mean-weighted','median-weighted')) {
                          .we <- tapply(weight[w],.ind,'mean')
                          .we <- .we ^ .power
                          .we <- .we / sum(.we)
                          .ind <- list(.ind,.we)
                        }
                        .ens[[paste0('ensemble_',.sp,'_',.m)]] <- .getEnsembleR(p[[.idw]],w=.ind,method = .m)
                      } else {
                        .ens <- .ens[[-which(names(.ens) == paste0('ensemble_',.sp,'_',.m))]]
                      }
                    }
                  } else {
                    .ind <- as.numeric(as.factor(mi$method))
                    if (.m %in% c('mean-weighted','median-weighted')) {
                      .we <- tapply(weight,.ind,'mean')
                      .we <- .we ^ .power
                      .we <- .we / sum(.we)
                      .ind <- list(.ind,.we)
                    }
                    
                    if (length(.id) > length(id)) {
                      .idw <- which(.id %in% id)
                      .ens[[paste0('ensemble_',.species,'_',.m)]] <- .getEnsembleR(p[[.idw]],w=.ind,method = .m)
                    } else {
                      .ens[[paste0('ensemble_',.species,'_',.m)]] <- .getEnsembleR(p,w=.ind,method = .m)
                    }
                  }
                }
              } else {
                .ind <- as.numeric(as.factor(mi$method))
                if (method2 %in% c('mean-weighted','median-weighted')) {
                  .we <- tapply(weight[w],.ind,'mean')
                  .we <- .we ^ .power
                  .we <- .we / sum(.we)
                  .ind <- list(.ind,.we)
                }
                if (length(.id) > length(id)) {
                  .idw <- which(.id %in% id)
                  if (filename == '') .ens <- .getEnsembleR(p[[.idw]],w=.ind,method = method2)
                  else .ens <- .getEnsembleR(p[[.idw]],w=.ind,method = method2,filename = filename,overwrite = overwrite,...)
                } else {
                  if (filename == '') .ens <- .getEnsembleR(p,w=.ind,method = method2)
                  else .ens <- .getEnsembleR(p,w=.ind,method = method2,filename = filename,overwrite = overwrite,...)
                }
                names(.ens) <- paste0('ensemble_',method2)
              }
            }
            
            if (filename != '' && nlyr(.ens) > 1) .ens <- writeRaster(.ens,filename,overwrite = overwrite,...)
            
            .ens
          }
)