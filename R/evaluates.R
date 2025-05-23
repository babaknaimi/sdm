# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  Jan. 2025
# Version 2.6
# Licence GPL v3
#--------



.statFix <- function(x) {
  if (is.numeric(x)) {
    x <- ifelse(x > 17,17,x)
    x <- ifelse(x < 1,1,x)
    x <- unique(x)
    x <- c('sensitivity','specificity','TSS','MCC','F1','Kappa','NMI','phi','ppv','npv','ccr','mcr','or','ommission','commission','predicted.prevalence')[x]
  } else {
    x <- tolower(x)
    for (i in seq_along(x)) {
      if (any(!is.na(pmatch(c("se"),x[i])))) x[i] <- 'sensitivity'
      else if (any(!is.na(pmatch(c("sp"),x[i])))) x[i] <- 'specificity'
      else if (any(!is.na(pmatch(c("ts"),x[i])))) x[i] <- 'TSS'
      else if (any(!is.na(pmatch(c("ka"),x[i])))) x[i] <- 'kappa'
      else if (any(!is.na(pmatch(c("nm"),x[i])))) x[i] <- 'NMI'
      else if (any(!is.na(pmatch(c("pp"),x[i])))) x[i] <- 'ppv'
      else if (any(!is.na(pmatch(c("np"),x[i])))) x[i] <- 'npv'
      else if (any(!is.na(pmatch(c("ccr"),x[i])))) x[i] <- 'ccr'
      else if (any(!is.na(pmatch(c("mcr"),x[i])))) x[i] <- 'mcr'
      else if (any(!is.na(pmatch(c("mcc"),x[i])))) x[i] <- 'MCC'
      else if (any(!is.na(pmatch(c("f1"),x[i])))) x[i] <- 'F1'
      else if (any(!is.na(pmatch(c("or"),x[i])))) x[i] <- 'or'
      else if (any(!is.na(pmatch(c("om"),x[i])))) x[i] <- 'ommission'
      else if (any(!is.na(pmatch(c("com"),x[i])))) x[i] <- 'commission'
      else if (any(!is.na(pmatch(c("pr"),x[i])))) x[i] <- 'predicted.prevalence'
      else if (any(!is.na(pmatch(c("ph"),x[i])))) x[i] <- 'phi'
    }
    x <- unique(x)
    w <- which(x %in% c('sensitivity','specificity','TSS','MCC','F1','Kappa','NMI','phi','ppv','npv','ccr','mcr','or','ommission','commission','predicted.prevalence'))
    x <- x[w]
  }
  x
}


.threshold <- function(o,p,th,stat=0) {
  if (missing(th)) th <- sort(unique(p))
  else th <- sort(unique(th))
  e <- matrix(nrow=length(th),ncol=18)
  colnames(e) <- c('threshold','sensitivity','specificity','TSS','MCC','F1','Kappa','NMI','phi','ppv','npv','ccr',
                   'mcr','or','ommission','commission','prevalence','obsPrevalence')
  
  e[,1] <- th
  for (i in seq_along(th)) {
    w <- which(p >= th[i])
    pt <- rep(0,length(p))
    pt[w] <- 1
    e[i,2:18] <- .evaluate.cmx(.cmx(o,pt))
  }
  
  w <- which(is.na(e[,"ppv"]) | is.na(e[,'npv']))
  if (length(w) > 0) e <- e[-w,,drop=FALSE]
  
  # 1: Se=SP
  w1 <- which.min(abs(e[,"sensitivity"] - e[,"specificity"]))
  # 2: Max(Se+Sp)
  w2 <- which.max(e[,"sensitivity"]+e[,"specificity"])
  # 3: Min.cost:
  w3 <- which.min((1 - e[,"sensitivity"])*e[,"prevalence"] + (1 - e[,"specificity"] )*(1 - e[,"prevalence"]))
  # 4: ROC:
  w4 <- which.min(((1 - e[,"sensitivity"])^2) + ((e[,"specificity"] -1)^2))
  # 5: Max(Kappa):
  w5 <- which.max(e[,"Kappa"])
  # 6: Max(npv+ppv) 
  w6 <- which.max(e[,"ppv"] + e[,"npv"])
  # 7: ppv=npv
  w7 <- which.min(abs(e[,"ppv"] - e[,"npv"]))
  # 8: Max(NMI):
  w8 <- which.max(e[,"NMI"])
  # 9: Max(CCR):
  w9 <- which.max(e[,"ccr"])
  # 10: PredictedPrevalence=ObservedPrevalence
  w10 <- which.min(abs(e[,"prevalence"] - e[,"obsPrevalence"]))
  # 11: Fixed_sensitivity:
  #w11 <- which.min(e[,"sensitivity"] > se)
  # Max(MCC):
  w11 <- which.max(e[,"MCC"])
  
  w <- which(o == 1)
  q <- quantile(p[w],probs=c(0.1,0.05,0.01,0),na.rm=TRUE)
  
  e2 <-matrix(nrow=length(q),ncol=18)
  colnames(e2) <- c('threshold','sensitivity','specificity','TSS','MCC','F1','Kappa','NMI','phi','ppv','npv','ccr',
                    'mcr','or','ommission','commission','prevalence','obsPrevalence')
  e2[,1] <- q
  for (i in seq_along(q)) {
    w <- which(p >= q[i])
    pt <- rep(0,length(p))
    pt[w] <- 1
    e2[i,2:18] <- .evaluate.cmx(.cmx(o,pt))
  }
  
  th.criteria <- c("sp=se","max(se+sp)","min(cost)","minROCdist","max(kappa)","max(ppv+npv)","ppv=npv","max(NMI)","max(ccr)","prevalence","max(MCC)","P10","P5","P1","P0")
  
  th <- e[c(w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11),unique(c(1,stat+1))]
  th2 <- e2[,unique(c(1,stat+1))]
  data.frame(criteria=th.criteria,round(rbind(th,th2),5))
}


.getEvalThresholds <- function(o,p,stat) {
  th <- sort(unique(p))
  e <- matrix(nrow=length(th),ncol=2)
  colnames(e) <- c('threshold',stat)
  
  e[,1] <- th
  for (i in seq_along(th)) {
    w <- which(p >= th[i])
    pt <- rep(0,length(p))
    pt[w] <- 1
    e[i,2] <- .evaluate.cmx(.cmx(o,pt))[stat]
  }
  if (stat == 'ppv' | stat == 'npv') e <- e[which(!is.na(e[,stat])),]
  e
}


#-----------------
#######################

.occurrence <- function(x) {
  if (is.logical(x)) {
    xx <- rep(0,length(x))
    xx[which(x)] <- 1
    return(xx)
  } else if (is.character(x) | is.factor(x)) {
    return(as.numeric(as.character(x)))
  } else return(x)
}
#--------
.deviance_binomial <- function(obs,pr,w=NULL) {
  pr <- ifelse(pr <= 0,0.001,pr)
  pr <- ifelse(pr >= 1,0.999,pr)
  o <- log( pr*obs + (1-pr)*(1-obs) )
  o <- o[o != -Inf]
  if (!is.null(w) && length(w) == length(o)) {
    w <- w / sum(w,na.rm=FALSE)
    o <- o * w
  }
  (-2 * sum(o)) / length(o)
}
#-------
.deviance_poisson <- function(obs,pr,w=NULL) {
  o <- ifelse(obs == 0, 0, (obs * log(obs/pr))) - (obs - pr)
  if (!is.null(w) && length(w) == length(o)) {
    w <- w / sum(w,na.rm=FALSE)
    o <- o * w
  }
  (-2 * sum(o)) / length(o)
}
#---------
.deviance_gaussian <- function(obs,pr,w=NULL) {
  o <- (obs - pr) * (obs - pr)
  if (!is.null(w) && length(w) == length(o)) {
    w <- w / sum(w,na.rm=FALSE)
    o <- o * w
  }
  sum(o) / length(o)
}
#----------

.deviance_laplace <- function(obs,pr,w=NULL) {
  o <- abs(obs - pr)
  if (!is.null(w) && length(w) == length(o)) {
    w <- w / sum(w,na.rm=FALSE)
    o <- o * w
  }
  sum(o) / length(o)
}
#------------
.cmx <- function(o,p) {
  cmx<-matrix(nrow=2,ncol=2)
  colnames(cmx) <- rownames(cmx) <- c('P','A')
  cmx[1,1] <- length(which(o == 1 & p == 1))
  cmx[2,2] <- length(which(o == 0 & p == 0))
  cmx[1,2] <- length(which(o == 0 & p == 1))
  cmx[2,1] <- length(which(o == 1 & p == 0))
  cmx[] <- as.numeric(cmx)
  cmx
}
#-------------

.evaluate.cmx <- function(cmx) {
  TP<-cmx[1,1];FP<-cmx[1,2];TN<-cmx[2,2];FN<-cmx[2,1]
  N <- sum(cmx)
  prev <- (TP+FN) / N
  pred.prev <- (TP + FP) / N
  ccr <- (TP + TN) / N
  sens <- TP / (TP + FN)
  spec <- TN / (FP + TN)
  ppv <- TP / (TP + FP)
  npv <- TN / (TN + FN)
  or <- (TP * TN) / (FP * FN)
  commission = FP/(FP + TN)
  ommission = FN/(TP + FN)
  mcr = (FP + FN) / N
  phi <- (TP*TN - FP*FN)/(sqrt((TP+FN)*(TN+FP)*(TP+FP)*(TN+FN)))
  
  Kappa <- ((TN + TP) - ((((TP + FN)*(TP + FP)) + ((FP + TN)*(FN + TN))) /N)) / 
    (N-((((TP+FN)*(TP+FP))+((FP+TN)*(FN+TN)))/N))
  
  mcc <- ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  f1s <- (2 * TP) / ((2*TP) + FP + FN)
  v <- as.vector(cmx) ; pp <- apply(cmx,1,sum) ; op <- apply(cmx,2,sum)
  v <- ifelse(v == 0,0.00001,v)
  pp <- ifelse(pp == 0,0.00001,pp)
  op <- ifelse(op == 0,0.00001,op)
  NMI <- 1-((-sum(v*log(v)) + sum(pp*log(pp))) / (N*log(N) - sum(op*log(op))))
  TSS <- sens+spec-1
  
  return(round(c(sensitivity=sens,specificity=spec,TSS=TSS,MCC=mcc,F1=f1s,Kappa=Kappa,NMI=NMI,phi=phi,ppv=ppv,npv=npv,ccr=ccr,
                 mcr=mcr,or=or,ommission=ommission,commission=commission,predicted.prevalence=pred.prev,prevalence=prev),4))
}
#---- 
#-----------
.auc <- function(o,p) {
  w1 <- which(o == 1)
  w0 <- which(o == 0)
  auc <- as.numeric(NA)
  if (length(w1) > 0 & length(w0) > 0) {
    auc <- as.vector(wilcox.test(p[w1], p[w0])$statistic)/(length(w1)*length(w0))
  }
  round(auc,3)
}
#-----

.cor <- function(o,p,method='pearson') {
  co <- try(cor(p, o, method=method), silent = TRUE)
  co.t <- try(cor.test(p, o, method=method), silent = TRUE)
  if (!inherits(co,"try-error")) {
    co <- round(co,3)
    if (!inherits(co.t,"try-error")) return(c(cor=co,p.value=co.t$p.value))
    else return(co)
  } else return(as.numeric(NA))
}
#-------
.rmse <- function(o,p) {
  o <- o - p
  sqrt(mean(o * o))
}

.mae <- function(o,p) {
  o <- o - p
  mean(abs(o))
}
#-----------


# computation part of CBI (continuous boyce index):
.cbi <- function(x,pr,wi,wh,method='spearman') {
  
  .nx <- length(x) # number of presence points
  .ny <- length(pr) # number of background points
  
  Fi <- rep(NA,length(wi))
  
  for (i in seq_along(wi)) {
    pi <- sum(as.numeric(x >= wi[i] & x <= wh[i])) / .nx
    ei <- sum(as.numeric(pr >= wi[i] & pr <= wh[i])) / .ny
    Fi[i] <- pi / ei
  }
  #-------
  .msuit <- (wi + wh) / 2 # mean suitability for each bin!
  
  
  .keep <- which(Fi != "NaN" & !is.na(Fi) & !is.infinite(Fi))
  #Fi <- Fi[.keep]
  .keep <- .keep[which(Fi[.keep] != c(Fi[.keep][-1],TRUE))]
  .stat <- cor(.msuit[.keep],Fi[.keep],method = method,use='complete.obs')
  
  list(CBI=.stat,Fi=Fi[.keep],suitability=.msuit)
}
#-----------
.boyce <- function(obs, pr, win=NULL,n=101,method='spearman') {
  # continuous boyce index
  # obs: presence-background data (background should represent the whole study area)
  # pr: predicted suitability corresponding to obs
  # win: windows width (default: NULL -> 10% of range)
  # n: number of bins (default: 101 )
  # correlation method (default: "spearman")
  #-------------
  
  
  if (length(obs) != length(pr)) stop('lengths of obs and pr are not the same!')
  
  if (length(which(pr > 0 & pr < 1)) == 0) stop('"pr" should be a vector of predicted probabilities!')
  
  
  
  .min <- min(pr,na.rm = TRUE)
  .max <- max(pr, na.rm = TRUE)
  
  if (is.null(win)) {
    win <- (.max - .min) / 10
  }
  #------
  if (length(win) > 1 || win > 0.9 || win < 0) {
    win <- (.max - .min) / 10
    warning(paste0('the selected value for "win" seems not reasonable and it is ignored; default (',round(win,4),') is replaced!'))
  }
  #--------
  wi <- seq(.min, .max - win, length.out = n)
  wh <- wi + win
  #----------------
  .x <- pr[obs == 1] # probability of presence locations
  
  .cbi(.x,pr,wi,wh,method=method)
}



#-----------
# 
# if (!isGeneric("sdmEvaluate")) {
#   setGeneric("sdmEvaluate", function(data,model,...)
#     standardGeneric("sdmEvaluate"))
# }  
# 
# 
# setMethod('sdmEvaluate', signature(model='sdmModel'),
#           function(data, model,...) {
#             #
#           }
# )


.evaluate.binomial <- function(o,p) {
  e <- new('sdmEvaluate')
  e@observed <- o
  e@predicted <- p
  e@statistics[['Prevalence']] <- round(length(which(o == 1)) / length(o),3)
  e@statistics[['AUC']] <- .auc(o,p)
  e@statistics[['COR']] <- .cor(o,p)
  e@statistics[['Deviance']] <- .deviance_binomial(o,p)
  e@threshold_based <- .threshold(o,p,stat = c(1:12,14:16))
  e
}
#-------

.evaluate.numerical <- function(o,p,distribution='gaussian') {
  e <- new('sdmEvaluate')
  e@observed <- o
  e@predicted <- p
  e@statistics[['COR']] <- .cor(o,p)
  if (distribution %in% c('gaussian','normal')) e@statistics[['Deviance']] <- .deviance_gaussian(o,p)
  else if (distribution == 'poisson') e@statistics[['Deviance']] <- .deviance_poisson(o,p)
  else if (distribution == 'laplase') e@statistics[['Deviance']] <- .deviance_laplace(o,p)
  e@statistics[['RMSE']] <- .rmse(o,p)
  e@statistics[['MAE']] <- .mae(o,p)
  e
}
#--------

#--------
if (!isGeneric("evaluates")) {
  setGeneric("evaluates", function(x,p,...)
    standardGeneric("evaluates"))
}  

setMethod('evaluates', signature(x='vector',p='vector'),
          function(x, p,distribution) {
            w <- which(!is.na(x) & !is.na(p))
            p <- p[w]; x <- x[w]
            if (length(x) != length(p)) stop('observed and predicted vectors should have the same length')
            if (missing(distribution) || is.null(distribution)) {
              if (.isBinomial(x)) distribution <- 'binomial'
              else {
                # guessing!!
                if (mean(round(x,0) - x) == 0) distribution <- 'poisson'
                else distribution <- 'gaussian'
              }
            }
            
            if (distribution %in% c('binomial','bernouli')) {
              x <- .occurrence(x)
              .evaluate.binomial(x,p)
            } else {
              .evaluate.numerical(x,p,distribution = distribution)
            }
          }
)


.extractEvaluation <- function(x,id,wtest=NULL,stat=NULL,opt=NULL) {
  mi <- x@run.info
  mi <- mi[mi$success,]
  if (missing(id) || is.null(id)) id <- mi$modelID
  if (is.null(wtest)) wtest <- colnames(mi)[9:7][which(as.matrix(mi[1,c(9,8,7)]))[1]]
  else {
    wtest <- .pmatch(wtest,c('training','test.dep','test.indep'))[1]
    if (is.na(wtest)) wtest <- colnames(mi)[9:7][which(as.matrix(mi[1,c(9,8,7)]))[1]]
  }
  
  .dist <- x@setting@distribution[1]
  
  
  if (.dist == 'binomial') {
    s1 <- c('AUC','COR','Deviance','obs.prevalence')
    
    s2 <- c('threshold','sensitivity','specificity','TSS','MCC','F1','Kappa','NMI','phi','ppv','npv','ccr','mcr','ommission', 'commission', 'prevalence')
    
  } else {
    s1 <- c('RMSE','COR','MAE','Deviance')
    
    s2 <- NULL
  }
  
  
  
  
  if (!is.null(stat)) {
    stat <- .pmatch(stat,c(s1,s2))
    stat <- stat[!is.na(stat)]
    if (length(stat) == 0) stat <- c('AUC','COR','Deviance','TSS','MCC','F1')
  } else {
    if (.dist == 'binomial') stat <- c('AUC','COR','Deviance','TSS','MCC','F1')
    else stat <- c('RMSE','COR','MAE','Deviance')
  }
  
  s1 <- stat[stat %in% s1]
  if (length(s1) == 0) s1 <- NULL
  else if ('obs.prevalence' %in% s1) s1[s1 == 'obs.prevalence'] <- 'Prevalence'
  
  if (!is.null(s2)) {
    s2 <- stat[stat %in% s2]
    if (length(s2) == 0) s2 <- NULL
  }
  
  th.criteria <- c("sp=se","max(se+sp)","min(cost)","minROCdist","max(kappa)","max(ppv+npv)","ppv=npv","max(NMI)","max(ccr)","prevalence","max(MCC)","P10","P5","P1","P0")
  if (!is.null(opt)) {
    if (is.numeric(opt)) {
      if (!opt %in% 1:15) {
        opt <- 2
        warning('opt (the criteria for optimum threshold) should be a number between 1:15; opt=2 is considered (i.e., max(se+sp))')
      }
    } else {
      opt <- .pmatch(opt,th.criteria)
      opt <- opt[!is.na(opt)]
      if (length(opt) == 0) {
        warning('opt (the criteria for optimum threshold) is not understood! max(se+sp) is considered...')
        opt <- 2
      } else {
        opt <- opt[1]
        opt <- which(th.criteria == opt)
      }
    } 
  } else opt <- 2
  
  mi <- mi[mi[,wtest],]
  
  w <- mi$modelID %in% id
  mi <- mi[w,]
  id <- as.character(mi$modelID)
  rownames(mi) <- id
  
  for (i in c(2,3,4)) mi[,i] <- as.character(mi[,i])
  
  o <- list()
  names(s1) <- s1
  names(s2) <- s2
  
  for (i in id) {
    sp <- mi[i,2]
    mo <- mi[i,3]
    r <- mi[i,4]
    if (!is.null(s1)) {
      o[[i]] <- lapply(s1,function(j) x@models[[sp]][[mo]][[i]]@evaluation[[wtest]]@statistics[[j]][[1]])
    }
    if (!is.null(s2)) {
      o[[i]] <- c(o[[i]],lapply(s2,function(j) x@models[[sp]][[mo]][[i]]@evaluation[[wtest]]@threshold_based[opt,j]))
    }
  }
  o
}

._getPerformance <- function(x,id,wtest=NULL,s1=NULL,s2=NULL,opt=2) {
  #internal function: given the exact parameters, the stat values are extracted from a model
  
  mi <- x@run.info
  w <- mi$modelID %in% id
  mi <- mi[w,]
  
  id <- as.character(mi$modelID)
  rownames(mi) <- id
  
  if (is.null(wtest)) wtest <- colnames(mi)[9:7][which(as.matrix(mi[1,c(9,8,7)]))[1]]
  else {
    wtest <- .pmatch(wtest,c('training','test.dep','test.indep'))[1]
    if (is.na(wtest)) wtest <- colnames(mi)[9:7][which(as.matrix(mi[1,c(9,8,7)]))[1]]
  }
  
  for (i in c(2,3,4)) mi[,i] <- as.character(mi[,i])
  
  o <- list()
  names(s1) <- s1
  names(s2) <- s2
  
  ws <- vector('list',length(c(s1,s2))) # for the case of NA in the evaluation
  names(ws) <- c(s1,s2)
  ws[] <- NA
  for (i in id) {
    sp <- mi[i,2]
    mo <- mi[i,3]
    if (!is.null(x@models[[sp]][[mo]][[i]]@evaluation[[wtest]])) {
      if (!is.null(s1)) {
        o[[i]] <- lapply(s1,function(j) x@models[[sp]][[mo]][[i]]@evaluation[[wtest]]@statistics[[j]][[1]])
      }
      if (!is.null(s2)) {
        o[[i]] <- c(o[[i]],lapply(s2,function(j) x@models[[sp]][[mo]][[i]]@evaluation[[wtest]]@threshold_based[opt,j]))
      }
    } else {
      o[[i]] <- c(o[[i]],ws)
    }
  }
  o
}
###############
.getEvalTable <- function(x,id,wtest,stat=NULL) {
  # This is to get evaluation table for a single model!
  # stat can be 1 (threshold-independent) OR 2 (threshold-dependent)
  mi <- x@run.info
  w <- which(mi$modelID == id)
  if (length(w) == 1) {
    
    if (missing(wtest) || is.null(wtest)) wtest <- colnames(mi)[9:7][which(as.matrix(mi[1,c(9,8,7)]))[1]]
    else {
      wtest <- .pmatch(wtest,c('training','test.dep','test.indep'))[1]
      if (is.na(wtest)) wtest <- colnames(mi)[9:7][which(as.matrix(mi[1,c(9,8,7)]))[1]]
    }
    
    sp <- as.character(mi$species)[w]
    mo <- as.character(mi$method)[w]
    i <- as.character(mi$modelID)[w]
    if (is.null(stat)) list(threhsold_independent=x@models[[sp]][[mo]][[i]]@evaluation[[wtest]]@statistics,
           threshold_based=x@models[[sp]][[mo]][[i]]@evaluation[[wtest]]@threshold_based)
    else if (stat == 1) x@models[[sp]][[mo]][[i]]@evaluation[[wtest]]@statistics
    else if (stat == 2) x@models[[sp]][[mo]][[i]]@evaluation[[wtest]]@threshold_based
    
  } 
  
}

#--------
if (!isGeneric("getEvaluation")) {
  setGeneric("getEvaluation", function(x,id,wtest,stat,opt,...)
    standardGeneric("getEvaluation"))
}  

setMethod('getEvaluation', signature(x='sdmModels'),
          function(x, id, wtest,stat,opt,...) {
            .dist <- x@setting@distribution[1]
            
            if (missing(id)) id <- NULL
            if (missing(wtest)) wtest <- NULL
            if (missing(stat)) {
              if (.dist == 'binomial') stat <-c('AUC','COR','Deviance','TSS','MCC','F1')
              else stat <- c('RMSE','COR','MAE','Deviance')
            }
            if (missing(opt)) opt <- 2
            
            if (!is.null(id) && length(id) == 1 && all(stat %in% 1:2)) {
              .getEvalTable(x,id = id,wtest=wtest,stat=stat)
            } else {
              e <- .extractEvaluation(x,id=id,wtest=wtest,stat=stat,opt=opt)
              if (length(e) > 0) {
                stat <- names(e[[1]])
                o <- data.frame(matrix(ncol=length(stat)+1,nrow=length(e)))
                colnames(o) <- c('modelID',stat)
                if (length(stat) > 1) o[,2:ncol(o)] <- t(sapply(e,function(x) unlist(x)))
                else o[,2] <- sapply(e,function(x) unlist(x))
                o[,1] <- as.numeric(names(e))
                o
              } else stop('No evaluation results for the specified models!')
            }
          }
)
#--------------------------------
# evaluation of predicted raster (e.g., ensemble) dataset:
setMethod('evaluates', signature(x='sdmdata',p='RasterLayer'),
          function(x, p,distribution,wtest=NULL,...) {
            if (is.null(x@info) || is.null(x@info@coords)) stop('sdmdata object does not contain spatial coordintes!')
            
            if (missing(distribution)) {
              if (x@species[[1]]@type %in% c('Abundance')) distribution <- 'poisson'
              else if (x@species[[1]]@type %in% c('Numerical')) distribution <- 'gaussian'
              else if (x@species[[1]]@type %in% c('Presence-Absence','Presence-Background')) distribution <- 'binomial'
              else distribution <- NULL
            }
            if (missing(wtest)) wtest <- NULL
            
            if ('test' %in% .getGroupNames(x,levels=TRUE)) {
              if (!is.null(wtest)) {
                if (wtest %in% c('test','test.indep','test.dep')) d <- .getDataFrame(x,grp='test',...)
                else if (wtest %in% c('train','training')) d <- .getDataFrame(x,grp='train',...)
                else {
                  warning('"test" is considered for wtest (it should be either "train" or "test")...!')
                  d <- .getDataFrame(x,grp='test',...)
                } 
              } else d <- .getDataFrame(x,grp='test',...)
            } else d <- .getDataFrame(x,...)
            #--------
            if (length(x@species.names) > 1) {
              w <- which(x@species.names %in% colnames(d))
              if (length(w) > 1) {
                stop('sdmdata object has multiple species. specify the name of species through the spn argument...')
              }
               spn <- x@species.names[w] 
            } else spn <- x@species.names
            
            o <- d[,spn]
            p <- extract(p,d[,colnames(coords(x))])
            
            evaluates(o, p, distribution = distribution)
          }
)
#----
setMethod('evaluates', signature(x='sdmdata',p='SpatRaster'),
          function(x, p,distribution,wtest=NULL,...) {
            if (is.null(x@info) || is.null(x@info@coords)) stop('sdmdata object does not contain spatial coordintes!')
            
            if (missing(distribution)) {
              if (x@species[[1]]@type %in% c('Abundance')) distribution <- 'poisson'
              else if (x@species[[1]]@type %in% c('Numerical')) distribution <- 'gaussian'
              else if (x@species[[1]]@type %in% c('Presence-Absence','Presence-Background')) distribution <- 'binomial'
              else distribution <- NULL
            }
            
            if (missing(wtest)) wtest <- NULL
            
            if ('test' %in% .getGroupNames(x,levels=TRUE)) {
              if (!is.null(wtest)) {
                if (wtest %in% c('test','test.indep','test.dep')) d <- .getDataFrame(x,grp='test',...)
                else if (wtest %in% c('train','training')) d <- .getDataFrame(x,grp='train',...)
                else {
                  warning('"test" is considered for wtest (it should be either "train" or "test")...!')
                  d <- .getDataFrame(x,grp='test',...)
                } 
              } else d <- .getDataFrame(x,grp='test',...)
            } else d <- .getDataFrame(x,...)
            #--------
            if (length(x@species.names) > 1) {
              w <- which(x@species.names %in% colnames(d))
              if (length(w) > 1) {
                stop('sdmdata object has multiple species. specify the name of species through the spn argument...')
              }
              spn <- x@species.names[w] 
            } else spn <- x@species.names
            
            o <- d[,spn]
            p <- extract(p,d[,colnames(coords(x))],ID=FALSE)[,1]
            
            evaluates(o, p, distribution = distribution)
          }
)
#------------
setMethod('evaluates', signature(x='sdmModels',p='SpatRaster'),
          function(x, p,distribution,wtest=NULL,...) {
            if (missing(distribution)) {
              if (x@data@species[[1]]@type %in% c('Abundance')) distribution <- 'poisson'
              else if (x@data@species[[1]]@type %in% c('Numerical')) distribution <- 'gaussian'
              else if (x@data@species[[1]]@type %in% c('Presence-Absence','Presence-Background')) distribution <- 'binomial'
              else distribution <- NULL
            }
            
            if (missing(wtest)) wtest <- NULL
            
            evaluates(x@data,p,distribution=distribution,wtest=wtest,...)
          }
)

#-------------
setMethod('evaluates', signature(x='sdmModels',p='missing'),
          function(x, p,distribution,wtest=NULL,species=NULL,...) {
            if (missing(distribution)) {
              if (x@data@species[[1]]@type %in% c('Abundance')) distribution <- 'poisson'
              else if (x@data@species[[1]]@type %in% c('Numerical')) distribution <- 'gaussian'
              else if (x@data@species[[1]]@type %in% c('Presence-Absence','Presence-Background')) distribution <- 'binomial'
              else distribution <- NULL
            }
            
            if (missing(wtest)) wtest <- NULL
            
            
            if (missing(species)) species <- NULL
            
            
            n <- x@data@species.names
            
            if (length(n) > 1) {
              if (is.null(species)) {
                warning('multiple species are available in the model and the species is not specified in the "species" argument. The first species is considered!')
                n <- n[1]
              } else {
                if (length(species) > 1) stop('only one species should be specified in the "species" argument!')
                else {
                  if (is.numeric(species)) {
                    if (length(n) < species) stop('The specified species is not available (species argument)!')
                    else n <- n[species]
                  } else if (is.character(species)) {
                    w <- which(n == species)
                    if (length(w) == 1) n <- n[w]
                    else stop('The specified species (in the "species" argument) is not available!')
                  } else stop('"species" should be a character or numeric!')
                }
              }
            }
            #---------
            en <- as.data.frame(x@data,sp=n)
            
            obs <- en[,n]
            
            en <- ensemble(x, en,...)
            if (length(x@data@species.names) > 1) {
              w <- which(grepl(n,colnames(en)))
              if (length(w) == 1) en <- en[,w]
              else {
                w <- which(x@data@species.names == n)
                en <- en[,w]
              }
            } else en <- en[,1]
            
            evaluates(obs,en, distribution = distribution)
            
          }
)

#-------------

if (!isGeneric("getReplication")) {
  setGeneric("getReplication", function(x,id,replication,species,run,index=FALSE,test=TRUE)
    standardGeneric("getReplication"))
}


setMethod('getReplication', signature(x='sdmModels'), 
          function(x,id,replication=NULL,species=NULL,run=NULL,index=FALSE,test=TRUE) {
            
            if (missing(species)) species <- NULL
            if (missing(replication)) replication <- NULL
            if (missing(id)) id <- NULL
            if (missing(run)) run <- NULL
            
            if (missing(index)) index <- FALSE
            if (missing(test)) test <- TRUE
            
            mi <- .getModel.info(x,w=id,species=species,method=NULL,replication=replication,run=run)
            
            if (nrow(mi) > 1) {
              if (all(is.na(mi$replication))) stop('There is no data partitioning (no replications...)!')
              
              if (max(c(length(unique(mi$replicationID)),length(unique(mi$replication)))) > 1) stop('Multiple replications are selected...!')
              
            }
            
            id <- unique(mi$replicationID)
            method <- unique(mi$replication)
            
            .r  <- x@replicates[[as.character(mi$species[1])]]
            
            
            if (length(unique(x@run.info$replication)) > 1) {
              .r <- .r[which(sapply(.r,function(x) x$method) == method)]
            }
            
            if (test) indx <- .r[[id]]$test
            else indx <- .r[[id]]$train
            
            if (index) {
              indx
            } else {
              .df <- as.data.frame(x@data,sp=as.character(mi$species[1]),grp='training')
              .df[.df$rID %in% indx,as.character(mi$species[1])]
            }
          }
)

