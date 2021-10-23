# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  Sep. 2021
# Version 1.5
# Licence GPL v3

.getFeature.linear <- function(x) {
  x
}


.getFeature.quad <- function(x) {
  x * x
}

.getFeature.cubic <- function(x) {
  x * x * x
}

.getFeature.poly <- function(x,degree=3,raw=TRUE) {
  d <- as.data.frame(poly(x,degree=degree,raw=raw))
  colnames(d) <- paste('poly',1:degree,sep='')
  d
}

#-------
.hinge <- function(x,th) {
  ifelse(x <= th,0,(x - th) / (max(x,na.rm=TRUE) - th))
}

.invhinge <- function(x,th) {
  ifelse(x >= th,0,1 - ((x - min(x,na.rm=TRUE)) / (th - min(x,na.rm=TRUE))))
}
#---------
.invthresh <- function(x,th) {
  ifelse(x >= th,0,1)
}

.thresh <- function(x,th) {
  ifelse(x <= th,0,1)
}

#-------
.getFeature.hinge <- function(x,knots=20) {
  .min <- min(x,na.rm=TRUE)
  .max <- max(x,na.rm=TRUE)
  k <- seq(.min,.max,length=knots)
  h1 <- as.data.frame(lapply(k[-length(k)],function(th,x,...) {
    .hinge(x,th)
  },x=x))
  colnames(h1) <- paste0('hi_',k[-length(k)])
  #----
  h2 <- as.data.frame(lapply(k[-1],function(th,x,...) {
    .invhinge(x,th)
  },x=x))
  colnames(h2) <- paste0('hd_',k[-1])
  cbind(h1,h2)
}
#-------
.getFeature.threshold <- function(x,knots=20) {
  .min <- min(x,na.rm=TRUE)
  .max <- max(x,na.rm=TRUE)
  k <- seq(.min,.max,length=knots)
  t1 <- as.data.frame(lapply(k[-length(k)],function(th,x,...) {
    .thresh(x,th)
  },x=x))
  colnames(t1) <- paste0('ti_',k[-length(k)])
  #----
  t2 <- as.data.frame(lapply(k[-1],function(th,x,...) {
    .invthresh(x,th)
  },x=x))
  colnames(t2) <- paste0('td_',k[-length(k)])
  cbind(t1,t2)
}


#------
.Jackard <- function(x,y) {
  sum(apply(data.frame(x,y),1,min)) / sum(apply(data.frame(x,y),1,max))
}
#----------
# 
# .getFeature.product <- function(data) {
#   apply(data,1,function(x) {
#     xx <- 1
#     for (i in seq_along(x))  xx <- xx * x[i]
#     xx
#   })
# }
#----------
.getFeature.product <- function(data) {
  x <- data[,1]
  for (i in 2:ncol(data)) {
    x <- x * data[,i]
  }
  x
}



.getSingleFeatureFrame <- function(n,cls,params=NULL,response=NULL) {
  o <- new('.featureFrame',var=n)
  if (cls == '.var') {
    o@feature.name <- n
    o@type <- 'linear'
  } else if (cls == '.quad') {
    o@feature.name <- paste0(n,'.quad')
    o@type <- 'quad'
  } else if (cls == '.cubic') {
    o@feature.name <- paste0(n,'.cubic')
    o@type <- 'cubic'
  } else if (cls == '.factor') {
    o@feature.name <- n
    o@type <- 'factor'
  } else if (cls == '.product') {
    o@feature.name <- paste(c('product_',n),collapse='.')
    o@type <- 'product'
  } else if (cls == '.log') {
    o@feature.name <- paste0(params[[1]],'.',n)
    o@type <- 'log'
    o@params <- params
  } else if (cls == '.poly') {
    d <- params$degree
    if (is.null(d)) d <- 3
    o@feature.name <- paste0(paste0(n,'.poly'),1:d)
    o@type <- 'poly'
    o@params <- params
  } else if (cls == '.hinge') {
    o@feature.name <- paste0(n,'.hinge')
    o@type <- 'hinge'
    o@params <- params
  } else if (cls == '.threshold') {
    o@feature.name <- paste0(n,'.threshold')
    o@type <- 'threshold'
    o@params <- params
  } else if (cls == '.func') {
    o@var <- character()
    o@feature.name <- deparse(params[[1]])
    o@type <- 'func'
    o@params <- params
  } else if (cls == '.simple.func') {
    o@var <- character()
    o@feature.name <- deparse(params[[1]])
    o@params <- params
    o@type <- 'simple.func'
  } else {
    o@feature.name <- paste0(cls,'.',n)
    o@type <- cls
    o@params <- params
  }
  o
}
#-----------

# need to be revised to support when there is multispecies with different indexes (records) for each species
.getFeaturetype <- function(d,f) {
  if (missing(f)) f <- d@sdmFormula
  if (inherits(f,'formula')) f <- .exFormula(f,as.data.frame(d))
  ff <- new('featuresFrame')
  o <- list()
  os <- om <- NULL
  for (m in f@model.terms) {
    fc <- class(m)
    if (fc %in% c('.var','.quad','.cubic','.factor')) o <- c(o,.getSingleFeatureFrame(slot(m,slotNames(m)[1]),fc))
    else if (fc == '.poly') o <- c(o,.getSingleFeatureFrame(m@x,fc,list(degree=m@degree,raw=m@raw)))
    else if (fc == '.product') o <- c(o,.getSingleFeatureFrame(m@x,fc))
    else if (fc == '.log') o <- c(o,.getSingleFeatureFrame(m@x,fc,list(as.character(m@term[[1]]))))
    else if (fc == '.func') {
      l <- as.list(m@x)
      if (as.character(l[[1]]) %in% c(':','*')) {
        o <- c(o,.getSingleFeatureFrame(as.character(.split.formula(m@x,l[[1]])),'.product'))
      } else if (as.character(l[[1]]) %in% c('^') && length(l) == 3) {
        if (l[[3]] == 2) o <- c(o,.getSingleFeatureFrame(as.character(l[[2]]),'.quad'))
        else if (l[[3]] == 3) o <- c(o,.getSingleFeatureFrame(as.character(l[[2]]),'.cubic'))
        else o <- c(o,.getSingleFeatureFrame('xxx',fc,list(m@term)))
      } else {
        o <- c(o,.getSingleFeatureFrame('xxx',fc,list(m@term)))
      }
      
    } else if (fc == '.simple.func') {
      o <- c(o,.getSingleFeatureFrame('xxx',fc,list(m@term)))
    } else if (fc == '.hinge') {
      o <- c(o,.getSingleFeatureFrame(m@x,fc,list(k=m@k)))
      
    } else if (fc == '.threshold') {
      o <- c(o,.getSingleFeatureFrame(m@x,fc,list(k=m@k)))
    } else if (fc == '.selectFrame') {
      # .getFeaturetype(d,f) # recursive
      # need to be checked and tested to select the final set!
    } else if (fc == '.nestedModel') {
      # need to be completed
    }
  }
  ff@vars <- unique(unlist(lapply(o,function(x) x@var)))
  ff@feature.types <- o
  
  if (!is.null(om)) {
    ff@vars <- unique(c(ff@vars,unlist(lapply(om,function(x) x@var))))
    ff@model.specific <- om
  }
  
  ff
}
#--------------



.getModelFrame <- function(x,data,response=NULL,dummy=FALSE) {
  o <- NULL
  if (all(x@vars %in% colnames(data))) {
    d <- data.frame(matrix(nrow=nrow(data),ncol=0))
    
    n <- lapply(x@feature.types,function(x) x@var)
    fn <- lapply(x@feature.types,function(x) x@feature.name)
    ft <- lapply(x@feature.types,function(x) x@type)
    
    un <- unique(n)
    u <- list()
    for (i in un) u[[i]] <- list()
    for (i in 1:length(n)) {
      u[[n[[i]]]][[length(u[[n[[i]]]])+1]] <-list(ft=ft[[i]],fn=fn[[i]],nr=i)
    }
    nrm <- c()
    for (i in 1:length(u)) {
      ft2 <- lapply(u[[i]],function(x) x$ft)
      fn2 <-lapply(u[[i]],function(x) x$fn)
      nr <- sapply(u[[i]],function(x) x$nr)
      
      if ('poly' %in% ft2) {
        w <- which(ft2 == 'poly')
        dg <- x@feature.types[[nr[w]]]@params$degree
        if (dg == 1) {
          if ('linear' %in% ft2) {
            nrm <- c(nrm,nr[w])
          } else {
            fn[[nr[w]]] <- n[[nr[w]]]
            ft[[nr[w]]] <- 'linear'
          }
        } else if (dg > 1) {
          if ('quad' %in% ft2) {
            w <- which(ft2 == 'quad')
            nrm <- c(nrm,nr[w])
          }
          
          if ('linear' %in% ft2) {
            w <- which(ft2 == 'linear')
            nrm <- c(nrm,nr[w])
          }
          if (dg > 2) {
            if ('cubic' %in% ft2) {
              w <- which(ft2 == 'cubic')
              nrm <- c(nrm,nr[w])
            }
          }
        }
      }
    }
    if (length(nrm) > 0) {
      ft <- ft[-nrm]
      fn <- fn[-nrm]
      n <- n[-nrm]
      x@feature.types <- x@feature.types[-nrm]
    }
    #---
    
    for (i in 1:length(n)) {
      if (ft[[i]] == 'linear') d[[fn[[i]]]] <- data[,n[[i]]]
      else if (ft[[i]] == 'factor') d[[fn[[i]]]] <- data[,n[[i]]]
      else if (ft[[i]] == 'quad') d[[fn[[i]]]] <- .getFeature.quad(data[,n[[i]]])
      else if (ft[[i]] == 'cubic') d[[fn[[i]]]] <- .getFeature.cubic(data[,n[[i]]])
      else if (ft[[i]] == 'poly') {
        temp <- do.call(".getFeature.poly",c(list(x=data[,n[[i]]]),x@feature.types[[i]]@params))
        colnames(temp) <- x@feature.types[[i]]@feature.name
        d <- cbind(d,temp)
      }
      else if (ft[[i]] == 'product') d[[fn[[i]]]] <- .getFeature.product(data[,x@feature.types[[i]]@var])
      else if (ft[[i]] == 'log') d[[fn[[i]]]] <- do.call(x@feature.types[[i]]@params[[1]],list(x=data[,n[[i]]]))
      else if (ft[[i]] %in% c('func','simple.func')) {
        temp <- model.frame(as.formula(paste('~',deparse(x@feature.types[[i]]@params[[1]]))),data=data)
        n <- colnames(temp)
        d[[n]] <- as.vector(temp[,1])
      } else if (ft[[i]] == 'threshold') d[[fn[[i]]]] <- .getFeature.threshold(data[,n[[i]]],knots = x@feature.types[[i]]@params$k)
      else if (ft[[i]] == 'hinge') d[[fn[[i]]]] <- .getFeature.hinge(data[,n[[i]]],knots = x@feature.types[[i]]@params$k)
    }
    
  } else stop('some of the specified variables in the formula do not exist in the data!')
  list(features=d,specis_specific=o) 
}
#-----------
.getFeatureNamesTypes <- function(x,merged=TRUE) {
  # get the name of features from a featureFrame onject (x)
  # If merged is FALSE, separately reported in a list
  n1 <- n2 <- data.frame(matrix(ncol=2,nrow = 0))
  colnames(n1) <- colnames(n2) <- c('name','type')
  for (i in seq_along(x@feature.types)) {
    n1 <-  rbind(n1,data.frame(name=x@feature.types[[i]]@feature.name,type=x@feature.types[[i]]@type))
  }
  
  for (i in seq_along(x@response.specific)) {
    n2 <-  rbind(n2,data.frame(name=x@response.specific[[i]]@feature.name,type=x@response.specific[[i]]@type))
  }
  if (nrow(n1) == 0) n1 <- NULL
  if (nrow(n2) == 0) n2 <- NULL
  
  if (merged) {
    n1 <- rbind(n1,n2)
    if (!is.null(n1)) {
      n1$name <- as.character(n1$name)
      n1$type <- as.character(n1$type)
    }
  } else {
    if (!is.null(n1)) {
      n1$name <- as.character(n1$name)
      n1$type <- as.character(n1$type)
    }
    if (!is.null(n2)) {
      n2$name <- as.character(n2$name)
      n2$type <- as.character(n2$type)
    }
    n1 <- list(features=n1,response.specific=n2)
  }
  n1
}
