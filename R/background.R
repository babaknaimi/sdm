# Author: Babak Naimi, naimi.b@gmail.com
# Date of last update :  May 2024
# Version 1.8
# Licence GPL v3
#--------------

# Functions for different methods of generating background:


.pseudo_gRandom.Raster <- function(preds,n=500,bias=NULL) {
  if (missing(bias)) bias <- NULL
  if (missing(n)) n <- 500
  
  if (is.null(bias)) {
    s <- sampleRandom(preds,n,cells=TRUE,xy=TRUE)
  } else {
    w <- Which(bias > 0,cells=TRUE)
    if (length(w) > 0) {
      if (length(w) < n) {
        warning(paste('the size of background (n) is less than the available pixels specified in the bias file, so n is changed to:'),length(w))
        n <- length(w)
      }
      s <- sample(w,n,prob = bias[w])
      .xy <- xyFromCell(bias,s)
      .c <- cellFromXY(preds,.xy)
      s <- data.frame(cell=.c,.xy,preds[.c])
    }
  }
  
  # if (!is.null(p) && ncol(p) == 2) {
  #   p.cells <- cellFromXY(preds,p)
  #   if (length(p.cells) > 0) {
  #     s <- s[which(!s[,'cell'] %in% p.cells),]
  #   }
  # }
  s[,-1]
}
#----
.pseudo_gRandom.Terra <- function(preds,n=500,bias=NULL) {
  if (missing(bias)) bias <- NULL
  if (missing(n)) n <- 500
  
  if (is.null(bias)) {
    s <- spatSample(preds,n,cells=TRUE,xy=TRUE,na.rm=TRUE)
  } else {
    bias <- ifel(bias > 0, bias,NA)
    
    w <- cells(bias)
    
    if (length(w) > 0) {
      if (length(w) < n) {
        warning(paste('the size of background (n) is less than the available pixels specified in the bias file, so n is changed to:'),length(w))
        n <- length(w)
      }
      s <- sample(w,n,prob = bias[w][,1])
      .xy <- xyFromCell(bias,s)
      .c <- cellFromXY(preds,.xy)
      s <- data.frame(cell=.c,.xy,preds[.c])
    }
  }
  
  s[,-1]
}
#----------
.pseudo_eRandom.Raster <- function(preds,n=1000,nclass=5,method='kmeans',factors=NULL) {
  
  if (missing(nclass) || is.null(nclass)) nclass <- 5
  
  if (missing(n)) n <- 500
  if (missing(factors)) factors <- NULL
  if (missing(method)) method <- 'kmeans'
  
  .n <- names(preds)
  
  if (!is.null(factors)) {
    if (is.numeric(factors)) {
      .nl <- c(1:nlayers(preds))
      w <- which(.nl %in% factors)
      if (length(w) > 0) {
        .nl <- .nl[-w]
        if (length(.nl) == 0) stop('all predictor layers are factor, but the method of eRandom for background generation works with continuous data; you may use gRandom instead!')
        else preds <- preds[[.nl]]
      }
      .n <- names(preds)
    } else if (is.character(factors)) {
      w <- which(.n %in% factors)
      if (length(w) > 0) {
        .n <- .n[-w]
        if (length(.n) == 0) stop('all predictor layers are factor, but the method of eRandom for background generation works with continuous data; you may use gRandom instead!')
        #else preds <- preds[[.n]]
      }
    } else warning('factors in bg setting should be either numeric or character; it is ignored!')
  }
  #---------
  
  .d <- as.data.frame(preds,na.rm=TRUE,xy=TRUE)
  
  if (nrow(.d) < n) {
    warning(paste('the size of background (n) is less than the available pixels, so n is changed to:'),nrow(.d))
    n <- nrow(.d)
  }
  
  k <- kmeans(.d[,.n,drop=FALSE],nclass, iter.max = 100)
  .nc <- floor(n / nclass)
  r <- raster(preds[[1]])
  .c <- cellFromXY(r,.d[,1:2])
  
  r[.c] <- k$cluster
  s <- sampleStratified(r,.nc,xy=TRUE,na.rm=TRUE)
  s <- data.frame(s[,2:3],preds[s[,1]])
  s
}
#----------
.pseudo_eRandom.Terra <- function(preds,n=1000,nclass=5,method='kmeans',factors=NULL) {
  
  if (missing(nclass) || is.null(nclass)) nclass <- 5
  
  if (is.null(method) || method %in% c('km','kmean')) method <- 'kmeans'
  
  if (method %in% c('binning','bin','categorizing','concat')) {
    if (.require('elsa')) method <- 'binning'
    else {
      method <- 'kmeans'
      warnings ('the selected algorithm for the eRandom method of background generation requires the package of elsa installed which is not available, therefore, the algorithm is set to "kmeans"...!')
    }
  }
  
  if (!method %in% c('binning','kmeans')) {
    method <- 'kmeans'
    warning('"kmeans" is set for the algorithm to generate background records using the eRandom method')
  }
  #--------
  .n <- names(preds)
  if (!is.null(factors)) {
    if (is.numeric(factors)) {
      .nl <- c(1:nlyr(preds))
      w <- which(.nl %in% factors)
      if (length(w) > 0) {
        .nl <- .nl[-w]
        if (length(.nl) == 0) stop('all predictor layers are factor, but the method of eRandom for background generation works with continuous data; you may use gRandom instead!')
        else preds <- preds[[.nl]]
      }
      .n <- names(preds)
    } else if (is.character(factors)) {
      w <- which(.n %in% factors)
      if (length(w) > 0) {
        .n <- .n[-w]
        if (length(.n) == 0) stop('all predictor layers are factor, but the "eRandom" background generation method works with continuous data; you may use gRandom instead!')
        #else preds <- preds[[.n]]
      }
    } else warning('factors in bg setting should be either numeric or character; it is ignored!')
  }
  #---------
  .n <- names(preds)
  if (method == 'kmeans') {
    if (.canProcessInMemory(preds,2)) {
      .d <- as.data.frame(preds,na.rm=TRUE,cells=TRUE)
      if (nrow(.d) < n) {
        warning(paste('the size of background (n) is less than the available pixels, so n is changed to:'),nrow(.d))
        n <- nrow(.d)
      }
    } else {
      .d <- spatSample(preds,min(c(20000,blocks(preds,n=2)[[2]][1]*ncol(preds) * 0.7)),na.rm=TRUE,cell=TRUE)
    }
    
    k <- kmeans(.d[,.n,drop=FALSE],nclass, iter.max = 100)
    .nc <- floor(n / nclass)
    .tb <- table(k$cluster)
    
    s <- c()
    if (any(.tb < .nc)) {
      .ntb <-as.numeric(names(.tb[.tb < .nc]))
      for (.cl in .ntb) {
        s <- c(s,sample(which(k$cluster == .cl),.nc,replace = TRUE))
      }
      #----
      if (any(.tb >= .nc)) {
        .ntb <-as.numeric(names(.tb[.tb >= .nc]))
        for (.cl in .ntb) {
          s <- c(s,sample(which(k$cluster == .cl),.nc,replace = FALSE))
        }
      }
      
    } else {
      .ntb <-as.numeric(names(.tb))
      for (.cl in .ntb) {
        s <- c(s,sample(which(k$cluster == .cl),.nc))
      }
    }
    #---
    
    
    
    #r <- rast(preds[[1]])
    #r[.d$cell] <- k$cluster
    #predict(preds,k)
    
    #s <- spatSample(r,.nc,"stratified",xy=TRUE,na.rm=TRUE)
    #.c <- cellFromXY(r,s[,1:2])
    xy <- xyFromCell(preds,.d$cell[s])
    s <- data.frame(xy,preds[.d$cell[s]])
    s
  } else {
    
    #--- The method of binning should be modified because of the method of concats
    # currently, this method only works based on the first two predictor layers!!!
    
    #if (nlyr(preds) > 2) warning('In the current version, the binning approach for eRandom background generation can work based on two predictors, so only the first 2 layers are considered!')
    
    if (nlyr(preds) > 2) preds <- pca(preds,TRUE)@data[[1:2]]
    
    prc <- .eval('categorize(preds,nc = nclass)',env = environment())
    
    prc <- as.factor(prc)
    x <- concats(prc[[1]],prc[[2]]) # it is only for two layers 
    xx <- as.numeric(x)
    .n <- freq(xx)
    .nx <- .n$value[.n$count == 1]
    #if (length(.nx) > 0) xx <- ifel(xx %in% .nx,xx@cpp@.xData$range_max,xx)
    if (length(.nx) > 0) xx <- ifel(xx %in% .nx,slot(xx,slotNames(xx)[1])$range_max,xx)
    
    
    s <- spatSample(xx,size=ceiling(n / nrow(.n)),method='stratified',na.rm=TRUE,cells=TRUE,replace=TRUE)
    .xy <- xyFromCell(preds,s$cell)
    data.frame(.xy,preds[s$cell])
    
  }
  
  
}
#----------
.pseudo_gDist.Terra <- function(preds,n=1000,p=NULL) {
  if (!is.null(p) && ncol(p) == 2) {
    p.cells <- cellFromXY(preds,p)
    if (length(p.cells) > 0) {
      r <- preds[[1]]
      r <- ifel(!is.na(r),1,0)
      r[p.cells] <- 0
      gd <- gridDist(r,target=0)
      gd <- ifel(gd <= 0, NA,gd)
      w <- cells(gd)
      
      if (length(w) > 0) {
        if (length(w) < n) {
          warning(paste('the size of background (n) is less than the available pixels, so n is changed to:'),length(w))
          n <- length(w)
        }
        s <- sample(w,n,prob = gd[w][,1])
        .xy <- xyFromCell(r,s)
        s <- data.frame(.xy,preds[s])
      }
    } else stop('dDist method is failed to generate background.... !')
  } else stop('dDist requires presence locations to generate backgroud....!')
  s
}
#----------
.pseudo_gDist.Raster <- function(preds,n=1000,p=NULL) {
  if (!is.null(p) && ncol(p) == 2) {
    p.cells <- cellFromXY(preds,p)
    if (length(p.cells) > 0) {
      r <- preds[[1]]
      r[!is.na(r)] <- 1
      
      r[p.cells] <- 0
      gd <- gridDistance(r,0)
      gd[is.na(r)] <- NA
      #w <- cellStats(gd,function(x,...) sum(x > 0,...))
      
      w <- Which(gd > 0, cells=TRUE)
      
      if (length(w) > 0) {
        if (length(w) < n) {
          warning(paste('the size of background (n) is less than the available pixels, so n is changed to:'),length(w))
          n <- length(w)
        }
        s <- sample(w,n,prob = gd[w])
        .xy <- xyFromCell(r,s)
        s <- data.frame(.xy,preds[s])
      }
    } else stop('gDist method is failed to generate background.... !')
  } else stop('dDist requires presence locations to generate backgroud....!')
  s
}
#----------
.pseudo_eDist.Terra <- function(preds,n=1000,p=NULL,method='mahal',power=1,factors=NULL) {
  
  if (missing(factors)) factors <- NULL
  
  if (!is.null(factors)) {
    if (is.numeric(factors)) {
      .nl <- c(1:nlyr(preds))
      w <- which(.nl %in% factors)
      if (length(w) > 0) {
        .nl <- .nl[-w]
        if (length(.nl) == 0) stop('all predictor layers are factor, but the method of eDist for background generation works with continuous data; you may use gRandom or gDist, instead!')
        else preds <- preds[[.nl]]
      }
    } else if (is.character(factors)) {
      .n <- names(preds)
      w <- which(.n %in% factors)
      if (length(w) > 0) {
        .n <- .n[-w]
        if (length(.n) == 0) stop('all predictor layers are factor, but the method of eDist for background generation works with continuous data; you may use gRandom or gDist, instead!')
        else preds <- preds[[.n]]
      }
    } else warning('factors in the bg setting should be either numeric or character; it is ignored!')
  }
  #---------
  if (missing(method) || is.null(method)) method <- 'mahal'
  if (missing(power) || is.null(power)) power <- 1
  
  
  #--------
  if (method != 'mahal') {
    # this is temporary until other methods of eDist (e.g., ENFA, DOMAIN, etc.) are implemented!
    method <- 'mahal'
    warning('The eDist method for background generation is changed to "mahal" (mahalanobis)...!')
  }
  #--------
  if (method == 'mahal') {
    .d <- as.data.frame(preds,na.rm=TRUE,cells=TRUE)
    s <- solve(cov(.d[,names(preds)]))
    .ex <- as.matrix(extract(preds,p))[,names(preds)]
    
    di <- mahalanobis(.d[,names(preds)],colMeans(.ex,na.rm = TRUE), cov=s,inverted = TRUE)
    #-----
    if (length(di) < n) {
      warning(paste('the size of background (n) is less than the available pixels, so n is changed to:'),length(di))
      n <- length(di)
    }
    
    s <- .d$cell[sample(length(di),n,prob = di^power)]
    .xy <- xyFromCell(preds,s)
    s <- data.frame(.xy,preds[s])
    s
  }
}
#----------
.pseudo_eDist.Raster <- function(preds,n=1000,p=NULL,method='mahal',power=1,factors=NULL) {
  
  if (missing(factors)) factors <- NULL
  
  if (!is.null(factors)) {
    if (is.numeric(factors)) {
      .nl <- c(1:nlayers(preds))
      w <- which(.nl %in% factors)
      if (length(w) > 0) {
        .nl <- .nl[-w]
        if (length(.nl) == 0) stop('all predictor layers are factor, but the method of eDist for background generation works with continuous data; you may use gRandom or gDist, instead!')
        else preds <- preds[[.nl]]
      }
    } else if (is.character(factors)) {
      .n <- names(preds)
      w <- which(.n %in% factors)
      if (length(w) > 0) {
        .n <- .n[-w]
        if (length(.n) == 0) stop('all predictor layers are factor, but the method of eDist for background generation works with continuous data; you may use gRandom or gDist, instead!')
        else preds <- preds[[.n]]
      }
    } else warning('factors in the bg setting should be either numeric or character; it is ignored!')
  }
  #---------
  if (missing(method) || is.null(method)) method <- 'mahal'
  if (missing(power) || is.null(power)) power <- 1
  
  
  
  if (method != 'mahal') {
    # this is temporary until other methods of eDist (e.g., ENFA, DOMAIN, etc.) are implemented!
    method <- 'mahal'
    warning('The eDist method for background generation is changed to "mahal" (mahalanobis)...!')
  }
  
  if (method == 'mahal') {
    .d <- as.data.frame(preds,na.rm=TRUE,xy=TRUE)
    .xy <- .d[,c(1:2)]
    .d <- data.frame(cell=cellFromXY(preds,.xy),.d[,-c(1:2)])
    s <- solve(cov(.d[,names(preds)]))
    .ex <- extract(preds,p)[,names(preds)]
    di <- mahalanobis(.d[,names(preds)],colMeans(.ex,na.rm = TRUE), cov=s,inverted = TRUE)
    #-----
    if (length(di) < n) {
      warning(paste('the size of background (n) is less than the available pixels, so n is changed to:'),length(di))
      n <- length(di)
    }
    
    s <- .d$cell[sample(length(di),n,prob = di^power)]
    .xy <- xyFromCell(preds,s)
    s <- data.frame(.xy,preds[s])
    s
  }
}
#----------
.pseudo.Raster <- function(preds,bg,p=NULL,factors=NULL) {
  
  if (missing(p)) p <- NULL
  
  if (missing(factors)) factors <- NULL
  
  if (missing(bg)) {
    bg <- list(method='gRandom',n=500,remove=FALSE)
  } else {
    names(bg) <- tolower(names(bg))
    nbg <- names(bg)
    nbg <- .pmatch(nbg,c('n','method','bias','remove','algorithm','power','nclass'))
    names(bg) <- nbg
    .bg <- list()
    if ('n' %in% nbg) .bg[['n']] <- bg[['n']]
    else .bg[['n']] <- 500
    
    if ('method' %in% nbg) {
      if (bg[['method']] %in% c('gRandom','random','rnd','grnd','georandom','gR','gRand','gr','geo')) {
        .bg[['method']] <- 'gRandom'
      } else if (bg[['method']] %in% c('eRandom','envrandom','ernd','eR','eRand','er')) {
        .bg[['method']] <- 'eRandom'
      } else if (bg[['method']] %in% c('gDistance','gd','gD','gDis','gdis','gDist','gdist','geoDist','geoD','geod')) {
        .bg[['method']] <- 'gDist'
      } else if (bg[['method']] %in% c('eDistance','environ','envDist','ed','eD','EDist','eDist','envD','edis','eDis')) {
        .bg[['method']] <- 'eDist'
      }
    } else .bg[['method']] <- 'gRandom'
    
    # if ('remove' %in% nbg && is.logical(bg[['remove']])) .bg[['remove']] <- bg[['remove']]
    # else .bg[['remove']] <- FALSE
    
    if ('bias' %in% nbg && inherits(bg[['bias']],'RasterLayer')) .bg[['bias']] <- bg[['bias']]
    
    if ('algorithm' %in% nbg) .bg[['algorithm']] <- bg[['algorithm']]
    
    if ('power' %in% nbg) .bg[['power']] <- bg[['power']]
    
    if ('nclass' %in% nbg) .bg[['nclass']] <- bg[['nclass']]
  }
  ########################################
  
  method <- .bg[['method']]
  
  
  if (method == 'gRandom') .pseudo_gRandom.Raster(preds,n=.bg[['n']],bias=.bg[['bias']])
  else if (method == 'gDist') .pseudo_gDist.Raster(preds,n=.bg[['n']],p=p)
  else if (method == 'eRandom') .pseudo_eRandom.Raster(preds,n=.bg[['n']],method= .bg[['algorithm']],nclass=.bg[['nclass']],factors=factors)
  else if (method == 'eDist') .pseudo_eDist.Raster(preds,n=.bg[['n']],p=p,method=.bg[['algorithm']],power=.bg[['power']],factors=factors)
  
}

#------

.pseudo.terra <- function(preds,bg,p=NULL,factors=NULL) {
  
  if (missing(p)) p <- NULL
  
  if (missing(factors)) factors <- NULL
  
  if (missing(bg)) {
    bg <- list(method='gRandom',n=500,remove=FALSE)
  } else {
    names(bg) <- tolower(names(bg))
    nbg <- names(bg)
    nbg <- .pmatch(nbg,c('n','method','bias','remove','algorithm','power','nclass'))
    names(bg) <- nbg
    .bg <- list()
    if ('n' %in% nbg) .bg[['n']] <- bg[['n']]
    else .bg[['n']] <- 500
    
    if ('method' %in% nbg) {
      if (bg[['method']] %in% c('gRandom','random','rnd','grnd','georandom','gR','gRand','gr','geo')) {
        .bg[['method']] <- 'gRandom'
      } else if (bg[['method']] %in% c('eRandom','envrandom','ernd','eR','eRand','er')) {
        .bg[['method']] <- 'eRandom'
      } else if (bg[['method']] %in% c('gDistance','gd','gD','gDis','gdis','gDist','gdist','geoDist','geoD','geod')) {
        .bg[['method']] <- 'gDist'
      } else if (bg[['method']] %in% c('eDistance','environ','envDist','ed','eD','EDist','eDist','envD','edis','eDis')) {
        .bg[['method']] <- 'eDist'
      }
    } else .bg[['method']] <- 'gRandom'
    
    # if ('remove' %in% nbg && is.logical(bg[['remove']])) .bg[['remove']] <- bg[['remove']]
    # else .bg[['remove']] <- FALSE
    
    if ('bias' %in% nbg && inherits(bg[['bias']],'SpatRaster')) .bg[['bias']] <- bg[['bias']]
    
    if ('algorithm' %in% nbg) .bg[['algorithm']] <- bg[['algorithm']]
    
    if ('power' %in% nbg) .bg[['power']] <- bg[['power']]
    
    if ('nclass' %in% nbg) .bg[['nclass']] <- bg[['nclass']]
  }
  ########################################
  
  method <- .bg[['method']]
  
  if (method == 'gRandom') .pseudo_gRandom.Terra(preds,n=.bg[['n']],bias=.bg[['bias']])
  else if (method == 'gDist') .pseudo_gDist.Terra(preds,n=.bg[['n']],p=p)
  else if (method == 'eRandom') .pseudo_eRandom.Terra(preds,n=.bg[['n']],method= .bg[['algorithm']],nclass=.bg[['nclass']],factors=factors)
  else if (method == 'eDist') .pseudo_eDist.Terra(preds,n=.bg[['n']],p=p,method=.bg[['algorithm']],power=.bg[['power']],factors=factors)
  
}
#----------

#--------
if (!isGeneric("background")) {
  setGeneric("background", function(x,n,method,bias,sp,setting)
    standardGeneric("background"))
}	




setMethod('background', signature(x='SpatRaster'), 
          function(x,n,method,bias,sp,setting) {
            if (missing(n)) stop('n (size) is not specified...!')
            if (missing(bias)) bias <- NULL
            else {
              if (!inherits(bias,'SpatRaster')) stop('the (optional) bias argument should be a SpatRaster!')
              
              if (nlyr(bias) > 1) stop('the (optional) bias argument should be a SpatRaster with a single layer (nlyr(x) is NOT 1)!')
            }
            if (missing(sp)) sp <- NULL
            
            bg <- list()
            
            if (!missing(setting)) {
              nbg <- names(setting)
              nbg <- .pmatch(nbg,c('algorithm','power','nclass','factors'))
              names(setting) <- nbg
              #----
              
              if ('algorithm' %in% nbg) bg[['algorithm']] <- setting[['algorithm']]
              
              if ('power' %in% nbg) bg[['power']] <- setting[['power']]
              
              if ('nclass' %in% nbg) bg[['nclass']] <- setting[['nclass']]
              
              if ('factors' %in% nbg) bg[['factors']] <- setting[['factors']]
            }
            
            
            
            if (method %in% c("gRandom","gr","grand","gR","random","rnd","geo","georandom","geoRandom")) {
              .pseudo_gRandom.Terra(x,n=n,bias=bias)
            } else if (method %in% c('eRandom','eR','er','erand','eRand','envR')) {
              .pseudo_eRandom.Terra(x,n=n,method=bg[['algorithm']],nclass=bg[['nclass']],factors=bg[['factors']])
            } else if (method %in% c('eDist','ed','eD','eDis','eDistance','envD','edistance')) {
              if (is.null(sp)) stop('to generate background based on eDist method, species presence locations are required (sp)!')
              
              if (inherits(sp,'SpatVector')) sp <- geom(sp)[,3:4]
              else if (inherits(sp,'data.frameORmatrix')) {
                if (ncol(sp) != 2) stop('Species presence locations should be provided either as Spatial points or a 2-columns data.frame/matrix')
                sp <- data.frame(sp)
              }
              
              
              .pseudo_eDist.Terra(x,n=n,p=sp,method=bg[['algorithm']],power=bg[['power']],factors=bg[['factors']])
            } else if (method %in% c('gDist','gd','gD','gDis','gDistance','geoD','gdistance')) {
              if (is.null(sp)) stop('to generate background based on gDist method, species presence locations are required (sp)!')
              
              if (inherits(sp,'SpatVector')) sp <- geom(sp)[,3:4]
              else if (inherits(sp,'data.frameORmatrix')) {
                if (ncol(sp) != 2) stop('Species presence locations should be provided either as Spatial points or a 2-columns data.frame/matrix')
                sp <- data.frame(sp)
              }
              
              .pseudo_gDist.Terra(x,n=n,p=sp)
            }
          }
)


setMethod('background', signature(x='Raster'), 
          function(x,n,method,bias,sp,setting) {
            if (missing(n)) stop('n (size) is not specified...!')
            if (missing(bias)) bias <- NULL
            else {
              if (!inherits(bias,'RasterLayer')) stop('the (optional) bias argument should be a RasterLayer!')
              
              if (nlayers(bias) > 1) stop('the (optional) bias argument should be a RasterLayer with a single layer (nlayer(x) is NOT 1)!')
            }
            if (missing(sp)) sp <- NULL
            
            bg <- list()
            
            if (!missing(setting)) {
              nbg <- names(setting)
              nbg <- .pmatch(nbg,c('algorithm','power','nclass','factors'))
              names(setting) <- nbg
              #----
              
              if ('algorithm' %in% nbg) bg[['algorithm']] <- setting[['algorithm']]
              
              if ('power' %in% nbg) bg[['power']] <- setting[['power']]
              
              if ('nclass' %in% nbg) bg[['nclass']] <- setting[['nclass']]
              
              if ('factors' %in% nbg) bg[['factors']] <- setting[['factors']]
            }
            
            
            
            if (method %in% c("gRandom","gr","grand","gR","random","rnd","geo","georandom","geoRandom")) {
              .pseudo_gRandom.Raster(x,n=n,bias=bias)
            } else if (method %in% c('eRandom','eR','er','erand','eRand','envR')) {
              .pseudo_eRandom.Raster(x,n=n,method=bg[['algorithm']],nclass=bg[['nclass']],factors=bg[['factors']])
            } else if (method %in% c('eDist','ed','eD','eDis','eDistance','envD','edistance')) {
              if (is.null(sp)) stop('to generate background based on eDist method, species presence locations are required (sp)!')
              
              if (inherits(sp,'SpatVector')) sp <- geom(sp)[,3:4]
              else if (inherits(sp,'data.frameORmatrix')) {
                if (ncol(sp) != 2) stop('Species presence locations should be provided either as Spatial points or a 2-columns data.frame/matrix')
                sp <- data.frame(sp)
              }
              
              
              .pseudo_eDist.Raster(x,n=n,p=sp,method=bg[['algorithm']],power=bg[['power']],factors=bg[['factors']])
            } else if (method %in% c('gDist','gd','gD','gDis','gDistance','geoD','gdistance')) {
              if (is.null(sp)) stop('to generate background based on gDist method, species presence locations are required (sp)!')
              
              if (inherits(sp,'SpatialPoints')) sp <- coordinates(sp)
              else if (inherits(sp,'data.frameORmatrix')) {
                if (ncol(sp) != 2) stop('Species presence locations should be provided either as Spatial points or a 2-columns data.frame/matrix')
                sp <- data.frame(sp)
              }
              
              .pseudo_gDist.Raster(x,n=n,p=sp)
            }
          }
)

