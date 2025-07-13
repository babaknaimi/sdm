# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Feb. 2019
# Last update: July 2025
# Version 1.3
# Licence GPL v3
#-----------------------------


.nicheSpace <- function(h,env1,env2,names=NULL) {
  # h: habitat (probability or PA); env1 & 2: two predictors (All numeric vectors!)
  if (length(h) != length(env1) || length(env1) != length(env2)) stop('The input arguments should have the same lengths!')
  
  w1 <-!is.na(h)
  w2 <- !is.na(env1)
  w3 <- !is.na(env2)
  w <- which(w1 & w2 & w3)
  env1 <- env1[w]
  env2 <- env2[w]
  h <- h[w]
  .min1 <- min(env1)
  .max1 <- max(env1)
  .min2 <- min(env2)
  .max2 <- max(env2)
  env1 <- ((env1 - .min1)/(.max1 - .min1)) * 99
  env2 <- ((env2 - .min2)/(.max2 - .min2)) * 99
  env1 <- env1 + 1
  env2 <- env2 + 1
  
  out <- raster(matrix(NA,nrow=100,ncol=100))
  rc <- data.frame(row=101 - env2,col=env1)
  
  cells <- cellFromRowCol(out,round(rc[,1]),round(rc[,2]))
  out[cells] <- h
  names(out) <- 'niche'
  .scp <- data.frame(c(.min1,.max1),c(.min2,.max2))
  colnames(.scp) <- names
  if (inherits(out,'Raster')) new('.nicheRaster',names=names,nicheRaster=out,scaleParams= .scp)
  else new('.nicheSpatRaster',names=names,nicheRaster=out,scaleParams= .scp)
}
#-------------------

.getNicheRaster <- function(x,y,h) {
  out <- raster(matrix(NA,nrow=100,ncol=100))
  rc <- data.frame(row=101 - y,col=x)
  
  cells <- cellFromRowCol(out,rc[,1],rc[,2])
  out[cells] <- h
  names(out) <- 'niche'
  out
}
#---
.getNicheSpatRaster <- function(x,y,h) {
  out <- rast(matrix(NA,nrow=100,ncol=100))
  ext(out) <- c(0,1,0,1)
  rc <- data.frame(row=101 - y,col=x)
  
  cells <- cellFromRowCol(out,rc[,1],rc[,2])
  out[cells] <- h
  names(out) <- 'niche'
  out
}

#-------------
.getEnvSpace <- function(x,limit=1e6) {
  # x is rasterstack/brick
  .n <- names(x)
  if (ncell(x) > (limit*1.3)) {
    df <- data.frame(sampleRandom(x,limit,xy=TRUE))
  } else {
    df <- as.data.frame(x,na.rm=TRUE,xy=TRUE)
  }
  #------
  
  .r <- data.frame(apply(df[,.n],2,range,na.rm=TRUE))
  
  .es <- list()
  for (v in .n) {
    .es[[v]] <- round(((df[,v] - .r[1,v]) / (.r[2,v] - .r[1,v]))*99 + 1)
  }
  new('.envSpace',names=.n,coords=df[,names(df)[!names(df) %in% .n]],scaledVariables=.es,scaleParams=.r)
}
#--------------
.getEnvSpace2 <- function(x,limit=1e6) {
  # x is SpatRaster
  .n <- names(x)
  if (ncell(x) > (limit*1.3)) {
    df <- data.frame(spatSample(x,limit,xy=TRUE,na.rm=TRUE))
  } else {
    df <- as.data.frame(x,na.rm=TRUE,xy=TRUE)
  }
  #------
  
  .r <- data.frame(apply(df[,.n],2,range,na.rm=TRUE))
  
  .es <- list()
  for (v in .n) {
    .es[[v]] <- round(((df[,v] - .r[1,v]) / (.r[2,v] - .r[1,v]))*99 + 1)
  }
  new('.envSpace',names=.n,coords=df[,names(df)[!names(df) %in% .n]],scaledVariables=.es,scaleParams=.r)
}

if (!isGeneric("niche")) {
  setGeneric("niche", function(x,h,n,.size,plot,out,...)
    standardGeneric("niche"))
}

setMethod('niche', signature(x='RasterStackBrick',h='RasterLayer'), 
          function(x,h,n,.size=1e6,plot,out,...) {
            
            if (missing(.size)) .size <- 1e6
            if (missing(plot) || !is.logical(plot)) plot <- TRUE
            if (missing(out) || !is.logical(out)) out <- FALSE
            
            if (nlayers(x) < 2) stop('The number of Raster layers in x is less than 2!')
            
            if (nlayers(x) > 2) {
              if (missing(n)) {
                n <- names(x)[1:2]
                warning('Since n is not specified, niche is generated based on the first two layers in x!')
              } else {
                if (is.character(n)) {
                  if (any(!n %in% names(x))) stop('The names specified in n do not exist in x!')
                  if (length(n) > 2) {
                    n <- n[1:2]
                    warning('The length of items in n should be 2; the first two items is considered!')
                  }
                } else {
                  if (is.numeric(n)) {
                    if (any(n > nlayers(x))) stop('n is not valid!')
                    if (length(n) > 2) {
                      n <- n[1:2]
                      warning('The length of items in n should be 2; the first two items is considered!')
                    }
                  }
                }
              }
              #----
              x <- x[[n]]
            }
            #-----------------
            .es <- .getEnvSpace(x,limit=.size)
            hv <- extract(h,.es@coords)
            nr <- .getNicheRaster(.es@scaledVariables[[1]],.es@scaledVariables[[2]],hv)
            .niche <- new('.nicheRaster',names=.es@names,nicheRaster=nr,scaleParams= .es@scaleParams)
            
            if (plot) {
              print(plot(.niche,...))
              if (out) return(.niche)
            } else  return(.niche)
          }
)
#--------


setMethod('niche', signature(x='RasterStackBrick',h='SpatialPoints'), 
          function(x,h,n,.size=1e6,plot,out,...) {
            .po <- !inherits(h,'SpatialPointsDataFrame') # if h is SpaialPoints -> PO
            
            if (missing(.size)) .size <- 1e6
            
            if (missing(plot) || !is.logical(plot)) plot <- TRUE
            if (missing(out) || !is.logical(out)) out <- FALSE
            
            if (nlayers(x) < 2) stop('The number of Raster layers in x is less than 2!')
            #==============================
            #----- n specifies the name of environmental variables (two first items)
            #----- and its third item the name of column in species data points
            #----- If it is not specified appropriately, it may be guessed, or the warning/error may be generated!
            if (missing(n)) {
              n <- names(x)[1:2]
              if (.po) {
                warning('Since n is not specified, niche is generated based on the first two layers in x!')
              } else {
                # to guess:
                if (ncol(h@data) > 1) {
                  w <- which(sapply(h@data,class) %in% c('numeric','integer'))
                  if (length(w) == 1 && all(range(h@data[,w],na.rm = TRUE) <= 1) && all(range(h@data[,w],na.rm = TRUE) >= 0)) {
                    .nh <- colnames(h@data)[w]
                    warning(paste0('The name of the column containing the species data in the SpatialPointsDataFrame (h) is not specified; So, ',.nh,' is used here!'))
                  } else if (length(w) == 0) {
                    .po <- TRUE
                    warning('Since the name of the species column is not specified in the argument n, it is assumed that the species data specified as the SpatialPointsDataFrame in the h argument is Presence-Only data!')
                  } else stop('The name of the column containing the species data in the SpatialPointsDataFrame (h) should be specified in the argument n; Example:... n=c("predictor1","predictor2","species")')
                } else {
                  if (all(range(h@data[,1],na.rm = TRUE) <= 1) && all(range(h@data[,1],na.rm = TRUE) >= 0)) {
                    .nh <- colnames(h@data)[1]
                    warning(paste0('The name of the column containing the species data in the SpatialPointsDataFrame (h) is not specified in n; So, ',.nh,' is used here!'))
                  } else {
                    .po <- TRUE
                    warning('Since the name of the species column is not specified in the argument n, it is assumed that the species data specified as the SpatialPointsDataFrame in the h argument is Presence-Only data!')
                  }
                }
              }
            } else {
              if (is.character(n)) {
                if (any(!n %in% c(names(x),names(h)))) stop('The names specified in n do not exist in x!')
                if (length(n) > 3) {
                  if (.po) {
                    n <- n[1:2]
                    warning('The first two names in n are used!')
                  } else {
                    .nh <- n[n %in% names(h)]
                    if (length(.nh) == 0) stop('No column name (species data) is specified in the n argument!')
                    else if ((length(.nh) > 1)) stop('Only one name (column name from species data) should be specified in the n argument (check the help)!')
                    
                    n <- n[n %in% names(x)]
                    if (length(n) > 2) {
                      n <- n[1:2]
                      warning(paste0('Only ',n[1],' and ',n[2],' are considered!'))
                    }
                  }
                } else if (length(n) == 3) {
                  if (.po) {
                    n <- n[1:2]
                    warning('The first two names in n are used!')
                  } else {
                    .nh <- n[n %in% names(h)]
                    if (length(.nh) == 0) stop('No column name (species data) is specified in the n argument!')
                    else if ((length(.nh) > 1)) stop('Only one name (column name from species data) should be specified in the n argument (check the help)!')
                    
                    n <- n[n %in% names(x)]
                  }
                } else if (length(n) == 2) {
                  if (!.po) {
                    .nh <- n[n %in% names(h)]
                    if (length(.nh) > 1) stop('Only one name (column name from species data) should be specified in the n argument (check the help)!')
                    else if (length(.nh) == 0) {
                      # to guess the species column (the column with the range between 0-1):
                      if (ncol(h@data) > 1) {
                        w <- which(sapply(h@data,class) %in% c('numeric','integer'))
                        if (length(w) == 1 && all(range(h@data[,w],na.rm = TRUE) <= 1) && all(range(h@data[,w],na.rm = TRUE) >= 0)) {
                          .nh <- colnames(h@data)[w]
                          warning(paste0('The name of the column containing the species data in the SpatialPointsDataFrame (h) is not specified; So, ',.nh,' is used here!'))
                        } else if (length(w) == 0) {
                          .po <- TRUE
                          warning('Since the name of the species column is not specified in the argument n, it is assumed that the species data specified as the SpatialPointsDataFrame in the h argument is Presence-Only data!')
                        } else stop('The name of the column containing the species data in the SpatialPointsDataFrame (h) should be specified in the argument n; Example:... n=c("predictor1","predictor2","species")')
                      } else {
                        if (all(range(h@data[,1],na.rm = TRUE) <= 1) && all(range(h@data[,1],na.rm = TRUE) >= 0)) {
                          .nh <- colnames(h@data)[1]
                          warning(paste0('The name of the column containing the species data in the SpatialPointsDataFrame (h) is not specified in n; So, ',.nh,' is used here!'))
                        } else {
                          .po <- TRUE
                          warning('Since the name of the species column is not specified in the argument n, it is assumed that the species data specified as the SpatialPointsDataFrame in the h argument is Presence-Only data!')
                        }
                      }
                    } else {
                      n <- n[n %in% names(x)]
                      if (length(n) < 2) {
                        n <- c(n,sample(.excludeVector(names(x),n),2-length(n)))
                        if (nlayers(x) > 2) warning(paste0('Predictors ',n[1],', and ',n[2],', are used!'))
                      }
                    }
                  }
                } else {
                  # n == 1 (the name of species may only be specified!)
                  if (!.po) {
                    .nh <- n[n %in% names(h)]
                    if (length(.nh) == 0) {
                      
                      n <- n[n %in% names(x)]
                      if (length(n) < 2) {
                        n <- c(n,sample(.excludeVector(names(x),n),2-length(n)))
                        if (nlayers(x) > 2) warning(paste0('Predictors ',n[1],', and ',n[2],', are used!'))
                      }
                      # to guess the species column (the column with the range between 0-1):
                      if (ncol(h@data) > 1) {
                        w <- which(sapply(h@data,class) %in% c('numeric','integer'))
                        if (length(w) == 1 && all(range(h@data[,w],na.rm = TRUE) <= 1) && all(range(h@data[,w],na.rm = TRUE) >= 0)) {
                          .nh <- colnames(h@data)[w]
                          warning(paste0('The name of the column containing the species data in the SpatialPointsDataFrame (h) is not specified; So, ',.nh,' is used here!'))
                        } else if (length(w) == 0) {
                          .po <- TRUE
                          warning('Since the name of the species column is not specified in the argument n, it is assumed that the species data specified as the SpatialPointsDataFrame in the h argument is Presence-Only data!')
                        } else stop('The name of the column containing the species data in the SpatialPointsDataFrame (h) should be specified in the argument n; Example:... n=c("predictor1","predictor2","species")')
                      } else {
                        if (all(range(h@data[,1],na.rm = TRUE) <= 1) && all(range(h@data[,1],na.rm = TRUE) >= 0)) {
                          .nh <- colnames(h@data)[1]
                          warning(paste0('The name of the column containing the species data in the SpatialPointsDataFrame (h) is not specified in n; So, ',.nh,' is used here!'))
                        } else {
                          .po <- TRUE
                          warning('Since the name of the species column is not specified in the argument n, it is assumed that the species data specified as the SpatialPointsDataFrame in the h argument is Presence-Only data!')
                        }
                      }
                    } else {
                      n <- names(x)[1:2]
                      if (nlayers(x) > 2) warning('Since the name of the two environmental variables are not specified in n, the niche is generated based on the first two layers in x!')
                    }
                  } else {
                    n <- n[n %in% names(x)]
                    if (length(n) < 2) {
                      n <- c(n,sample(.excludeVector(names(x),n),2-length(n)))
                      if (nlayers(x) > 2) warning(paste0('Predictors ',n[1],', and ',n[2],', are used!'))
                    }
                  }
                }
              } else stop('The n argument should be a character vector containing 3 items, the first two, specifying the name of the predictors, and the third specifies the column in species points data (not needed if it is Presence-Only)!')
            }
            #=================
            x <- x[[n]]
            
            if (.po) {
              .es <- .getEnvSpace(x,limit=.size)
              .cp <- cellFromXY(x,h)
              .ca <- cellFromXY(x,.es@coords)
              .ca <- .ca[!.ca %in% .cp]
              hv <- data.frame(x[c(.ca,.cp)])
              .niche <- .nicheSpace(h=c(rep(0,length(.ca)),rep(1,length(.cp))),env1 = hv[,n[1]],env2 = hv[,n[2]],names=n)
            } else {
              hv <- data.frame(extract(x,h))
              .niche <- .nicheSpace(h=h@data[,.nh],env1 = hv[,n[1]],env2 = hv[,n[2]],names=n)
            }
            
            if (plot) {
              print(plot(.niche,...))
              if (out) return(.niche)
            } else  return(.niche)
          }
)
#--------

setMethod('niche', signature(x='sdmdata'), 
          function(x,h,n,.size=1e6,plot,out,...) {
            if (missing(.size)) .size <- 1e6
            if (missing(plot) || !is.logical(plot)) plot <- TRUE
            if (missing(out) || !is.logical(out)) out <- FALSE
            
            if ((length(x@features.name) - length(x@factors)) < 2) stop('The number of predictor variables in x is less than 2!')
            nf <- .excludeVector(x@features.name,x@factors)
            #----
            if (missing(n)) {
              n <- nf[1:2]
              if (length(nf) > 2) warning('Since n is not specified, niche is generated based on the first two predictors in the sdmdata object!')
              .nh <- x@species.names[1]
              if (length(x@species.names) > 1) warning('Since the species name is not specified in n argument, the first species in the sdmdata object is considered!')
            } else {
              if (any(n %in% x@species.names)) {
                .nh <- n[n %in% x@species.names]
                if (length(.nh) > 1) {
                  .nh <- .nh[1]
                  warning('More than one species name is specified in n; the first one is considered!')
                }
              }
              n <- n[n %in% nf]
              if (length(n) > 2) {
                n <- n[1:2]
                warning('More than two predictors is specified in n; the first two are considered!')
              } else if (length(n) == 1) {
                n <- c(n,sample(.excludeVector(nf,n),2-length(n)))
                if (length(nf) > 2) warning(paste0('Predictors ',n[1],', and ',n[2],', are used!'))
              } else if (length(n) == 0) {
                n <- nf[1:2]
                if (length(nf) > 2) warning(paste0('It seems that the predictor names specified in n do not exist; Predictors ',n[1],', and ',n[2],', are used!'))
              }
            }
            #------------
            df <- as.data.frame(x)
            .niche <- .nicheSpace(h=df[,.nh],env1 = df[,n[1]],env2 = df[,n[2]],names=n)
            
            if (plot) {
              print(plot(.niche,...))
              if (out) return(.niche)
            } else  return(.niche)
          }
)
#----------

setMethod('niche', signature(x='RasterStackBrick',h='sdmdata'), 
          function(x,h,n,.size=1e6,plot,out,...) {
            if (is.null(coords(h))) stop('sdmdata object does not have spatial coordinates!')
            else {
              df <- h[1:nrow(h@features),]
              nxy <- colnames(coords(h))
            }
            
            nf <- .excludeVector(h@features.name,h@factors)
            #----
            if (missing(.size)) .size <- 1e6
            if (missing(plot) || !is.logical(plot)) plot <- TRUE
            if (missing(out) || !is.logical(out)) out <- FALSE
            #-----
            if (missing(n)) {
              n <- nf[1:2]
              if (length(nf) > 2) warning('Since n is not specified, niche is generated based on the first two predictors in the sdmdata object!')
              .nh <- h@species.names[1]
              if (length(h@species.names) > 1) warning('Since the species name is not specified in n argument, the first species in the sdmdata object is considered!')
            } else {
              if (any(n %in% h@species.names)) {
                .nh <- n[n %in% h@species.names]
                if (length(.nh) > 1) {
                  .nh <- .nh[1]
                  warning('More than one species name is specified in n; the first one is considered!')
                }
              } else .nh <- h@species.names[1]
              n <- n[n %in% nf]
              if (length(n) > 2) {
                n <- n[1:2]
                warning('More than two predictors is specified in n; the first two are considered!')
              } else if (length(n) == 1) {
                n <- c(n,sample(.excludeVector(nf,n),2-length(n)))
                if (length(nf) > 2) warning(paste0('Predictors ',n[1],', and ',n[2],', are used!'))
              } else if (length(n) == 0) {
                n <- nf[1:2]
                if (length(nf) > 2) warning(paste0('It seems that the predictor names specified in n do not exist; Predictors ',n[1],', and ',n[2],', are used!'))
              }
            }
            #------------
            if (h@species[[.nh]]@type == "Presence-Only") {
              df <- df[,nxy]
              coordinates(df) <- as.formula(paste('~',paste(nxy,collapse = '+')))
            } else {
              df <- df[,c(nxy,.nh)]
              coordinates(df) <- as.formula(paste('~',paste(nxy,collapse = '+')))
            }
            #----
            n <- c(n,.nh)
            niche(x=x,h=df,n=n,.size = .size,plot = plot,out=out,...)
          }
)

#################
##################

setMethod('niche', signature(x='SpatRaster',h='SpatRaster'), 
          function(x,h,n,.size=1e6,plot,out,...) {
            
            if (missing(.size)) .size <- 1e6
            if (missing(plot) || !is.logical(plot)) plot <- TRUE
            if (missing(out) || !is.logical(out)) out <- FALSE
            
            if (nlyr(x) < 2) stop('The number of Raster layers in x is less than 2!')
            
            if (nlyr(x) > 2) {
              if (missing(n)) {
                n <- names(x)[1:2]
                warning('Since n is not specified, niche is generated based on the first two layers in x!')
              } else {
                if (is.character(n)) {
                  if (any(!n %in% names(x))) stop('The names specified in n do not exist in x!')
                  if (length(n) > 2) {
                    n <- n[1:2]
                    warning('The length of items in n should be 2; the first two items is considered!')
                  }
                } else {
                  if (is.numeric(n)) {
                    if (any(n > nlyr(x))) stop('n is not valid!')
                    if (length(n) > 2) {
                      n <- n[1:2]
                      warning('The length of items in n should be 2; the first two items is considered!')
                    }
                  }
                }
              }
              #----
              x <- x[[n]]
            }
            #-----------------
            .es <- .getEnvSpace2(x,limit=.size)
            hv <- extract(h,.es@coords,ID=FALSE)
            nr <- .getNicheSpatRaster(.es@scaledVariables[[1]],.es@scaledVariables[[2]],hv)
            .niche <- new('.nicheSpatRaster',names=.es@names,nicheRaster=nr,scaleParams= .es@scaleParams)
            
            if (plot) {
              print(plot(.niche,...))
              if (out) return(.niche)
            } else  return(.niche)
          }
)
#--------


setMethod('niche', signature(x='SpatRaster',h='SpatVector'), 
          function(x,h,n,.size=1e6,plot,out,...) {
            .po <- dim(h)[2] == 0 # if h is SpaialPoints -> PO
            
            if (missing(.size)) .size <- 1e6
            
            if (missing(plot) || !is.logical(plot)) plot <- TRUE
            if (missing(out) || !is.logical(out)) out <- FALSE
            
            if (nlyr(x) < 2) stop('The number of Raster layers in x is less than 2!')
            #==============================
            #----- n specifies the name of environmental variables (two first items)
            #----- and its third item the name of column in species data points
            #----- If it is not specified appropriately, it may be guessed, or the warning/error may be generated!
            if (missing(n)) {
              n <- names(x)[1:2]
              .nh <- n[n %in% names(h)]
              if (.po) {
                warning('Since n is not specified, niche is generated based on the first two layers in x!')
              } else {
                # to guess:
                if (ncol(h) > 1) {
                  w <- which(sapply(head(h),class) %in% c('numeric','integer'))
                  if (length(w) == 1 && all(range(h[,w],na.rm = TRUE) <= 1) && all(range(h[,w],na.rm = TRUE) >= 0)) {
                    .nh <- names(h)[w]
                    warning(paste0('The name of the column containing the species data in the SpatialPointsDataFrame (h) is not specified; So, ',.nh,' is used here!'))
                  } else if (length(w) == 0) {
                    .po <- TRUE
                    warning('Since the name of the species column is not specified in the argument n, it is assumed that the species data specified as the SpatialPointsDataFrame in the h argument is Presence-Only data!')
                  } else stop('The name of the column containing the species data in the SpatialPointsDataFrame (h) should be specified in the argument n; Example:... n=c("predictor1","predictor2","species")')
                } else {
                  if (all(range(h[,1],na.rm = TRUE) <= 1) && all(range(h[,1],na.rm = TRUE) >= 0)) {
                    .nh <- names(h)[1]
                    warning(paste0('The name of the column containing the species data in the SpatialPointsDataFrame (h) is not specified in n; So, ',.nh,' is used here!'))
                  } else {
                    .po <- TRUE
                    warning('Since the name of the species column is not specified in the argument n, it is assumed that the species data specified as the SpatialPointsDataFrame in the h argument is Presence-Only data!')
                  }
                }
              }
            } else {
              if (is.character(n)) {
                if (any(!n %in% c(names(x),names(h)))) stop('The names specified in n do not exist in x!')
                if (length(n) > 3) {
                  if (.po) {
                    n <- n[1:2]
                    warning('The first two names in n are used!')
                  } else {
                    .nh <- n[n %in% names(h)]
                    if (length(.nh) == 0) stop('No column name (species data) is specified in the n argument!')
                    else if ((length(.nh) > 1)) stop('Only one name (column name from species data) should be specified in the n argument (check the help)!')
                    
                    n <- n[n %in% names(x)]
                    if (length(n) > 2) {
                      n <- n[1:2]
                      warning(paste0('Only ',n[1],' and ',n[2],' are considered!'))
                    }
                  }
                } else if (length(n) == 3) {
                  if (.po) {
                    n <- n[1:2]
                    .nh <- n[n %in% names(h)]
                    warning('The first two names in n are used!')
                  } else {
                    .nh <- n[n %in% names(h)]
                    if (length(.nh) == 0) stop('No column name (species data) is specified in the n argument!')
                    else if ((length(.nh) > 1)) stop('Only one name (column name from species data) should be specified in the n argument (check the help)!')
                    
                    n <- n[n %in% names(x)]
                  }
                } else if (length(n) == 2) {
                  if (!.po) {
                    .nh <- n[n %in% names(h)]
                    if (length(.nh) > 1) stop('Only one name (column name from species data) should be specified in the n argument (check the help)!')
                    else if (length(.nh) == 0) {
                      # to guess the species column (the column with the range between 0-1):
                      if (ncol(h) > 1) {
                        w <- which(sapply(head(h),class) %in% c('numeric','integer'))
                        if (length(w) == 1 && all(range(h[,w],na.rm = TRUE) <= 1) && all(range(h[,w],na.rm = TRUE) >= 0)) {
                          .nh <- names(h)[w]
                          warning(paste0('The name of the column containing the species data in the SpatialPointsDataFrame (h) is not specified; So, ',.nh,' is used here!'))
                        } else if (length(w) == 0) {
                          .po <- TRUE
                          warning('Since the name of the species column is not specified in the argument n, it is assumed that the species data specified as the SpatialPointsDataFrame in the h argument is Presence-Only data!')
                        } else stop('The name of the column containing the species data in the SpatialPointsDataFrame (h) should be specified in the argument n; Example:... n=c("predictor1","predictor2","species")')
                      } else {
                        if (all(range(h[,1],na.rm = TRUE) <= 1) && all(range(h[,1],na.rm = TRUE) >= 0)) {
                          .nh <- names(h)[1]
                          warning(paste0('The name of the column containing the species data in the SpatialPointsDataFrame (h) is not specified in n; So, ',.nh,' is used here!'))
                        } else {
                          .po <- TRUE
                          warning('Since the name of the species column is not specified in the argument n, it is assumed that the species data specified as the SpatialPointsDataFrame in the h argument is Presence-Only data!')
                        }
                      }
                    } else {
                      n <- n[n %in% names(x)]
                      if (length(n) < 2) {
                        n <- c(n,sample(.excludeVector(names(x),n),2-length(n)))
                        if (nlyr(x) > 2) warning(paste0('Predictors ',n[1],', and ',n[2],', are used!'))
                      }
                    }
                  }
                } else {
                  # n == 1 (the name of species may only be specified!)
                  if (!.po) {
                    .nh <- n[n %in% names(h)]
                    if (length(.nh) == 0) {
                      
                      n <- n[n %in% names(x)]
                      if (length(n) < 2) {
                        n <- c(n,sample(.excludeVector(names(x),n),2-length(n)))
                        if (nlyr(x) > 2) warning(paste0('Predictors ',n[1],', and ',n[2],', are used!'))
                      }
                      # to guess the species column (the column with the range between 0-1):
                      if (ncol(h) > 1) {
                        w <- which(sapply(head(h),class) %in% c('numeric','integer'))
                        if (length(w) == 1 && all(range(h[,w],na.rm = TRUE) <= 1) && all(range(h[,w],na.rm = TRUE) >= 0)) {
                          .nh <- names(h)[w]
                          warning(paste0('The name of the column containing the species data in the SpatialPointsDataFrame (h) is not specified; So, ',.nh,' is used here!'))
                        } else if (length(w) == 0) {
                          .po <- TRUE
                          warning('Since the name of the species column is not specified in the argument n, it is assumed that the species data specified as the SpatialPointsDataFrame in the h argument is Presence-Only data!')
                        } else stop('The name of the column containing the species data in the SpatialPointsDataFrame (h) should be specified in the argument n; Example:... n=c("predictor1","predictor2","species")')
                      } else {
                        if (all(range(h[,1],na.rm = TRUE) <= 1) && all(range(h[,1],na.rm = TRUE) >= 0)) {
                          .nh <- names(h)[1]
                          warning(paste0('The name of the column containing the species data in the SpatialPointsDataFrame (h) is not specified in n; So, ',.nh,' is used here!'))
                        } else {
                          .po <- TRUE
                          warning('Since the name of the species column is not specified in the argument n, it is assumed that the species data specified as the SpatialPointsDataFrame in the h argument is Presence-Only data!')
                        }
                      }
                    } else {
                      n <- names(x)[1:2]
                      if (nlyr(x) > 2) warning('Since the name of the two environmental variables are not specified in n, the niche is generated based on the first two layers in x!')
                    }
                  } else {
                    n <- n[n %in% names(x)]
                    if (length(n) < 2) {
                      n <- c(n,sample(.excludeVector(names(x),n),2-length(n)))
                      if (nlyr(x) > 2) warning(paste0('Predictors ',n[1],', and ',n[2],', are used!'))
                    }
                  }
                }
              } else stop('The n argument should be a character vector containing 3 items, the first two, specifying the name of the predictors, and the third specifies the column in species points data (not needed if it is Presence-Only)!')
            }
            #=================
            x <- x[[n]]
            
            if (.po) {
              .es <- .getEnvSpace2(x,limit=.size)
              .cp <- cellFromXY(x,h)
              .ca <- cellFromXY(x,.es@coords)
              .ca <- .ca[!.ca %in% .cp]
              hv <- data.frame(x[c(.ca,.cp)])
              .niche <- .nicheSpace(h=c(rep(0,length(.ca)),rep(1,length(.cp))),env1 = hv[,n[1]],env2 = hv[,n[2]],names=n)
            } else {
              hv <- data.frame(extract(x,h,ID=FALSE))
              .niche <- .nicheSpace(h=data.frame(h[, names(h)[1]])[,1], env1 = hv[,n[1]], env2 = hv[,n[2]], names=n)
            }
            
            if (plot) {
              print(plot(.niche,...))
              if (out) return(.niche)
            } else  return(.niche)
          }
)
#--------

setMethod('niche', signature(x='SpatRaster',h='sdmdata'), 
          function(x,h,n,.size=1e6,plot,out,...) {
            if (is.null(coords(h))) stop('sdmdata object does not have spatial coordinates!')
            else {
              df <- h[1:nrow(h@features),]
              nxy <- colnames(coords(h))
            }
            
            nf <- .excludeVector(h@features.name,h@factors)
            #----
            if (missing(.size)) .size <- 1e6
            if (missing(plot) || !is.logical(plot)) plot <- TRUE
            if (missing(out) || !is.logical(out)) out <- FALSE
            #-----
            if (missing(n)) {
              n <- nf[1:2]
              if (length(nf) > 2) warning('Since n is not specified, niche is generated based on the first two predictors in the sdmdata object!')
              .nh <- h@species.names[1]
              if (length(h@species.names) > 1) warning('Since the species name is not specified in n argument, the first species in the sdmdata object is considered!')
            } else {
              if (any(n %in% h@species.names)) {
                .nh <- n[n %in% h@species.names]
                if (length(.nh) > 1) {
                  .nh <- .nh[1]
                  warning('More than one species name is specified in n; the first one is considered!')
                }
              } else {
                .nh <- h@species.names[1]
              }
              n <- n[n %in% nf]
              if (length(n) > 2) {
                n <- n[1:2]
                warning('More than two predictors is specified in n; the first two are considered!')
              } else if (length(n) == 1) {
                n <- c(n,sample(.excludeVector(nf,n),2-length(n)))
                if (length(nf) > 2) warning(paste0('Predictors ',n[1],', and ',n[2],', are used!'))
              } else if (length(n) == 0) {
                n <- nf[1:2]
                if (length(nf) > 2) warning(paste0('It seems that the predictor names specified in n do not exist; Predictors ',n[1],', and ',n[2],', are used!'))
              }
            }
            #------------
            if (h@species[[.nh]]@type == "Presence-Only") {
              df <- df[,nxy]
              df <- vect(df,nxy)
            } else {
              df <- data.frame(df[,c(nxy,.nh)]) # Wrapped this to a dataframe
              df <- vect(df,nxy)
            }
            #----
            n <- c(n,.nh)
            niche(x=x,h=df,n=n,.size = .size,plot = plot,out=out,...)
          }
)