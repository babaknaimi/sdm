# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Feb. 2024
# Last Update :  Oct. 2024
# Version 1.1
# Licence GPL v3
#-------


.fcor <- function(x,y,...) {
  cor(x,y,method='spearman',...)
}
#----
.nicheSimilarity <- function(x,y) {
  
  # two input rasters will be compared globally
  
  .s <- global(c(x,y),'sum',na.rm=TRUE)
  
  x <- x / .s[1,1]
  y <- y / .s[2,1]
  
  D <- 1- (0.5 * global(abs(x - y),'sum',na.rm=TRUE)[1,1])
  
  sq <- (sqrt(x)-sqrt(y))^2
  
  .s <- global(sq,'sum',na.rm=TRUE)[1,1]
  
  Imod <- 1 - (0.5 * sqrt(.s))
  
  Icor <- 1 - (0.5 * (sqrt(.s)^2))
  
  R <- (global((x+y)*log(x+y),'sum',na.rm=TRUE)[1,1] - 
          global(x*log(x),'sum',na.rm=TRUE)[1,1] - 
          global(y * log(y),'sum',na.rm=TRUE)[1,1]) / (2*log(2))
  
  O <- global(x * y, 'sum', na.rm=TRUE)[1,1] / sqrt(global(x ^ 2, 'sum', na.rm=TRUE)[1,1] * global(y ^ 2, 'sum', na.rm=TRUE)[1,1])
  
  
  BC <- global(min(c(x,y)),'sum',na.rm=TRUE)[1,1] / global(x+y,'sum',na.rm=TRUE)[1,1]
  
  COR <- layerCor(c(x,y),.fcor,use="pairwise.complete.obs")[1,2]
  
  return (c(D=D,Imod=Imod,Icor=Icor,R=R,O=O,BC=BC,COR=COR))
}

#--------------



.nicheSimilarityV <- function(x,y) {
  # two input vectors will be compared globally
  x <- x / sum(x,na.rm=TRUE)
  y <- y / sum(y,na.rm=TRUE)
  
  D <- 1 - (0.5*sum(abs(x-y),na.rm=TRUE))
  
  sq <- (sqrt(x)-sqrt(y)) ^ 2
  
  Imod <- 1 - (0.5 * sqrt(sum(sq,na.rm=TRUE)))
  
  Icor <- 1 - (0.5 * (sqrt(sum(sq,na.rm=TRUE))^2))
  
  lg1 <- (x+y)*log(x+y); lg2 <- x*log(x); lg3 <- y * log(y)
  
  R <- (sum(lg1,na.rm=TRUE) - sum(lg2,na.rm=TRUE) - sum(lg3,na.rm=TRUE)) / (2*log(2))
  o1 <- x*y; o2 <- x^2; o3 <- y^2
  O <- sum(o1,na.rm=TRUE)/sqrt(sum(o2,na.rm=TRUE) * sum(o3,na.rm=TRUE))
  bc1 <- apply(data.frame(x=x,y=y),1,min)
  bc2 <- x + y
  BC <- sum(bc1,na.rm=TRUE) / sum(bc2,na.rm=TRUE)
  COR <- cor(x,y,method='spearman',use="pairwise.complete.obs")
  return (c(D=D,Imod=Imod,Icor=Icor,R=R,O=O,BC=BC,COR=COR))
}



.nicheSimilarityP <- function(x,y,s=3) {
  # two input rasters will be compared partially
  # s specifies the number of breaks in partitioning the probability map
  
  cells <- 1:ncell(x)
  x <- x[][,1]; y <- y[][,1]
  w <- which(!is.na(x))
  x <- x[w]; y <- y[w]
  w <- which(!is.na(y))
  x <- x[w]; y <- y[w]
  
  M <- min(x)
  brk <- (max(x) - M) / s
  lbl <- rep(NA,length(x))
  lbl[which(x >= M & x <= brk)] <- 1
  for(i in 2:s) lbl[which(x > ((i-1)*brk) & x <= (i*brk))] <- i
  out <- data.frame(matrix(ncol=7,nrow=0))
  colnames(out) <- c("D","Imod","Icor","R","O","BC","COR")
  for (i in 1:s) {
    w <- which(lbl == i)
    xx <- x[w]
    yy <- y[w]
    xx <- xx / sum(xx)
    yy <- yy / sum(yy)
    D <- 1 - (0.5*sum(abs(xx-yy)))
    
    sq <- (sqrt(xx)-sqrt(yy))^2
    Imod <- 1 - (0.5 * sqrt(sum(sq,na.rm=TRUE)))
    Icor <- 1 - (0.5 * (sqrt(sum(sq,na.rm=TRUE))^2))
    
    lg1 <- (xx+yy)*log(xx+yy); lg2 <- xx*log(xx); lg3 <- yy * log(yy)
    R <- (sum(lg1,na.rm=TRUE) - sum(lg2,na.rm=TRUE) - sum(lg3,na.rm=TRUE)) / (2*log(2))
    o1 <- xx*yy; o2 <- xx^2; o3 <- yy^2
    O <- sum(o1,na.rm=TRUE)/sqrt(sum(o2,na.rm=TRUE) * sum(o3,na.rm=TRUE))
    bc1 <- apply(data.frame(x=xx,y=yy),1,min)
    bc2 <- xx + yy
    BC <- sum(bc1,na.rm=TRUE) / sum(bc2,na.rm=TRUE)
    COR <- cor(xx,yy,method='spearman',use="pairwise.complete.obs")
    out <- rbind(out,data.frame(D=D,Imod=Imod,Icor=Icor,R=R,O=O,BC=BC,COR=COR))
  }
  
  return (out)
}


#----------
if (!isGeneric("nicheSimilarity")) {
  setGeneric("nicheSimilarity", function(x,y,stat=NULL,w=NULL,...)
    standardGeneric("nicheSimilarity"))
}


setMethod("nicheSimilarity", signature(x='SpatRaster',y='SpatRaster'),
          function(x,y,stat=NULL,w=NULL,...) {
            
            if (nlyr(x) > 1 | nlyr(y) > 1) stop('x and y should have a single raster layer (probability of occurrance)!')
            
            if (ncell(x) != ncell(y)) stop('the number of cells in x and y are not the same!')
            
            if (missing(w)) w <- NULL
            
            if (missing(stat) || is.null(stat)) stat <- NULL
            else {
              if (!is.character(stat)) warning('stat should be character; it is ignored; (stat="all" is used)!')
              else {
                w <- which(!tolower(stat) %in% c('imod','icor','d','r','o','bc','cor','all'))
                if (length(w) > 0) {
                  if (!any(tolower(stat) %in% c('imod','icor','d','r','o','bc','cor','all'))) {
                    warning('stat should be selected from c("Imod","Icor","D","O","BC","R","COR"), or be "all" for all stats (stat="all" is used)! ')
                  } else {
                    if ('all' %in% tolower(stat)) stat <- NULL
                    else {
                      stat <- stat[tolower(stat) %in% c('imod','icor','d','r','o','bc','cor')]
                      stat <- c("Imod","Icor","D","O","BC","R","COR")[tolower(c("Imod","Icor","D","O","BC","R","COR")) %in% tolower(stat)]
                      warning('some of the specified stats are unknown and they are ignored!')
                    }
                  }
                } else {
                  if ('all' %in% tolower(stat)) stat <- NULL
                  else stat <- c("Imod","Icor","D","O","BC","R","COR")[tolower(c("Imod","Icor","D","O","BC","R","COR")) %in% tolower(stat)]
                }
              }
            }
            #################
            
            if (!is.null(w)) {
              if (is.numeric(w)) {
                if (length(w) == 1) {
                  .ns <- .nicheSimilarityP(x,y,s=w)
                  if (!is.null(stat)) .ns <- .ns[,stat]
                } else {
                  if (length(w) < 30) stop('the number of cells specified in w are not enough!')
                  
                  if (all(w %in% 1:ncell(x))) {
                    x <- x[w][,1]
                    y <- y[w][,1]
                    .ns <- .nicheSimilarityV(x,y)
                    if (!is.null(stat)) .ns <- .ns[stat]
                  } else {
                    .w <- length(w)
                    w <- w[w %in% c(1:ncell(x))]
                    if (length(w) > 30 && length(w) > (.w / 3)) {
                      warning('Some of specified cells in "w" are beyound existing cells in x, so they are excluded!')
                      x <- x[w][,1]
                      y <- y[w][,1]
                      .ns <- .nicheSimilarityV(x,y)
                      if (!is.null(stat)) .ns <- .ns[stat]
                    } else {
                      stop('Did you specify cell numbers in "w"?! They are not matched with the cell numbers in the input layers')
                    }
                    
                  }
                }
              } else stop('w is not numeric!')
            } else {
              .ns <- .nicheSimilarity(x, y)
              if (!is.null(stat)) .ns <- .ns[stat]
            }
            
            .ns
            
          }
)



setMethod("nicheSimilarity", signature(x='SpatRaster',y='missing'),
          function(x,y,stat=NULL,w=NULL,...) {
            
            if (nlyr(x) != 2) stop('since y is missing, x should have two raster layers (probability of occurrance)!')
            
            if (missing(w)) w <- NULL
            
            y <- x[[2]]
            
            x <- x[[1]]
            
            if (missing(stat) || is.null(stat)) stat <- NULL
            else {
              if (!is.character(stat)) warning('stat should be character; it is ignored; (stat="all" is used)!')
              else {
                w <- which(!tolower(stat) %in% c('imod','icor','d','r','o','bc','cor','all'))
                if (length(w) > 0) {
                  if (!any(tolower(stat) %in% c('imod','icor','d','r','o','bc','cor','all'))) {
                    warning('stat should be selected from c("Imod","Icor","D","O","BC","R","COR"), or be "all" for all stats (stat="all" is used)! ')
                  } else {
                    if ('all' %in% tolower(stat)) stat <- NULL
                    else {
                      stat <- stat[tolower(stat) %in% c('imod','icor','d','r','o','bc','cor')]
                      stat <- c("Imod","Icor","D","O","BC","R","COR")[tolower(c("Imod","Icor","D","O","BC","R","COR")) %in% tolower(stat)]
                      warning('some of the specified stats are unknown and they are ignored!')
                    }
                  }
                } else {
                  if ('all' %in% tolower(stat)) stat <- NULL
                  else stat <- c("Imod","Icor","D","O","BC","R","COR")[tolower(c("Imod","Icor","D","O","BC","R","COR")) %in% tolower(stat)]
                }
              }
            }
            #################
            if (!is.null(w)) {
              if (is.numeric(w)) {
                if (length(w) == 1) {
                  .ns <- .nicheSimilarityP(x,y,s=w)
                  if (!is.null(stat)) .ns <- .ns[,stat]
                } else {
                  if (length(w) < 30) stop('the number of cells specified in w are not enough!')
                  
                  if (all(w %in% 1:ncell(x))) {
                    x <- x[w][,1]
                    y <- y[w][,1]
                    .ns <- .nicheSimilarityV(x,y)
                    if (!is.null(stat)) .ns <- .ns[stat]
                  } else {
                    .w <- length(w)
                    w <- w[w %in% c(1:ncell(x))]
                    if (length(w) > 30 && length(w) > (.w / 3)) {
                      warning('Some of specified cells in "w" are beyound existing cells in x, so they are excluded!')
                      x <- x[w][,1]
                      y <- y[w][,1]
                      .ns <- .nicheSimilarityV(x,y)
                      if (!is.null(stat)) .ns <- .ns[stat]
                    } else {
                      stop('Did you specify cell numbers in "w"?! They are not matched with the cell numbers in the input layers')
                    }
                    
                  }
                }
              } else stop('w is not numeric!')
            } else {
              .ns <- .nicheSimilarity(x, y)
              if (!is.null(stat)) .ns <- .ns[stat]
            }
            
            
            
            .ns
            
          }
)

#---------



setMethod("nicheSimilarity", signature(x='.nicheSpatRaster',y='.nicheSpatRaster'),
          function(x,y,stat=NULL,...) {
            
            y <- y@nicheRaster
            
            x <- x@nicheRaster
            
            if (missing(stat) || is.null(stat)) stat <- NULL
            else {
              if (!is.character(stat)) warning('stat should be character; it is ignored; (stat="all" is used)!')
              else {
                w <- which(!tolower(stat) %in% c('imod','icor','d','r','o','bc','cor','all'))
                if (length(w) > 0) {
                  if (!any(tolower(stat) %in% c('imod','icor','d','r','o','bc','cor','all'))) {
                    warning('stat should be selected from c("Imod","Icor","D","O","BC","R"), or be "all" for all stats (stat="all" is used)! ')
                  } else {
                    if ('all' %in% tolower(stat)) stat <- NULL
                    else {
                      stat <- stat[tolower(stat) %in% c('imod','icor','d','r','o','bc','cor')]
                      stat <- c("Imod","Icor","D","O","BC","R","COR")[tolower(c("Imod","Icor","D","O","BC","R","COR")) %in% tolower(stat)]
                      warning('some of the specified stats are unknown and they are ignored!')
                    }
                  }
                } else {
                  if ('all' %in% tolower(stat)) stat <- NULL
                  else stat <- c("Imod","Icor","D","O","BC","R","COR")[tolower(c("Imod","Icor","D","O","BC","R","COR")) %in% tolower(stat)]
                }
              }
            }
            #################
            
            .ns <- .nicheSimilarity(x, y)
            
            if (is.null(stat) ) .ns
            else .ns[stat]
            
          }
)


