# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Feb. 2024
# Last Update :  Feb. 2024
# Version 1.0
# Licence GPL v3
#-------




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
  
  return (c(D=D,Imod=Imod,Icor=Icor,R=R,O=O,BC=BC))
}

#--------------


#----------
if (!isGeneric("nicheSimilarity")) {
  setGeneric("nicheSimilarity", function(x,y,stat=NULL,...)
    standardGeneric("nicheSimilarity"))
}


setMethod("nicheSimilarity", signature(x='SpatRaster',y='SpatRaster'),
          function(x,y,stat=NULL,...) {
            
            if (nlyr(x) > 1 | nlyr(y) > 1) stop('x and y should have a single raster layer (probability of occurrance)!')
            
            
            if (missing(stat) || is.null(stat)) stat <- NULL
            else {
              if (!is.character(stat)) warning('stat should be character; it is ignored; (stat="all" is used)!')
              else {
                w <- which(!tolower(stat) %in% c('imod','icor','d','r','o','bc','all'))
                if (length(w) > 0) {
                  if (!any(tolower(stat) %in% c('imod','icor','d','r','o','bc','all'))) {
                    warning('stat should be selected from c("Imod","Icor","D","O","BC","R"), or be "all" for all stats (stat="all" is used)! ')
                  } else {
                    if ('all' %in% tolower(stat)) stat <- NULL
                    else {
                      stat <- stat[tolower(stat) %in% c('imod','icor','d','r','o','bc')]
                      stat <- c("Imod","Icor","D","O","BC","R")[tolower(c("Imod","Icor","D","O","BC","R")) %in% tolower(stat)]
                      warning('some of the specified stats are unknown and they are ignored!')
                    }
                  }
                } else {
                  if ('all' %in% tolower(stat)) stat <- NULL
                  else stat <- c("Imod","Icor","D","O","BC","R")[tolower(c("Imod","Icor","D","O","BC","R")) %in% tolower(stat)]
                }
              }
            }
            #################
            
            .ns <- .nicheSimilarity(x, y)
            
            if (is.null(stat) ) .ns
            else .ns[stat]
            
          }
)



setMethod("nicheSimilarity", signature(x='SpatRaster',y='missing'),
          function(x,y,stat=NULL,...) {
            
            if (nlyr(x) != 2) stop('since y is missing, x should have two raster layers (probability of occurrance)!')
            
            y <- x[[2]]
            
            x <- x[[1]]
            
            if (missing(stat) || is.null(stat)) stat <- NULL
            else {
              if (!is.character(stat)) warning('stat should be character; it is ignored; (stat="all" is used)!')
              else {
                w <- which(!tolower(stat) %in% c('imod','icor','d','r','o','bc','all'))
                if (length(w) > 0) {
                  if (!any(tolower(stat) %in% c('imod','icor','d','r','o','bc','all'))) {
                    warning('stat should be selected from c("Imod","Icor","D","O","BC","R"), or be "all" for all stats (stat="all" is used)! ')
                  } else {
                    if ('all' %in% tolower(stat)) stat <- NULL
                    else {
                      stat <- stat[tolower(stat) %in% c('imod','icor','d','r','o','bc')]
                      stat <- c("Imod","Icor","D","O","BC","R")[tolower(c("Imod","Icor","D","O","BC","R")) %in% tolower(stat)]
                      warning('some of the specified stats are unknown and they are ignored!')
                    }
                  }
                } else {
                  if ('all' %in% tolower(stat)) stat <- NULL
                  else stat <- c("Imod","Icor","D","O","BC","R")[tolower(c("Imod","Icor","D","O","BC","R")) %in% tolower(stat)]
                }
              }
            }
            #################
            
            .ns <- .nicheSimilarity(x, y)
            
            if (is.null(stat) ) .ns
            else .ns[stat]
            
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
                w <- which(!tolower(stat) %in% c('imod','icor','d','r','o','bc','all'))
                if (length(w) > 0) {
                  if (!any(tolower(stat) %in% c('imod','icor','d','r','o','bc','all'))) {
                    warning('stat should be selected from c("Imod","Icor","D","O","BC","R"), or be "all" for all stats (stat="all" is used)! ')
                  } else {
                    if ('all' %in% tolower(stat)) stat <- NULL
                    else {
                      stat <- stat[tolower(stat) %in% c('imod','icor','d','r','o','bc')]
                      stat <- c("Imod","Icor","D","O","BC","R")[tolower(c("Imod","Icor","D","O","BC","R")) %in% tolower(stat)]
                      warning('some of the specified stats are unknown and they are ignored!')
                    }
                  }
                } else {
                  if ('all' %in% tolower(stat)) stat <- NULL
                  else stat <- c("Imod","Icor","D","O","BC","R")[tolower(c("Imod","Icor","D","O","BC","R")) %in% tolower(stat)]
                }
              }
            }
            #################
            
            .ns <- .nicheSimilarity(x, y)
            
            if (is.null(stat) ) .ns
            else .ns[stat]
            
          }
)


