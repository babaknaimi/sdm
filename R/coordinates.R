# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March 2016
# Last Update :  Jan. 2024
# Version 1.4
# Licence GPL v3

#-------------
if (!isGeneric('coords')) {
  setGeneric('coords', function(obj,...)
    standardGeneric('coords'))
}

setMethod('coords', signature(obj='sdmdata'),
          function(obj,...) {
            if (!is.null(obj@info) && !is.null(obj@info@coords)) obj@info@coords[,2:3]
          }
)

setMethod('coords', signature(obj='sdmModels'),
          function(obj,...) {
            if (!is.null(obj@data@info) && !is.null(obj@data@info@coords)) obj@data@info@coords[,2:3]
          }
)


if (!isGeneric("coords<-")) {
  setGeneric("coords<-", function(object,value)
    standardGeneric("coords<-"))
}


setReplaceMethod('coords', signature(object='sdmdata'),
                 function(object,value) {
                   if (inherits(value,'matrix')) {
                     if (nrow(value) != nrow(object@features)) stop('number of rows in the provided coordinates is not the same as the number of records in the sdmdata object!')
                     if (is.null(object@info)) object@info <- new('.info')
                     
                     if (ncol(value) == 2) {
                       colnames(value) <- tolower(colnames(value))
                       value <- data.frame(rID=object@features$rID,value)
                       if (all(c('x','y') %in% colnames(value))) {
                         object@info@coords <- as.matrix(value[,c('rID','x','y')])
                       } else if (all(c('lon','lat') %in% colnames(value))) {
                         object@info@coords <- as.matrix(value[,c('rID','lon','lat')])
                       } else object@info@coords <- as.matrix(value)
                     } else if (ncol(value) == 3) {
                       if (!'rID' %in% colnames(value)) stop('The coordinates should have two columns (x and y)')
                       if (nrow(value) != nrow(object@features)) stop('number of rows in the provided coordinates is not the same as the number of records in the sdmdata object!')
                       value <- data.frame(value)
                       if (!all(value$rID %in% object@features$rID)) stop('The record IDs (rID) in the coordinate matrix is not the same as the rIDs in the sdmdata object')
                       
                       if (all(c('x','y') %in% colnames(value))) {
                         object@info@coords <- as.matrix(value[,c('rID','x','y')])
                       } else if (all(c('lon','lat') %in% colnames(value))) {
                         object@info@coords <- as.matrix(value[,c('rID','lon','lat')])
                       } else object@info@coords <- as.matrix(value)
                       
                     } else stop('The coordinates should have two columns (x and y)')
                   } else if (inherits(value,'data.frame')) {
                     if (nrow(value) != nrow(object@features)) stop('number of rows in the provided coordinates is not the same as the number of records in the sdmdata object!')
                     if (is.null(object@info)) object@info <- new('.info')
                     
                     if (ncol(value) == 2) {
                       colnames(value) <- tolower(colnames(value))
                       value <- data.frame(rID=object@features$rID,value)
                       if (all(c('x','y') %in% colnames(value))) {
                         object@info@coords <- as.matrix(value[,c('rID','x','y')])
                       } else if (all(c('lon','lat') %in% colnames(value))) {
                         object@info@coords <- as.matrix(value[,c('rID','lon','lat')])
                       } else object@info@coords <- as.matrix(value)
                     } else if (ncol(value) == 3) {
                       if (!'rID' %in% colnames(value)) stop('The coordinates should have two columns (x and y)')
                       if (nrow(value) != nrow(object@features)) stop('number of rows in the provided coordinates is not the same as the number of records in the sdmdata object!')
                       if (!all(value$rID %in% object@features$rID)) stop('The record IDs (rID) in the coordinate data.frame is not the same as the rIDs in the sdmdata object')
                       
                       if (all(c('x','y') %in% colnames(value))) {
                         object@info@coords <- as.matrix(value[,c('rID','x','y')])
                       } else if (all(c('lon','lat') %in% colnames(value))) {
                         object@info@coords <- as.matrix(value[,c('rID','lon','lat')])
                       } else object@info@coords <- as.matrix(value)
                       
                     } else stop('The coordinates should have two columns (x and y)')
                   } else if (inherits(value,'character')) {
                     if (length(value) != 2) stop('value should provide either the coordinate matrix or to change the names of coordinate columns in sdmdata (character with the length of 2)!')
                     if (!is.null(object@info) && !is.null(object@info@coords)) {
                       colnames(object@info@coords)[2:3] <- value
                       .eval('cat("\nNames of coordinates in the sdmdata is changed...")',environment())
                     } else stop('There is no coordinates in the sdmdata!')
                   }
                   
                   object
                 }
)