# Author: Babak Naimi, naimi.b@gmail.com
# Date of last update :  Oct. 2023
# Version 3.8
# Licence GPL v3

#------
# create a .group class # .newgroup('training','test',list(test=30:40))
.newgroup <- function(name,values,index) {
  g <- new('.group',name=name)
  if (!is.null(values)) {
    if ((length(values) != length(index)) || !is.list(index)) stop('group values and index list are not match')
    g@values <- data.frame(indexID=1:length(values),values=values)
    g@indices <- index
  } else if (!is.null(index) && is.list(index) && !is.null(names(index))) {
    g@values <- data.frame(indexID=1:length(index),values=names(index))
    g@indices <- index
  }
  g
}
#-------
.getGroupNames <- function(d,levels=FALSE) {
  if (levels) {
    if (!is.null(names(d@groups))) {
      nn <- c()
      for (n in names(d@groups)) nn <- c(nn,as.character(d@groups[[n]]@values[,2]))
      nn
    } else NULL
  } else names(d@groups)
}
#-----
# .newgroup creates a new group class, and .newGroup call it and add it to the data object!
.newGroup <- function(d,name,values=NULL,index=NULL) {
  d@groups[[name]] <- .newgroup(name,values,index)
  d
}
#-----
.getGroupLevels <- function(d,g=NULL) {
  if (!is.null(g)) {
    g <- .pmatch(g,.getGroupNames(d))
    g <- g[!is.na(g)][1]
    if (is.na(g)) stop(paste('group',g,'does not exist!'))
    as.character(d@groups[[g]]@values[,2])
  }
}
#-----

###################

# tbl <- c(122,'test','train') (example1 : tbl is vector)
# 
# tbl <- data.frame(matrix(NA,ncol=3,nrow=3)) (example 2: tbl is data.frame)
# tbl[,1] <- c(4,5,6)
# tbl[,2] <- rep('test',3)
# tbl[,3] <- rep('train',3)
# x: sdmdata
# the group for the specified ids will be changed from the first group to the second
.updateGroup <- function(x,tbl,sp=NULL) {
  if (is.vector(tbl)) {
    id <- as.numeric(tbl[1])
    g1 <- tbl[2]
    g2 <- tbl[3]
  } else {
    id <- as.numeric(tbl[,1])
    g1 <- tbl[,2]
    g2 <- tbl[,3]
  }
  gr <- .getGroupNames(x)
  for (g in gr) {
    for (i in seq_along(id)) {
      if (g1[i] %in% names(x@groups[[g]]@indices)) {
        if (id[i] %in% .getGroupIndex(x,g1[i])) {
          x@groups[[g]]@indices[[g1[i]]] <- x@groups[[g]]@indices[[g1[i]]][-which(x@groups[[g]]@indices[[g1[i]]] == id[i])]
          x@groups[[g]]@indices[[g2[i]]] <- c(x@groups[[g]]@indices[[g2[i]]] , id[i])
        }
      }
    }
  }
  x
}
#-----------
#######################
.getSpeciesNames <- function(d,n=NULL) {
  if (is.null(n)) d@species.names
  else d@species.names[d@species.names %in% n]
}

#-----
# convert the sdmData to a data.frame (if list is TRUE, for multi-species, a list is returned)
.getSpeciesDF <- function(d,sp,id,list=TRUE) {
  if (missing(sp)) sp <- d@species.names
  if (missing(id)) id <- d@features$rID
  o <- list()
  for (s in sp) {
    o[[s]] <- data.frame(rID=id,value=NA)
    if (d@species[[s]]@type == 'Presence-Absence') {
      o[[s]][id %in% d@species[[s]]@presence,2] <- 1
      o[[s]][id %in% d@species[[s]]@absence,2] <- 0
    } else if (d@species[[s]]@type == 'Presence-Only') {
      o[[s]][id %in% d@species[[s]]@presence,2] <- 1
    } else if (d@species[[s]]@type == 'Presence-Background') {
      o[[s]][id %in% d@species[[s]]@presence,2] <- 1
      o[[s]][id %in% d@species[[s]]@background,2] <- 0
    } else if (d@species[[s]]@type == 'Abundance') {
      o[[s]][id %in% d@species[[s]]@abundance$rID,2] <- d@species[[s]]@abundance[d@species[[s]]@abundance$rID %in% id,2]
    } else if (d@species[[s]]@type == 'Absence-Only!') {
      o[[s]][id %in% d@species[[s]]@absence,2] <- 0
    } else if (d@species[[s]]@type == 'Abundance-constant!') {
      o[[s]][id %in% d@species[[s]]@abundance$rID,2] <- d@species[[s]]@abundance[1,2]
    }
  }
  
  for (i in 1:length(o)) {
    o[[i]] <- o[[i]][!is.na(o[[i]]$value),]
  }

  if (list) {
    o
  } else {
    .o <- o[[1]]
    if ('species' %in% colnames(.o)) .o[['species____name']] <- names(o)[1]
    else .o[['species']] <- names(o)[1]
    
    if (length(o) > 1) {
      for (i in 2:length(o)) {
        .o2 <- o[[i]]
        if ('species' %in% colnames(.o2)) .o[['species____name']] <- names(o)[i]
        else .o2[['species']] <- names(o)[i]
        
        .o <- rbind(.o,.o2)
      }
    }
    .o
  }
  
}

#-----
.getSpeciesIndex <- function(d,sp=NULL) {
  id <- c()
  sp <- .getSpeciesNames(d,sp)
  for (i in seq_along(sp)) {
    id <- c(id,d@species[[sp[i]]]@presence,d@species[[sp[i]]]@absence,d@species[[sp[i]]]@abundance[,1],d@species[[sp[i]]]@background,d@species[[sp[i]]]@Multinomial[,1])
  }
  sort(unique(id))
}
#-----
.getGroupIndex=function(d,g=NULL) {
  if (!is.null(g)) {
    id <- c()
    for (gg in g) {
      gl <- unlist(lapply(unlist(strsplit(gg,':')),.trim))
      if (length(gl) > 1) {
        if (!gl[1] %in% .getGroupNames(d)) stop(paste('group',gl[1],'does not exist!'))
        if (!gl[2] %in% .getGroupNames(d,TRUE)) stop(paste('group level',gl[2],'does not exist!'))
        id <- c(id,d@groups[[gl[1]]]@indices[[gl[2]]])
      } else {
        if (!gl %in% c(.getGroupNames(d),.getGroupNames(d,TRUE))) stop(paste(gl,' is neither a group nor a group level!'))
        if (gl %in% .getGroupNames(d)) {
          for (gv in d@groups[[gl]]@values[,2]) id <- c(id,d@groups[[gl]]@indices[[gv]])
        } else {
          for (gn in .getGroupNames(d)) {
            if (gl %in% d@groups[[gn]]@values[,2]) id <- c(id,d@groups[[gn]]@indices[[gl]])
          }
        }
      }
    }
    unique(id)
  }
}
#-----
.getTimeIndex <- function(d,.t=NULL) {
  if (!is.null(d@info) && !is.null(d@info@time)) {
    id <- NULL
    if (is.null(.t)) {
      id <- d@info@time
    } else {
      if (is.numeric(.t)) {
        w <- which(d@info@time[,1] %in% .t)
        if (length(w) > 0) id <- d@info@time[w,]
      } else {
        w <- which(d@info@time[,2] %in% .t)
        if (length(w) > 0) id <- d@info@time[w,]
      }
    }
    
    
    id
  }
}
#-----
.getIndex <- function(d,sp=NULL,groups=NULL,time=NULL) {
  id1 <- id2 <- id3 <- NULL
  if (is.null(sp)) id1 <- d@features$rID
  else {
    id1 <- .getSpeciesIndex(d,sp)
  }
  
  if (!is.null(groups)) {
    id2 <- .getGroupIndex(d,groups)
  }
  
  if (!is.null(time)) {
    id3 <- .getTimeIndex(d,time)
  }
  
  if (is.null(id2)) id <- id1
  else id <- id1[id1 %in% id2]
  
  if (is.null(id3)) id
  else id[id %in% id3]
}

#-----------

if (!isGeneric('.addLog<-')) {
  setGeneric('.addLog<-', function(d,value)
    standardGeneric('.addLog<-'))
}
#---
setReplaceMethod('.addLog','sdmdata', 
                 function(d,value) {
                   d@errorLog <- c(d@errorLog,value)
                   d
                 }
)
#----------

# detect the type of species data (i.e., pa, po, ab); when all are 0 returns ao, and when variance is 0, returns ab_constant 
.speciesType <- function(x) {
  u <- unique(x)
  if (is.numeric(x)) {
    if (length(u) > 2) {
      if (all(c(x - round(x)) == 0)) return('Abundance')
      else return('Numeric')
    } else if (length(u) == 2) {
      if ((all(sort(u) == c(0,1)) || all(sort(u) == c(-1,1)))) return('Presence-Absence')
      else {
        if (all(c(x - round(x)) == 0)) return('Abundance')
        else return('Numeric')
      }
    } else {
      if (u == 1) return('Presence-Only')
      else if (u == 0 || u == -1) return('Absence-Only!') # ao is absence only! ONLY to detect and handle this kind of records
      else return('Constant (single numeric value)!') # ab_constant to detect the species data when the variance is 0
    } 
  } else if (is.logical(x)) {
    if (length(u) == 2) return('Presence-Absence')
    else {
      if (u) return('Presence-Only')
      else return('Absence-Only!')
    }
  } else {
    x <- as.character(x)
    u <- unique(x)
    if (length(u) > 2) return('Presence-Only')
    else if (length(u) == 2) {
      if (all(sort(u) == c('0','1')) || all(sort(u) == c('-1','1'))) return('Presence-Absence')
      else return('Presence-Only')
    } else return('Presence-Only')
  }
}

#----------


# get species data.frame from the input data:
.getSpecies <- function(data,nsp,bg=FALSE,id.train=NULL,id.test=NULL) {
  species <- list()
  for (n in nsp) {
    if (!is.null(id.test)) {
      typ1 <- .speciesType(data[id.train,n])
      typ2 <- .speciesType(data[id.test,n])
      if (typ1 == typ2) typ <- typ1
      else stop('train and test data seem different (for example, one is presence-only and the other is presence-absence...)!')
    } else typ <- .speciesType(data[,n])
    
    if (typ == 'Presence-Absence') {
      species[[n]] <- new('.species.data')
      species[[n]]@name <- n
      w <- as.numeric(data[,n])
      species[[n]]@presence <- data$rID[which(w == 1)]
      if (bg) {
        species[[n]]@background <- data$rID[which(w %in% c(0,-1))]
        species[[n]]@type <- 'Presence-Background'
      } else {
        species[[n]]@absence <- data$rID[which(w %in% c(0,-1))]
        species[[n]]@type <- typ
      }
    } else if (typ == 'Presence-Only') {
      if (is.numeric(data[,n]) || is.logical(data[,n])) {
        species[[n]] <- new('.species.data')
        species[[n]]@name <- n
        species[[n]]@presence <- data$rID
        species[[n]]@type <- typ
      } else {
        w <- as.character(data[,n])
        u <- unique(w)
        if ('0' %in% u) {
          u <- u[u != '0']
          bg <- data$rID[which(w == '0')]
        } else bg <- NULL
        for (uu in u) {
          species[[uu]] <- new('.species.data')
          species[[uu]]@name <- uu
          species[[uu]]@presence <- data$rID[which(w == uu)]
          species[[uu]]@type <- typ
        }
        if (!is.null(bg)) {
          for (uu in u) {
            species[[uu]]@background <- bg
            species[[uu]]@type <- 'Presence-Background'
          }
        }
      }
    } else if (typ == 'Abundance') {
      species[[n]] <- new('.species.data')
      species[[n]]@name <- n
      species[[n]]@abundance <- data.frame(rID=data$rID,abundance=data[,n])
      species[[n]]@type <- typ
    } else if (typ == 'Numeric') {
      species[[n]] <- new('.species.data')
      species[[n]]@name <- n
      species[[n]]@numerical <- data.frame(rID=data$rID,value=data[,n])
      species[[n]]@type <- typ
    } else if (typ == 'Abundance_constant!') {
      species[[n]] <- new('.species.data')
      species[[n]]@name <- n
      species[[n]]@abundance <- data.frame(rID=data$rID,abundance=data[,n])
      species[[n]]@type <- typ
      warning(paste('for species',n,', the variance in abundance data is ZERO!'))
    } else if (typ == 'Multinomial') {
      species[[n]] <- new('.species.data')
      species[[n]]@name <- n
      species[[n]]@multinomial <- data.frame(rID=data$rID,name=data[,n])
      species[[n]]@type <- typ
    } else if (typ == 'Absence-Only!') {
      species[[n]] <- new('.species.data')
      species[[n]]@name <- n
      if (bg) {
        species[[n]]@background <- data$rID
        species[[n]]@type <- 'Background'
      } else {
        species[[n]]@absence <- data$rID
        species[[n]]@type <- typ
      }
      warning(paste('for species',n,', there is no presence record (all values are 0 or absence)!'))
    }
  }
  species
}

#----------
# remove duplicate records, and the rows that species columns contain NA OR all (or any) columns contain NA
.dataClean <- function(x,nsp,ANY=TRUE) {
  # ANY: whether a record should be removed if only one/some predictor variable is/are NA, or all should be NA to be removed?
  rm.na <-0; rm.duplicate <- 0
  w <- nrow(x)
  x <- unique(x)
  if (nrow(x) < w) rm.duplicate <- w - nrow(x)
  
  w <- nrow(x)
  if (!missing(nsp)) {
    if (length(nsp) > 1) ww <- which(apply(x[,which(colnames(x) %in% nsp)],1,function(x){any(is.na(x))}))
    else ww <- which(is.na(x[,which(colnames(x) %in% nsp)]))
    if (length(ww) > 0) {
      x <- x[-ww,]
      rm.na <- w - nrow(x)
    }
  } else {
    if (ANY) ww <- which(apply(x,1,function(x){any(is.na(x))}))
    else ww <- which(apply(x,1,function(x){all(is.na(x))}))
    if (length(ww) > 0) {
      x <- x[-ww,]
      rm.na <- w - nrow(x)
    }
  }
  list(x,c(na=rm.na,duplicate=rm.duplicate))
}
#----------


.speciesDetect <- function(data) {
  # to detect species columns with presence/absence data
  # also detect factors, and lon/lat coordinate columns
  nsp <-  nxy <- nFact <- nf <- nt <- NULL
  w <- which(unlist(lapply(data,.isBinomial)))
  if (length(w) == 0) {
    stop ('No species variable is detected, spcify the species variable in the formula or use an appropriate data structure...')
  } else {
    varNames <- colnames(data)
    nsp <- varNames[w]
    if (length(varNames) > length(nsp)) {
      nf <- varNames[-w]
      w <- which(unlist(lapply(data[,nf],function(x) class(x) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr"))))
      if (length(w) > 0) {
        nt <- nf[w]
        nf <- .excludeVector(nf,nt)
      }
      w <- which(unlist(lapply(data[,nf],function(x) class(x) %in% c('character','factor'))))
      if (length(w) > 0) {
        nFact <- nf[w]
        nf <- .excludeVector(nf,nFact)
      }
    }
    w <- .which.is.coords(nf)
    
    if (length(w == 2)) {
      nxy <- w
      nf <- .excludeVector(nf,nxy)
    }
  }
  list(nsp=nsp,nf=nf,nFact=nFact,nxy=nxy,nt=nt)
}
#--------

.selectData <- function(d) {
  
  nf <- d@sdmFormula@vars@numeric$names
  
  # the Formula in select term should be in the form of either all variables - those to keep (e.g., ~. - v1 - v2)
  # or list of variables to be checked for selection (e.g., ~ v1 + v2 + v3 ; OR: ~ .  [all numeric variables])
  if (!is.null(d@sdmFormula@data.terms)) {
    .cls <- sapply(d@sdmFormula@data.terms,class)
    #function(x) inherits(x,'.selectFrame')
    if ('.selectFrame' %in% .cls) {
      .slf <- d@sdmFormula@data.terms[[which(.cls == '.selectFrame')]]
      if ('.' %in% .slf@vars || length(.slf@vars) == 0 || is.null(.slf@vars)) .svar <- nf
      else .svar <- .slf@vars
      #----
      .svar <- .svar[.svar %in% nf]
      if (all(nf %in% .svar)) {
        if (!is.null(.slf@keep) || length(.slf@keep) != 0) {
          .slf@keep <- .slf@keep[.slf@keep %in% nf]
          if (length(.slf@keep) == 0) {
            .slf@keep <- NULL
            warning('The variables specified in keep in the select argument in Formula are ignored because either they do not exist or they are not numeric variables!')
            .addLog(d) <- 'The variables specified in keep in the select argument in Formula are ignored because either they do not exist or they are not numeric variables!'
          }
        }
      } else {
        .slf@keep <- nf[!nf %in% .svar]
        .svar <- nf # among all nf variables check collinearity and keep those specified in keep!
      }
      #----
      if (length(.svar) < 2) {
        .slf@method <- ""
        warning('select is ignored: The number of variables specified for select in Formula should be at least 2 variables!')
        .addLog(d) <- 'select is ignored: The number of variables specified for select in Formula should be at least 2 variables!'
      } else {
        if (!tolower(.slf@method) %in% c('vifstep','vifcor','vif')) {
          .w <- .pmatch(tolower(.slf@method),c('vifstep','vifcor','vif'))
          if (is.na(.w)) {
            .addLog(d) <- 'The method in the select argument in formula is not identified (should be either vifstep or vifcor) so it is ignored!'
            warning('The method in the select argument in formula is not identified (should be either vifstep or vifcor) so it is ignored!')
          } else {
            if (.w == 'vif') {
              .w <- 'vifstep'
              .addLog(d) <- '"vifstep" is assumed for the method in the select argument in formula!'
            }
            d@sdmFormula@data.terms[[which(.cls == '.selectFrame')]]@method <- .slf@method <- .w
          }
        } else {
          if (tolower(.slf@method) == 'vif') {
            d@sdmFormula@data.terms[[which(.cls == '.selectFrame')]]@method <- .slf@method <- 'vifstep'
            .addLog(d) <- '"vifstep" is assumed for the method in the select argument in formula!'
          } else d@sdmFormula@data.terms[[which(.cls == '.selectFrame')]]@method <- .slf@method <- tolower(.slf@method)
        }
      }
      #------
      
      
      if (.slf@method %in% c('vifstep','vifcor')) {
        if (is.null(.slf@keep) || length(.slf@keep) == 0) {
          if (.require('usdm')) {
            if (.slf@method == 'vifstep') {
              if (is.null(.slf@th) || length(.slf@th) == 0) {
                d@sdmFormula@data.terms[[which(.cls == '.selectFrame')]]@th <- .slf@th <- 10
              } else {
                if (.slf@th < 1) {
                  warning('Are you sure the right threshold is specified for the vifstep method in the select argument in the Formula?! (for vifstep, a commonly used threshold is between 5 to 10)')
                  .addLog(d) <-'Are you sure the right threshold is specified for the vifstep method in the select argument in the Formula?! (for vifstep, a commonly used threshold is between 5 to 10)'
                }
              }
              #----
              
              # vifstep
              .v <- .eval("vifstep(d@features[,.svar],th=.slf@th,size=nrow(d@features)+1)",env=environment())
              # exclude based on vifstep:
              if (!inherits(.v,'try-error')) {
                if (length(.v@excluded) > 0) {
                  if (length(.v@excluded) == 1) .addLog(d) <- paste0('One variable had collinearity problem and so, it is excluded. Its name is: ',.v@excluded)
                  else .addLog(d) <- paste0(length(.v@excluded),' variables had collinearity problems and they are excluded. They include: ',paste(.v@excluded,collapse = ', '))
                  #----
                  .v <- .eval('usdm::exclude(d@features,.v)',env=environment())
                  if (!inherits(.v,'try-error')) {
                    d@features <- .v
                  } else {
                    .addLog(d) <- .v
                  }
                  
                } else {
                  .addLog(d) <- 'No input variables from those specified in the select argument in Formula had collinearity problem so none of them are excluded! '
                }
              } else {
                .addLog(d) <- .v
              }
              #-----
            } else {
              if (is.null(.slf@th) || length(.slf@th) == 0) {
                d@sdmFormula@data.terms[[which(.cls == '.selectFrame')]]@th <- .slf@th <- 0.9
              } else {
                if (.slf@th >= 1) {
                  warning('A wrong threshold is specified for the vifcor method in the select argument in the Formula. It is changed to 0.9!')
                  .addLog(d) <-'A wrong threshold is specified for the vifcor method in the select argument in the Formula. It is changed to 0.9!'
                }
              }
              #----
              # vifcor:
              .v <- .eval("vifcor(d@features[,.svar],th=.slf@th,size=nrow(d@features)+1)",env=environment())
              # exclude based on vifcor:
              if (!inherits(.v,'try-error')) {
                if (length(.v@excluded) > 0) {
                  if (length(.v@excluded) == 1) .addLog(d) <- paste0('One variable had collinearity problem and so, it is excluded. Its name is: ',.v@excluded)
                  else .addLog(d) <- paste0(length(.v@excluded),' variables had collinearity problems and they are excluded. They include: ',paste(.v@excluded,collapse = ', '))
                  #----
                  .v <- .eval('usdm::exclude(d@features,.v)',env=environment())
                  if (!inherits(.v,'try-error')) {
                    d@features <- .v
                  } else {
                    .addLog(d) <- .v
                  }
                  
                } else {
                  .addLog(d) <- 'No input variables from those specified in the select argument in Formula had collinearity problem so none of them are excluded! '
                }
              } else {
                .addLog(d) <- .v
              }
            }
            
          } else {
            warning('usdm package is required for multicollinearity check but it is not installed!')
          }
        } else { # when keep is specified in select(...)
          #-------
          
          if (.require('usdm')) {
            if (.slf@method == 'vifstep') {
              if (is.null(.slf@th) || length(.slf@th) == 0) {
                d@sdmFormula@data.terms[[which(.cls == '.selectFrame')]]@th <- .slf@th <- 10
              } else {
                if (.slf@th < 1) {
                  warning('Are you sure the right threshold is specified for the vifstep method in the select argument in the Formula?! (for vifstep, a commonly used threshold is between 5 to 10)')
                  .addLog(d) <-'Are you sure the right threshold is specified for the vifstep method in the select argument in the Formula?! (for vifstep, a commonly used threshold is between 5 to 10)'
                }
              }
              #----
              
              # vifstep
              .v <- .eval("vifstep(d@features[,.svar],th=.slf@th,keep=.slf@keep,size=nrow(d@features)+1)",env=environment())
              # exclude based on vifstep:
              if (!inherits(.v,'try-error')) {
                if (length(.v@excluded) > 0) {
                  if (length(.v@excluded) == 1) .addLog(d) <- paste0('One variable had collinearity problem and so, it is excluded. Its name is: ',.v@excluded)
                  else .addLog(d) <- paste0(length(.v@excluded),' variables had collinearity problems and they are excluded. They include: ',paste(.v@excluded,collapse = ', '))
                  #----
                  .v <- .eval('usdm::exclude(d@features,.v)',env=environment())
                  if (!inherits(.v,'try-error')) {
                    d@features <- .v
                  } else {
                    .addLog(d) <- .v
                  }
                  
                } else {
                  .addLog(d) <- 'No input variables from those specified in the select argument in Formula had collinearity problem so none of them are excluded! '
                }
              } else {
                .addLog(d) <- .v
              }
              #-----
            } else {
              if (is.null(.slf@th) || length(.slf@th) == 0) {
                d@sdmFormula@data.terms[[which(.cls == '.selectFrame')]]@th <- .slf@th <- 0.9
              } else {
                if (.slf@th >= 1) {
                  warning('A wrong threshold is specified for the vifcor method in the select argument in the Formula. It is changed to 0.9!')
                  .addLog(d) <-'A wrong threshold is specified for the vifcor method in the select argument in the Formula. It is changed to 0.9!'
                }
              }
              #----
              # vifcor:
              .v <- .eval("vifcor(d@features[,.svar],th=.slf@th,keep=.slf@keep,size=nrow(d@features)+1)",env=environment())
              # exclude based on vifcor:
              if (!inherits(.v,'try-error')) {
                if (length(.v@excluded) > 0) {
                  if (length(.v@excluded) == 1) .addLog(d) <- paste0('One variable had collinearity problem and so, it is excluded. Its name is: ',.v@excluded)
                  else .addLog(d) <- paste0(length(.v@excluded),' variables had collinearity problems and they are excluded. They include: ',paste(.v@excluded,collapse = ', '))
                  #----
                  .v <- .eval('usdm::exclude(d@features,.v)',env=environment())
                  if (!inherits(.v,'try-error')) {
                    d@features <- .v
                  } else {
                    .addLog(d) <- .v
                  }
                  
                } else {
                  .addLog(d) <- 'No input variables from those specified in the select argument in Formula had collinearity problem so none of them are excluded! '
                }
              } else {
                .addLog(d) <- .v
              }
            }
            
          } else {
            warning('usdm package is required for multicollinearity check but it is not installed!')
          }
          
          
          #--------
        }
        
      }
      
      #*****
    } 
    #++++
  }
  d
}
#----

# create sdmdata object:
.createSdmdata <- function(train,formula=NULL,test,bg=NULL,crs=NULL,author=NULL,website=NULL,citation=NULL,help=NULL,description=NULL,date=NULL,license=NULL) {
  if (missing(test)) test <- NULL
  
  nFact <- nf <- nxy <- nsp <- ng <- nt <- ni <- NULL
  
  if (is.null(formula)) {
    nnFact <- nnf <- nnxy <- nnt <- NULL
    w <- .speciesDetect(train)
    if (!is.null(w$nFact)) {
      nFact <- w$nFact
      nnFact <- paste(paste('f(',w$nFact,')',sep=''),collapse='+')
    }
    if (!is.null(w$nf)) {
      nf <- w$nf
      nnf <- paste(w$nf,collapse='+')
    }
    if (!is.null(w$nxy)) {
      nxy <- w$nxy
      nnxy <- paste(paste('coords(',paste(w$nxy,collapse='+'),')',sep=''),collapse='+')
    }
    if (!is.null(w$nt)) {
      nt <- w$nt
      nnt <- paste(paste('time(',w$nt,')',sep=''),collapse='+')
    }
    formula <- as.formula(paste(paste(w$nsp,collapse="+"),'~',paste(c(nnf,nnFact,nnxy,nnt),collapse='+')),env = parent.frame())
  }
  
  exf <- .exFormula(formula,train)
  nall <- c(exf@vars@names)
  
  
  if (!.varExist(train,nall)) stop('one (or more) variable(s) specified in the formula does not exist in the train data...!')
  
  nsp <- exf@vars@species
  
  d <- new('sdmdata')
  
  d@sdmFormula <- exf
  
  w <- .dataClean(train)
  if (any(w[[2]] > 0)) {
    train <- w[[1]]
    ww <- c()
    if (w[[2]][1] > 0) .addLog(d) <- paste(w[[2]][1],'records with NA from the train data are removed')
    if (w[[2]][2] > 0) .addLog(d) <- paste(w[[2]][2],'duplicarted records from the train data are removed')
  }
  
  train$rID <- 1:nrow(train)
  
  train <- .int.to.numeric(train)
  
  
  if (!is.null(bg)) {
    w <- which(nsp %in% colnames(bg))
    if (length(w) > 0) {
      for (nnsp in nsp[w]) bg[,nnsp] <- 0
    }
    
    w <- which(!nsp %in% colnames(bg))
    if (length(w) > 0) {
      nnsp <- nsp[w]
      nnsp <- matrix(0,nrow=nrow(bg),ncol=length(nnsp))
      colnames(nnsp) <- nsp[w]
      bg <- cbind(bg,nnsp)
    }
    bg <- as.data.frame(bg)
    
    if (!.varExist(data.frame(bg),nall)) stop('One (or more) predictor variable(s) does not exist in the background data...!')
    
    w <- .dataClean(bg)
    if (any(w[[2]] > 0)) {
      bg <- w[[1]]
      ww <- c()
      if (w[[2]][1] > 0) .addLog(d) <- paste(w[[2]][1],'records with NA from the background data are removed')
      if (w[[2]][2] > 0) .addLog(d) <- paste(w[[2]][2],'duplicarted records from the background data are removed')
    }
    
    
    for (n in nsp) {
      if (is.factor(train[,n]) || is.character(train[,n])) {
        train[,n] <- as.character(train[,n])
        bg[,n] <- as.character(bg[,n])
      }
    }
    
    bg$rID <- (nrow(train)+1):(nrow(bg)+nrow(train))
    train <- rbind(train[,c('rID',nall)],bg[,c('rID',nall)])
    bg <- bg$rID
    exf <- .exFormula(formula,train)
    d@sdmFormula <- exf
  }
  
  if (!is.null(test)) {
    if (!.varExist(test,nall)) stop('One (or more) variable(s) specified in the formula does not exist in the test data...!')
    w <- .dataClean(test)
    if (any(w[[2]] > 0)) {
      test <- w[[1]]
      ww <- c()
      if (w[[2]][1] > 0) .addLog(d) <- paste(w[[2]][1],'records with NA from the test data are removed')
      if (w[[2]][2] > 0) .addLog(d) <- paste(w[[2]][2],'duplicarted records from the test data are removed')
    }
    
    test$rID <- (nrow(train)+1):(nrow(test)+nrow(train))
    
    if (!is.null(bg)) {
      w <- unlist(lapply(nsp,function(x) .speciesType(test[,x])))
      w <- unique(w)
      if (length(w) > 1) stop(paste('The independent test dataset has different types of records including',paste(w,collapse=', ')))
      else if (w == "Presence-Only") {
        d <- .newGroup(d,'training',index=list(train=train$rID,test=c(test$rID,bg)))
        cat('WARNING:\n The independent test dataset contains only presence records, thus, the background (pseuso-absences) is used as "Absence" in the dataset...!\n')
      } else d <- .newGroup(d,'training',index=list(train=train$rID,test=test$rID))
    } else d <- .newGroup(d,'training',index=list(train=train$rID,test=test$rID))
    test <- .int.to.numeric(test)
  } else {
    d <- .newGroup(d,'training',index=list(train=train$rID))
  }
  
  #-------
  
  if (is.null(nf) & is.null(nFact)) {
    nf <- exf@vars@numeric$names
    nFact <- names(exf@vars@categorical)
    
    nf <- .excludeVector(nf,nFact)
    if (!is.null(exf@vars@others)) {
      if ("x_coordinate" %in% exf@vars@others$type) nxy <- exf@vars@others$name[exf@vars@others$type %in% c("x_coordinate","y_coordinate")]
      if ("Date/Time" %in% exf@vars@others$type) nt <- exf@vars@others$name[exf@vars@others$type %in% c("Date/Time")]
      if ("group" %in% exf@vars@others$type) ng <- exf@vars@others$name[exf@vars@others$type %in% c("group")]
      if ("Info" %in% exf@vars@others$type) ni <- exf@vars@others$name[exf@vars@others$type %in% c("Info")]
    } else {
      w <- !colnames(train) %in% c(nall,'rID')
      if (any(w)) {
        nxy <- .which.is.coords(colnames(train)[w])
        if (!is.null(test) && !.varExist(test,nxy)) nxy <- NULL
        ww <- unlist(lapply(which(w),function(x) class(train[,x]))) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr")
        if (any(ww)) {
          nt <- colnames(train)[which(w)[which(ww)]]
          nf <- .excludeVector(nf,nt)
          nFact <- .excludeVector(nFact,nt)
        }
      }
    }
    
    nf <- .excludeVector(nf,nxy)
    nFact <- .excludeVector(nFact,nxy)
  }
  #-------
  nall <- c(nsp,nf,nFact,nxy,ng,ni,nt,'rID')
  if (!is.null(test)) {
    train <- rbind(train[,nall],test[,nall])
    rm(test)
    species <- .getSpecies(train,nsp,bg=!is.null(bg),id.train = d@groups$training@indices$train,id.test = d@groups$training@indices$test)
  } else {
    train <- train[,nall]
    species <- .getSpecies(train,nsp,bg=!is.null(bg))
  }
  #----
  
  if (!is.null(ng)) {
    for (n in ng) {
      ww <- as.character(train[,n])
      u <- unique(ww)
      if (length(u) == 1) warning(paste('the grouping variable',n,'is ignored; it is constant!'))
      else {
        w <- list()
        for (uu in u) {
          w[[uu]] <- train$rID[which(ww == uu)]
        }
        d <- .newGroup(d,n,index=w)
      }
    }
  }
  #------
  if (!is.null(c(nxy,ni,nt,website,help,description,date,license)) || !is.null(c(citation,author))) {
    d@info <- new('.info')
    if (!is.null(nxy)) {
      d@info@coords <- as.matrix(train[,c('rID',nxy)])
      if (!is.null(crs) && is.character(crs)) d@info@crs <- crs
    }
    if (!is.null(ni)) d@info@info <- train[,c('rID',ni)]
    
    if (!is.null(nt)) {
      dt <- data.frame(rID=train$rID)
      for (n in nt) {
        if ((class(train[,n]) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr"))[1]) {
          dt <- cbind(dt,train[,n])
          colnames(dt)[ncol(dt)] <- n
        } else {
          w <- unlist(lapply(exf@data.terms,class))
          w <- exf@data.terms[which(w == ".time")]
          w <- w[[which(unlist(lapply(w,function(x) x@terms[[1]] == n)))]]
          if (length(w@terms) > 1 && is.null(names(w@terms)) && length(names(w@terms)) == length(w@terms)) {
            w <- do.call(.char2time,c(list(d=train[,n]),w@terms[2:length(w@terms)]))
            if ((class(w) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr"))[1]) dt[,n] <- w
            else warning(paste('a time-based format is not detected for variable',n,", so it is IGNORED; it must have a detectable character format, or being one of time-based classes including: 'POSIXct', 'POSIXt', 'Date', 'yearmon','yearqtr'"))
          } else {
            w <- .char2time(train[,n])
            if ((class(w) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr"))[1]) dt[,n] <- w
            else warning(paste('a time-based format is not detected for variable',n,", so it is IGNORED; it must have a detectable character format, or being one of time-based classes including: 'POSIXct', 'POSIXt', 'Date', 'yearmon','yearqtr'"))
          }
        }
      }
      if (ncol(dt) > 1) d@info@time <- dt
    }
    
    if (!is.null(c(website,help,description,date,license)) || !is.null(c(citation,author))) {
      d@info@metadata <- .newMetadata(authors=author,web=website,cit=citation,desc=description,date=date,license=license,help=help)
    }
  }
  #----------
  
  d@features.name <- c(nf,nFact)
  if (!is.null(nFact)) d@factors <- nFact
  d@species <- species
  d@species.names <- names(species)
  for (n in nFact) train[,n] <- factor(train[,n])
  if (!is.null(d@features.name)) d@features <- train[,c('rID',nf,nFact)]
  #---
  d <- .selectData(d)
  
  if (!is.null(d@features.name) && !all(d@features.name %in% colnames(d@features))) {
    d@features.name <- d@features.name[d@features.name %in% colnames(d@features)]
    .rem <- !(d@sdmFormula@vars@numeric$names %in% colnames(d@features))
    if (any(.rem)) {
      d@sdmFormula@vars@names <- .excludeVector(d@sdmFormula@vars@names,d@sdmFormula@vars@numeric$names[.rem])
      d@sdmFormula@vars@numeric <- d@sdmFormula@vars@numeric[!.rem,]
    }
  }
  d
}
#-------

.Extract <- function(x,cells,factors) {
  n <- names(x)
  
  if (length(factors) == 1) {
    x2 <- values(x[[factors]])[cells]
  } else x2 <- values(x[[factors]])[cells,]
  
  if (length(n) > length(factors)) {
    x1 <- x[[-factors]][cells]
    d <- data.frame(x1,x2)
    colnames(d) <- c(n[-factors],n[factors])
  } else {
    d <- data.frame(x2)
    colnames(d) <- n[factors]
  }
  
  for (.n in n[factors]) d[,.n] <- factor(d[,.n])
  return(d )
}
#---------


if (!isGeneric("sdmData")) {
  setGeneric("sdmData", function(formula, train, predictors, test,bg, filename, crs,impute,metadata,...)
    standardGeneric("sdmData"))
}


setMethod('sdmData', signature(train='data.frame',predictors='missing'), 
          function(formula,train,predictors,test=NULL,bg=NULL,filename=NULL,crs=NULL,impute,metadata,...) {
            if(missing(test)) test <- NULL
            if(missing(filename)) filename <- NULL
            if(missing(crs)) crs <- NULL
            if(missing(formula)) formula <- NULL
            if(missing(bg)) bg <- NULL
            if(missing(metadata)) metadata <- NULL
            #---
            if(missing(impute)) impute <- "neighbor"
            else if (is.null(impute)) impute <- "none"
            else if (is.logical(impute)) {
              if (impute) impute <- "neighbor"
              else impute <- "none"
            } else if (is.character(impute)) {
              if (impute %in% c('neigh','neighbor','neighbour','neighCell')) impute <- 'neighbor'
              else {
                warning('Currently, the only available imputation method is "neighbor", so the impute is changed to it...')
                impute <- 'neighbor'
              }
            }
            
            #---
            if (!.sdmOptions$getOption('sdmLoaded')) suppressMessages(.addMethods())
            #---
            if (is.null(metadata)) metadata <- list()
            if (!is.list(metadata)) {
              warning('metadata is ignored (it should be a list)!')
              metadata <- list()
            }
            
            n <- tolower(names(metadata))
            for (i in seq_along(n)) {
              if (any(!is.na(pmatch(c("aut"),n[i])))) n[i] <- 'author'
              else if (any(!is.na(pmatch(c("web"),n[i])))) n[i] <- 'website'
              else if (any(!is.na(pmatch(c("cit"),n[i])))) n[i] <- 'citation'
              else if (any(!is.na(pmatch(c("hel"),n[i])))) n[i] <- 'help'
              else if (any(!is.na(pmatch(c("des"),n[i])))) n[i] <- 'description'
              else if (any(!is.na(pmatch(c("dat"),n[i])))) n[i] <- 'date'
              else if (any(!is.na(pmatch(c("lic"),n[i])))) n[i] <- 'license'
            }
            names(metadata) <- n
            author <- metadata[['author']]
            website <- metadata[['website']]
            citation <- metadata[['citation']]
            help <- metadata[['help']]
            description <- metadata[['description']]
            date <- metadata[['date']]
            license <- metadata[['license']]
            #-----------
            ############
            
            if (inherits(train,'sf')) {
              
              if (!.require('sf')) stop('To handle the train data, the package "sf" is required, but it is not installed!')
              
              .xy <- .eval('st_coordinates(train)',env = environment())
              nxy <- colnames(.xy)
              
              if (!is.null(test)) {
                if (inherits(test,'sf')) {
                  .xyt <- .eval('st_coordinates(test)',env = environment())
                  nxyt <- colnames(.xyt)
                  if (ncol(test) == 1) test <- data.frame(SPECIES=rep(1,nrow(test)),.xyt)
                  else test <- data.frame(as(test,'data.frame')[,-which(colnames(test) == "geometry")],.xyt)
                  if (all(nxy != nxyt)) {
                    colnames(test)[unlist(lapply(nxyt,function(x) which(colnames(test) == x)))] <- nxy
                  } 
                }
              }
              
              if (!.eval("is.na(st_crs(train))",env = environment())) crs <- .eval("st_crs(train)$wkt",env = environment())
              
              if (ncol(train) == 1) test <- data.frame(SPECIES=rep(1,nrow(train)),.xy)
              else train <- data.frame(as(train,'data.frame')[,-which(colnames(train) == "geometry")],.xy)
              
              
              if (!missing(formula)) {
                if (!all(nxy %in% all.vars(formula))) {
                  if ('.' %in% all.vars(formula)) {
                    ww <- .exFormula(formula,train)
                    nw <- .excludeVector(colnames(train),c(ww@vars@species,nxy))
                    if (length(nw) > 0) {
                      w <- as.character(.exFormula(formula,train,FALSE)@vars@species)
                      if (length(w) == 0) formula[[2]] <- terms.formula(formula,data=train[,nw])[[2]]
                      else {
                        formula[[3]] <- terms.formula(formula,data=train[,nw],simplify = TRUE)[[3]]
                        if (colnames(train)[1] == 'SPECIES') {
                          colnames(train)[1] <- w[1]
                          if (!is.null(test) && colnames(test)[1] == 'SPECIES') colnames(test)[1] <- w[1]
                        }
                      }
                      
                    }
                  }
                  
                  formula <- update(formula,as.formula(paste('~ . + coords(',paste(nxy,collapse='+'),')',sep='')))
                }
              } else formula <- as.formula(paste('~ . + coords(',paste(nxy,collapse='+'),')',sep=''))
              
            }
            #----
            .createSdmdata(train = train, formula = formula, test = test,bg=bg,crs = crs,author=author,website=website,citation=citation,help=help,description=description,date=date,license=license)
          }
)

#-------  

setMethod('sdmData', signature(formula='data.frame',train='formula',predictors='missing'), 
          function(formula,train,predictors,test=NULL,bg=NULL,filename=NULL,crs=NULL,impute,metadata,...) {
            # to make it user friendly
            if(missing(test)) test <- NULL
            if(missing(filename)) filename <- NULL
            if(missing(crs)) crs <- NULL
            if(missing(bg)) bg <- NULL
            if(missing(metadata)) metadata <- NULL
            
            if(missing(impute)) impute <- "neighbor"
            else if (is.null(impute)) impute <- "none"
            else if (is.logical(impute)) {
              if (impute) impute <- "neighbor"
              else impute <- "none"
            } else if (is.character(impute)) {
              if (impute %in% c('neigh','neighbor','neighbour','neighCell')) impute <- 'neighbor'
              else {
                warning('Currently, the only available imputation method is "neighbor", so the impute is changed to it...')
                impute <- 'neighbor'
              }
            }
            
            #---
            
            sdmData(formula=train, train=formula,test=test,bg=bg,filename=filename,crs=crs,metadata=metadata,...)
          }
)

setMethod('sdmData', signature(formula='data.frame',train='missing',predictors='missing'), 
          function(formula,train,predictors,test=NULL,bg=NULL,filename=NULL,crs=NULL,impute,metadata,...) {
            # to make it user friendly
            if(missing(test)) test <- NULL
            if(missing(filename)) filename <- NULL
            if(missing(crs)) crs <- NULL
            if(missing(bg)) bg <- NULL
            if(missing(metadata)) metadata <- NULL
            
            if(missing(impute)) impute <- "neighbor"
            else if (is.null(impute)) impute <- "none"
            else if (is.logical(impute)) {
              if (impute) impute <- "neighbor"
              else impute <- "none"
            } else if (is.character(impute)) {
              if (impute %in% c('neigh','neighbor','neighbour','neighCell')) impute <- 'neighbor'
              else {
                warning('Currently, the only available imputation method is "neighbor", so the impute is changed to it...')
                impute <- 'neighbor'
              }
            }
            #---
            
            sdmData(train=formula,test=test,bg=bg,filename=filename,crs=crs,metadata=metadata,...)
          }
)

setMethod('sdmData', signature(train='SpatialPoints',predictors='missing'), 
          function(formula,train,predictors,test=NULL,bg=NULL,filename=NULL,crs=NULL,impute=TRUE,metadata=NULL,...) {
            
            if(missing(test)) test <- NULL
            if(missing(filename)) filename <- NULL
            if(missing(crs)) crs <- NULL
            else {
              if (is.character(crs)) {
                if (crs(train,asText=TRUE) == "" || is.na(crs(train,asText=TRUE))) crs(train) <- crs
                #---
                if (!is.null(test) && inherits(test,'Spatial') && is.na(crs(test,asText=TRUE))) crs(test) <- crs
              }
            }
            
            if(missing(bg)) bg <- NULL
            if(missing(metadata)) metadata <- NULL
            
            #---
            if(missing(impute)) impute <- "neighbor"
            else if (is.null(impute)) impute <- "none"
            else if (is.logical(impute)) {
              if (impute) impute <- "neighbor"
              else impute <- "none"
            } else if (is.character(impute)) {
              if (impute %in% c('neigh','neighbor','neighbour','neighCell')) impute <- 'neighbor'
              else {
                warning('Currently, the only available imputation method is "neighbor", so the impute is changed to it...')
                impute <- 'neighbor'
              }
            }
            
            #---
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            #---
            if (!is.na(projection(train))) crs <- projection(train)
            else {
              if (!is.null(crs) && is.character(crs)) projection(train) <- crs
            }
            #-----
            
            nxy <- coordnames(train)
            
            .cls <- as.character(class(train))
            
            if (!is.null(test)) {
              .cls0 <- as.character(class(test))
              
              if (.cls == .cls0) {
                nxyt <- coordnames(test)
                if (.cls0 == 'SpatialPoints') test <- data.frame(SPECIES=rep(1,length(test)),as(test,'data.frame'))
                else test <- as(test,'data.frame')
                if (all(nxy != nxyt)) {
                  colnames(test)[unlist(lapply(nxyt,function(x) which(colnames(test) ==x)))] <- nxy
                } 
              }
            }
            
            
            if (.cls == 'SpatialPoints') train <- data.frame(SPECIES=rep(1,length(train)),as(train,'data.frame'))
            else train <- as(train,'data.frame')
            
            if (!missing(formula)) {
              if (!all(nxy %in% all.vars(formula))) {
                if ('.' %in% all.vars(formula)) {
                  ww <- .exFormula(formula,train)
                  nw <- .excludeVector(colnames(train),c(ww@vars@species,nxy))
                  if (length(nw) > 0) {
                    w <- as.character(.exFormula(formula,train,FALSE)@vars@species)
                    if (length(w) == 0) formula[[2]] <- terms.formula(formula,data=train[,nw])[[2]]
                    else {
                      formula[[3]] <- terms.formula(formula,data=train[,nw],simplify = TRUE)[[3]]
                      if (colnames(train)[1] == 'SPECIES') {
                        colnames(train)[1] <- w[1]
                        if (!is.null(test) && colnames(test)[1] == 'SPECIES') colnames(test)[1] <- w[1]
                      }
                    }
                     
                  }
                }
                
                formula <- update(formula,as.formula(paste('~ . + coords(',paste(nxy,collapse='+'),')',sep='')))
              }
            } else formula <- as.formula(paste('~ . + coords(',paste(nxy,collapse='+'),')',sep=''))
            
            #------
            if (is.null(metadata)) metadata <- list()
            if (!is.list(metadata)) {
              warning('metadata is ignored (it should be a list)!')
              metadata <- list()
            }
            
            n <- tolower(names(metadata))
            for (i in seq_along(n)) {
              if (any(!is.na(pmatch(c("aut"),n[i])))) n[i] <- 'author'
              else if (any(!is.na(pmatch(c("web"),n[i])))) n[i] <- 'website'
              else if (any(!is.na(pmatch(c("cit"),n[i])))) n[i] <- 'citation'
              else if (any(!is.na(pmatch(c("hel"),n[i])))) n[i] <- 'help'
              else if (any(!is.na(pmatch(c("des"),n[i])))) n[i] <- 'description'
              else if (any(!is.na(pmatch(c("dat"),n[i])))) n[i] <- 'date'
              else if (any(!is.na(pmatch(c("lic"),n[i])))) n[i] <- 'license'
            }
            names(metadata) <- n
            author <- metadata[['author']]
            website <- metadata[['website']]
            citation <- metadata[['citation']]
            help <- metadata[['help']]
            description <- metadata[['description']]
            date <- metadata[['date']]
            license <- metadata[['license']]
            #-----------
            
            
            .createSdmdata(train = train, formula = formula, test = test,bg=bg,crs = crs,author=author,website=website,citation=citation,help=help,description=description,date=date,license=license) 
            
          }
)


setMethod('sdmData', signature(train='SpatialPoints',predictors='Raster'), 
          function(formula,train,predictors,test=NULL,bg=NULL,filename=NULL,crs=NULL,impute,metadata=NULL,...) {
            if(missing(test)) test <- NULL
            if(missing(filename)) filename <- NULL
            if(missing(crs)) crs <- NULL
            else {
              if (is.character(crs)) {
                if (crs(train,asText=TRUE) == "" || is.na(crs(train,asText=TRUE))) crs(train) <- crs
                #---
                if (!is.null(test) && inherits(test,'Spatial') && is.na(crs(test,asText=TRUE))) crs(test) <- crs
              }
            }
            
            
            if(missing(bg)) bg <- NULL
            if(missing(metadata)) metadata <- NULL
            #---
            if(missing(impute)) impute <- "neighbor"
            else if (is.null(impute)) impute <- "none"
            else if (is.logical(impute)) {
              if (impute) impute <- "neighbor"
              else impute <- "none"
            } else if (is.character(impute)) {
              if (impute %in% c('neigh','neighbor','neighbour','neighCell')) impute <- 'neighbor'
              else {
                warning('Currently, the only available imputation method is "neighbor", so the impute is changed to it...')
                impute <- 'neighbor'
              }
            }
            
            #---
            if (!.sdmOptions$getOption('sdmLoaded')) suppressMessages(.addMethods())
            #---
            #-----------------------
            if (!is.na(projection(train))) crs <- projection(train)
            else {
              if (!is.null(crs) && is.character(crs)) projection(train) <- crs
              else if (projection(predictors) != "" && !is.na(projection(predictors))) {
                projection(train) <- projection(predictors)
              }
            }
            #-----
            
            
            #----
            testc <- as.character(class(test))
            trainc <- as.character(class(train))
            
            errLog <- list()
            
            nxy <- coordnames(train)[1:2]
            wF <- is.factor(predictors)
           
            if (!is.null(test)) {
              if (inherits(test,'SpatialPoints')) {
                
                if (is.na(projection(test))) {
                  if (!is.na(projection(train))) projection(test) <- projection(train)
                  else if (!is.na(projection(predictors))) projection(test) <- projection(predictors)
                }
                
                
                nxyt <- coordnames(test)
                if (testc == 'SpatialPointsDataFrame') test <- as(test,'data.frame')
                else test <- data.frame(coordinates(test))
                if (all(nxy != nxyt)) {
                  colnames(test)[unlist(lapply(nxyt,function(x) which(colnames(test) ==x)))] <- nxy
                } 
              } else if (!is.data.frame(test)) stop('test data should be a data.frame or in the same class as the train data!')
              
              if (nxy[1] %in% colnames(test) & nxy[2] %in% colnames(test)) {
                cells <- cellFromXY(predictors,test[,nxy])
                cNA <- is.na(cells)
                if (any(cNA)) {
                  if (all(cNA)) stop('Test dataset has no overlap with the predictors...!')
                  wNA <- which(cNA)
                  test <- test[-wNA,]
                  cells <- cells[-wNA]
                  errLog <- c(errLog,paste(length(wNA),'records were removed from the test dataset because of no overlap with the predictors.'))
                  rm(wNA)
                }
                rm(cNA)
                if (!any(wF)) test.p <- data.frame(predictors[cells])
                else test.p <- .Extract(predictors,cells,which(wF))
                rm(cells)
                colnames(test.p) <- names(predictors)
                
                w <- colnames(test) %in% colnames(test.p)
                if (any(w)) {
                  test <- test[,colnames(test)[!w]]
                  errLog <- c(errLog,paste('WARNING: The variables',colnames(test)[w],'were removed from the test dataset as they exist in the predictors as well.'))
                }
              } else stop('the coordinate names in the train and test datasets do not match!')
            }
            
            
            
            if (trainc == 'SpatialPointsDataFrame') train <- as(train,'data.frame')
            else train <- coordinates(train)
            
            cells <- cellFromXY(predictors,train[,nxy])
            cNA <- is.na(cells)
            if (any(cNA)) {
              if (all(cNA)) stop('Train data has no overlap with the predictors...!')
              wNA <- which(cNA)
              train <- train[-wNA,]
              cells <- cells[-wNA]
              errLog <- c(errLog,paste(length(wNA),'records were removed from the train dataset as they had no overlap with the predictors.'))
            }
            rm(cNA)
            
            if (!any(wF)) train.p <- data.frame(predictors[cells])
            else train.p <- .Extract(predictors,cells,which(wF))
            rm(cells)
            colnames(train.p) <- names(predictors)
            w <- colnames(train) %in% colnames(train.p)
            if (any(w)) {
              train <- train[,colnames(train)[!w]]
              errLog <- c(errLog,paste('WARNING: The variables',colnames(train)[w],'were removed from the train dataset as they exist in the predictors as well.'))
            }
            
            if (trainc == 'SpatialPointsDataFrame') {
              train <- data.frame(train,train.p)
              rm(train.p)
              if (!is.null(test)) {
                test <- data.frame(test,test.p)
                rm(test.p)
              }
              
              if (!missing(formula)) {
                if (!all(nxy %in% all.vars(formula))) {
                  if ('.' %in% all.vars(formula)) {
                    ww <- .exFormula(formula,train)
                    nw <- .excludeVector(colnames(train),c(ww@vars@species,nxy))
                    w <- as.character(.exFormula(formula,train,FALSE)@vars@species)
                    if (length(w) == 0) formula[[2]] <- terms.formula(formula,data=train[,nw])[[2]]
                    else formula[[3]] <- terms.formula(formula,data=train[,nw],simplify = TRUE)[[3]]
                  }
                  formula <- update(formula,as.formula(paste('~ . + coords(',paste(nxy,collapse='+'),')',sep='')))
                }
              } else formula <- as.formula(paste('~ . + coords(',paste(nxy,collapse='+'),')',sep=''))
            } else {
              if (!missing(formula)) {
                ww <- as.character(.exFormula(formula,train.p,FALSE)@vars@species)
                if (length(ww) > 1) {
                  warning('While SpatialPoints can be used for only 1 species, more names are defined in the formula! The first name is considered!')
                  errLog <- c(errLog,'WARNING: While SpatialPoints can be used for only 1 species, more names are defined in the formula! The first name is considered!')
                  w <- ww[1]
                } else if (length(ww) == 0) w <- 'SPECIES'
                
                if (!all(nxy %in% all.vars(formula))) {
                  
                  if ('.' %in% all.vars(formula)) {
                    if (length(ww) == 0) formula[[2]] <- terms.formula(formula,data=train.p)[[2]]
                    else formula[[3]] <- terms.formula(formula,data=train.p,simplify = TRUE)[[3]]
                  }
                  formula <- update(formula,as.formula(paste('~ . + coords(',paste(nxy,collapse='+'),')',sep='')))
                }
              } else {
                formula <- as.formula(paste('SPECIES ~ . + coords(',paste(nxy,collapse='+'),')',sep=''))
                w <- 'SPECIES'
              }
              
              train <- data.frame(SPECIES=rep(1,nrow(train)),train)
              colnames(train)[1] <- w
              train <- data.frame(train,train.p)
              rm(train.p)
              
              if (!is.null(test)) {
                test <- data.frame(SPECIES=rep(1,nrow(test)),test,test.p)
                colnames(test)[1] <- w
                rm(test.p)
              }
            } 
            ###################################
            
            if (!is.null(bg)) {
              if (inherits(bg,'list')) {
                ww <- .exFormula(formula,train)
                #--
                #--- .p (presence locations; only needed for gDist method)
                nsp <- ww@vars@species # name of species
                if (length(nsp) > 0) {
                  w <- c()
                  for (n in nsp) {
                    w <- c(w,which(train[,n] == 1))
                  }
                  w <- unique(w)
                  .p <- train[w,nxy]
                } else .p <- train[,nxy] # assuming all records are presence-only!
                #----
                #.p 
                bg <- .pseudo.Raster(predictors,bg = bg,p = .p,factors = names(ww@vars@categorical))
                colnames(bg)[1:2] <- nxy
              } else if (is.numeric(bg)) {
                bg <- .pseudo_gRandom.Raster(predictors,n=bg[1])
                colnames(bg)[1:2] <- nxy
              } else if (is.data.frame(bg) || is.matrix(bg)) bg <- data.frame(bg)
              else {
                bg <- NULL
                warning('bg is not recognised so (ignored)! It should be either a list (settings of generating background), or a single number (number of background), or a data.frame (background records)')
              }
            }
            
            #################################
            
            #-----------
            if (is.null(metadata)) metadata <- list()
            if (!is.list(metadata)) {
              warning('metadata is ignored (it should be a list)!')
              metadata <- list()
            }
            
            n <- tolower(names(metadata))
            for (i in seq_along(n)) {
              if (any(!is.na(pmatch(c("aut"),n[i])))) n[i] <- 'author'
              else if (any(!is.na(pmatch(c("web"),n[i])))) n[i] <- 'website'
              else if (any(!is.na(pmatch(c("cit"),n[i])))) n[i] <- 'citation'
              else if (any(!is.na(pmatch(c("hel"),n[i])))) n[i] <- 'help'
              else if (any(!is.na(pmatch(c("des"),n[i])))) n[i] <- 'description'
              else if (any(!is.na(pmatch(c("dat"),n[i])))) n[i] <- 'date'
              else if (any(!is.na(pmatch(c("lic"),n[i])))) n[i] <- 'license'
            }
            names(metadata) <- n
            author <- metadata[['author']]
            website <- metadata[['website']]
            citation <- metadata[['citation']]
            help <- metadata[['help']]
            description <- metadata[['description']]
            date <- metadata[['date']]
            license <- metadata[['license']]
            #-----------
            
            d <- .createSdmdata(train = train, formula = formula, test = test,bg=bg,crs = crs,author=author,website=website,citation=citation,help=help,description=description,date=date,license=license)
            
            if (length(errLog) > 0) {
              for (i in seq_along(errLog)) .addLog(d) <- errLog[[i]]
            }
            d
          }
)
#-----------
setMethod('sdmData', signature(train='SpatVector',predictors='SpatRaster'), 
          function(formula,train,predictors,test=NULL,bg=NULL,filename=NULL,crs=NULL,impute=TRUE,metadata=NULL,...) {
            if(missing(test)) test <- NULL
            if(missing(filename)) filename <- NULL
            
            if(missing(crs)) crs <- NULL
            else {
              if (is.character(crs)) {
                if (crs(train,proj=TRUE) == "" || is.na(crs(train,proj=TRUE))) crs(train) <- crs
                #---
                if (!is.null(test) && inherits(test,'SpatVector') && crs(test,proj=TRUE) == "") crs(test) <- crs
              }
            }
            
            if(missing(bg)) bg <- NULL
            if(missing(metadata)) metadata <- NULL
            
            # the current imputation method is "neighbor"!
            if(missing(impute)) impute <- "neighbor"
            else if (is.null(impute)) impute <- "none"
            else if (is.logical(impute)) {
              if (impute) impute <- "neighbor"
              else impute <- "none"
            } else if (is.character(impute)) {
              if (impute %in% c('neigh','neighbor','neighbour','neighCell')) impute <- 'neighbor'
              else {
                warning('Currently, the only available imputation method is "neighbor", so the impute is changed to it...')
                impute <- 'neighbor'
              }
            }
            
            #---
            if (!.sdmOptions$getOption('sdmLoaded')) suppressMessages(.addMethods())
            #---
            errLog <- list()
            #----
            #-----------------------
            if (crs(train) != "" && !is.na(crs(train))) crs <- crs(train)
            else {
              if (!is.null(crs) && is.character(crs)) crs(train) <- crs
              else if (crs(predictors) != "" && !is.na(crs(predictors))) {
                crs(train) <- crs(predictors)
              }
            }
            #-----
            
            nxy <- colnames(crds(train))
            wF <- is.factor(predictors)
            
            if (!is.null(test)) {
              if (inherits(test,'SpatVector')) {
                
                if (crs(test) == "" || is.na(crs(test))) {
                  if (crs(train) != "" && !is.na(crs(train))) crs(test) <- crs(train)
                  else if (crs(predictors) != "" && !is.na(crs(predictors))) crs(test) <- crs(predictors)
                }
                
                testXY <- crds(test)
                test.df <- as.data.frame(test)
                
                if (nrow(test.df) == 0) stop('The independent "test" dataset has NO species information (has ONLY spatial coordinates)...!')
                
                if (all(nxy %in% colnames(test.df))) {
                  w <- which(colnames(test.df) %in% nxy)
                  test.df[,-w]
                }
                
                test.df <- data.frame(test.df,testXY)
                
                rm(testXY)
                
              } else if (is.data.frame(test)) {
                if (!all(nxy %in% colnames(test))) {
                  #warning('The records in test data have no coordinates!')
                  errLog <- c(errLog,'The test dataset has no coordinates!')
                }
                test.df <- test
              } else stop('test data should be either a data.frame or with the same class as the train data!')
              
              if (nxy[1] %in% colnames(test.df) & nxy[2] %in% colnames(test.df)) {
                cells <- cellFromXY(predictors,test.df[,nxy])
                cNA <- is.na(cells)
                
                if (any(cNA)) {
                  
                  if (all(cNA)) {
                    if (inherits(test,'SpatVector')) {
                      if ((crs(predictors,proj=TRUE) != "" && crs(test,proj=TRUE) != "") && (crs(predictors,proj=TRUE) != crs(test,proj=TRUE))) {
                        test <- project(test,predictors)
                        cells <- cellFromXY(predictors,geom(test)[,3:4])
                        cNA <- is.na(cells)
                        if (all(cNA)) stop('No spatial overlap between records in the "test" dataset and the predictors...!')
                        else {
                          test.df <- as.data.frame(test)
                          if (all(nxy %in% colnames(test.df))) {
                            w <- which(colnames(test.df) %in% nxy)
                            test.df[,-w]
                          }
                          
                          test.df <- data.frame(test.df,geom(test)[,3:4])
                        }
                      } else stop('No spatial overlap between records in the "test" dataset and the predictors...!')
                    } else stop('No spatial overlap between records in the "test" dataset and the predictors...!')
                  }
                  
                  wNA <- which(cNA)
                  test <- test[-wNA,]
                  test.df <- test.df[-wNA,]
                  cells <- cells[-wNA]
                  errLog <- c(errLog,paste(length(wNA),'records were removed from the test dataset because of they were located outside of the bounding box (extent) of the predictors.'))
                  rm(wNA)
                }
                rm(cNA)
                #---------
                
                
                test.p <- data.frame(predictors[cells])
                #if (!any(wF)) test.p <- data.frame(predictors[cells])
                #else test.p <- .Extract(predictors,cells,which(wF))
                rm(cells)
                #colnames(test.p) <- names(predictors)
                
                wNA <- which(apply(test.p,1,function(x) any(is.na(x)))) 
                
                if (length(wNA) > 0) {
                  #---- imputation method: "neighbor"
                  if (impute == "neighbor") {
                    for (w in wNA) {
                      .cxy <- cellFromXY(predictors,test.df[w,nxy])
                      if (!is.na(.cxy)) {
                        .nec <- .getNeighCellRandom(.cxy,predictors)
                        if (!is.na(.nec)) {
                          x <- predictors[.nec]
                          if (all(!is.na(x[1,]))) {
                            test.p[w,names(x)] <- x
                          }
                        }
                      }
                    }
                    #-------
                    wNA <- which(apply(test.p,1,function(x) any(is.na(x)))) 
                    if (length(wNA) > 0) {
                      test <- test[-wNA,]
                      test.p <- test.p[-wNA,]
                      test.df <- test.df[-wNA,]
                      errLog <- c(errLog,paste(length(wNA),'records were removed from the test dataset because they were located on cells with missing values (NA) in the predictors.'))
                    }
                  } else {
                    test <- test[-wNA,]
                    test.p <- test.p[-wNA,]
                    test.df <- test.df[-wNA,]
                    errLog <- c(errLog,paste(length(wNA),'records were removed from the test dataset because they were located on cells with missing values (NA) in the predictors.'))
                  }
                }
                #----
                
                
                w <- colnames(test.df) %in% colnames(test.p)
                if (any(w)) {
                  test.df <- test.df[,colnames(test.df)[!w]]
                  test <- test[,names(test)[!w]]
                  errLog <- c(errLog,paste('WARNING: The variables',colnames(test.df)[w],'were removed from the test dataset as they exist in the predictors as well.'))
                }
              } else stop('the coordinate names in the train and test datasets do not match!')
              
            }
            #-----------------------
            #-----------------------
            
            
            .PO <- FALSE # defined to separate between SpatVectors without and With attributes (case of SpatialPoints vs SpatialPointsDataFrame)
            
            trainXY <- crds(train)
            
            if (dim(train)[2] == 0) {
              train$SPECIES <- 1
              .PO <- TRUE
            }
            
            train.df <- as(train,'data.frame')
            train.df <- data.frame(train.df,trainXY)
            
            cells <- cellFromXY(predictors,trainXY)
            cNA <- is.na(cells)
            
            if (any(cNA)) {
              if (all(cNA)) {
                
                if ((crs(predictors,proj=TRUE) != "" && crs(train,proj=TRUE) != "") && (crs(predictors,proj=TRUE) != crs(train,proj=TRUE))) {
                  train <- project(train,predictors)
                  cells <- cellFromXY(predictors,crds(train))
                  cNA <- is.na(cells)
                  if (all(cNA)) stop('No spatial overlap between records in the "train" dataset and the predictors...!')
                  else {
                    warning('The species dataset has a different CRS than the predictors dataset; so it is projected to the CRS of "predictors"...!')
                    train.df <- as.data.frame(train)
                    if (all(nxy %in% colnames(train.df))) {
                      w <- which(colnames(train.df) %in% nxy)
                      train.df[,-w]
                    }
                    
                    train.df <- data.frame(train.df,crds(train))
                  }
                } else stop('No spatial overlap between records in the "train" dataset and the predictors...!')
                
              }
              wNA <- which(cNA)
              train <- train[-wNA,]
              train.df <- train.df[-wNA,]
              cells <- cells[-wNA]
              errLog <- c(errLog,paste(length(wNA),'records were removed from the "train" dataset because of they were located outside of the bounding box (extent) of the predictors.'))
            }
            rm(cNA)
            
            
            train.p <- predictors[cells]
            
            rm(cells)
            
            
            #---
            wNA <- which(apply(train.p,1,function(x) any(is.na(x)))) 
            
            if (length(wNA) > 0) {
              #---- imputation method: "neighbor"
              if (impute == "neighbor") {
                for (w in wNA) {
                  .cxy <- cellFromXY(predictors,train.df[w,nxy])
                  if (!is.na(.cxy)) {
                    .nec <- .getNeighCellRandom(.cxy,predictors)
                    if (!is.na(.nec)) {
                      x <- predictors[.nec]
                      if (all(!is.na(x[1,]))) {
                        train.p[w,names(x)] <- x
                      }
                    }
                  }
                }
                #-------
                wNA <- which(apply(train.p,1,function(x) any(is.na(x)))) 
                if (length(wNA) > 0) {
                  train <- train[-wNA,]
                  train.p <- train.p[-wNA,]
                  train.df <- train.df[-wNA,]
                  errLog <- c(errLog,paste(length(wNA),'records were removed from the "train" dataset because they were located on cells with missing values (NA) in the predictors.'))
                }
              } else {
                train <- train[-wNA,]
                train.p <- train.p[-wNA,]
                train.df <- train.df[-wNA,]
                errLog <- c(errLog,paste(length(wNA),'records were removed from the "train" dataset because they were located on cells with missing values (NA) in the predictors.'))
              }
            }
            #----
            #---------
            
            
            
            w <- colnames(train.df) %in% colnames(train.p)
            if (any(w)) {
              train.df <- train.df[,colnames(train.df)[!w]]
              train <- train[,names(train)[!w]]
              errLog <- c(errLog,paste('WARNING: The variables',colnames(train)[w],'were removed from the train dataset as they exist in the predictors as well.'))
            }
            #######################################################
            if (!.PO) {
              train <- data.frame(train.df,train.p)
              #rm(train.p,train.df)
              if (!is.null(test)) {
                test <- data.frame(test.df,test.p)
                rm(test.p,test.df)
              }
              
              if (!missing(formula)) {
                #if ((length(formula) == 3 && formula[[2]] == '.') || length(formula) == 2) nsp <- as.character(.exFormula(formula,train.df)@vars@species)
                if (!all(nxy %in% all.vars(formula))) {
                  if ('.' %in% all.vars(formula)) {
                    ww <- .exFormula(formula,train)
                    nw <- .excludeVector(colnames(train),c(ww@vars@species,nxy))
                    w <- as.character(.exFormula(formula,train,FALSE)@vars@species)
                    nw <- nw[nw %in% colnames(train.p)]
                    if (length(w) == 0) formula[[2]] <- terms.formula(formula,data=train[,nw])[[2]]
                    else formula[[3]] <- terms.formula(formula,data=train[,nw],simplify = TRUE)[[3]]
                  }
                  formula <- update(formula,as.formula(paste('~ . + coords(',paste(nxy,collapse='+'),')',sep='')))
                }
              } else formula <- as.formula(paste('~ . + coords(',paste(nxy,collapse='+'),')',sep=''))
            } else {
              if (!missing(formula)) {
                ww <- as.character(.exFormula(formula,train.p,FALSE)@vars@species)
                if (length(ww) > 1) {
                  warning('While SpatialPoints can be used for only 1 species, more names are defined in the formula! The first name is considered!')
                  errLog <- c(errLog,'WARNING: While SpatialPoints can be used for only 1 species, more names are defined in the formula! The first name is considered!')
                  w <- ww[1]
                } else if (length(ww) == 0) w <- 'SPECIES'
                
                names(train) <- w
                
                if (!all(nxy %in% all.vars(formula))) {
                  
                  if ('.' %in% all.vars(formula)) {
                    if (length(ww) == 0) formula[[2]] <- terms.formula(formula,data=train.p)[[2]]
                    else formula[[3]] <- terms.formula(formula,data=train.p,simplify = TRUE)[[3]]
                  }
                  formula <- update(formula,as.formula(paste('~ . + coords(',paste(nxy,collapse='+'),')',sep='')))
                }
              } else {
                formula <- as.formula(paste('SPECIES ~ . + coords(',paste(nxy,collapse='+'),')',sep=''))
                w <- 'SPECIES'
              }
              
              train <- data.frame(train)
              #colnames(train)[1] <- w
              train <- data.frame(train,train.p)
              
              
              if (!is.null(test)) {
                test <- data.frame(test.df,test.p)
                rm(test.p,test.df)
                #test <- data.frame(SPECIES=rep(1,nrow(test)),test,test.p)
                #colnames(test)[1] <- w
                rm(test.p)
              }
            } 
            #--------------------------
            rm(train.p,train.df)
            
            if (!is.null(bg)) {
              if (inherits(bg,'list')) {
                ww <- .exFormula(formula,train)
                #--
                #--- .p (presence locations; only needed for gDist method)
                nsp <- ww@vars@species # name of species
                if (length(nsp) > 0) {
                  w <- c()
                  for (n in nsp) {
                    w <- c(w,which(train[,n] == 1))
                  }
                  w <- unique(w)
                  .p <- train[w,nxy]
                } else .p <- train[,nxy] # assuming all records are presence-only!
                #----
                #.p 
                bg <- .pseudo.terra(predictors,bg = bg,p = .p,factors = names(ww@vars@categorical))
                colnames(bg)[1:2] <- nxy
              } else if (is.numeric(bg)) {
                bg <- .pseudo_gRandom.Terra(predictors,n=bg[1])
                colnames(bg)[1:2] <- nxy
              } else if (is.data.frame(bg) || is.matrix(bg)) bg <- data.frame(bg)
              else {
                bg <- NULL
                warning('bg is not recognised so (ignored)! It should be either a list of settings used to generate background, or a data.frame of background records')
              }
            }
            ########################
            if (is.null(metadata)) metadata <- list()
            if (!is.list(metadata)) {
              warning('metadata is ignored (it should be a list)!')
              metadata <- list()
            }
            
            n <- tolower(names(metadata))
            for (i in seq_along(n)) {
              if (any(!is.na(pmatch(c("aut"),n[i])))) n[i] <- 'author'
              else if (any(!is.na(pmatch(c("web"),n[i])))) n[i] <- 'website'
              else if (any(!is.na(pmatch(c("cit"),n[i])))) n[i] <- 'citation'
              else if (any(!is.na(pmatch(c("hel"),n[i])))) n[i] <- 'help'
              else if (any(!is.na(pmatch(c("des"),n[i])))) n[i] <- 'description'
              else if (any(!is.na(pmatch(c("dat"),n[i])))) n[i] <- 'date'
              else if (any(!is.na(pmatch(c("lic"),n[i])))) n[i] <- 'license'
            }
            names(metadata) <- n
            author <- metadata[['author']]
            website <- metadata[['website']]
            citation <- metadata[['citation']]
            help <- metadata[['help']]
            description <- metadata[['description']]
            date <- metadata[['date']]
            license <- metadata[['license']]
            #-----------
            
            d <- .createSdmdata(train = train, formula = formula, test = test,bg=bg,crs = crs,author=author,website=website,citation=citation,help=help,description=description,date=date,license=license)
            
            if (length(errLog) > 0) {
              for (i in seq_along(errLog)) .addLog(d) <- errLog[[i]]
            }
            d
          }
)
#---------

setMethod('sdmData', signature(train='data.frame',predictors='SpatRaster'), 
          function(formula,train,predictors,test=NULL,bg=NULL,filename=NULL,crs=NULL,impute=TRUE,metadata=NULL,...) {
            if(missing(test)) test <- NULL
            if(missing(filename)) filename <- NULL
            if(missing(crs) || !is.character(crs)) crs <- ""
            
            if(missing(bg)) bg <- NULL
            if(missing(metadata)) metadata <- NULL
            
            # the current imputation method is "neighbor"!
            if(missing(impute)) impute <- "neighbor"
            else if (is.null(impute)) impute <- "none"
            else if (is.logical(impute)) {
              if (impute) impute <- "neighbor"
              else impute <- "none"
            } else if (is.character(impute)) {
              if (impute %in% c('neigh','neighbor','neighbour','neighCell')) impute <- 'neighbor'
              else {
                warning('Currently, the only available imputation method is "neighbor", so the impute is changed to it...')
                impute <- 'neighbor'
              }
            }
            #---
            if (inherits(train,'sf')) {
              
              train <- vect(train)
              
              if (!is.null(test)) {
                if (inherits(test,'sf')) {
                  test <- vect(test)
                } else if (inherits(train,'data.frame')) {
                  nxy <- .which.is.coords(colnames(test))
                  test <- vect(test,geom=nxy)
                } else stop('(Independent) test data should be either data.frame (with coordinates) or a spatial points object (e.g., SpatVector, sf, SpatialPoints)')
              }
              
            } else {
              if (!missing(formula)) {
                .tmp <- cbind(train[1:5,],spatSample(predictors,size=5,na.rm=TRUE))
                ww <- .exFormula(formula,.tmp)
                if (!is.null(ww@vars@others) && "x_coordinate" %in% ww@vars@others$type) {
                  nxy <- ww@vars@others$name[c(which(ww@vars@others$type == "x_coordinate"),which(ww@vars@others$type == "y_coordinate"))]
                  train <- vect(train,geom=nxy,crs=crs)
                  
                  if (!is.null(test)) {
                    if (inherits(test,'sf') || inherits(test,'SpatialPoints')) {
                      test <- vect(test)
                    } else if (inherits(test,'data.frame')) {
                      if (!all(nxy %in% colnames(test))) stop('coordinates columns do not exist in test data!')
                      test <- vect(test,geom=nxy,crs=crs)
                    } else {
                      if (!inherits(test,'SpatVector')) stop('test data should be either data.frame (with coordinates) or a spatial points object (e.g., SpatVector, sf, SpatialPoints)')
                    }
                  }
                  formula <- .rmCoordsInFormula(formula)
                } else stop("You need to have coordinates (which should be introduced in the formula using 'coords(x+y)') in the species data.frame (train) when predictors is a Raster object!")
            } else stop('Since train data is a data.frame, you need to specify coordinates columns in the formula argument...!')
            #----
            }
            
            sdmData(formula=formula, train=train,predictors=predictors,test=test,bg=bg,filename=filename,crs=crs,metadata=metadata,impute=impute,...)
          }
)
