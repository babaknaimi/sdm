# Author: Babak Naimi, naimi.b@gmail.com
# Date of last update :  May 2024
# Version 1.6
# Licence GPL v3

#------



#----
.where <- function(f, x) {
  vapply(x, f, logical(1))
}
#------
.int.to.numeric <- function(data) {
  w <- which(unlist(lapply(data,is.integer)))
  if (length(w) > 0) {
    for (i in w) data[,i] <- as.numeric(data[[i]])
  }
  data
}
#----------
.is.Date <- function(x) {
  inherits(x, c("Date", "POSIXt"))
}
#----

#--------------------

.which.is.coords <- function(n) {
  nxy <- NULL
  .n <- tolower(n)
  w <- which(.n %in% c('lon','long','longitude'))
  if (length(w) == 1 & any(.n %in% c('lat','latitude'))) {
    nxy <- c(n[w],n[.n %in% c('lat','latitude')])
    return(nxy)
  }
  
  w <- which(.n %in% c('coords.x','coords.x1'))
  if (length(w) == 1 & any(.n %in% c('coords.y','coords.x2'))) {
    nxy <- c(n[w],n[.n %in% c('coords.y','coords.x2')])
    return(nxy)
  }
  
  w <- which(.n == 'x')
  if (length(w) == 1 & any(.n == c('y'))) {
    nxy <- c(n[w],n[.n %in% c('y')])
    return(nxy)
  }
  
}
#---------------

.normalize <- function(x,except=NULL) {
  w <- !.where(is.factor,x)
  if (!is.null(except)) {
    w[except] <- FALSE
  }
  if (any(w)) {
    xx <- x
    for (i in seq_along(w)) {
      if (w[i]) {
        xx[,i] <- xx[,i] - mean(xx[,i],na.rm=TRUE)
        if (sd(x[,i],na.rm=TRUE) != 0) xx[,i] <- xx[,i] / sd(x[,i],na.rm=TRUE)
      }
    }
  }
  xx
}
#-----------
# given a vector of colnames, their corresponding col numbers are returned
.colNumber <- function(d,n) {
  unlist(lapply(n,function(x) which(colnames(d) == x)))
}
#----------
.char2time <- function(d,...) {
  if (length(list(...)) > 0) {
    tst <- try(as.POSIXct(d[1],...),silent=TRUE)
    if (!inherits(tst, "try-error") && !is.na(tst)) return(as.POSIXct(d,...))
    else {
      tst <- try(as.POSIXct(d[1]),silent=TRUE)
      if (!inherits(tst, "try-error") && !is.na(tst)) return(as.POSIXct(d))
      else {
        tst <- try(as.Date(d[1],...),silent=TRUE)
        if (!inherits(tst, "try-error") && !is.na(tst)) return(as.Date(d,...))
        else {
          tst <- try(as.Date(d[1]),silent=TRUE)
          if (!inherits(tst, "try-error") && !is.na(tst)) return(as.Date(d))
          else return(NA)
        }
      }
    }
  } else {
    tst <- try(as.POSIXct(d[1]),silent=TRUE)
    if (!inherits(tst, "try-error") && !is.na(tst)) return(as.POSIXct(d))
    else {
      tst <- try(as.Date(d[1]),silent=TRUE)
      if (!inherits(tst, "try-error") && !is.na(tst)) return(as.Date(d))
      else return(NA)
    }
  }
}
#-----------
# check whether the names (vars) do exist in data (data.frame)
.varExist <-function(data,vars) {
  all(vars %in% names(data))
}
#----------

.isBinomial <- function(x) {
  if (is.numeric(x)) {
    u <- unique(x)
    if (length(u) > 2) return(FALSE)
    else if (length(u) == 2) return(all(sort(u) == c(0,1)) | all(sort(u) == c(-1,1)))
    else return(u == 1)
  } else if (is.logical(x)) return(TRUE)
  else {
    x <- as.character(x)
    u <- unique(x)
    if (length(u) > 2) return(FALSE)
    else if (length(u) == 2) return(all(sort(u) == c('0','1')) | all(sort(u) == c('-1','1')))
    else return(u == '1')
  }
}
#----------

#---- Levenshtein distance (similarity of two strings):
.LD <- function(s,t) {
  sl <- unlist(strsplit(s,''))
  tl <- unlist(strsplit(t,''))
  if (s == t) return(0)
  else if (length(sl) == 0) return(length(tl))
  else if (length(tl) == 0) return(length(sl))
  v0 <- 0:length(tl)
  v1 <- rep(NA,length(tl)+1)
  for (i in seq_along(sl)) {
    v1[1] <- i
    for (j in seq_along(tl)) {
      if (sl[i] == tl[j]) cost <- 0
      else cost <- 1
      v1[j+1] <- min(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost)
    }
    for (j in seq_along(v0)) {
      v0[j] <- v1[j]
    }
  }
  return(v1[length(tl)+1])
}
#-----

.LD2 <- function(s,t) {
  sl <- unlist(strsplit(s,''))
  tl <- unlist(strsplit(t,''))
  xx <- c(length(sl),length(tl))
  xx <- min(xx)/max(xx)
  if (s == t) return(0)
  else if (length(sl) == 0) return(length(tl))
  else if (length(tl) == 0) return(length(sl))
  v0 <- 0:length(tl)
  v1 <- rep(NA,length(tl)+1)
  for (i in seq_along(sl)) {
    v1[1] <- i
    for (j in seq_along(tl)) {
      if (sl[i] == tl[j]) cost <- 0
      else cost <- 1
      v1[j+1] <- min(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost)
    }
    for (j in seq_along(v0)) {
      v0[j] <- v1[j]
    }
  }
  return(v1[length(tl)+1]*xx)
}
#------

.agrep <- function(n,choices, r=seq(0,0.3,0.05)) {
  # r is a range can be used for max distance
  for (i in r) {
    w <- agrep(n,choices,ignore.case = TRUE,max.distance = list(all=i))
    if (length(w) > 0) break
  }
  if (length(w) > 1) {
    d <- unlist(lapply(choices[w],function(x) .LD(n,x)))
    w2 <- which(d == min(d))
    if (length(w2) == 1) choices[w][w2]
    else NA
  } else if (length(w) == 1) choices[w]
  else NA
}

.pmatch <- function(n,choices) {
  for (i in seq_along(n)) {
    if (n[i] != '') {
      if (!n[i] %in% choices) {
        u <- try(match.arg(tolower(n[i]),tolower(choices)),silent=TRUE)
        if (!inherits(u,"try-error")) {
          n[i] <- choices[which(tolower(choices) == u)]
        } else {
          n[i] <- .agrep(n[i],choices)
        }
      }
    } else n[i] <- NA
  }
  n
}

##########
# old version:
.pmatch2 <- function(n,choices) {
  for (i in seq_along(n)) {
    if (n[i] != '') {
      if (!n[i] %in% choices) {
        u <- try(match.arg(tolower(n[i]),tolower(choices)),silent=TRUE)
        if (!inherits(u,"try-error")) {
          n[i] <- choices[which(tolower(choices) == u)]
        } else {
          u <- unlist(strsplit(n[i],''))
          w1 <- which(unlist(lapply(choices,function(x) tolower(strsplit(x,'')[[1]][1]) == tolower(u[1]))))
          w2 <- unlist(lapply(choices,function(x)  length(which(tolower(u) %in% tolower(strsplit(x,'')[[1]])))/length(u)))
          w4 <- unlist(lapply(choices,function(x)  .LD2(n[i],x)))
          w3 <- which(w2 > 0.5)
          if (length(w1) > 0) {
            if (length(w3) > 0) {
              w <- w1[w1 %in% w3]
              if (length(w) > 1) {
                w <- w[which(w2[w] == max(w2))]
                if (length(w) == 1) n[i] <- choices[w]
                else if (length(w1) == 1 && w2[w1] > 0.2) n[i] <- choices[w1]
                else n[i] <- NA
              } else if (length(w) == 1) {
                n[i] <- choices[w]
              } else {
                if (length(w1) == 1 && w2[w1] > 0.2) n[i] <- choices[w1]
                #else stop(paste('no match is found for',n[i]))
                else n[i] <- NA
              }
            } else {
              if (length(w1) == 1 && w2[w1] > 0.2) n[i] <- choices[w1]
              #else stop(paste('no match is found for',n[i]))
              else n[i] <- NA
            }
          } else {
            if (length(which(w2 > 0.7)) > 0) {
              w <- which(w2 > 0.7)[which(w2 > 0.7) %in% which(w4 < 3)]
              if (length(w) == 1) n[i] <- choices[w]
              #else stop(paste('no match is found for',n[i]))
              else n[i] <- NA
            } else n[i] <- NA
          } 
        }
      }
    }
  }
  
  if ('' %in% n) {
    w <- which(n == '')
    for (i in w) if (!choices[i] %in% n) n[i] <- choices[i]
    w <- which(n == '')
    if (length(w) == 1) {
      ww <- which(!choices %in% n)
      if (length(ww) == 1) n[w] <- choices[ww]
    }
  }
  #if (length(unique(n)) < length(n)) stop('repeated arguments!')
  n
}
#-------
.canProcessInMemory <- function(x,n=1) {
  # copied partially from mem_info in the terra package!
  opt <- .eval("terra:::spatOptions()",env=environment())
  opt$ncopies = n
  v <- slot(x,slotNames(x)[1])$mem_needs(opt)
  return(round(v[5]) != 0)
}
#----
.getCells <- function(.nc,r,nr) {
  .s <- (r-1)*.nc
  c((.s+1):(.s + (nr*.nc)))
}
#---
#--------
.trim <- function(x) {
  x <- strsplit(x,'')[[1]]
  if (x[1] == ' ') {
    i <- 1
    while(x[i] == ' ') {
      i <- i+1
    }
    x <- x[i:length(x)]
  }
  #-----
  if (x[length(x)] == ' ') {
    i <- length(x)
    while(x[i] == ' ') {
      i <- i-1
    }
    x <- x[1:i]
  }
  #----
  paste(x,collapse='')
}
#-------
