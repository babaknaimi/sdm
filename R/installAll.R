# Author: Babak Naimi, naimi.b@gmail.com
# Last Update :  Jan 2024
# Version 2.0
# Licence GPL v3
#--------

.eval <- function(x,env) {
  eval(parse(text=x),envir=env)
}
#------
.is.installed <- function(n) {
  names(n) <- n
  sapply(n, function(x) length(unlist(lapply(.libPaths(), function(lib) find.package(x, lib, quiet=TRUE, verbose=FALSE)))) > 0)
}
#---------

.require <-function(x) {
  x <- as.character(x)
  xx <- unlist(lapply(.libPaths(), function(lib) find.package(x, lib, quiet=TRUE, verbose=FALSE)))
  if (length(xx) > 0) {
    .loaded <- eval(parse(text=paste0('require(',x,')')))
    return (.loaded)
  } else FALSE
}
#----------
.loadLib <- function(pkgs) {
  options(warn=-1)
  return(unlist(lapply(pkgs,function(x) {
    all(unlist(lapply(x,function(p) {.require(p)})))
  })))
  options(warn=0)
}
#---------

.getPackageList <- function() {
  methodInfo <- NULL
  pkgs <- c()
  lst <- list.files(system.file("methods/sdm", package="sdm"),pattern='R$',full.names = TRUE)
  for (l in lst) {
    source(l,local=TRUE)
    p <- methodInfo$packages
    p <- p[!p == '.tmp']
    pkgs <- c(pkgs,p)
  }
  p <- c('shiny','raster','shinyBS','leaflet','usdm','devtools','mmap','ggplot2','gridExtra')
  unique(c(pkgs,p))
}
#------------- 
# List of the packages that should be installed from GitHUB:
.getPackageGitHubList <- function() {
  c('babaknaimi/mraster')
}
#--------------

if (!isGeneric("installAll")) {
  setGeneric("installAll", function(pkgs,update,...)
    standardGeneric("installAll"))
}


setMethod('installAll', signature(pkgs='ANY'),
          function(pkgs,update=FALSE,...) {
            if (missing(update)) update <- FALSE
            pl <- .getPackageList()
            plG <- .getPackageGitHubList()
            plGr <- sapply(strsplit(plG,'/'),function(x) x[2]) # name of the github packages
            
            if (!update) {
              p <- pl[!.is.installed(pl)]
              pG <- plG[!.is.installed(plGr)]
              if (length(c(p,pG)) > 0) {
                s <- rep(TRUE,length(c(p,pG)))
                if (length(p) > 0) {
                  #s <- rep(TRUE,length(p))
                  for (i in seq_along(p)) {
                    pi <- try(install.packages(p[i],...),silent = TRUE)
                    if (inherits(pi, "try-error")) s[i] <- FALSE
                  }
                }
                #---
                if (length(pG) > 0) {
                  s[(length(p)+1):length(c(p,pG))] <- .eval("devtools::install_github(pG,quiet=TRUE)",env=environment())
                }
                
                if (any(!s)) {
                  if (any(s)) {
                    cat(paste('\n',length(c(p,pG)[s]),' packages are successfully installed...\n'))
                    cat(paste('The following packages could not be installed:\n.... ',paste(c(p,pG)[!s],collapse=', '),'\n'))
                  } 
                } else cat(paste('\n ',length(c(p,pG)[s]),' packages are successfully installed...\n'))
              } else cat(paste('\n All required packages have been already installed!\n'))
              
            } else {
              p <- pl[!pl %in% c('stats','utils','parallel','base','grDevice','tools','methods','graphics','compiler','datasets','profile','grid')]
              pG <- plG[!plGr %in% c('stats','utils','parallel','base','grDevice','tools','methods','graphics','compiler','datasets','profile','grid')]
              if (length(c(p,pG)) > 0) {
                s <- rep(TRUE,length(c(p,pG)))
                if (length(p) > 0) {
                  .detachPackage(p)
                  pi <- p[.is.installed(p)]
                  if (length(pi) > 0) pi <- try(remove.packages(pi),silent = TRUE)
                  
                  
                  for (i in seq_along(p)) {
                    pi <- try(install.packages(p[i],...),silent = TRUE)
                    if (inherits(pi, "try-error")) s[i] <- FALSE
                  }
                }
                
                if (length(pG) > 0) {
                  pGr <- sapply(strsplit(pG,'/'),function(x) x[2])
                  .detachPackage(pGr)
                  pGi <- pGr[.is.installed(pGr)]
                  if (length(pGi) > 0) {
                    pGi <- try(remove.packages(pGi),silent = TRUE)
                    .xx <- try(.eval("devtools::install_github(pG,quiet=TRUE)"),silent = TRUE)
                    if (inherits(.xx,'try-error')) .xx <- FALSE
                    else .xx <- TRUE
                    s[(length(p)+1):length(c(p,pG))] <- .xx
                  }
                }
                if (any(!s)) {
                  if (any(s)) {
                    cat(paste('\n',length(c(p,pG)[s]),' packages are successfully installed or updated...\n'))
                    cat(paste('The following packages could not be installed:\n.... ',paste(c(p,pG)[!s],collapse=', '),'\n'))
                  }
                } else cat(paste('\n ',length(c(p,pG)[s]),' packages are successfully installed or updated...\n'))
              } else cat(paste('\n There is no package to install!\n'))
            }
            .addMethods()
          }
)

