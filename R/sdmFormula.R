# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  Feb. 2024
# Version 2.4
# Licence GPL v3

.newFormulaFunction <- function(cls,name,args,setFrame=NULL,getFeature=NULL) {
  new('.formulaFunction',cls=cls,name=name,args=args,setFrame=setFrame,getFeature=getFeature)
}

# adding the classes of formula functions into the corresponding container:
# for each function, a specific class is defined in which the name of feature, 
# and its arguments as well as the function to generate the feature from dataset is specified

.sdmFormulaFuncs <- new('.formulaFunctions')

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.auto',
                                                            representation(x='character',
                                                                           features='character',
                                                                           stat='characterORnull',
                                                                           term='call'
                                                            ))),
                                         name=c('auto','Auto'),args=c('x','features','stat')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.hinge',
                                                            representation(x='character',
                                                                           k='numericORnull',
                                                                           thresholds='numericORnull',
                                                                           feature.name='characterORnull',
                                                                           term='call'
                                                            ),
                                                            prototype(
                                                              k=20
                                                            ))),
                                         name=c('hinge','h','H','Hinge','hing','Hing'),
                                         args=c('x','k'),
                                         setFrame = function(x,param) {
                                           x@thresholds <- seq(param$min,param$max,length=x@k)
                                           x@feature.name <- c(paste0(x@x,'_hi_',x@thresholds[-length(x@thresholds)]),paste0(x@x,'_hd_',x@thresholds[-1]))
                                           x
                                         },
                                         getFeature = function(x,.var) {
                                           k <- x@thresholds
                                           h1 <- as.data.frame(lapply(k[-length(k)],function(th,x,...) {
                                             .hinge(.var,th)
                                           },x=.var))
                                           #----
                                           h2 <- as.data.frame(lapply(k[-1],function(th,x,...) {
                                             .invhinge(x,th)
                                           },x=.var))
                                           
                                           d <- cbind(h1,h2)
                                           colnames(d) <- x@feature.name
                                           d
                                         }))


.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.quad',
                                                            representation(x='character',
                                                                           feature.name='characterORnull',
                                                                           term='call'
                                                            ))),
                                         name=c('quad','q','Q','Quad','quadratic','Quadratic'),args=c('x'),
                                         getFeature = function(x) {
                                           x * x
                                         }))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.cubic',
                                                            representation(x='character',
                                                                           feature.name='characterORnull',
                                                                           term='call'
                                                            ))),
                                         name=c('cubic','c','C','Cubic'),args=c('x')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.factor',
                                                            representation(x='character',
                                                                           levels='character',
                                                                           feature.name='characterORnull',
                                                                           term='call'
                                                            ))),
                                         name=c('factor','f','F','Factor','fact','Fact'),args=c('x')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.threshold',
                                                            representation(x='character',
                                                                           k='numeric',
                                                                           thresholds='numericORnull',
                                                                           feature.name='characterORnull',
                                                                           term='call'
                                                            ),
                                                            prototype(
                                                              k=20
                                                            ))),
                                         name=c('threshold','th','Th','thereshold','thresh','Thresh'),args=c('x','k')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.poly',
                                                            representation(x='character',
                                                                           degree='numeric',
                                                                           raw='logical',
                                                                           feature.name='characterORnull',
                                                                           term='call'
                                                            ),
                                                            prototype(
                                                              degree=3,
                                                              raw=TRUE
                                                            ))),
                                         name=c('poly','Po','po','Poly'),args=c('x','degree','raw')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.interaction',
                                                            representation(x='character',
                                                                           depth='numeric',
                                                                           feature.name='characterORnull',
                                                                           term='call'
                                                            ),
                                                            prototype(
                                                              depth=1
                                                            ))),
                                         name=c('int','interaction','in','In','Int','INT'),args=c('x','depth')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.product',
                                                            representation(x='character',
                                                                           feature.name='characterORnull',
                                                                           term='call'
                                                            ))),
                                         name=c('product','p','P','Product','prod','Prod','pr'),args=c('x')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.func',
                                                            representation(x='call',
                                                                           varName='characterORnull',
                                                                           feature.name='characterORnull',
                                                                           term='call'
                                                            ))),
                                         name=c('I'),args=c('x')))



.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.simple.func',
                                                            representation(x='name',
                                                                           varName='characterORnull',
                                                                           feature.name='characterORnull',
                                                                           term='call'
                                                            ))),
                                         name=c('FUN'),args=c('x')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.log',
                                                            representation(x='character',
                                                                           feature.name='characterORnull',
                                                                           term='call'
                                                            ))),
                                         name=c('log'),args=c('x')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.log10',
                                                            representation(x='character',
                                                                           feature.name='characterORnull',
                                                                           term='call'
                                                            ))),
                                         name=c('log10'),args=c('x')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.exp',
                                                            representation(x='character',
                                                                           feature.name='characterORnull',
                                                                           term='call'
                                                            ))),
                                         name=c('exp'),args=c('x')))



.sdmFormulaFuncs$setClasses() # set the classes
# 

########################################################

#------ split a formula (f) based on the sep
.split.formula <- function(f,sep='~') {
  # based on split_formula in the package Formula by Achim Zeileis
  o <- list()
  if (length(f) > 2) {
    while (length(f) > 2 && f[[1]] == sep) {
      o <- c(f[[3]],o)
      f <- f[[2]]
    }
    c(f,o)
  } else  if (length(f) == 2) {
    o[[1]] <- f[[2]]
    o
  } else {
    o[[1]] <- f
    o
  }
}
#-------
.split.formula.for.select <- function(x) {
  x <- .split.formula(x)[[1]]
  .n <- x[[1]]
  .x <- .split.formula(x,.n)
  
  if (length(.x) > 1 && length(.x[[1]]) > 1) {
    o <-.x[-1]
    n <- rep(as.character(.n),length(.x)-1)
  } else {
    o <- .x
    n <- rep(as.character(.n),length(.x))
  }
  
  while (length(.x) > 1) {
    if (length(.x[[1]]) > 1) {
      .n <- .x[[1]][[1]]
      .x <- .split.formula(.x[[1]],.n)
      
      if (length(.x) > 1 && length(.x[[1]]) > 1) {
        o <- c(.x[-1],o)
        n <- c(rep(as.character(.n),length(.x)-1),n)
      } else {
        if (length(.x) > 1) {
          o <- c(.x,o)
          n <- c(rep(as.character(.n),length(.x)),n)
        } else {
          o <- c(.x,o)
          n <- c(as.character(.n),n)
        }
      }
      
    } else {
      break
    } 
  }
  
  if (any(n == '-')) {
    .k <- as.character(o[which(n != '-')])
    if (length(.k) == 0) .k <- as.character(o[1])
    return(list(select=.k,keep = as.character(o[which(n[-1] == '-')+1])))
  } else {
    return(list(select=as.character(o)))
  }
}
#----

.fixFormula <- function(formula) {
  if (length(formula) == 3) {
    nsp <- trim(unlist(strsplit(as.character(formula[2]),'[+|]')))
    nsp <- nsp[nsp != '']
    if ('.' %in% nsp) {
      nsp <- nsp[nsp != '.']
      if (length(nsp) == 0) formula <- as.formula(paste(formula[1],formula[3]),env = parent.frame())
      else formula <- as.formula(paste(paste(nsp,collapse='+'),formula[1],formula[3]),env = parent.frame())
    }
  }
  formula
}
#------- excludes items in y from x
.excludeVector <- function(x,y) {
  w <- unlist(lapply(y,function(z){which(x == z)}))
  if (length(w) > 0) x <- x[-w]
  if (length(x) == 0) x <- NULL
  x
}
#-------

.getRhsFromFormula <- function(f,env=parent.frame()) {
  if (length(f) == 3) as.formula(paste('~',deparse(f[[3]])),env=env)
  else if (length(f) == 2) f
  else as.formula(paste('~',deparse(f)),env=env)
}
#------

# remove coords term from formula!
.rmCoordsInFormula <- function(f) {
  
  sf <- .split.formula(f)
  lhs <- rhs <- NULL
  if (length(sf) > 2) stop('in the right hand side of the formula, the `~` can only be in m(...)')
  
  if (length(sf) == 1) rhs <- sf[[1]]
  else {
    lhs <- sf[[1]]
    rhs <- sf[[2]]
  }
  
  if (length(rhs) == 2) rhsi <- list(rhs)
  else rhsi <- .split.formula(rhs,'+')
  #-----
  temp <- unlist(lapply(rhsi,function(x) as.character(x)[[1]] %in% c('coords','coord','coordinates','coordinate','xy','geom','coods','crds','crd','cords','cord')))
  if (any(temp)) {
    #.nxy <- as.character(.split.formula(rhsi[[which(temp)]][[2]],'+'))
    if (is.null(lhs)) as.formula(paste('~',paste(rhsi[-which(temp)],collapse = '+')))
    else {
      
      if (length(lhs) > 1 && lhs[[1]] == '+') {
        
        as.formula(paste(paste(.split.formula(lhs,'+'),collapse = '+'),'~',paste(rhsi[-which(temp)],collapse = '+')),env = parent.frame())
      } else {
        as.formula(paste(as.character(lhs),'~',paste(rhsi[-which(temp)],collapse = '+')),env = parent.frame())
      }
      
    }
  } else f
  
}


#-----------

.getDataParams <- function(data,nf=NULL,nFact=NULL,id=NULL) {
  if (!is.null(id)) data <- data[id,,drop=FALSE]
  
  if (all(c(is.null(nf),is.null(nFact)))) {
    nf <- colnames(data)
    nFact <- nf[.where(is.factor,data)]
    nFact <- c(nFact,nf[.where(is.character,data)])
    
    if (length(nFact) == 0) nFact <- NULL
    
    nf <- .excludeVector(nf,nFact)
    
    nt <- nf[.where(.is.Date,data[,nf])]
    
    if (length(nt) > 0 && all(is.na(nt))) nt <- NULL
    
    nf <- .excludeVector(nf,nt)
  }
  
  if (!is.null(nf)) {
    df<- data[,nf,drop=FALSE]
    df <- data.frame(names=nf,min=round(apply(df,2,min,na.rm=TRUE),5),
                     max=round(apply(df,2,max,na.rm=TRUE),5),
                     mean=round(apply(df,2,mean,na.rm=TRUE),5),
                     sd=round(apply(df,2,sd,na.rm=TRUE),5),
                     type=apply(df,2,function(x) {
                       x <- x[!is.na(x)]
                       if (all((x - round(x,0)) == 0)) "integer" 
                       else "numeric"
                     }),
                     unique_count=apply(df,2,function(x) length(unique(x))))
  }
  
  if (!is.null(nFact)) {
    dFact <- list()
    for (n in nFact) {
      .tab <- table(data[,n])
      dFact[[n]] <- data.frame(levels=names(.tab),count=as.vector(.tab))
    }
  }
  
  if (any(c(is.null(nf),is.null(nFact)))) {
    if (is.null(nf)) dFact
    else df
  } else {
    list(continuous = df, categorical = dFact)
  }
  
}
#------


#################---- detect the terms in the nested formula (model) inside the main formula:
# .nested_terms <- function(x) {
#   .nm <- new('.nestedModel',term=x)
#   #--------
#   a <- c('formula','method','setting','output')
#   s <- list()
#   n <- names(x[-1])
#   if (length(n) > 0) {
#     for (i in seq_along(n)) {
#       if (n[i] != '') {
#         .w <- .pmatch(tolower(n[i]),a)
#         if (!is.na(.w)) n[i] <- .w
#       }
#     }
#   } else n <- rep('',length(x[-1]))
#   
#   if (length(n[n != '']) > 0 && !all(n[n != ''] %in% a)) stop('some arguments in the formula:model (m) function is unknown!')
#   if (length(x) > 5) stop('the arguments in the formula:model (m) function are not match!')
#   for (i in 1:length(n)) {
#     if (n[i] != '') s[[n[i]]] <- x[[i+1]]
#     else s[[a[i]]] <- x[[i+1]]
#   }
#   x <- s[['formula']]
#   
#   if (length(x) > 1) {
#     if (x[[1]] == '~') {
#       if (length(x) == 3) {
#         x <- .fixFormula(x)
#         l <- .split.formula(x[[3]],'+')
#         if (length(.split.formula(x[[2]],'+')) > 1) stop('nested formula in the rhs, cannot be multi-response!')
#         .nm@response <- as.character(x[[2]])
#       } else l <- .split.formula(x[[2]],'+')
#     } else l <- .split.formula(x,'+')
#   } else l <- list(x)
#   #----------
#   
#   .nm@terms <- lapply(l,.term)
#   x <- s[['method']]
#   if (!is.null(x)) {
#     w <- unlist(.sdmMethods$getMethodNames(alt = TRUE))
#     names(w) <- NULL
#     
#     if (length(x) > 1) {
#       if (x[[1]] == 'c') {
#         x <- tolower(as.character(x[-1]))
#         if (!all(x %in% w)) {
#           if (!any(x %in% w)) .nm@method <- NULL
#           else .nm@method <- x[x %in% w]
#         } else .nm@method <- x
#       }
#     } else {
#       x <- tolower(as.character(x))
#       if (!x %in% w) .nm@method <- NULL
#       else .nm@method <- x
#     }
#   }
#   
#   #----
#   x <- s[['setting']]
#   if (!is.null(x)) {
#     if (length(x) > 1) {
#       if (x[[1]] == 'c') x[1] <- call('list')
#       
#       .nm@setting <- eval(x)
#     }
#   }
#   #-------
#   
#   x <- s[['output']]
#   if (!is.null(x)) {
#     if (length(x) == 1) {
#       
#       x <- tolower(as.character(x))
#       x <- .pmatch(x,c('prediction','residual'))
#       if (!is.na(x)) .nm@output <- x
#     } else .nm@output <- 'prediction'
#   } else .nm@output <- 'prediction'
#   #-------
#   .nm
# }
#----
.nested_terms <- function(x) {
  .nm <- new('.nestedModel',term=x)
  #--------
  a <- c('formula','method','setting','output')
  s <- list()
  n <- names(x[-1])
  if (length(n) > 0) {
    for (i in seq_along(n)) {
      if (n[i] != '') {
        .w <- .pmatch(tolower(n[i]),a)
        if (!is.na(.w)) n[i] <- .w
      }
    }
  } else n <- rep('',length(x[-1]))
  
  if (length(n[n != '']) > 0 && !all(n[n != ''] %in% a)) stop('some arguments in the formula:model (m) function is unknown!')
  if (length(x) > 5) stop('the arguments in the formula:model (m) function are not match!')
  for (i in 1:length(n)) {
    if (n[i] != '') s[[n[i]]] <- x[[i+1]]
    else s[[a[i]]] <- x[[i+1]]
  }
  x <- s[['formula']]
  
  if (length(x) > 1) {
    if (x[[1]] == '~') {
      if (length(x) == 3) {
        x <- .fixFormula(x)
        l <- all.vars(x[[3]])
        #l <- .split.formula(x[[3]],'+')
        if (length(.split.formula(x[[2]],'+')) > 1) stop('nested formula in the rhs, cannot be multi-response!')
        .nm@response <- as.character(x[[2]])
      } else l <- all.vars(x[[2]])
    } else l <- all.vars(x)
  } else l <- list(x)
  #----------
  
  .nm@predictors <- l
  x <- s[['method']]
  if (!is.null(x)) {
    w <- unlist(.sdmMethods$getMethodNames(alt = TRUE))
    names(w) <- NULL
    
    if (length(x) > 1) {
      if (x[[1]] == 'c') {
        x <- tolower(as.character(x[-1]))
        if (!all(x %in% w)) {
          if (!any(x %in% w)) .nm@method <- NULL
          else .nm@method <- x[x %in% w]
        } else .nm@method <- x
      }
    } else {
      x <- tolower(as.character(x))
      if (!x %in% w) .nm@method <- NULL
      else .nm@method <- x
    }
  }
  
  #----
  x <- s[['setting']]
  if (!is.null(x)) {
    if (length(x) > 1) {
      if (x[[1]] == 'c') x[1] <- call('list')
      
      .nm@setting <- eval(x)
    }
  }
  #-------
  
  x <- s[['output']]
  if (!is.null(x)) {
    if (length(x) == 1) {
      
      x <- tolower(as.character(x))
      x <- .pmatch(x,c('prediction','residual'))
      if (!is.na(x)) .nm@output <- x
    } else .nm@output <- 'prediction'
  } else .nm@output <- 'prediction'
  #-------
  .nm
}
#----------
#----------
.exPCA <- function(x) {
  a <- c('formula','n')
  #-----
  s <- list()
  if (length(x) == 1) {
    .pc <- new('.pcaSetting',vars='.',n=3,term=x)
  } else {
    .n <- names(x)
    if (length(.n) > 0) {
      for (i in seq_along(.n)) {
        if (.n[i] != '') {
          .w <- .pmatch(tolower(.n[i]),a)
          if (!is.na(.w)) .n[i] <- .w
        }
      }
    } else .n <- rep('',length(x))
    
    if (length(.n[.n != '']) > 0 && !all(.n[.n != ''] %in% a)) stop('some arguments in pca function is unknown!')
    
    if (length(x) > 3) stop('the arguments in the pca function are not match!')
    
    for (i in 2:length(.n)) {
      if (.n[i] != '') s[[.n[i]]] <- x[[i]]
      else s[[a[i-1]]] <- x[[i]]
    }
    
    if (length(s[['formula']]) > 1) {
      if (s[['formula']][[1]] != '+') stop('something in `pca` is wrong; pca(var1+var2+var3,n=3); pca(.,n="90%"); pca(.,n="auto"); pca(.,n=3)')
      else {
        l <- .split.formula(s[['formula']],'+')
        if (any(unlist(lapply(l,function(x) {length(x) > 1})))) stop('something wrong with pca; example: pca(var1+var2+var3,n=3); pca(.,n="90%"); pca(.,n="auto"); pca(.,n=3)')
      }
    } else l <- list(s[['formula']])
    
    .pc <- new('.pcaSetting',term=x)
    
    if (!is.null(s[['n']])) .pc@n <- s[['n']]
    else .pc@n <- 3
    
    if (any(as.character(l) %in% c("NULL","."))) .pc@vars <- '.'
    else .pc@vars <- as.character(l)
  }
  .pc
}
#--------
.exScale <- function(x) {
  a <- c('formula','method')
  s <- list()
  if (length(x) == 1) {
    n <- new('.scaleSetting',vars='.',method="minmax")
  } else {
    n <- names(x)
    if (length(n) > 0) {
      for (i in seq_along(n)) {
        if (n[i] != '') {
          .w <- .pmatch(tolower(n[i]),a)
          if (!is.na(.w)) n[i] <- .w
        }
      }
    } else n <- rep('',length(x))
    
    if (length(n[n != '']) > 0 && !all(n[n != ''] %in% a)) stop('some arguments in the scale function is unknown!')
    if (length(x) > 3) stop('the arguments in the scale function are not correct!')
    for (i in 2:length(n)) {
      if (n[i] != '') s[[n[i]]] <- x[[i]]
      else s[[a[i-1]]] <- x[[i]]
    }
    
    if (length(s[['formula']]) > 1) {
      if (s[['formula']][[1]] != '+') stop('something in `scale` is wrong; example: scale(var1+var2+var3,method="minmax")...')
      else {
        l <- .split.formula(s[['formula']],'+')
        if (any(unlist(lapply(l,function(x) {length(x) > 1})))) stop('something is wrong in the scale function in formula; example: scale(var1+var2+var3,method="minmax")')
      }
    } else l <- list(s[['formula']])
    
    n <- new('.scaleSetting',term=x)
    
    if (!is.null(s[['method']])) n@method <- s[['method']]
    else n@method <- 'minmax'
    if (any(as.character(l) %in% c("NULL","."))) n@vars <- '.'
    else n@vars <- as.character(l)
  }
  n
}
#----

.exInteraction <- function(x) {
  a <- .sdmFormulaFuncs$getFuncs('int')$int@args
  s <- list()
  
  n <- names(x)
  if (length(n) > 0) {
    for (i in seq_along(n)) {
      if (n[i] != '') {
        .w <- .pmatch(tolower(n[i]),a)
        if (!is.na(.w)) n[i] <- .w
      }
    }
  } else n <- rep('',length(x))
  
  if (length(n[n != '']) > 0 && !all(n[n != ''] %in% a)) stop('Some arguments in the interaction function in the formula is unknown!')
  if (length(x) > 3) stop('Arguments in the interaction function in the formula are not correct!')
  for (i in 2:length(n)) {
    if (n[i] != '') s[[n[i]]] <- x[[i]]
    else s[[a[i-1]]] <- x[[i]]
  }
  
  if (length(s[['x']]) > 1) {
    if (s[['x']][[1]] != '+') stop('something in `interaction` is wrong: Examples: int(., depth=2); int(x1 + x2 + x3, depth=3) !')
    else {
      l <- .split.formula(s[['x']],'+')
      if (any(unlist(lapply(l,function(x) {length(x) > 1})))) stop('something wrong with scale; example: scale(var1+var2+var3,method="minmax")')
    }
  } else l <- list(s[['x']])
  
  n <- new('.interaction',term=x)
  
  if (!is.null(s[['depth']])) {
    if (is.numeric(s[['depth']]) && s[['depth']] > 0) n@depth <- ceiling(s[['depth']])
  } else n@depth <- 1
  
  if (any(as.character(l) %in% c("NULL","."))) n@x <- '.'
  else n@x <- as.character(l)
  n
}

#--------- detect the class of the term in the formula
.term <- function(x) {
  if (length(x) == 1) {
    if (x == '.') return(new('.all.vars',names=as.character(x)))
    else return(new('.var',name=as.character(x)))
  } else if (length(x) > 1) {
    if (as.character(x[[1]]) %in% c('m','M','mo','model','mod','MODEL','Model')) {
      .nested_terms(x)
    } else if (as.character(x[[1]]) %in% c('select','se','sel','SEL','SELECT','SE')) {
      .select.terms(x)
    } else if (as.character(x[[1]]) %in% c('coords','coord','coordinates','coordinate','xy','geom','coods','crds','crd','cords','cord') || any(!is.na(pmatch(c("co"),tolower(as.character(x[[1]])))))) {
      .exCoords(x)
    } else if (as.character(x[[1]]) %in% c('g','G','gr','group','Group','GROUP','GR','gro','grop','grp')) {
      .exGroup(x)
    } else if (as.character(x[[1]]) %in% c('time','Time','tim','TIME','Tim')) {
      .exTime(x)
    } else if (as.character(x[[1]]) %in% c('info','Info','inf','INFO')) {
      .exInfo(x)
    } else if (as.character(x[[1]]) %in% c('*',':','product','p','P','Product','prod','Prod','pro','Pro')) {
      .exProduct(x)
    } else if (as.character(x[[1]]) %in% c('scale','SCALE','sc')) {
      .exScale(x)
    } else if (as.character(x[[1]]) %in% c('int','INT','interaction','Int','in','In')) {
      .exInteraction(x)
    } else if (as.character(x[[1]]) %in% c('pca','PCA','prcomp','princomp')) {
      .exPCA(x)
    } else .exFunc(x)
  }
}
#-------
.exGroup <- function(x) {
  new('.grouping',group.var=as.character(x[[2]]),term=x)
}
#------
.exTime <- function(x) {
  xx <- as.character(x[[1]])
  s <- list()
  n <- names(x)
  if (!is.null(n)) {
    for (i in 2:length(n)) {
      if (n[i] != '') s[[n[i]]] <- x[[i]]
      else s[[(i-1)]] <- x[[i]]
    }
  } else {
    for (i in 2:length(x)) s[[(i-1)]] <- x[[i]]
  }
  new('.time',name=xx,terms=s,term=x)
}
#--------
.exInfo <- function(x) {
  n <- NULL
  if (length(x[[2]]) > 1 && as.character(x[[2]][[1]]) %in% c('|','+')) {
    if (as.character(x[[2]][[1]]) == '+') n <- as.character(.split.formula(x[[2]],'+'))
    else n <- as.character(.split.formula(x[[2]],'|'))
  } else {
    n <- as.character(x[[2]])
  } 
  if (!is.null(n)) new('.Info',names=n)
}


#---------
.exFunc <- function(x) {
  xx <- as.character(x[[1]])
  mn <- .sdmFormulaFuncs$funcNames
  names(mn) <- mn
  if (xx %in% mn) xx <- names(mn)[mn == xx]
  else {
    mnlist <- lapply(mn,function(x) .sdmFormulaFuncs$funcs[[x]]@name)
    u <- unlist(lapply(mnlist,function(x) xx %in% x))
    if (any(u)) xx <- names(u)[which(u)]
    else xx <- NULL
  }
  if (!is.null(xx)) {
    ss <- .sdmFormulaFuncs$getFuncs(xx)
    a <- ss[[xx]]@args
    cls <- new(ss[[xx]]@cls[[2]])
    n <- names(x)
    n <- n[2:length(n)]
    
    if (!is.null(n)) {
      if (length(n[n != ''] > 0) && !all(n[n != ''] %in% a)) stop(paste0('some arguments in function ',xx, ' is unknown!'))
      for (i in 1:length(n)) {
        if (n[i] != '') {
          if (inherits(x[[i+1]],'name')) {
            slot(cls,n[i]) <- as.character(x[[i+1]])
          } else slot(cls,n[i]) <- x[[i+1]]
        } else {
          if (inherits(x[[i+1]],'name')) slot(cls,a[i]) <- as.character(x[[i+1]])
          else slot(cls,a[i]) <- x[[i+1]]
        }
      }
    } else {
      for (i in 2:length(x)) {
        if (inherits(x[[i]],'name')) slot(cls,a[i-1]) <- as.character(x[[i]])
        else slot(cls,a[i-1]) <- x[[i]]
      }
    }
  } else {
    if (exists(as.character(x[[1]]),mode='function')) {
      cls <- new('.simple.func')
      cls@x <- x[[2]]
      cls@varName <- all.vars(x)
    } else stop(paste(as.character(x[[1]]),'is not a known function!'))
  }
  
  cls@term <- x
  cls
}
#--------
.exProduct <- function(x) {
  cls <- new('.product')
  if (as.character(x[[1]]) %in% c('*',':')) {
    cls@x <- unique(as.character(.split.formula(x,as.character(x[[1]]))))
  } else {
    if (length(x) == 2) {
      xx <- x[[2]]
      cls@x <- unique(as.character(.split.formula(xx,as.character(xx[[1]]))))
    } else {
      xx <- c()
      for (i in 2:length(x)) xx <- c(xx,as.character(x[[i]]))
      cls@x <- unique(xx)
    }
  }
  
  cls@term <- x
  cls
}
#--------
.exCoords <- function(x) {
  if (length(x[[2]]) > 1 && as.character(x[[2]][[1]]) %in% c('|','+')) {
    xy <- as.character(x[[2]])[2:3]
  } else if (length(x) == 3) {
    xy <- c(as.character(x[[2]]),as.character(x[[3]]))
  } else stop('in formula, coordinates are not properly defined; Example: ~...+coords(x+y)+...')
  new('.coord.vars',xy=xy)
}


#-------------

.select.terms <- function(x) {
  a <- c('formula','method','th')
  s <- list()
  n <- names(x)
  if (length(n) > 0) {
    for (i in seq_along(n)) {
      if (n[i] != '') {
        .w <- .pmatch(tolower(n[i]),a)
        if (!is.na(.w)) n[i] <- .w
      }
    }
  } else n <- rep('',length(x))
  
  if (length(n[n != '']) > 0 && !all(n[n != ''] %in% a)) stop('some arguments in select function is unknown!')
  if (length(x) > 5) stop('the arguments in select function are not match!')
  for (i in 2:length(n)) {
    if (n[i] != '') s[[n[i]]] <- x[[i]]
    else s[[a[i-1]]] <- x[[i]]
  }
  
  if (length(s[['formula']]) > 1) {
    l <- .split.formula.for.select(s[['formula']])
    
  } else l <- list(select=as.character(s[['formula']]))
  
  n <- new('.selectFrame',term=x)
  
  if (!is.null(s[['th']])) n@th <- s[['th']]
  if (!is.null(s[['method']])) n@method <- s[['method']]
  
  n@vars <- l$select
  n@keep <- l$keep
  n
}
#----
#--------
.trim <- function(x) {
  x <- strsplit(x,'')[[1]]
  paste(x[x != ' '],collapse='')
}
#-------

########

# .exFormula extract terms in formula and detect what each term is. it may be a model.term (including a
# variable, a function, a nested model, etc.) or a data.term (including coordinates, select function, group, etc.)
.exFormula <- function(f,data,detect=TRUE) {
  
  f <- .fixFormula(f)
  v <- colnames(data)
  
  nFact <- nf <- ng <- ni <- nt <- nxy <- nsp <- NULL
  
  nall <- n <- all.vars(f)
  
  if ('.' %in% n) {
    n <- n[-which(n == '.')]
    nall <- unique(c(v,n))
  }
  
  nFact <- v[.where(is.factor,data)]
  nFact <- c(nFact,v[.where(is.character,data)])
  
  if (length(nFact) == 0) nFact <- NULL
  else {
    if (any(nFact %in% nall)) nFact <- nFact[nFact %in% nall]
    else nFact <- NULL
  }
  
  sf <- .split.formula(f)
  lhs <- rhs <- NULL
  if (length(sf) > 2) stop('in the right hand side of the formula, the `~` can only be in m(...)')
  f <- new('sdmFormula',formula=f)
  if (length(sf) == 1) rhs <- sf[[1]]
  else {
    lhs <- sf[[1]]
    rhs <- sf[[2]]
  }
  
  if (!is.null(lhs)) {
    lhs <- .split.formula(lhs,'+')
    nsp <- as.character(lhs)
    nall <- .excludeVector(nall,nsp)
    n <- .excludeVector(n,nsp)
    nFact <- .excludeVector(nFact,nsp)
  } else {
    if (detect) {
      w <- which(unlist(lapply(data,.isBinomial)))
      if (length(w) > 0) {
        nsp <- v[w]
        nsp <- nsp[!nsp %in% n]
        nall <- .excludeVector(nall,nsp)
        nFact <- .excludeVector(nFact,nsp)
        lhs <- as.list(nsp)
      } else nsp <- NULL
    } else nsp <- NULL
    
  }
  #f@vars <- nall
  
  
  #f@species <- as.character(lhs)
  
  if (length(rhs) == 2) rhsi <- list(rhs)
  else rhsi <- .split.formula(rhs,'+')
  
  
  temp <- unlist(lapply(rhsi,function(x) as.character(x)[[1]] %in% c('coords','coord','coordinates','coordinate','xy','geom','coods','crds','crd','cords','cord')))
  if (any(temp)) nxy <- as.character(.split.formula(rhsi[[which(temp)]][[2]],'+'))
  
  vars <- .excludeVector(nall,c(n,nxy,nFact))
  
  w <- unlist(lapply(rhsi,function(x) x == '.'))
  if (any(w)) {
    if (!is.null(vars)) rhsi <- c(rhsi[!w],.split.formula(as.formula(paste('~',paste(vars,collapse='+')))[[2]],'+')) 
    else rhsi <- rhsi[!w]
  }
  
  if (!is.null(nFact)) {
    for (i in seq_along(rhsi)) {
      if (length(as.character(rhsi[[i]])) == 1) {
        if (any(nFact == as.character(rhsi[[i]]))) {
          rhsi[[i]] <- as.formula(paste('~f(',nFact[nFact == as.character(rhsi[[i]])],')'))[[2]]
        }
      }
    }
  }
  
  nf <- .excludeVector(nall,c(nxy,nFact))
  
  func.cls <- unlist(lapply(.sdmFormulaFuncs$funcNames,function(x) .sdmFormulaFuncs$funcs[[x]]@cls[[2]]))
  temp <- lapply(rhsi,.term)
  .fc <- unlist(lapply(temp,class)) # classes of the items in formula
  #---
  if ('.selectFrame' %in% .fc) {
    if (!is.null(temp[[which(.fc == '.selectFrame')]]@keep)) {
      w <- FALSE
      for (.k in temp[[which(.fc == '.selectFrame')]]@keep) {
        wt <- sapply(rhsi,function(x) x == as.name(.k))
        if (!any(wt)) {
          rhsi <- c(rhsi,as.name(.k))
          w <- TRUE
        }
      }
      if (w) {
        temp <- lapply(rhsi,.term)
        .fc <- unlist(lapply(temp,class)) # classes of the items in formula
      }
    }
  }
  #----
  wt <- which(.fc %in% c('.var','.nestedModel',func.cls))
  if (length(wt) > 0) f@model.terms <- temp[wt]
  wt <- which(.fc %in% c('.coord.vars','.grouping','.Info','.time','.scaleSetting','.selectFrame','.pcaSetting'))
  
  if (length(wt) > 0) {
    f@data.terms <- c(f@data.terms,temp[wt])
    w <- unlist(lapply(f@data.terms,class))
    
    if (".grouping" %in% w) {
      wt <- f@data.terms[which(w == ".grouping")]
      ng <- sapply(wt,function(x) x@group.var)
      nf <- .excludeVector(nf,ng)
      nFact <- .excludeVector(nFact,ng)
    }
    #---
    if ('.Info' %in% w) {
      wt <- f@data.terms[which(w == ".Info")]
      ni <- sapply(wt,function(x) as.character(x@names))
      nf <- .excludeVector(nf,ni)
      nFact <- .excludeVector(nFact,ni)
    }
    #----
    if ('.time' %in% w) {
      wt <- f@data.terms[which(w == ".time")]
      nt <- sapply(wt,function(x) as.character(x@terms[1]))
      nf <- .excludeVector(nf,nt)
      nFact <- .excludeVector(nFact,nt)
    } else {
      if (!is.null(nf)) nt <- nf[.where(.is.Date,data[,nf])]
      if (length(nt) > 0) {
        nf <- .excludeVector(nf,nt)
        nt <- nt[1]
        for (i in 1:length(f@model.terms)) {
          if (inherits(f@model.terms[[i]], '.var') && f@model.terms[[i]]@name == nt ) {
            f@model.terms <- f@model.terms[-i]
            f@data.terms <- c(f@data.terms,.exTime(.split.formula(as.formula(paste('~time(',nt,')')))[[1]]))
            break
          }
        }
      } else nt <- NULL
    }
    #----
  } else {
    if (!is.null(nf)) nt <- nf[.where(.is.Date,data[,nf])]
    if (length(nt) > 0) {
      nf <- .excludeVector(nf,nt)
      nt <- nt[1]
      for (i in 1:length(f@model.terms)) {
        if (inherits(f@model.terms[[i]], '.var') && f@model.terms[[i]]@name == nt ) {
          f@model.terms <- f@model.terms[-i]
          f@data.terms <- c(f@data.terms,.exTime(.split.formula(as.formula(paste('~time(',nt,')')))[[1]]))
          break
        }
      }
    } else nt <- NULL
  }
  #-----
  if (!is.null(f@model.terms)) {
    w <- unlist(lapply(f@model.terms,class))
    if ('.factor' %in% w) {
      wi <- which(w == '.factor')
      for (i in wi) {
        w <- as.character(f@model.terms[[i]]@x)
        w <- .excludeVector(w,'+')
        nFact <- unique(c(nFact,w))
        if (is.factor(data[,w])) f@model.terms[[i]]@levels <- levels(data[,w])
        else f@model.terms[[i]]@levels <- sort(unique(as.character(data[,w])))
      }
    }
  }
  
  nf <- .excludeVector(nf,nFact)
  
  if (!is.null(nf)) w1 <- .getDataParams(data,nf=nf)
  else w1 <- NULL
  if (!is.null(nFact)) w2 <- .getDataParams(data,nFact=nFact)
  else w2 <- NULL
  
  if (any(!is.null(c(nxy,ng,ni,nt)))) {
    #w3 <- data.frame(name='xxx',type='xxx',details='xxx')[0,]
    w3 <- data.frame(name='xxx',type='xxx')[0,]
    if (!is.null(nxy)) {
      #w3 <- rbind(w3,data.frame(name=nxy,type=c('x_coordinate','y_coordinate'),details=c(paste('mean:',round(mean(data[,nxy[1]],na.rm=T),2)),paste('mean:',round(mean(data[,nxy[2]],na.rm=T),2)))))
      w3 <- rbind(w3,data.frame(name=nxy,type=c('x_coordinate','y_coordinate')))
    }
    
    if (!is.null(nt)) {
      #w3 <- rbind(w3,data.frame(name=nt,type='Date/Time',details=paste0('range: ',paste(range(data[,nt]),collapse=':'))))
      w3 <- rbind(w3,data.frame(name=nt,type='Date/Time'))
    }
    
    if (!is.null(ng)) {
      w3 <- rbind(w3,data.frame(name=ng,type=rep('group',length(ng))))
    }
    
    if (!is.null(ni)) {
      w3 <- rbind(w3,data.frame(name=ni,type=rep('Info',length(ni))))
    }
  } else w3 <- NULL
  
  f@vars <- new('.variables',names=c(nsp,nf,nFact,nxy,nt,ng,ni),
                species=nsp,numeric=w1,categorical=w2,others=w3)
  if (!is.null(w2)) {
    nFact <- names(w2)
    w <- which(sapply(f@model.terms,class) == '.factor')
    if (length(w) > 0) {
      w <- which(!nFact %in% sapply(f@model.terms[w],function(x) x@x))
      if (length(w) > 0) {
        for (i in w) {
          f@model.terms <- c(f@model.terms,
                             new('.factor',x=nFact[i],levels=w2[[nFact[i]]]$levels))
        }
      }
    } else {
      for (i in 1:length(nFact)) {
        n <- nFact[i]
        f@model.terms <- c(f@model.terms,
                           new('.factor',x=n,levels=w2[[n]]$levels,term=.eval(paste0("call('f',substitute(",n,'))'),environment())))
      }
    }
  }
  
  f
}
#-----------
###################

