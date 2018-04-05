\name{sdmSetting}

\alias{sdmSetting}
\alias{sdmSetting,ANY,sdmdata,character-method}
\alias{show,sdmSetting-method}

\title{creating sdmSetting object}

\description{
Creates sdmSetting object that holds settings to fit and evaluate the models. It can be used to reproduce a study.
}

\usage{
sdmSetting(formula,data,methods,interaction.depth=1,n=1,replication=NULL,cv.folds=NULL,
     test.percent=NULL,bg=NULL,bg.n=NULL,var.importance=NULL,response.curve=TRUE,
     var.selection=FALSE,ncore=1L,modelSettings=NULL,seed=NULL,parallelSettings=NULL,...)


}

\arguments{
  \item{formula}{specify the structure of the model}
  \item{data}{sdm data object or data.frame including species and feature data}
  \item{methods}{character, name of the algorithms}
  \item{interaction.depth}{level of interactions between predictors}
  \item{n}{number of replicates (run)}
  \item{replication}{replication method (e.g., 'subsampling', 'bootstrapping', 'cv')}
  \item{cv.folds}{number of folds if cv (cross-validation) is in the selected replication methods}
  \item{test.percent}{test percentage if subsampling is in the selected replication methods}
  \item{bg}{method to generate background}
  \item{bg.n}{number of background records}
  \item{var.importance}{logical, whether variable importance should be calculated}
  \item{response.curve}{method to calculate variable importance}
  \item{var.selection}{logical, whether variable selection should be considered}
  \item{ncore}{number of cores to parallelize processing}
  \item{modelSettings}{optional list; settings for modelling methods can be specified by users}
  \item{seed}{default is NULL; either logical specify whether a seed for random number generator should be considered, or a numerical to specify the exact seed number}
  \item{parallelSettings}{default is NULL; a list include settings items for parallel processing. The parallel setting items include ncore, method, type, hosts, doParallel, and fork; see details for more information.}
  \item{...}{additional arguments}
  }
\details{
using sdmSetting, the feature types, interaction.depth and all settings of the model can be defined. This function generate a sdmSetting object that can be specifically helpful for reproducibility. The object can be shared by a user that may be used for other studies.

If a user aims to reproduce the same results for every time the code is running with the same data and settings, a seed number should be specified. Through the \code{seed} argument, a user can specify \code{NULL}, means a seed should not be set (if a random sampling is incorporated in the modelling procedure, for different runs the results would be different); \code{TRUE}, means a seed should be set (the seed number is randomly selected and used everytime the same setting is incorporated); \code{a number}, means the seed will be set to the number specified by the user. 

For parallel processing, a list of items can be passed to \code{parallelSettings}, include:

\code{ncore}: defines the number of cores (it can also be specified outside of this list, but will be removed in future)

\code{method}: defines the platform/set of functions to run the parallelisation. Currently, two options of 'parallel', and 'foreach' is implemented. default is 'parallel'

\code{doParallel}: Optional, definition to register for a backend for parallel processing (currently when method='foreach'). It should be provided as an R expression.

\code{cluster}: Optional, if a cluster is already created and started, it can be introduced through this item to be used as the parallel processing platform (currently when method='parallel')

\code{hosts}: A list of addresses for the accessible hosts (remote clusters) to be registered and used in parallel processing (may not work appropriately as it is still under development!)

\code{fork}: Logical, Available for non-windows operating system and specifies whether a fork solution should be used for the parallelisation. Default is TRUE.


}

\value{
 an object of class \code{sdmSettings}
}


\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}

\url{http://biogeoinformatics.org}
}

\examples{
\dontrun{
file <- system.file("external/pa_df.csv", package="sdm")

df <- read.csv(file)

head(df) 

d <- sdmData(sp~b15+NDVI,train=df)

# generate sdmSettings object:
s <- sdmSetting(sp~., methods=c('glm','gam','brt','svm','rf'),
          replication='sub',test.percent=30,n=10,modelSettings=list(brt=list(n.trees=500)))

s


}
}

\keyword{spatial,species,data}
