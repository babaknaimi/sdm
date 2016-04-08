\name{sdm}

\alias{sdm}
\alias{sdm,formula,sdmdata,character-method}


\title{Fit and evaluate species distribution models}

\description{
Fits sdm for single or multiple species using single or multiple methods specified by a user in \code{methods} argument, and evaluates their performance.
}

\usage{
sdm(formula, data, methods,...)
}

\arguments{
  \item{formula}{Specifies the structure of the model, types of features, etc.}
  \item{data}{a \code{sdmdata} object created using \code{\link{sdmData}} function}
  \item{methods}{Character. Specifies the methods, used to fit the models}
  \item{...}{additional arguments}
  }
\details{
sdm fits multiple models and can be used to generate multiple runs (replicates) of each method through partitioning (using one or several partitioning methods including: \code{subsampling}, \code{cross-validation}, and \code{bootstrapping}.

Each model is evaluated against training data, and if available, splitted data (through partitioning; called dependent test data as well, i.e., "dep.test") and/or indipendent test data ("indep.test").

User should make sure that the methods are available and the required packages for them are installed before putting their names in the function, otherwise, the methods that cannot be run for any reason, are excluded by the function. It is a good practice to call \code{\link{installAll}} function (just one time when the sdm is installed), that tries to install all the packages that may be needed somewhere in the \code{sdm} package.

A new method can be adopted and added to the package by a user using \code{\link{add}} function. It is also possible to get an instance of an existing method, override the setting and definition, and then add it with a new name (e.g., my.glm).

The output would be a single object (\code{sdmModels}) that can be read/reproduced everywhere (e.g., on a new machine). A setting object can also be taken (exported) out of the output \code{sdmModels} object, that can be used to reproduce the same practice but given new conditions (i.e., new dataset, area. etc.)
}

\value{
 an object of class \code{sdmModels}
}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://biogeoinformatics.org}
}

\seealso{# }

\examples{
\dontrun{
file <- system.file("external/pa_df.csv", package="sdm")

df <- read.csv(file)

head(df) 

d <- sdmData(sp~b15+NDVI,train=df)

d
#----
# Example 1: fit using 3 models, and no evaluation (evaluation based on training dataset):

m <- sdm(sp~b15+NDVI,data=d,methods=c('glm','gam','gbm'))

m

# Example 3: fit using 5 models, and 
# evaluates using 10 runs of subsampling replications taking 30 percent as test:

m <- sdm(sp~b15+NDVI,data=d,methods=c('glm','gam','gbm','svm','rf'),
          replication='sub',test.percent=30,n=10)

m


# Example 3: fits using 5 models, and 
# evaluates using 10 runs of both 5-folds cross-validation and bootsrapping replication methods

m <- sdm(sp~.,data=d,methods=c('gbm','tree','mars','mda','fda'),
          replication=c('cv','boot'),cv.folds=5,n=10)

m

}

}

\keyword{spatial,species,sdm, model}
