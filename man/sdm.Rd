\name{sdm}

\alias{sdm}
\alias{sdm,ANY,sdmdata,character-method}
\alias{sdm,sdmdata,.sdmCorSetting,ANY-method}
\alias{sdm,ANY,sdmdata,.sdmCorSetting-method}


\title{Fit and evaluate species distribution models}

\description{
Fits sdm for single or multiple species using single or multiple methods specified by a user in \code{methods} argument, and evaluates their performance.
}

\usage{
sdm(formula, data, methods,...)
}

\arguments{
  \item{formula}{specifies the structure of the model, types of features, etc.}
  \item{data}{a \code{sdmdata} object created using \code{\link{sdmData}} function}
  \item{methods}{a character, specifies the methods, used to fit the models}
  \item{...}{additional arguments}
  }
\details{
sdm fits multiple models and can be used to generate multiple runs (replicates) of each method through partitioning (using one or several partitioning methods including: \code{subsampling}, \code{cross-validation}, and \code{bootstrapping}.

Each model is evaluated against training data, and if available, splitted data (through partitioning; called dependent test data as well, i.e., "dep.test") and/or indipendent test data ("indep.test").

User should make sure that the methods are available and the required packages for them are installed before putting their names in the function, otherwise, the methods that cannot be run for any reason, are excluded by the function. It is a good practice to call \code{\link{installAll}} function (just one time when the sdm is installed), that tries to install all the packages that may be needed somewhere in the \code{sdm} package.

A new method can be adopted and added to the package by a user using \code{\link{add}} function. It is also possible to get an instance of an existing method, override the setting and definition, and then add it with a new name (e.g., my.glm).

The output would be a single object (\code{sdmModels}) that can be read/reproduced everywhere (e.g., on a new machine). A setting object can also be taken (exported) out of the output \code{sdmModels} object, that can be used to reproduce the same practice but given new conditions (i.e., new dataset, area. etc.)

To speed up the model fitting, you may use parallel processing (a high-performance computing solution) by providing a list of items that can be passed to \code{parallelSetting} argument. The items in the list includes:

\code{ncore}: defines the number of cores (it can also be specified outside of this list)

\code{method}: a character (default: "parallel"), defines the parallelising engine. Currently, three options are available including 'parallel', 'foreach', and 'future'. 

\code{doParallel}: optional, definition to register for a backend for parallel processing (needed when method = 'foreach'). It should be provided as an R expression like the following example:

expression(registerDoParallel(parallelSetting@cl))

The above example is based on the function available in doParallel package. Other packages can also be used to provide and register backend technologies (e.g., doMC)

\code{cluster}: optional, in case a cluster is created and available (e.g., using cl <- parallel::makeCluster(2)), the cluster object can be introduced here to be used as the parallel processing engine, otherwise, it is handled by the sdm package.

\code{hosts}: optional, to use remote machines/clusters in the parallel processing, a character vector with the addresses (names or IPs) of the accessible (on the network) remote clusters can be provided here to be registered and used in parallel processing (still under development so it may not work appropriately!)

\code{fork}: logical, Available for non-windows operating systems and specifies whether a fork solution should be used for the parallelisation. Default is TRUE for non-windows OS and FALSE for windows.

NOTE: Only use parallelSetting when you deal with a big dataset or large number of models otherwise, it makes the procedure slower rather than faster if the procedure is quick without parallelising!


}

\value{
 an object of class \code{sdmModels}
}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org/}
}



\examples{
\dontrun{
file <- system.file("external/pa_df.csv", package = "sdm")

df <- read.csv(file = file)

head(df) 

d <- sdmData(formula = sp ~ b15 + NDVI, train = df)

d
#----
# Example 1: fit using 3 models, and no evaluation (evaluation based on training dataset):

m <- sdm(formula = sp ~ b15 + NDVI, data = d, methods = c('glm', 'gam', 'gbm'))

m

# Example 3: fit using 5 models, and 
# evaluates using 10 runs of subsampling replication taking 30 percent as test:

m <- sdm(formula = sp ~ b15 + NDVI, data = d, methods = c('glm', 'gam', 'gbm', 'svm', 'rf'),
          replication = 'sub', test.percent = 30, n = 10)

m


# Example 3: fit using 5 models, and 
# evaluates using 10 runs of both 5-folds cross-validation and bootsrapping replication methods

m <- sdm(formula = sp ~., data = d, methods = c('gbm', 'tree', 'mars', 'mda', 'fda'),
          replication = c('cv', 'boot'), cv.folds = 5, n = 10)

m

# Example 4: fit using 3 models; evaluate the models using subsampling, 
# and override the default settings for the method brt:

m <- sdm(formula = sp ~ b15 + NDVI, data = d, methods = c('glm', 'gam', 'brt'), test.percent = 30,
          modelSettings = list(brt = list(n.trees = 500, train.fraction = 0.8)))

m

}

}

\keyword{spatial}
\keyword{model}
\keyword{data}
\keyword{species}
