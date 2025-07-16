\name{sdmSetting}

\alias{sdmSetting}
\alias{sdmSetting,ANY,sdmdata,character-method}
\alias{show,sdmSetting-method}

\title{creating sdmSetting object}

\description{
Creates sdmSetting object that holds settings to fit and evaluate the models. It can be used to reproduce a study.
}

\usage{
sdmSetting(formula, data, methods, interaction.depth = 1, n = 1, replication = NULL, cv.folds = NULL,
     test.percent = NULL, bg = NULL, bg.n = NULL, var.importance = NULL, response.curve = TRUE,
     var.selection = FALSE, modelSettings = NULL, seed = NULL, parallelSetting = NULL,...)


}

\arguments{
  \item{formula}{specifies the structure of the model}
  \item{data}{sdm data object or data.frame including species and feature data}
  \item{methods}{character vector, name(s) of the algorithms}
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
  \item{modelSettings}{optional list, settings for modelling methods, can be specified by users}
  \item{seed}{default is NULL, either logical, specifies whether a seed for random number generator should be considered, or a numeric to specify the exact seed number}
  \item{parallelSetting}{default is NULL, a list, includes setting items for parallel processing. The items in parallel setting include: ncore, method, type, hosts, doParallel, and fork. See details for more information.}
  \item{...}{additional arguments}
  }
\details{
using sdmSetting, the feature types, interaction.depth and all settings of the model can be defined. This function generates a sdmSetting object that can be specifically helpful for reproducibility. The object can be shared with other users or may be used for other studies.

If a user aims to reproduce the same results for every time the code is running with the same data and settings, a seed number should be specified. Through the \code{seed} argument, a user can specify \code{NULL}, means a seed should not be set (if a random sampling is incorporated in the modelling procedure, for different runs the results would be different); \code{TRUE}, means a seed should be set (the seed number is randomly selected and used everytime the same setting is incorporated); \code{a number}, means the seed will be set to the number specified by the user. 

For parallel processing, a list of items can be passed to \code{parallelSetting}, including:

\code{ncore}: defines the number of cores (it can also be specified outside of this list)

\code{method}: character (default: "parallel"), defines the parallelising engine. Currently, three options are available including 'parallel', 'foreach', and 'future'. 

\code{doParallel}: optional, definition to register for a backend for parallel processing (needed when method = 'foreach'). It should be provided as an R expression like the following example:

expression(registerDoParallel(parallelSetting@cl))

The above example is based on the function available in doParallel package. Other packages can also be used to provide and register backend technologies (e.g., doMC)

\code{cluster}: optional, in case a cluster is created and available (e.g., using cl <- parallel::makeCluster(2)), the cluster object can be introduced here to be used as the parallel processing engine, otherwise, it is handled by the sdm package.

\code{hosts}: optional, to use remote machines/clusters in the parallel processing, a character vector with the addresses (names or IPs) of the accessible (on the network) remote clusters can be provided here to be registered and used in parallel processing (still under development so it may not work appropriately!)

\code{fork}: logical, available for non-windows operating systems and specifies whether a fork solution should be used for the parallelisation. Default is TRUE for non-windows OS and FALSE for windows.


NOTE: Only use parallelSetting when you deal with a big dataset or large number of models otherwise, it makes the procedure slower rather than faster if the procedure is quick without parallelising!

}

\value{
 an object of class \code{sdmSettings}
}


\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org/}
}

\examples{
\dontrun{
file <- system.file("external/pa_df.csv", package = "sdm")

df <- read.csv(file = file)

head(x = df) 

d <- sdmData(formula = sp ~ b15 + NDVI, train = df)

# generate sdmSettings object:
s <- sdmSetting(formula = sp ~., methods = c('glm', 'gam', 'brt', 'svm', 'rf'),
          replication = 'sub', test.percent = 30, n = 10, modelSettings = list(brt = list(n.trees = 500)))

s


}
}

\keyword{spatial}
\keyword{model}
\keyword{data}

