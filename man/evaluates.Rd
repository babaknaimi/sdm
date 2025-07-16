\name{evaluates}

\alias{evaluates}
\alias{getEvaluation}
\alias{getReplication}
\alias{evaluates,vector,vector-method}
\alias{evaluates,sdmdata,RasterLayer-method}
\alias{evaluates,sdmdata,SpatRaster-method}
\alias{evaluates,sdmModels,SpatRaster-method}
\alias{evaluates,sdmModels,missing-method}
\alias{getEvaluation,sdmModels-method}
\alias{getReplication,sdmModels-method}
\alias{show,sdmEvaluate-method}


\title{evaluate for accuracy}

\description{
evaluates for accuracy 
}

\usage{
evaluates(x, p,...)

getEvaluation(x, id, wtest, stat, opt,...)

getReplication(x, id, replication, species, run, index, test)
}

\arguments{
  \item{x}{a numeric vector or a \code{sdmdata} object including the observed values; a \code{sdmModels} object in \code{getEvaluation}}
  \item{p}{a numeric vector or a \code{RasterLayer} including the predicted values}
  \item{id}{a single numeric value, indicates the modelID}
  \item{wtest}{character, which test data should be used: "training", "test.dep", or "test.indep"?}
  \item{stat}{character, statistics that should be extracted from the \code{sdmEvaluate} object}
  \item{opt}{a numeric value, indicates which threshold optimisation criterion should be considered if a threshold-based statistic is selected in stat}
  \item{species}{optional (default: NULL); a character vector, specifies the name of species for which the replication is returned}
  \item{replication}{a character, specifies the name of the replication method}
  \item{run}{a single numeric value, specifies the replication ID}
  \item{index}{logical (default: FALSE); specifies whether the index or species data of drawn records should be returned}
  \item{test}{logical (default: TRUE); specifies whether the test partition should be returned or training partition}
  \item{...}{additional arguments (see details)}
  }
\details{
Evaluates the preformance (accuracy) given the observed values, and the predicted values. As additional argument, the distribution of data can be specified (through \code{distribution}), that can be either of \code{'binomial'}, \code{'gaussian'}, \code{'laplase'}, or \code{'poisson'}. If not specified, it will be guessed by the function!

\code{getEvaluation} can be used to get the evaluation results from a fitted model (\code{sdmModels} object that is output of the \code{sdm} function). Each model in \code{sdmModels} has a modelID, that can be specified in \code{w} argument. If \code{w} is not specified or more than a modelID is specified, then a data.frame is generated that contains the statistics specified in \code{stat}. For a single model (if the length of \code{w} is 1), \code{stat} can be 1 (threhold_independent statistics), or 2 (threshold_based statistics) or NULL (both groups). If more than a model is specified (\code{w} is either NULL or has a length greater than 1), stat can be the name of statistics such as \code{'AUC', 'COR', 'Deviance', 'obs.prevalence', 'threshold', 'sensitivity', 'specificity', 'TSS','MCC', 'Kappa', 'NMI', 'phi', 'ppv', 'npv', 'ccr', 'prevalence'}. 
If either of the threshold_based stats is selected, \code{opt} can also be specified to select one of the criteria for optimising the threshold. The possible value can be between 1 to 15 for \code{"sp=se", "max(se+sp)", "min(cost)", "minROCdist", "max(kappa)", "max(ppv+npv)", "ppv=npv", "max(NMI)", "max(ccr)", "prevalence", "max(MCC)", "P10", "P5", "P1", "P0"} criteria, respectively. P10, P5, P1 refer to 10th, 5th, and 1st percentiles of presence records in the evaluation dataset, respectively for which the suitability value is used as the threshold. By choosing P0, the minimum suitability value across presence records is selected as the threshold.

\code{getReplication} returns portion of records randomly selected through data partitioning using one of the replication methods (e.g., 'cv', 'boot', 'sub').
}

\value{
 an object of class \code{sdmEvaluate} from \code{evaluates} function
 
 a list or data.frame from \code{getEvaluation} function
}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org/}

}

\seealso{# }

\examples{
\dontrun{
file <- system.file("external/model.sdm", package = "sdm")

m <- read.sdm(filename = file) # a sdmModels Object (fitted using sdm function)

getModelInfo(x = m)

# there are 4 models in the sdmModels object

# so let's take a look  at all the results for the model with a modelID of 1

# evaluation using training data (both threshold_independent and threshold_based groups):

getEvaluation(m, w = 1, wtest = 'training') 

getEvaluation(m, w = 1, wtest = 'training', stat = 1) # stat = 1 (threshold_independent)

getEvaluation(m, w = 1, wtest = 'test.dep', stat = 2) # stat = 2 (threshold_based)

getEvaluation(m, w = 1:3, wtest = 'test.dep', stat = c('AUC','TSS'), opt = 2) 

getEvaluation(m, opt = 1) # all models

getEvaluation(m, stat = c('TSS', 'Kappa', 'AUC'), opt = 1) # all models


############

#example for evaluation:

evaluates(x = c(1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0),
          p = c(0.69, 0.04, 0.05, 0.95, 0.04, 0.65, 0.09, 0.61, 0.75, 0.84, 0.15))

##############
# Example for getReplication:


df <- read.csv(file) # load a csv file

head(df)

d <- sdmData(formula = sp ~ b15 + NDVI, train = df) # sdmdata object

d
#----
# fit SDMs using 2 methods and a subsampling replication method with 2 replications:

m <- sdm(formula = sp ~ b15 + NDVI, data = d, methods = c('glmpoly', 'gbm'), replication = 'sub', test = 30, n = 2)

m


# randomly drawn species records for test data in the second replication (run) of subsampling:
getReplication(formula = m, replication = 'sub', run = 2) 

getReplication(formula = m, replication = 'sub', run = 2, test = FALSE) # drawn record in the training partition

ind <- getReplication(x = m, replication = 'sub', run = 2, index = TRUE) # index of the selected test record 

head(ind)

.df <- as.data.frame(m@data) # convert sdmdata object in the model to data.frame

head(.df)

.df <- .df[.df$rID \%in\% ind, ] # the full test dataset drawn (second replication)

pr <- predict(m,.df) # predictions of all the methods for the test dataset 



pr <- predict(m, .df) # predictions of all the methods for the test dataset 


head(pr) 

e <- evaluates(x = .df$sp, p = pr[,1]) # evaluates for the first method using the selected test data

e@statistics

e@threshold_based


}
}

\keyword{spatial}
\keyword{data}
\keyword{accuracy}