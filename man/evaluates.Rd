\name{evaluates}

\alias{evaluates}
\alias{getEvaluation}
\alias{evaluates,vector,vector-method}
\alias{getEvaluation,sdmModels-method}

\title{evaluate for accuracy}

\description{
evaluates for accuracy 
}

\usage{
evaluates(x,p,...)

getEvaluation(x,w,wtest,stat,opt,...)
}

\arguments{
  \item{x}{a numeric vector including the observed values; or a \code{sdmModels} object}
  \item{p}{a numeric vector including the predicted values}
  \item{w}{a numeric vector indicates model IDs}
  \item{wtest}{which test, training, dep.test, or indep.test?}
  \item{stat}{statistics that should be extracted from the \code{sdmEvaluate} object}
  \item{opt}{a numeric value indicates which threshold optimisation criteria should be considered if a threshold-based statistic is selected in stat}
  \item{...}{additional arguments (see details)}
  }
\details{
Evaluates the preformance (accuracy) given the obsetved values, and the predicted values. As additional argument, the distribution of data can be specified (through \code{distribution}), that can be either of \code{'binomial'}, \code{'gaussian'}, \code{'laplase'}, or \code{'poisson'}. If not specified, it will be guessed by the function!

\code{getEvaluation} can be used to get the evaluation results from a fitted model (\code{sdmModels} object that is output of the \code{sdm} function). Each model in \code{sdmModels} has a modelID, that can be specified in \code{w} argument. If \code{w} is not specified or more than a modelID is specified, then a data.frame is generated that contains the statistics specified in \code{stat}. For a single model (if length \code{w} is 1), \code{stat} can be 1 (threhold_independent statistics), or 2 (threshold_based statistics) or NULL (both groups). If more than a model is specified (\code{w} is either NULL or has a length greater than 1), stat can be the name of statistics such as \code{'AUC', 'COR', 'Deviance', 'obs.prevalence', 'threshold', 'sensitivity', 'specificity', 'TSS', 'Kappa', 'NMI', 'phi', 'ppv', 'npv', 'ccr', 'prevalence'}. 
If either of the thershold_based stats are selected, \code{opt} can be also specified to select one of the criteria for optimising the threshold. The possible value can be between 1 to 10 for \code{"sp=se", "max(se+sp)", "min(cost)", "minROCdist", "max(kappa)", "max(ppv+npv)", "ppv=npv", "max(NMI)", "max(ccr)", "prevalence"} criteria, respectively.


}

\value{
 an object of class \code{sdmEvaluate} from \code{evaluates} function
 
 a list or data.frame from \code{getEvaluation} function
}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}

\url{http://biogeoinformatics.org}

}

\seealso{# }

\examples{
\dontrun{
file <- system.file("external/model.sdm", package="sdm")

m <- read.sdm(file) # a sdmModels Object (fitted using sdm function)

getModelInfo(m)

# there are 4 models in the sdmModels objects

# so let's take a look  at all the results for the model with modelID 1

# evaluation using training data (both threshod_independent and threshold_based groups):

getEvaluation(m,w=1,wtest='training') 

getEvaluation(m,w=1,wtest='training',stat=1) # stat=1 (threshold_independent)

getEvaluation(m,w=1,wtest='test.dep',stat=2) # stat=2 (threshold_based)

getEvaluation(m,w=1:3,wtest='test.dep',stat=c('AUC','TSS'),opt=2) 

getEvaluation(m,opt=1) # all models

getEvaluation(m,stat=c('TSS','Kappa','AUC'),opt=1) # all models


############

example for evaluation:

evaluates(x=c(1,1,0,1,0,0,0,1,1,1,0),
          p=c(0.69,0.04,0.05,0.95,0.04,0.65,0.09,0.61,0.75,0.84,0.15))

}
}

\keyword{spatial,species,data}