\name{ensemble}


\docType{methods}

\alias{ensemble}
\alias{ensemble,sdmModels-method}
\alias{ensemble,sdmModels,Raster-method}
\alias{ensemble,sdmModels,data.frame-method}

\title{Ensemble Forecasting of SDMs}

\description{
Make a Raster object with a weighted averaging over all predictions from several fitted model in a sdmModel object.
}
\usage{
\S4method{ensemble}{sdmModels}(x, newdata, filename="",setting,...)
}


\arguments{
  \item{x}{a sdmModels object}
  \item{newdata}{Raster* object or data.frame, can be either predictors or the results of the predict function}
  \item{filename}{character, output file name}
  \item{setting}{list, contains the parameters that are used in the ensemble procedure; see details}
  \item{...}{additional arguments passed to the predict function}
  
  }
\details{
ensemble function uses the fitted models in an \code{sdmModels} object to generate an ensemble/consensus of predictions by individual models. Several methods do exist for this procedure, that are (or will be) implemented in this function, and can be defined in the method argument.
A list can be introduced in the \code{setting} argument in which several parameters can be set including:

- \code{method}: specify which ensemble method should be used. Several methods are implemented including:
 
 -- 'unweighted': unweighted averaging/mean.
 
 -- 'weighted': weighted averaging.
 
 -- 'median': median.
 
 -- 'pa': mean of predicted presence-absence values (predicted probability of occurrences are first converted to presence-absence using a threshold, then they are averaged).
 
 -- 'mean-weighted': A two step mean that is when several replications are fitted for each modelling methods (e.g., through bootstrapping or cross-validation), using this method an unweighted mean is taken over the predicted values of different replications of each method (i.e., within model averaging), then a weighted mean is used to combine them into final ensemble values (i.e., between models averaging).

-- 'mean-unweighted': Same as the previous one, but an unweighted mean is also used for the second step (instead of weighted mean).

-- 'median-weighted': Same as the 'mean-weighted, but the median is used instead of unweighted mean.

-- 'median-unweighted': another two-step method, median is used for the first step and unweighted mean is used for the second step.

-- 'uncertainty' or 'entropy': this method generates the uncertainty among the models' predictions that can be interpreted as model-based uncertainty or inconsistency among different models. It ranges between 0 and 1, 0 means all the models predicted the same value (either presence or absence), and 1 referes to maximum uncertainy, e.g., half of the models predicted presence (or absence) and the other half predicted the oposite value.


- \code{stat}: if the method='weighted' is used, this specify which evaluation statistics can be used as weight in the weighted averaging procedure. Alternatively, one may directly introduce weights (see the next argument)

- \code{weights}: an optional numeric vector (with a length equal to the models that are successfully fitted) to specify the weights for weighted averaging procedure (if the method='weighted' is specified)

- \code{id}: specify the model IDs that should be considered in the ensemble procedure. If missing, all the models that are successfully fitted are considered.

- \code{wtest}: specify which test dataset ("training","test.dep","test.indep") should be used to extract the statistic (stat) values as weights (if a relevant method is specified)

- \code{opt}: If either of the thershold_based stats are selected, \code{opt} can be also specified to select one of the criteria for optimising the threshold. The possible value can be between 1 to 10 for \code{"sp=se", "max(se+sp)", "min(cost)", "minROCdist", "max(kappa)", "max(ppv+npv)", "ppv=npv", "max(NMI)", "max(ccr)", "prevalence"} criteria, respectively.

- \code{power}: default: 1, a numeric value to which the weights are raised. Greater value than 1 affects weighting scheme (for the methods e.g., "weighted") to increase the weights for the models with greater weight. For example, if weights are c(0.2,0.2,0.2,0.4), raising them to power 2 would be resulted to new weights as c(0.1428571,0.1428571, 0.1428571, 0.5714286) that causes greater influence of the models with  greater performances to the ensemble output.



}

\value{
 - a Raster object if \code{predictors} is a Raster object
 
 - a numeric vector if \code{predictors} is a data.frame object
}

\references{
#
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org}
}

\seealso{# }

\examples{
\dontrun{


file <- system.file("external/species.shp", package="sdm") # get the location of the species data

species <- shapefile(file) # read the shapefile

path <- system.file("external", package="sdm") # path to the folder contains the data

lst <- list.files(path=path,pattern='asc$',full.names = T) # list the name of the raster files 


# stack is a function in the raster package, to read/create a multi-layers raster dataset
preds <- stack(lst) # making a raster object

d <- sdmData(formula=Occurrence~., train=species, predictors=preds)

d

# fit the models (5 methods, and 10 replications using bootstrapping procedure):
m <- sdm(Occurrence~.,data=d,methods=c('rf','tree','fda','mars','svm'),
          replicatin='boot',n=10)
    
# ensemble using weighted averaging based on AUC statistic:    
p1 <- ensemble(m, newdata=preds, filename='ens.img',setting=list(method='weighted',stat='AUC'))
plot(p1)

# ensemble using weighted averaging based on TSS statistic
# and optimum threshold critesion 2 (i.e., Max(spe+sen)) :    
p2 <- ensemble(m, newdata=preds, filename='ens2.img',setting=list(method='weighted',
                                                                  stat='TSS',opt=2))
plot(p2)

}


}

\keyword{spatial}
\keyword{sdm}
\keyword{modelling}
\keyword{learning}