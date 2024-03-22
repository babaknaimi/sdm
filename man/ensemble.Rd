\name{ensemble}


\docType{methods}

\alias{ensemble}
\alias{ensemble,sdmModels-method}
\alias{ensemble,sdmModels,Raster-method}
\alias{ensemble,sdmModels,SpatRaster-method}
\alias{ensemble,sdmModels,data.frame-method}

\title{Ensemble Forecasting of SDMs}

\description{
Make a Raster object with a weighted averaging over all predictions from several fitted model in a sdmModel object.
}
\usage{
\S4method{ensemble}{sdmModels}(x, newdata, filename="",setting,overwrite=FALSE,pFilename="",...)
}


\arguments{
  \item{x}{a sdmModels object}
  \item{newdata}{raster object or data.frame, can be either predictors or the results of the \code{predict} function}
  \item{filename}{optional character, output file name (if newdata is raster object)}
  \item{setting}{list, contains the parameters that are used in the ensemble procedure; see details}
  \item{overwrite}{logical, whether existing filename is overwritten (if exists and filename is given) }
  \item{pFilename}{it is ignored if newdata is the output of \code{predict}, otherwise, since the ensemble first call \code{predict}, it specifies the filename to write the output of \code{predict (if newdata is raster) }}
  \item{...}{additional arguments pass to the \code{writeRaster} function (if used)}
}
\details{
ensemble function uses the fitted models in an \code{sdmModels} object to generate an ensemble/consensus of predictions by multiple individual models. Several \code{ensemble} methods are available and can be defined in the setting argument.

A list of settings can be introduced in the \code{setting} argument including:

- \code{method}: a character vector specifies which ensemble method(s) should be employed (multiple choice is possible). The details about the available methods are provided at the end of this page.

- \code{stat}: if the - \code{method='weighted'} is used, it specifies which evaluation metrics can be used as weight in the weighted averaging procedure. Alternatively, one may directly introduce weights (see the next argument).

- \code{weights}: an optional numeric vector (with a length equal to the models that are successfully fitted) to specify the weights for weighted averaging procedure (if the method='weighted' is specified).

- \code{id}: specifies the model IDs that should be considered in the ensemble procedure. If missing, all the models that are successfully fitted are considered.

- \code{expr}: A character or an expression specifies a condition to select models for the ensemble procedure. For example: \code{expr='auc > 0.7'} only use models with AUC accuracy greater than 0.7. OR  \code{expr='auc > 0.7 & tss > 0.5'} subsets models based on both AUC and TSS metrics.

- \code{wtest}: specifies which test dataset ("training","test.dep","test.indep") should be used to extract the statistic (stat) values as weights (if a relevant method is specified)

- \code{opt}: if a thershold_based metric is used in  is selected in \code{stat} or in \code{expr}, \code{opt} specifies the threshold selection criterion. The possible value can be between 1 to 14 for \code{"sp=se", "max(se+sp)", "min(cost)", "minROCdist", "max(kappa)", "max(ppv+npv)", "ppv=npv", "max(NMI)", "max(ccr)", "prevalence", "P10", "P5", "P1", "P0"} criteria, respectively.

- \code{power}: default: 1, a numeric value to which the weights are raised. Greater value than 1 affects weighting scheme (for the methods e.g., "weighted") to increase the weights for the models with greater weight. For example, if weights are c(0.2,0.2,0.2,0.4), raising them to power 2 would be resulted to new weights as c(0.1428571,0.1428571, 0.1428571, 0.5714286) that causes greater contribution of the models with  greater performances to the ensemble output.


---> The available ensemble methods (to be specified in \code{method}) include:
 
 -- 'unweighted': unweighted averaging/mean.
 
 -- 'weighted': weighted averaging.
 
 -- 'median': median.
 
 -- 'pa': mean of predicted presence-absence values (predicted probabilities are first converted to presence-absence given a threshold (\code{opt} defines which threshold optimisation strategy should be used), then they are averaged).
 
 -- 'mean-weighted': A two step averaging, that can be used when several replications are available for each modelling methods (e.g., fitted through bootstrapping or cross-validation resampling); it first takes an unweighted mean over the predicted values of multiple replications for each method (within model averaging), then a weighted mean is employed to combine the probabilities of different methods (between models averaging).

-- 'mean-unweighted': Same as the previous one, but an unweighted mean is also used for the second step (instead of weighted mean).

-- 'median-weighted': Same as the 'mean-weighted, but the median is used in the first step.

-- 'median-unweighted': another two-step method, median is used for the first step and unweighted mean is used for the second step.

----> in addition to tne ensemble methods, some other methods are available to generate some outputs that can represent uncertainty:

-- 'uncertainty' or 'entropy': this method generates the uncertainty among the models' predictions that can be interpreted as model-based uncertainty or inconsistency among different models. It ranges between 0 and 1, 0 means all the models predicted the same value (either presence or absence), and 1 referes to maximum uncertainy, e.g., half of the models predicted presence (or absence) and the other half predicted the oposite value.

-- 'cv': Coefficient of variation of probabilities generated from multiple models

-- 'stdev': Standard deviation of probabilities generated from multiple models

-- 'ci': This generates confidence interval length (marginal error) which assigns the difference between upper and lower limits of confidence interval to each pixel (upper - lower). The default level of confidence interval is 95\% (i.e., \code{alpha = 0.05}), unless a different \code{alpha} is defined in \code{setting}. In case two separate upper and lower rasters are needed, by using the following codes, the upper and lower limits can be calculated:

\code{en <- ensemble(x, newdata, setting=list(method=c('mean','ci')))} # taking unweighted averaging and ci

# en[[1]] is the mean of all probabilities and en[[2]] is the ci
\code{ci.upper <- en[[1]] + en[[2]] / 2} # adding marginal error (half of the generated ci) to mean
\code{ci.lower <- en[[1]] - en[[2]] / 2} # subtracting marginal error from mean

\code{plot(ci.upper,main='Upper limit of Confidence Interval - alpha = 0.05')}

\code{plot(ci.lower,main='Lower limit of Confidence Interval - alpha = 0.05')}
}

\value{
 - a Raster object if \code{predictors} is a Raster object
 
 - a numeric vector (or a data.frame) if \code{predictors} is a data.frame object
}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org/}
}

\seealso{# }

\examples{
\dontrun{


file <- system.file("external/species.shp", package="sdm") # get the location of the species data

species <- vect(file) # read the shapefile

path <- system.file("external", package="sdm") # path to the folder contains the data

lst <- list.files(path=path,pattern='asc$',full.names = T) # list the name of the raster files 


# stack is a function in the raster package, to read/create a multi-layers raster dataset
preds <- rast(lst) # making a raster object

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