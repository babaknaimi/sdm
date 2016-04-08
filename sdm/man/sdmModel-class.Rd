\name{sdmModels-classes}

\docType{class}

\alias{sdmModels-class}
\alias{sdmEvaluate-class}
\alias{sdmFormula-class}

\alias{listORcharacter-class}
\alias{characterORnull-class}

\alias{CRSorNULL-class}
\alias{formulaORnull-class}
\alias{numericORnull-class}
\alias{characterORmissing-class}
\alias{listORnull-class}
\alias{functionORnull-class}
\alias{matrixORnull-class}
\alias{data.frameORnull-class}
\alias{functionORcharacter-class}
\alias{environmentORnull-class}


\alias{show,sdmModels-method}

\title{sdmModels classes}

\description{

An S4 class to keep all the information of fitted models as well as their evaluations.

}

\section{Slots}{
Slots for sdmModels objects:
	\describe{
\item{\code{data}}{a sdmdata object}

\item{\code{recordIDs}}{Contains the species data}

\item{\code{setting}}{A data.frame containing predictor variables}

\item{\code{run.info}}{a data.frame containing info on runs}

\item{\code{replicates}}{The names of categorical variables (if any)}

\item{\code{models}}{a list contains all fitted objects and relevant information (e.g., evaluation) }



\bold{ Slots for sdmEvaluate objects:}
\item{\code{observed}}{a numeric vector of observed values}

\item{\code{predicted}}{a numeric vector of predicted values}

\item{\code{statistics}}{a list of threshold-independent statistics}
\item{\code{threshold_based}}{a data.frame of threshold-based statistics}

\bold{ Slots for sdmFormula objects:}
\item{\code{formula}}{input formula}

\item{\code{vars}}{character, name of variables}

\item{\code{model.terms}}{the formula terms used in model fitting}
\item{\code{data.terms}}{the formula terms used to manipulate data}
  
  }
  

}

\author{Babak Naimi 

\email{naimi.b@gmail.com}

\url{http://biogeoinformatics.org}
}



\keyword{classes}
\keyword{spatial}
