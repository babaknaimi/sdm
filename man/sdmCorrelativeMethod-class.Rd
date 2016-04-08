\docType{class}
\name{sdmCorrelativeMethod-class}
\alias{sdmCorrelativeMethod-class}
\title{sdmCorrelativeMethod class}
\description{
An S4 class representing sdm dataset
}
\section{Slots}{

\describe{
\item{\code{name}}{Modelling method name}

\item{\code{aliases}}{Alternative names for the method}

\item{\code{dataArgument.names}}{A list keeps the name of data agruments in both fit and predict functions}

\item{\code{packages}}{The required external package by the method}

\item{\code{modelTypes}}{Specifies whether the model is presence-absence, presence-only, abundance, or multinomial}

\item{\code{fitParams}}{a list of parameters needed by the method}

\item{\code{fitSettings}}{a list of setting parameters for the method}

\item{\code{settingRules}}{a function that adjust the setting parameters according to data}

\item{\code{fitFunction}}{The main function use for fitting the model}

\item{\code{tuneParams}}{a list of parameters to be tuned before the final fitting}

\item{\code{predictParams}}{a list of parameters needed by predict function}

\item{\code{predictSettings}}{a list of setting parameters for prediction}

\item{\code{predictFunction}}{The main predict function}

\item{\code{metadata}}{a metadata object containing the information about who creates the object, date, etc.}

\item{\code{.temp.env}}{an environment object containing the functions defined by a user that is not from a package}
}}

