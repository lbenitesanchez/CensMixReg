\name{initial.values.fmr.smn.cr}
\alias{initial.values.fmr.smn.cr}
\docType{data}
\title{Initial values for the FMR-SMN-CR}
\description{Obtained the initial values for the parameter beta_j,sigma_j, p_j and nu the SMN distribution.
}
\usage{
initial.values.fmr.smn.cr(cc, y, x, g=2, algorithm="k-medoids", family="T",
lower=1, upper=20, space=0.1, plotLog = TRUE, searchNU=TRUE, printNU=TRUE,
saveFigure = FALSE)
}

\arguments{
  \item{cc}{Vector of censoring indicators. For each observation: 0 if non-censored, 1 if censored.}

  \item{y}{Vector of responses in case of right censoring.}

  \item{x}{Matrix or vector of covariates for each component}

  \item{g}{Numbers of components}

  \item{algorithm}{It indicates the algorithm: "trim-kmeans","MinMax_kmeans", "k-means" and "k-medoids"}

  \item{family}{"T", "Normal", "Slash" or "NormalC"}

  \item{lower}{lower value of the search of the nu or gamma parameter}

  \item{upper}{upper value of the search of the nu or gamma parameter}

  \item{space}{number: increment of the sequence}

  \item{plotLog}{TRUE or FALSE the figure of profile}

  \item{searchNU}{TRUE or FALSE the search the nu or gamma parameter}

  \item{printNU}{TRUE or FALSE the print the nu or gamma parameter}

  \item{saveFigure}{TRUE or FALSE save the figure of profile}
}



\examples{
  #See examples for the wage.rates function linked above.
}


\keyword{datasets}



