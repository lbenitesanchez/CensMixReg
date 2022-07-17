\name{fmr.smn.cr}
\alias{fmr.smn.cr}
\docType{package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Censored mixture regression models based in the Scale Mixture of Normal (SMN) distribution
}
\description{
  Performs a Finite Mixture Regression (FMR) with censored based in the SMN  using EM-type algorithm for iteratively computing maximum likelihood estimates of the parameters.
}
\usage{
  fmr.smn.cr(cc, y, x, Abetas = NULL, sigma2 = NULL, pii = NULL, nu=NULL, g = NULL,
               family = "Normal", error = 0.00001, iter.max = 100)
}

\arguments{
  \item{cc}{Vector of censoring indicators. For each observation: 0 if non-censored, 1 if censored.}

  \item{y}{Vector of responses in case of right censoring.}

  \item{x}{Matrix or vector of covariates for each component}

  \item{Abetas}{Parameters of vector regression dimension \eqn{(p_j + 1)} include or not intercept, j=1,...,G}

  \item{sigma2}{Initial value for the EM algorithm. Each of them must be a vector of length g.(the algorithm considers the number of components to be adjusted based on the size of these vectors)}

  \item{pii}{Initial value for the EM algorithm. Each of them must be a vector of length g.(the algorithm considers the number of components to be adjusted based on the size of these vectors)}

  \item{nu}{Initial value for the EM algorithm, nu it's degrees of freedom. Value of one size 1 (If Student's t or Slash) or size 2 (if Contaminated Normal)}

  \item{g}{Numbers of components}

  \item{family}{"T":  fits a t-student regression mixture for censured data or "Normal":  fits a Normal regression mixture censored data or "Slash": fits a Slash regression mixture censored data or "NormalC": fits a Contaminated Normal regression mixture censored data}

  \item{error}{define the stopping criterion of the algorithm}

  \item{iter.max}{the maximum number of iterations of the EM algorithm}
}

%\details{


  %}

\author{Luis Benites Sanchez \email{lbenitesanchez@gmail.com}, Victor Hugo Lachos \email{hlachos@ime.unicamp.br}
}

\references{
Zeller, C. B., Cabral, C. R. B. and Lachos, V. H. (2016). Robust mixture regression modeling based on scale mixtures of skew-normal distributions. Test, 25, 375-396.
}

\keyword{package}

\seealso{
  \code{\link{fmr.smn.cr}},\code{\link{wage.rates}}
}

\examples{
  #See examples for the fmr.smn.cr function linked above.
}
\keyword{Mixture Regression}
