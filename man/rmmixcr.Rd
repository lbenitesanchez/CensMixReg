\name{rmmixcr}
\alias{rmmixcr}
\title{Random multivariate Finite Mixture Censored generator }
\description{
Random generator of multivariate FM-SMSN distributions.}
\usage{
rmmixcr(n, pii, mu, Sigma, shape, nu, percCensu, family)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of observations}
  \item{pii}{a vector of weights for the mixture (dimension of the number \code{g} of clusters). Must sum to one!}
  \item{mu}{a list of \code{g} lists with each list containing the necessary parameters of the selected \code{family} }
  \item{Sigma}{a list of \code{g} arguments of matrices of initial values (dimension pxp) for the scale parameters.}
  \item{shape}{a list of \code{g} arguments of vectors of initial values (dimension p)for the skewness parameters.}
  \item{nu}{nu it's degrees of freedom. }
  \item{percCensu}{Matrix of censoring indicators. For each observation: 0 if non-censored, 1 if censored.}
  \item{family}{ distribution family to be used in fitting ("t", "Skew.t", "Skew.slash", "Skew.normal", "Normal") }
}

\examples{
mu          <- Sigma <- shape <- family <- list()
mu[[1]]     <- c(-3,-1)
mu[[2]]     <- c(2,3)
Sigma[[1]]  <- matrix(c(8,1,1,5.5), 2,2)
Sigma[[2]]  <- matrix(c(2,1,1,2.5), 2,2)
family[[1]] <- "Skew.t"
family[[2]] <- "Skew.normal"
shape[[1]]  <- c(4,5)
shape[[2]]  <- c(1,1)
nu          <- c(3,3)
pii         <- c(0.65,0.35)
percCensu   <- c(0.15,0.15)

data        <- rmmixcr(100, pii , mu , Sigma, shape, nu, percCensu, family)
}

\author{Luis Benites Sanchez \email{lbenitesanchez@gmail.com}, Victor Hugo Lachos \email{hlachos@ime.unicamp.br}, Edgar J. Lopez Moreno \email{edgar.javier.lopez.moreno@gmail.com}
}

\seealso{ \code{\link{smsn.mmix}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

