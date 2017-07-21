\name{CensMmix}
\alias{CensMmix}
\docType{package}
\title{
Censored multivariate finite mixture model
}
\description{
Performs a Finite Mixture Censored multivariate (FM-MC) Student-t and Normal distribution using using EM-type algorithm for iteratively computing maximum likelihood estimates of the parameters.
}
\usage{
CensMmix(cc, y, nu=3, mu=NULL, Sigma = NULL, pii = NULL, g = NULL, get.init = TRUE,
         criteria = TRUE, group = FALSE, family = "Normal", error = 0.0001,
         iter.max = 300, uni.Sigma = FALSE, obs.prob= FALSE, kmeans.param = NULL)
}

\arguments{
\item{cc}{Vector of censoring indicators. For each observation: 0 if non-censored, 1 if censored.}

\item{y}{Vector of responses in case of right censoring.}

\item{nu}{Initial value for the EM algorithm, nu it's degrees of freedom. Value of one size 1 (If Student's t)}

\item{mu}{Initial value for the EM algorithm. Each of them must be a vector of length g.(the algorithm considers the number of components to be fitted based on the size of these vectors)}

\item{Sigma}{a list of \code{g} arguments of matrices of initial values (dimension pxp) for the scale parameters}

\item{pii}{Initial value for the EM algorithm. The vector of initial values (dimension g) for the weights for each cluster. Must sum one!}

\item{g}{Numbers of components}

\item{get.init}{TRUE or FALSE. It indicates if the program (TRUE) is get the initial values or if the user (FALSE) entered manually the initial values.}

\item{criteria}{It indicates if are calculated the criterion selection methods (AIC, BIC, EDC and ICL)}

\item{group}{TRUE or FALSE.}

\item{family}{"t":  fits a t-student regression mixture for censured data or "Normal":  fits a Normal regression mixture censored data}

\item{error}{define the stopping criterion of the algorithm}

\item{iter.max}{the maximum number of iterations of the EM algorithm}

\item{uni.Sigma}{TRUE: if the covariance matrix are equals or FALSE if are not equal}

\item{obs.prob}{TRUE or FALSE.}

\item{kmeans.param}{Parameters for the k-means clustering algorithm}

}

%\details{


%}

\author{Luis Benites Sanchez \email{lbenitesanchez@gmail.com},
Victor Hugo Lachos \email{hlachos@ime.unicamp.br},
Edgar J. Lopez Moreno \email{edgar.javier.lopez.moreno@gmail.com}
}
\encoding{latin1}
\references{
Arellano-Valle, R. B., Castro, L., Gonzalez-Farias, G. & Munos Gajardo, K. (2012). Student-t censored regression model: properties and inference. Statistical Methods and Applications, 21(4), 453-473.

Dempster, A., Laird, N. & Rubin, D. (1977). Maximum likelihood from incomplete data via the EM algorithm. Journal of the Royal Statistical Society, Series B,39, 1-38.

Peel, D. & McLachlan, G. J. (2000). Robust mixture modelling using the t distribution. Statistics and Computing,10(4), 339-348.

Karlsson, M. & Laitila, T. (2014). Finite mixture modeling of censored regression models. Statistical Papers,55(3), 627-642.

Basso,R.M.,Lachos,V.H.,Cabral,C.R.B. & Ghosh,P. (2010). Robust mixture modeling based on scale mixtures of skew-normal distributions. Computational Statistics & Data Analysis, 54(12), 2926-2941.

Basford, K., Greenway, D.,McLachlan,G. & Peel,D. (1997). Standard errors of fitted component means of normal mixtures. Computational Statistics,12, 1-18.
}

\seealso{
\code{\link{CensMmix}},\code{\link{concentration}}
}

\examples{
#See examples for the CensMmix function linked above.
}
\keyword{Finite Mixture}
