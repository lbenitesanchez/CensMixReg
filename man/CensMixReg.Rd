\name{CensMixReg}
\alias{CensMixReg}
\docType{package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Censored linear mixture regression models
}
\description{
Performs a Finite Mixture Censored (FM-CR) using using EM-type algorithm for iteratively computing maximum likelihood estimates of the parameters.
}
\usage{
CensMixReg(cc, y, x1, Abetas = NULL, medj= NULL, sigma2 = NULL, nu=NULL, pii = NULL,
g = NULL, get.init = TRUE, criteria = TRUE, group = FALSE, family = "T",
error = 0.00001, iter.max = 100, obs.prob= FALSE, kmeans.param = NULL)
}

\arguments{
\item{cc}{Vector of censoring indicators. For each observation: 0 if non-censored, 1 if censored.}

\item{y}{Vector of responses in case of right censoring.}

\item{x1}{Matrix or vector of covariates.}

\item{Abetas}{Parameters of vector regression dimension \eqn{(p + 1)} include intercept}

\item{medj}{Initial value for the EM algorithm. Each of them must be a vector of length g.(the algorithm considers the number of components to be fitted based on the size of these vectors)}

\item{sigma2}{Initial value for the EM algorithm. Each of them must be a vector of length g.(the algorithm considers the number of components to be adjusted based on the size of these vectors)}

\item{pii}{Initial value for the EM algorithm. Each of them must be a vector of length g.(the algorithm considers the number of components to be adjusted based on the size of these vectors)}

\item{nu}{Initial value for the EM algorithm, nu it's degrees of freedom. Value of one size 1 (If Student's t)}

\item{g}{Numbers of components}

\item{get.init}{TRUE or FALSE. It indicates if the program (TRUE) is get the initial values or if the user (FALSE) entered manually the initial values.}

\item{criteria}{It indicates if are calculated the criterion selection methods (AIC, BIC, EDC and ICL)}

\item{group}{TRUE or FALSE.}

\item{family}{"T":  fits a t-student regression mixture for censured data or "Normal":  fits a Normal regression mixture censored data}

\item{error}{define the stopping criterion of the algorithm}

\item{iter.max}{the maximum number of iterations of the EM algorithm}

\item{obs.prob}{TRUE or FALSE.}

\item{kmeans.param}{Parameters for the k-means clustering algorithm}

}

%\details{


%}

\author{Luis Benites Sanchez \email{lbenitesanchez@gmail.com}, Victor Hugo Lachos \email{hlachos@ime.unicamp.br}, Edgar J. Lopez Moreno \email{edgar.javier.lopez.moreno@gmail.com}
}

\references{
Benites, L., Lachos, V.H., Cabral, C.R.B. (2015). Robust Regression Modeling for Censored Data Based on Mixtures of Student-t Distributions. Technical Report 5, Universidade Estadual de Campinas. \url{http://www.ime.unicamp.br/sites/default/files/rp05-15.pdf}

Karlsson, M. & Laitila, T. (2014). Finite mixture modeling of censored regression models. Statistical papers, 55(3), 627-642.

Massuia, M. B., Cabral, C. R. B., Matos, L. A. & Lachos, V. H. (2014). Influence diagnostics for student-t censored linear regression models. Statistics, (ahead-of-print), 1-21.

Arellano-Valle, R., Castro, L., Gonzalez-Farias, G. & Munoz-Gajardo, K. (2012). Student-t censored regression model: properties and inference. Statistical Methods & Applications, 21, 453-473.

Garay, A. M., Lachos, V. H., Bolfarine, H. & Cabral, C. R. (2015). Linear censored regression models with scale mixtures of normal distributions. Statistical Papers, pages 1-32.

Arellano-Valle, R. B., Castro, L., Gonzalez-Farias, G. & Munos Gajardo, K. (2012). Student-t censored regression model: properties and inference. Statistical Methods and Applications, 21(4), 453-473.

Dempster, A., Laird, N. & Rubin, D. (1977). Maximum likelihood from incomplete data via the EM algorithm. Journal of the Royal Statistical Society, Series B,39, 1-38.

Peel, D. & McLachlan, G. J. (2000). Robust mixture modelling using the t distribution. Statistics and Computing,10(4), 339-348.

Karlsson, M. & Laitila, T. (2014). Finite mixture modeling of censored regression models. Statistical Papers,55(3), 627-642.

Basso,R.M.,Lachos,V.H.,Cabral,C.R.B. & Ghosh,P. (2010). Robust mixture modeling based on scale mixtures of skew-normal distributions. Computational Statistics & Data Analysis, 54(12), 2926-2941.

Basford, K., Greenway, D.,McLachlan,G. & Peel,D. (1997). Standard errors of fitted component means of normal mixtures. Computational Statistics,12, 1-18.
}

\keyword{package}

\seealso{
\code{\link{CensMixReg}},\code{\link{wage.rates}}
}

\examples{
#See examples for the CensMixReg function linked above.
}
\keyword{Mixture Regression}
