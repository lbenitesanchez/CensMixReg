\name{imm.fm.smn.cr}
\alias{imm.fm.smn.cr}
\title{Information matrix}

\description{
Calculate the information matrix of returned analysis based on the
model \code{family} choice (univariate case, p=1).
}
\usage{
imm.fm.smn.cr(cc, y,x1,model)
}
\arguments{
\item{y}{Vector of responses in case of right censoring.}
\item{x1}{Matrix or vector of covariates.}
\item{cc}{Vector of censoring indicators. For each observation: 0 if non-censored, 1 if censored.}
\item{model}{ a variable returned by \code{\link{fm.smn.cr}} }
}

\value{
  Estimate the Information Matrix of the parameters.
}

\examples{ ## see \code{\link{wage.rates}}
}
\author{Luis Benites Sanchez \email{lbenitesanchez@gmail.com}, Victor Hugo Lachos \email{hlachos@ime.unicamp.br}, Edgar J. Lopez Moreno \email{edgar.javier.lopez.moreno@gmail.com}
}

\seealso{ \code{\link{fm.smn.cr}}}

\keyword{models}
