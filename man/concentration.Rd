\name{concentration}
\alias{concentration}
\docType{data}
\title{
Concentration levels
}
\description{
The dataset corresponds to concentration levels of certain dissolved trace metals in freshwater streams across the Commonwealth of Virginia. This dataset consist of the concentration levels of the dissolved trace metals copper (Cu), lead (Pb), zinc (Zn), calcium (Ca) and magnesium (Mg) from 184 independent randomly selected sites in freshwater streams across Virginia. The Cu, Pb, and Zn concentrations are reported in ug/L of water, whereas Ca and Mg concentration are suitably reported in mg/L of water. Since the measurements are taken at different times, the presence of multiple limit of detection values are possible for each trace metal (VDEQ (2003)). The limit of detection but Cu and Pb is the 0.1ug/L,
1.0mg/L for Zn while Ca and Mg has limit of 0.5mg/L and 1.0mg/L.
The percentage of left-censored values of 2.7\% for (Ca), 4.9\% for (Cu), 9.8\% for (Mg) are
small in comparison to 78.3\% for (Pb) and 38.6\% for (Zn), also note that 17.9\% of the streams
had 0 nondetected trace metals, 39.1\% had 1, 37.0\% had 2, 3.8\% had 3, 1.1\% had 4 and 1.1\%
had 5.
}

\usage{data(concentration)}
\format{\code{concentration} is a data frame with 184 cases (rows) with 5 variables (columns).
}
\details{
For a complete description of data concentration levels see VDEQ (2003).}
\encoding{latin1}
\source{
Hoffman, H. & Johnson, R. (2014). Pseudo-likelihood estimation of multivariate normal parameters in the presence of left-censored data. Journal of Agricultural, Biological, and Environmental Statistics, pages 156-171.

VDEQ (2003). The quality of virginia non-tidal streams: First year report. VDEQ Technical Bulletin WQA/2002-2001, Office of Water Quality and Assessments, Virginia Department of Environmental Quality, pages 13-16.
}

\examples{
\dontrun{
data(concentration)
p  <- 5
y  <- as.matrix(concentration[,1:p])
cc <- as.matrix(concentration[,(p+1):(2*p)])

fitT <- CensMmix(cc, y, nu=3, mu=NULL, Sigma = NULL,  pii = NULL, g = 2,
                  get.init = TRUE, criteria = TRUE, group = TRUE,
                  family = "t", error = 0.00001, iter.max = 350,
                  uni.Sigma = TRUE, obs.prob= FALSE, kmeans.param = NULL)

fitN <- CensMmix(cc, y, nu=3, mu=NULL, Sigma = NULL,  pii = NULL, g = 2,
                  get.init = TRUE, criteria = TRUE, group = TRUE,
                  family = "Normal", error = 0.00001, iter.max = 350,
                  uni.Sigma = TRUE, obs.prob= FALSE, kmeans.param = NULL)

#Standard error (SE)
fitT$res$SE
fitN$res$SE

}
}
\keyword{datasets}
