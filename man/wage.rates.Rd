\name{wage.rates}
\alias{wage.rates}
\docType{data}
\title{Wage Rates of 753 Women}
\description{Wage rates of 753 married white women with left censoring.
}
\usage{data(wage.rates)}
\format{
A data frame with 753 observed wage rates of married white women in 1975. Some wage rates are set equal to zero, this means that those wives did not work in 1975, therefore, these observations are considered left censored at zero.

 \describe{
    \item{inlf}{
  =1 if in labor force, 1975
    }
    \item{hours}{
 hours worked, 1975
    }
    \item{kidslt6}{
 # kids < 6 years
    }
    \item{kidsge6}{
# kids 6-18
    }
    \item{age}{
 woman's age in yrs
    }
    \item{educ}{
 years of schooling
    }
    \item{wage}{
 estimated wage from earns., hours
    }
    \item{huseduc}{
husband's years of schooling
    }
    \item{faminc}{
  family income, 1975
    }
    \item{motheduc}{
 mother's years of schooling
    }
    \item{fatheduc}{
 father's years of schooling
    }
    \item{unem}{
  unem. rate in county of resid.
    }
    \item{city}{
  =1 if live in SMSA
    }
    \item{exper}{
  actual labor mkt exper
    }
    \item{nwifeinc}{
  (faminc - wage*hours)/1000
    }
    \item{expersq}{
 exper^2
    }
 }

}

\encoding{latin1}

\source{
Mroz, T.A. 1987. "The sensitivity of an empirical model of married women's hours of work to economic and statistical assumptions".
}


\examples{
\dontrun{
##Load the data
 library(CensRegMod)
 data(wage.rates)

##Set the response y and covariate x

#Set the response y and covariate x
 y    <- wage.rates$wage
 x1   <- cbind(1,wage.rates$age,wage.rates$educ,wage.rates$hours/1000)
 cc   <- c(rep(0,428),rep(1,325))

##Fits a left mixture censored Student-t model to the data
 fitT <- CensMixReg(cc, -y, -x1, Abetas = NULL, medj= NULL, sigma2 = NULL, nu=NULL, pii = NULL,
 g = 2, get.init = TRUE, criteria = TRUE,  group = FALSE, family = "T", error = 0.00001,
 iter.max = 500, obs.prob= FALSE, kmeans.param = NULL)

 #########################################
 #A simple output example

 ---------------------------------------------------------
           Censored linear mixture regression model
---------------------------------------------------------

Observations = 753
-----------
Estimates
-----------

        Estimate
beta0    5.69541
beta1   -0.00958
beta2    0.41823
beta3    2.87485
mu1      0.18454
mu2    -22.24469
sigma1   3.80318
sigma2   8.51759
pii1     0.99177
nu       3.17652
------------------------
Model selection criteria
------------------------

         Loglik      AIC      BIC      EDC
Value -1239.528 2499.057 2545.298 2533.939
-------
Details
-------

Iterations = 482
Processing time = 1.009154 mins



##Fits a left mixture censored Normal model to the data
 fitN <- CensMixReg(cc, -y, -x1, Abetas = NULL, medj= NULL, sigma2 = NULL, nu=NULL, pii = NULL,
 g = 2, get.init = TRUE, criteria = TRUE, group = FALSE, family = "Normal", error = 0.00001,
 iter.max = 500, obs.prob= FALSE, kmeans.param = NULL)

#Standard error (SE) when random errors follow a finite mixture of Student-t distributions.
 IFmixRegT <- imm.fmcr(y,x1,cc, fitT)
 SEmixT    <- sqrt(diag(solve(IFmixRegT$IM)))

#Standard error when random errors follow a finite mixture of normal distributions.
 IFmixRegN <- imm.fmcr(y,x1,cc, fitN)
 SEmixN    <- sqrt(diag(solve(IFmixRegN$IM)))
 }
}


\keyword{datasets}



