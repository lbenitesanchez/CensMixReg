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

#Load the data
 data(wage.rates)

#Set the response y and covariate x
 y    <- wage.rates$wage
 x1   <- cbind(1,wage.rates$age,wage.rates$educ,wage.rates$hours/1000)
 cc   <- c(rep(0,428),rep(1,325))

########################################################################
#Example for regression modelling of censored data based on
#Mixtures of Scale Mixtures of Normal (SMN) distributions
########################################################################
#Obtain the initial values
 initial    <- initial.values.fm.smn.cr(cc, y,x1,g=2,algorithm="k-means"
 ,family="T",lower=1,upper=20,space=1,plotLog = TRUE,searchNU=TRUE,
 printNU=FALSE, saveFigure = FALSE)

##Fits a left mixture censored Student-t model to the data
 fitT <- fm.smn.cr(cc, y, x1, Abetas = initial$Abetas, medj = initial$medj
 , sigma2 = initial$sigma2, pii = initial$pii, nu=initial$nu, g = 2,
 family = "T", error = 0.0001, iter.max = 500)

##Fits a left mixture censored Normal model to the data
 fitN <- fm.smn.cr(cc, y, x1, Abetas = initial$Abetas, medj = initial$medj
 , sigma2 = initial$sigma2, pii = initial$pii, nu=initial$nu, g = 2,
 family = "Normal", error = 0.0001, iter.max = 500)

######################################################################
#Example for finite mixture of regression models for censored data
#based on scale mixtures of modelling of censored data based on
#Mixtures of SMN distributions
######################################################################
#Obtain the initial values
 initial    <- initial.values.fmr.smn.cr(cc, y,x1,g=2,algorithm="k-means"
 ,family="T",lower=1,upper=20,space=1,plotLog = TRUE,searchNU=TRUE,
 printNU=FALSE, saveFigure = FALSE)

##Fits a left mixture censored Student-t model to the data
 fitT <- fmr.smn.cr( cc, y, x1, Abetas = initial$Abetas,
 sigma2 = initial$sigma2, pii = initial$pii, nu=initial$nu,
 g = 2, family = "T", error = 10^-4, iter.max = 500)

##Fits a left mixture censored Normal model to the data
 fitN <- fmr.smn.cr(cc, y, x1, Abetas = initial$Abetas,
 sigma2 = initial$sigma2, pii = initial$pii, nu=initial$nu,
 g = 2, family = "Normal", error = 10^-4, iter.max = 500)

}




\keyword{datasets}



