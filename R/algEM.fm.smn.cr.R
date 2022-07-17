CensMixRegEM <- function(cc, y, x1, Abetas = NULL, medj= NULL, sigma2 = NULL, pii = NULL, nu=NULL, g = NULL, family = "Normal", error = 0.00001, iter.max = 100, aitken = TRUE)
{

  if(ncol(as.matrix(y)) > 1) stop("This function is only for univariate response y!")
  if((family != "T") && (family != "Normal") && (family != "CN") && (family != "Slash")) stop(paste("Family",family,"not recognized.\n",sep=" "))
  if((length(g) == 0) && ((length(sigma2)==0) ||  (length(pii)==0) || (ncol(Abetas)==0) ))  stop("The model is not specified correctly.\n")
  if((length(g)!= 0) && (g < 1)) stop("g must be greater than 0.\n")

  p    <- ncol(x1)
  n    <- length(y)

  Lim1 <- y #Esta sendo considerado o censura tipo 2
  Lim2 <- rep(Inf,n)

  ################################################################################
  ###                                     Normal
  ################################################################################
  if (family == "Normal")
  {
    start.time  <- Sys.time()

    n            <- length(y)
    p            <- ncol(x1)-1
    x            <- as.matrix(x1[,2:(p+1)])
    beta0        <- Abetas[1] ## intercepto
    betas        <- as.matrix(Abetas[2:(p+1)])   ### parameters of regression with dimension "p"
    varphi       <- rep(0,g)  ### alphas
    mu           <- mu1 <- matrix(0,n,g)

    for (k in 1:g)
    {
      varphi[k]  <- beta0+medj[k]
      mu1[,k]    <- x%*%betas
      mu[,k]     <- mu1[,k]+varphi[k]
    }

    teta         <- c(Abetas, medj, sigma2, pii);#print(Abetas)

    criterio     <- 1
    count        <- 0


    lk = lk1 = lk2    <- sum(log(d.mixedN(cc, y, pii, mu, sigma2))) #log-likelihood
    teta1 <- teta2    <- c(Abetas, medj, sigma2, pii)

    while((criterio > error) && (count <= iter.max))
    {
      count  <- count + 1
      #print(count)
      tal     <- matrix(0, n, g)
      soma1   <- matrix(0, p,1)
      soma2   <- matrix(0, p, p)

      for (j in 1:g)
      {
        ### E-step: calculando os momentos

        NCensEUY  <- NCensurEsperUY(y,mu[,j],sigma2[j],nu=NULL,0,type="Normal")
        u0        <- NCensEUY$EUY0
        u1        <- NCensEUY$EUY1
        u2        <- NCensEUY$EUY2
        #aux1<-MomN(mu[,j],sigma2[j],y)

        CensEUY   <- CensEsperUY1(mu[cc==1,j],sigma2=sigma2[j],nu=0,delta=0,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type="Normal", cens="2")
        u0[cc==1] <- CensEUY$EUY0
        u1[cc==1] <- CensEUY$EUY1
        u2[cc==1] <- CensEUY$EUY2

        d1         <- dNormal(cc, y, mu[,j], sigma2[j])
        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin

        d2         <- d.mixedN(cc, y, pii, mu, sigma2)
        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

        tal[,j]    <- d1*pii[j]/d2
        ### M-step: atualizar os parametros ###

        pii[j]     <- (1/n)*sum(tal[,j])
        sigma2[j]  <- sum(tal[,j]*(u2+u0*varphi[j]^2+u0*mu1[,j]^2-2*u1*varphi[j]-2*u1*mu1[,j]+2*u0*varphi[j]*mu1[,j])) / sum(tal[,j])
        varphi[j]  <- sum(tal[,j]*(u1-u0*mu1[,j]))/sum(u0*tal[,j])
        soma1      <- soma1 + t(x)%*%diag(tal[,j]/sigma2[j])%*%(u1-u0*varphi[j])
        soma2      <- soma2 + t(x)%*%diag(tal[,j]*c(u0)/sigma2[j])%*%x
      }

      betas         <- solve(soma2)%*%soma1

      pii[g]        <- 1 - (sum(pii) - pii[g])

      zero.pos      <- NULL
      zero.pos      <- which(pii == 0)
      if(length(zero.pos) != 0)
      {
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }

      #if (pii[1]< 0.5 && g==2)
      #{
      #  medj        <- as.vector(c(medj[2],medj[1]))
      #  varphi      <- as.vector(c(varphi[2], varphi[1]))
      #  pii         <- as.vector(c(pii[2], pii[1]))
      #  sigma2      <- as.vector(c(sigma2[2], sigma2[1]))
      #}

      beta0         <- sum(pii*varphi)

      for (k in 1:g)
      {
        mu1[,k]     <- x%*%betas
        mu[,k]      <- mu1[,k]+varphi[k]
        medj[k]     <- varphi[k]-beta0
      }

      Abetas        <- c(beta0,betas)#;print(betas)

      auxlog        <- d.mixedN(cc, y, pii, mu, sigma2)

      if(length(which(auxlog == 0)) > 0) auxlog[which(auxlog == 0)] <- .Machine$double.xmin
      ###########
      teta3 <- c(Abetas, medj, sigma2, pii)
      #if(count>2) r = abs(norm((teta3 - teta2), type="2")/norm((teta2 - teta1), type="2"))
      teta1 <- teta2
      teta2 <- teta3
      ###########
      if(aitken==TRUE)
      {
        lk3         <- sum(log(auxlog))

        if(count<2){ criterio <- abs(lk2 - lk3)/abs(lk3)
        }else {
          c         <- (lk3 - lk2)/(lk2 - lk1)
          tmp2      <- lk2 + (lk3 - lk2)/(1-c)#;print(tmp)
          criterio  <- abs(tmp2 - lk3)
        }
        lk2         <- lk3
      }else{
        lk1        <- sum(log(auxlog))
        criterio   <- abs(lk1/lk2-1)
        lk2        <- lk1
      }
      #print(criterio)
    } #End while and the estimation process
    lk         <- lk2

    end.time        <- Sys.time()
    time.taken      <- end.time - start.time
  }

  ################################################################################
  ###                       Mixture Student-t                                  ###
  ################################################################################

  if (family == "T")
  {
    start.time  <- Sys.time() #Begin Time
    #Parameters of function optimize for estimate nu (degrees of freedom)
    ERRO        <- 1e-6
    TOLERANCIA  <- 1e-6
    MAX_NU      <- 20
    MIN_NU      <- 1.01

    n           <- length(y)
    p           <- ncol(x1)-1
    x           <- as.matrix(x1[,2:(p+1)])
    beta0       <- Abetas[1] ## intercepto
    betas       <- as.matrix(Abetas[2:(p+1)])   ### parameters of regression with dimension "p"
    varphi      <- rep(0,g)  ### alphas
    mu          <- mu1 <- matrix(0,n,g)

    for (k in 1:g)
    {
      varphi[k]  <- beta0+medj[k]
      mu1[,k]    <- x%*%betas
      mu[,k]     <- mu1[,k]+varphi[k]
    }

    teta         <- c(Abetas, medj, sigma2, pii,nu);#print(Abetas)

    criterio     <- 1
    count        <- 0


    lk = lk1 = lk2     <- sum(log(d.mixedT(cc, y, pii, mu, sigma2,nu))) #log-likelihood
    teta1 <- teta2     <- c(Abetas, medj, sigma2, pii,nu)

    while((criterio > error) && (count <= iter.max))
    {#Begin while
      count        <- count + 1
      #print(count)

      tal          <- matrix(0, n, g)
      soma1        <- matrix(0, p,1)
      soma2        <- matrix(0, p, p)

      for (j in 1:g) #loop for each components
      {#Beginf FOR
        ### E-step: computing the moments
        NCensEUY    <- NCensurEsperUY(y,mu[,j],sigma2[j],nu,0,type="T")
        u0          <- NCensEUY$EUY0
        u1          <- NCensEUY$EUY1
        u2          <- NCensEUY$EUY2

        CensEUY     <- CensEsperUY1(mu[cc==1,j],sigma2=sigma2[j],nu=nu,delta=0,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type="T", cens="2")
        u0[cc==1]   <- CensEUY$EUY0
        u1[cc==1]   <- CensEUY$EUY1
        u2[cc==1]   <- CensEUY$EUY2

        d1          <- dT(cc, y, mu[,j], sigma2[j],nu)
        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
        d2          <- d.mixedT(cc, y, pii, mu, sigma2,nu)

        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin
        tal[,j]     <- d1*pii[j]/d2

        ### M-step: update the parameters ###
        pii[j]     <- (1/n)*sum(tal[,j])
        sigma2[j]  <- sum(tal[,j]*(u2+u0*varphi[j]^2+u0*mu1[,j]^2-2*u1*varphi[j]-2*u1*mu1[,j]+2*u0*varphi[j]*mu1[,j])) / sum(tal[,j])
        varphi[j]  <- sum(tal[,j]*(u1-u0*mu1[,j]))/sum(u0*tal[,j])
        soma1      <- soma1 + t(x)%*%diag(tal[,j]/sigma2[j])%*%(u1-u0*varphi[j])
        soma2      <- soma2 + t(x)%*%diag(tal[,j]*c(u0)/sigma2[j])%*%x
      }#End FOR

      betas         <- solve(soma2)%*%soma1

      pii[g]        <- 1 - (sum(pii) - pii[g])

      zero.pos      <- NULL
      zero.pos      <- which(pii == 0)
      if(length(zero.pos) != 0)
      {
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }

      #if (pii[1]< 0.5 && g==2)
      #{
      #  medj        <- as.vector(c(medj[2],medj[1]))
      #  varphi      <- as.vector(c(varphi[2], varphi[1]))
      #  pii         <- as.vector(c(pii[2], pii[1]))
      #  sigma2      <- as.vector(c(sigma2[2], sigma2[1]))
      #}

      beta0         <- sum(pii*varphi)

      for (k in 1:g)
      {
        mu1[,k]     <- x%*%betas
        mu[,k]      <- mu1[,k]+varphi[k]
        medj[k]     <- varphi[k]-beta0
      }

      Abetas        <- c(beta0,betas)#;print(betas)

      ft            <- function(nu) sum(log(d.mixedT(cc, y, pii, mu, sigma2,nu)))
      #nu            <- optimize(f=ft, interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,maximum=TRUE,tol=TOLERANCIA)$maximum
      nu            <- optim(nu, ft, control = list(fnscale = -1), method = "L-BFGS-B", lower = rep(MIN_NU, 1), upper = rep(MAX_NU,1))$par

      teta          <- c(Abetas, medj, sigma2, pii,nu)
      auxlog        <- d.mixedT(cc, y, pii, mu, sigma2,nu)

      if(length(which(auxlog == 0)) > 0) auxlog[which(auxlog == 0)] <- .Machine$double.xmin
      ###########
      teta3 <- c(Abetas, medj, sigma2, pii,nu)
      #if(count>2) r = abs(norm((teta3 - teta2), type="2")/norm((teta2 - teta1), type="2"))
      teta1 <- teta2
      teta2 <- teta3
      ###########
      if(aitken==TRUE)
      {
        lk3         <- sum(log(auxlog))

        if(count<2){ criterio <- abs(lk2 - lk3)/abs(lk3)
        }else {
          c         <- (lk3 - lk2)/(lk2 - lk1)
          tmp2      <- lk2 + (lk3 - lk2)/(1-c)#;print(tmp)
          criterio  <- abs(tmp2 - lk3)
        }
        lk2         <- lk3
      }else{
        lk1         <- sum(log(auxlog))
        criterio    <- abs(lk1/lk2-1)
        lk2         <- lk1
      }
      #print(criterio)
    } #End while and the estimation process
    lk              <- lk2
    end.time        <- Sys.time()
    time.taken      <- end.time - start.time
  }


  ################################################################################
  ###                                     Slash                                ###
  ################################################################################

  if (family == "Slash")
  {
    start.time  <- Sys.time()

    ERRO        <- 1e-8
    TOLERANCIA  <- 1e-8
    MAX_NU      <- 20
    MIN_NU      <- 1.01

    n           <- length(y)
    p           <- ncol(x1)-1
    x           <- as.matrix(x1[,2:(p+1)])
    beta0       <- Abetas[1] ## intercepto
    betas       <- as.matrix(Abetas[2:(p+1)])   ### parameters of regression with dimension "p"
    varphi      <- rep(0,g)  ### alphas
    mu          <- mu1 <- matrix(0,n,g)

    for (k in 1:g)
    {
      varphi[k]  <- beta0+medj[k]
      mu1[,k]    <- x%*%betas
      mu[,k]     <- mu1[,k]+varphi[k]
    }

    teta        <- c(Abetas, medj, sigma2, pii,nu)

    criterio    <- 1
    count       <- 0

    lk = lk1 = lk2     <- sum(log(d.mixedSL(cc, y, pii, mu, sigma2,nu))) #log-likelihood
    teta1 <- teta2 <-  c(Abetas, medj, sigma2, pii,nu)

    while((criterio > error) && (count <= iter.max))
    {
      count       <- count + 1
      #print(count)

      tal         <- matrix(0, n, g)
      soma1       <- matrix(0, p,1)
      soma2       <- matrix(0, p, p)

      for (j in 1:g)
      {
        ### E-step: calculando os momentos
        NCensEUY   <- NCensurEsperUY(y,mu[,j],sigma2[j],nu,0,type="Slash")
        u0         <- NCensEUY$EUY0
        u1         <- NCensEUY$EUY1
        u2         <- NCensEUY$EUY2

        CensEUY    <- CensEsperUY1(mu[cc==1,j],sigma2=sigma2[j],nu=nu,delta=0,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type="Slash", cens="2")
        u0[cc==1]  <- CensEUY$EUY0
        u1[cc==1]  <- CensEUY$EUY1
        u2[cc==1]  <- CensEUY$EUY2

        d1         <- dSL(cc, y, mu[,j], sigma2[j],nu)

        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
        d2         <- d.mixedSL(cc, y, pii, mu, sigma2,nu)
        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

        tal[,j]    <- d1*pii[j]/d2

        ### M-step: atualizar os parametros ###

        pii[j]     <- (1/n)*sum(tal[,j])
        sigma2[j]  <- sum(tal[,j]*(u2+u0*varphi[j]^2+u0*mu1[,j]^2-2*u1*varphi[j]-2*u1*mu1[,j]+2*u0*varphi[j]*mu1[,j])) / sum(tal[,j])
        varphi[j]  <- sum(tal[,j]*(u1-u0*mu1[,j]))/sum(u0*tal[,j])
        soma1      <- soma1 + t(x)%*%diag(tal[,j]/sigma2[j])%*%(u1-u0*varphi[j])
        soma2      <- soma2 + t(x)%*%diag(tal[,j]*c(u0)/sigma2[j])%*%x
        #        Abetas[,j] <- solve(soma2)%*%soma1
        #        mu[,j]<- x%*%Abetas[,j]
      }

      betas         <- solve(soma2)%*%soma1

      pii[g]        <- 1 - (sum(pii) - pii[g])

      zero.pos      <- NULL
      zero.pos      <- which(pii == 0)
      if(length(zero.pos) != 0)
      {
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }

      #if (pii[1]< 0.5 && g==2)
      #{
      #  medj        <- as.vector(c(medj[2],medj[1]))
      #  varphi      <- as.vector(c(varphi[2], varphi[1]))
      #  pii         <- as.vector(c(pii[2], pii[1]))
      #  sigma2      <- as.vector(c(sigma2[2], sigma2[1]))
      #  # nu         <- as.vector(c(nu[2], nu[1]))
      #}

      beta0         <- sum(pii*varphi)

      for (k in 1:g)
      {
        mu1[,k]     <- x%*%betas
        mu[,k]      <- mu1[,k]+varphi[k]
        medj[k]     <- varphi[k]-beta0
      }

      Abetas        <- c(beta0,betas)#;print(betas)
      ft            <- function(nu)sum(log(d.mixedSL(cc, y, pii, mu, sigma2,nu)))
      #nu            <- optimize(f=ft, interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,maximum=TRUE,tol=TOLERANCIA)$maximum
      nu            <- optim(nu, ft, control = list(fnscale = -1), method = "L-BFGS-B", lower = rep(MIN_NU, 1), upper = rep(MAX_NU,1))$par

      #    nu  <- nlminb(start = nu,  ft,
      #           gradient = NULL,
      #           hessian  = NULL,
      #           scale    = 1,
      #           control  = list(),
      #           lower    = rep(MIN_NU, g),
      #           upper    = rep(MAX_NU,g))$par

      # print(nu)

      teta          <- c(Abetas, medj, sigma2, pii,nu)
      auxlog        <- d.mixedSL(cc, y, pii, mu, sigma2,nu)

      if(length(which(auxlog == 0)) > 0) {auxlog[which(auxlog == 0)] <- .Machine$double.xmin}

      ###########
      teta3 <- c(Abetas, medj, sigma2, pii,nu)
      #if(count>2) r = abs(norm((teta3 - teta2), type="2")/norm((teta2 - teta1), type="2"))
      teta1 <- teta2
      teta2 <- teta3
      ###########
      if(aitken==TRUE)
      {
        lk3         <- sum(log(auxlog))

        if(count<2){ criterio <- abs(lk2 - lk3)/abs(lk3)
        }else {
          c         <- (lk3 - lk2)/(lk2 - lk1)
          tmp2      <- lk2 + (lk3 - lk2)/(1-c)
          criterio  <- abs(tmp2 - lk3)#;print(criterio)
        }
        lk2         <- lk3
      }else{
        lk1        <- sum(log(auxlog))
        criterio   <- abs(lk1/lk2-1)
        lk2         <- lk1
      }
      #print(criterio)
    } #End while and the estimation process
    lk         <- lk2

    end.time        <- Sys.time()
    time.taken      <- end.time - start.time
  }

  ################################################################################
  ###                                     Slash                                ###
  ################################################################################

  if (family == "NormalC")
  {
    start.time  <- Sys.time()

    ERRO        <- 1e-8
    TOLERANCIA  <- 1e-8
    MAX_NU      <- 20
    MIN_NU      <- 1.01

    n           <- length(y)
    p           <- ncol(x1)-1
    x           <- as.matrix(x1[,2:(p+1)])
    beta0       <- Abetas[1] ## intercepto
    betas       <- as.matrix(Abetas[2:(p+1)])   ### parameters of regression with dimension "p"
    varphi      <- rep(0,g)  ### alphas
    mu          <- mu1 <- matrix(0,n,g)

    for (k in 1:g)
    {
      varphi[k] <- beta0+medj[k]
      mu1[,k]   <- x%*%betas
      mu[,k]    <- mu1[,k]+varphi[k]
    }

    teta        <- c(Abetas, medj, sigma2, pii,nu)

    criterio    <- 1
    count       <- 0

    lk = lk1 = lk2     <- sum(log(d.mixedCN(cc, y, pii, mu, sigma2,nu))) #log-likelihood
    teta1 <- teta2 <-  c(Abetas, medj, sigma2, pii,nu)

    while((criterio > error) && (count <= iter.max))
    {
      count        <- count + 1
      #print(count)

      tal          <- matrix(0, n, g)
      soma1        <- matrix(0, p,1)
      soma2        <- matrix(0, p, p)

      for (j in 1:g)
      {
        ### E-step: calculando os momentos
        NCensEUY   <- NCensurEsperUY(y,mu[,j],sigma2[j],nu,0,type="NormalC")
        u0         <- NCensEUY$EUY0
        u1         <- NCensEUY$EUY1
        u2         <- NCensEUY$EUY2

        CensEUY    <- CensEsperUY1(mu[cc==1,j],sigma2=sigma2[j],nu=nu,delta=0,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type="NormalC", cens="2")
        u0[cc==1]  <- CensEUY$EUY0
        u1[cc==1]  <- CensEUY$EUY1
        u2[cc==1]  <- CensEUY$EUY2

        d1         <- dCN(cc, y, mu[,j], sigma2[j],nu)

        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
        d2         <- d.mixedCN(cc, y, pii, mu, sigma2,nu)
        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

        tal[,j]    <- d1*pii[j]/d2

        ### M-step: atualizar os parametros ###

        pii[j]     <- (1/n)*sum(tal[,j])
        sigma2[j]  <- sum(tal[,j]*(u2+u0*varphi[j]^2+u0*mu1[,j]^2-2*u1*varphi[j]-2*u1*mu1[,j]+2*u0*varphi[j]*mu1[,j])) / sum(tal[,j])
        varphi[j]  <- sum(tal[,j]*(u1-u0*mu1[,j]))/sum(u0*tal[,j])
        soma1      <- soma1 + t(x)%*%diag(tal[,j]/sigma2[j])%*%(u1-u0*varphi[j])
        soma2      <- soma2 + t(x)%*%diag(tal[,j]*c(u0)/sigma2[j])%*%x
      }

      betas        <- solve(soma2)%*%soma1

      pii[g]       <- 1 - (sum(pii) - pii[g])

      zero.pos     <- NULL
      zero.pos     <- which(pii == 0)
      if(length(zero.pos) != 0)
      {
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }

      #if (pii[1]< 0.5 && g==2)
      #{
      #  medj        <- as.vector(c(medj[2],medj[1]))
      #  varphi      <- as.vector(c(varphi[2], varphi[1]))
      #  pii         <- as.vector(c(pii[2], pii[1]))
      #  sigma2      <- as.vector(c(sigma2[2], sigma2[1]))
      #  # nu         <- as.vector(c(nu[2], nu[1]))
      #}

      beta0         <- sum(pii*varphi)

      for (k in 1:g)
      {
        mu1[,k]     <- x%*%betas
        mu[,k]      <- mu1[,k]+varphi[k]
        medj[k]     <- varphi[k]-beta0
      }

      Abetas        <- c(beta0,betas)#;print(betas)
      ft            <- function(nu)sum(log(d.mixedCN(cc, y, pii, mu, sigma2,nu)))
      #nu            <- optimize(f=ft, interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,maximum=TRUE,tol=TOLERANCIA)$maximum
      nu            <- optim(nu, ft, control = list(fnscale = -1), method = "L-BFGS-B", lower = rep(MIN_NU, 1), upper = rep(MAX_NU,1))$par

      #    nu  <- nlminb(start = nu,  ft,
      #           gradient = NULL,
      #           hessian  = NULL,
      #           scale    = 1,
      #           control  = list(),
      #           lower    = rep(MIN_NU, g),
      #           upper    = rep(MAX_NU,g))$par

      # print(nu)

      teta          <- c(Abetas, medj, sigma2, pii,nu)
      auxlog        <- d.mixedCN(cc, y, pii, mu, sigma2,nu)

      if(length(which(auxlog == 0)) > 0) {auxlog[which(auxlog == 0)] <- .Machine$double.xmin}

      ###########
      teta3 <- c(Abetas, medj, sigma2, pii,nu)
      #if(count>2) r = abs(norm((teta3 - teta2), type="2")/norm((teta2 - teta1), type="2"))
      teta1 <- teta2
      teta2 <- teta3
      ###########
      if(aitken==TRUE)
      {
        lk3         <- sum(log(auxlog))

        if(count<2){ criterio <- abs(lk2 - lk3)/abs(lk3)
        }else {
          c         <- (lk3 - lk2)/(lk2 - lk1)
          tmp2      <- lk2 + (lk3 - lk2)/(1-c)
          criterio  <- abs(tmp2 - lk3)#;print(criterio)
        }
        lk2         <- lk3
      }else{
        lk1        <- sum(log(auxlog))
        criterio   <- abs(lk1/lk2-1)
        lk2         <- lk1
      }
      #print(criterio)
    } #End while and the estimation process
    lk         <- lk2

    end.time        <- Sys.time()
    time.taken      <- end.time - start.time
  }


  ################################################################################
  if(family=="Normal")
  {
    #Row names of ttable
    namesrowAbetas <- c(); for(i in 1:length(Abetas)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
    namesrowMu     <- c(); for(i in 1:g){namesrowMu[i]<- paste("mu",i,sep="")}
    namesrowSigmas <- c(); for(i in 1:g){namesrowSigmas[i]<- paste("sigma",i,sep="")}
    namesrowPii    <- c(); for(i in 1:(g-1)){namesrowPii[i]<- paste("p",i,sep="")}
    ttable         <- data.frame(c(Abetas, medj, sigma2, pii[1:(g-1)]))
    rownames(ttable) <- c(namesrowAbetas,namesrowMu,namesrowSigmas,namesrowPii)
    colnames(ttable) <- c("Estimate")
  }

  if(family=="T" || family=="Slash")
  {
    #Row names of ttable
    namesrowAbetas <- c(); for(i in 1:length(Abetas)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
    namesrowMu     <- c(); for(i in 1:g){namesrowMu[i]<- paste("mu",i,sep="")}
    namesrowSigmas <- c(); for(i in 1:g){namesrowSigmas[i]<- paste("sigma",i,sep="")}
    namesrowPii    <- c(); for(i in 1:(g-1)){namesrowPii[i]<- paste("p",i,sep="")}
    ttable         <- data.frame(c(Abetas, medj, sigma2, pii[1:(g-1)],nu))

    rownames(ttable) <- c(namesrowAbetas,namesrowMu,namesrowSigmas,namesrowPii,"nu")
    colnames(ttable) <- c("Estimate")
  }

  if(family=="NormalC")
  {
    #Row names of ttable
    namesrowAbetas <- c(); for(i in 1:length(Abetas)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
    namesrowMu     <- c(); for(i in 1:g){namesrowMu[i]<- paste("mu",i,sep="")}
    namesrowSigmas <- c(); for(i in 1:g){namesrowSigmas[i]<- paste("sigma",i,sep="")}
    namesrowPii    <- c(); for(i in 1:(g-1)){namesrowPii[i]<- paste("p",i,sep="")}
    ttable         <- data.frame(c(Abetas, medj, sigma2, pii[1:(g-1)],nu[1],nu[2]))

    rownames(ttable) <- c(namesrowAbetas,namesrowMu,namesrowSigmas,namesrowPii,"nu","gamma")
    colnames(ttable) <- c("Estimate")
  }


  if(family == "Normal") d <- g*p + g + (g-1) #Abeta + Sigma + pi
  if(family == "T")      d <- g*p + g + (g-1) + 1 #Abeta + sigma + pi + nu
  if(family == "Slash")  d <- g*p + g + (g-1) + 1 #Abeta + sigma + pi + nu
  if(family == "NormalC")d <- g*p + g + (g-1) + 1  + 1 #Abeta + sigma + pi + nu + gamma

  aic   <- -2*lk + 2*d
  bic   <- -2*lk + log(n)*d
  edc   <- -2*lk + 0.2*sqrt(n)*d
  aic_c <- -2*lk + 2*n*d/(n-d-1)
  abic  <- -2*lk + d*log((n+2)/24)
  #icl <- -2*icl + log(n)*d
  obj.out <- list(time = time.taken, group = apply(tal, 1, which.max), Abetas = Abetas, medj=medj, sigma2 = sigma2, pii = pii, nu=nu, ttable=ttable, aic = aic, bic = bic, edc = edc,  aic_c=aic_c, abic=abic, loglik=lk, iter = count, n = length(y), u0=u0,u1=u1,u2=u2,convergence=criterio < error)

  class(obj.out) <- family
  obj.out
}

#CensMixRegEM(cc, y, x1, Abetas = NULL, medj= NULL, sigma2 = NULL, pii = NULL, nu=NULL, g = NULL, family = "Normal", error = 0.00001, iter.max = 100, aitken = TRUE)

