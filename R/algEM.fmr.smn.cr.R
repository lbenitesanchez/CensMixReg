#-----------------------------------------------------------------------------------------------------------------------#
#Nesta terceira funcao esta se considerando xij ao inves do xi e considerando nu1=nu2=....=nuG
algEM.fmr.smn.cr <- function(cc, y, x, Abetas = NULL, sigma2 = NULL, pii = NULL, nu=NULL, g = NULL, family = "Normal", error = 0.00001, iter.max = 100)
{
  if(ncol(as.matrix(y)) > 1) stop("This function is only for univariate response y!")
  if((family != "T") && (family != "Normal") && (family != "NormalC") && (family != "Slash")) stop(paste("Family",family,"not recognized.\n",sep=" "))
  if((length(g) == 0) && ((length(sigma2)==0) ||  (length(pii)==0) || (ncol(Abetas)==0) ))  stop("The model is not specified correctly.\n")
  if((length(g)!= 0) && (g < 1)) stop("g must be greater than 0.\n")

  p    <- ncol(x)
  n    <- length(y)

  Lim1 <- y
  Lim2 <- rep(Inf,n)

  ################################################################################
  ###                                     Normal
  ################################################################################
  initial_values   <- c(Abetas,sigma2,pii)
  if (family == "Normal")
  {
    start.time     <- Sys.time()
    mu             <- matrix(0,n,g)
    for(k in 1:g){mu[,k]<- x%*%Abetas[,k]}

    criterio       <- 1
    count          <- 0
    lk = lk1 = lk2 <- sum(log(d.mixedN(cc, y, pii, mu, sigma2)))## log-likelihood

    while((criterio > error) && (count <= iter.max))
    {
      count  <- count + 1
      #print(count)

      tal          <- matrix(0, n, g)
      soma1        <- matrix(0, p,1)
      soma2        <- matrix(0, p, p)

      for (j in 1:g)
      {
        ### E-step: calculando os momentos
        NCensEUY   <- NCensurEsperUY(y,mu[,j],sigma2[j],nu=NULL,0,type="Normal")

        u0         <- NCensEUY$EUY0
        u1         <- NCensEUY$EUY1
        u2         <- NCensEUY$EUY2

        #aux1<-MomN(mu[,j],sigma2[j],y)
        CensEUY    <- CensEsperUY1(mu[cc==1,j],sigma2=sigma2[j],nu=0,delta=0,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type="Normal", cens="2")
        u0[cc==1]  <- CensEUY$EUY0
        u1[cc==1]  <- CensEUY$EUY1
        u2[cc==1]  <- CensEUY$EUY2

        d1         <- dNormal(cc, y, mu[,j], sigma2[j])

        if(length(which(d1 == 0)) > 0) {d1[which(d1 == 0)] <- .Machine$double.xmin}
        d2         <- d.mixedN(cc, y, pii, mu, sigma2)
        if(length(which(d2 == 0)) > 0) {d2[which(d2 == 0)] <- .Machine$double.xmin}

        tal[,j]    <- d1*pii[j]/d2

        ### M-step: atualizar os parametros ###
        pii[j]     <- (1/n)*sum(tal[,j])
        sigma2[j]  <- sum(tal[,j]*(u2-2*u1*mu[,j]+u0*mu[,j]^2))/sum(tal[,j])
        soma1      <- t(x)%*%diag(tal[,j])%*%u1
        soma2      <- t(x)%*%diag(tal[,j]*u0)%*%x
        Abetas[,j] <- solve(soma2)%*%soma1
        mu[,j]     <- x%*%Abetas[,j]
      } #End ---> for (j in 1:g)

      pii[g]       <- 1 - (sum(pii) - pii[g])
      zero.pos     <- NULL
      zero.pos     <- which(pii == 0)
      if(length(zero.pos) != 0)
      {
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }

      if (pii[1]< 0.5 & g==2)
      {
        mu         <- cbind(mu[,2],mu[,1])
        Abetas     <- cbind(Abetas[,2], Abetas[,1])
        pii        <- as.vector(c(pii[2], pii[1]))
        sigma2     <- as.vector(c(sigma2[2], sigma2[1]))
      }

      auxlog       <- d.mixedN(cc, y, pii, mu, sigma2)
      if(length(which(auxlog == 0)) > 0) auxlog[which(auxlog == 0)] <- .Machine$double.xmin

      lk3          <- sum(log(auxlog))

      if(count<2){criterio <- abs(lk2 - lk3)/abs(lk3)
      }else{
        tmp        <- (lk3 - lk2)/(lk2 - lk1)
        tmp2       <- lk2 + (lk3 - lk2)/(1-tmp)
        criterio   <- abs(tmp2 - lk3)#; print(criterio)
      }

      lk2          <- lk3
    } #End ---> while((criterio > error) && (count <= iter.max))
    lk   <- lk2
    #if (criteria == TRUE)
    #{
      cl <- apply(tal, 1, which.max)
      #    icl <- 0
      #    for (j in 1:g) icl<-icl+sum(log(pii[j]*dNormal(cc, y, mu[,j], sigma2[j])))
    #}
    end.time       <- Sys.time()
    time.taken     <- end.time - start.time
    EP              <- im.fmr.smn.cr(cc, y,x,Abetas,sigma2,pii,nu,family)$EP
    print(Abetas)
    print(class(Abetas))
    print(pii)
    print(sigma2)
    parameters      <- cbind(c(do.call(c, as.list(Abetas) ),pii[1:(g-1)],sigma2))
    table           <- data.frame(parameters,EP)
    colnames(table) <- c("Estimate","Std. Error")

    k <- 0
    namesrowAbetas   <- c()
    for(j in 1:g)
      for(i in 1:length(Abetas[,j]))
      {
        k <- k + 1
        namesrowAbetas[k] <- paste("beta",i-1,j,sep="")
      }
    namesrowSigmas   <- c(); for(i in 1:g){namesrowSigmas[i] <- paste("sigma",i,sep="")}
    namesrowPii      <- c(); for(i in 1:(g-1)){namesrowPii[i]  <- paste("pii",i,sep="")}
    rownames(table) <- c(namesrowAbetas,namesrowPii,namesrowSigmas)
  }


  ################################################################################
  ###                                     Student-t
  ################################################################################

  if (family == "T")
  {
    start.time     <- Sys.time()
    ERRO           <- 1e-6
    TOLERANCIA     <- 1e-6
    MAX_NU         <- 5
    MIN_NU         <- 1.01
    mu             <- matrix(0,n,g)

    for (k in 1:g){mu[,k]<- x%*%Abetas[,k]}

    criterio       <- 1
    count          <- 0

    lk = lk1 = lk2 <- sum(log(d.mixedT(cc, y, pii, mu, sigma2,nu)))## log-likelihood

    while((criterio > error) && (count <= iter.max))
    {
      count        <- count + 1
      #print(count)

      tal          <- matrix(0, n, g)
      soma1        <- matrix(0, p,1)
      soma2        <- matrix(0, p, p)

      for(j in 1:g)
      {
        ### E-step: calculando os momentos
        NCensEUY   <- NCensurEsperUY(y,mu[,j],sigma2[j],nu,0,type="T")
        u0         <- NCensEUY$EUY0
        u1         <- NCensEUY$EUY1
        u2         <- NCensEUY$EUY2

        CensEUY    <- CensEsperUY1(mu[cc==1,j],sigma2=sigma2[j],nu=nu,delta=0,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type="T", cens="2")
        u0[cc==1]  <- CensEUY$EUY0
        u1[cc==1]  <- CensEUY$EUY1
        u2[cc==1]  <- CensEUY$EUY2

        d1         <- dT(cc, y, mu[,j], sigma2[j],nu)
        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
        d2         <- d.mixedT(cc, y, pii, mu, sigma2,nu)

        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

        tal[,j]    <- d1*pii[j]/d2

        ### M-step: atualizar os parametros ###
        pii[j]     <- (1/n)*sum(tal[,j])
        sigma2[j]  <- sum(tal[,j]*(u2-2*u1*mu[,j]+u0*mu[,j]^2)) /sum(tal[,j])
        soma1      <- t(x)%*%diag(tal[,j])%*%u1
        soma2      <- t(x)%*%diag(c(u0)*tal[,j])%*%x

        Abetas[,j] <- solve(soma2)%*%soma1
        mu[,j]     <- x%*%Abetas[,j]
      }#End ---> for (j in 1:g)

      pii[g]       <- 1 - (sum(pii) - pii[g])
      zero.pos     <- NULL
      zero.pos     <- which(pii == 0)
      if(length(zero.pos) != 0)
      {
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }


      if (pii[1]< 0.5 & g==2)
      {
        mu         <- cbind(mu[,2],mu[,1])
        Abetas     <- cbind(Abetas[,2], Abetas[,1])
        pii        <- as.vector(c(pii[2], pii[1]))
        sigma2     <- as.vector(c(sigma2[2], sigma2[1]))
      }

      ft           <- function(nu)sum(log(d.mixedT(cc, y, pii, mu, sigma2,nu)))
      nu           <- optimize(f=ft, interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,maximum=TRUE,tol=TOLERANCIA)$maximum

      auxlog       <- d.mixedT(cc, y, pii, mu, sigma2,nu)
      if(length(which(auxlog == 0)) > 0) auxlog[which(auxlog == 0)] <- .Machine$double.xmin

      lk3          <- sum(log(auxlog))

      if(count<2){criterio <- abs(lk2 - lk3)/abs(lk3)
      }else{
        tmp        <- (lk3 - lk2)/(lk2 - lk1)
        tmp2       <- lk2 + (lk3 - lk2)/(1-tmp)
        criterio   <- abs(tmp2 - lk3)#; print(criterio)
      }

      lk2          <- lk3
    }#End ---> while((criterio > error) && (count <= iter.max))
    lk             <- lk2
    #if (criteria == TRUE)
    #{
      cl <- apply(tal, 1, which.max)
      #   icl <- 0
      #   for (j in 1:g) icl<-icl+sum(log(pii[j]*dT(cc, y, mu[,j], sigma2[j],nu)))
    #}
    end.time       <- Sys.time()
    time.taken     <- end.time - start.time

    EP              <- im.fmr.smn.cr(cc, y,x,Abetas,sigma2,pii,nu,family)$EP
    parameters      <- cbind(c(do.call(c, as.list(Abetas)),pii[1:(g-1)],sigma2,nu))
    table           <- data.frame(parameters,c(EP,0))
    colnames(table) <- c("Estimate","Std. Error")

    k <- 0
    namesrowAbetas   <- c()
    for(j in 1:g)
      for(i in 1:length(Abetas[,j]))
      {
        k <- k + 1
        namesrowAbetas[k] <- paste("beta",i-1,j,sep="")
      }
    namesrowSigmas   <- c(); for(i in 1:g){namesrowSigmas[i] <- paste("sigma",i,sep="")}
    namesrowPii      <- c(); for(i in 1:(g-1)){namesrowPii[i]  <- paste("pii",i,sep="")}
    namesrowNu       <- paste("nu",sep="")
    rownames(table)  <- c(namesrowAbetas,namesrowPii,namesrowSigmas,namesrowNu)
  }


  ################################################################################
  ###                                     Slash
  ################################################################################

  if (family == "Slash")
  {
    start.time     <- Sys.time()
    ERRO           <- 1e-6
    TOLERANCIA     <- 1e-6
    MAX_NU         <- 5
    MIN_NU         <- 1.01
    mu             <- matrix(0,n,g)

    for (k in 1:g){mu[,k]<- x%*%Abetas[,k]}

    criterio       <- 1
    count          <- 0

    lk = lk1 = lk2 <- sum(log(d.mixedSL(cc, y, pii, mu, sigma2,nu)))## log-likelihood

    while((criterio > error) && (count <= iter.max))
    {
      count <- count + 1
      #print(count)

      tal          <- matrix(0, n, g)
      soma1        <- matrix(0, p,1)
      soma2        <- matrix(0, p, p)

      for (j in 1:g)
      {
        j  ### E-step: calculando os momentos
        NCensEUY   <- NCensurEsperUY(as.numeric(y),mu[,j],sigma2[j],nu,0,type="Slash")
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
        sigma2[j]  <- sum(tal[,j]*(u2-2*u1*mu[,j]+u0*mu[,j]^2)) /sum(tal[,j])
        soma1      <- t(x)%*%diag(tal[,j])%*%u1
        soma2      <- t(x)%*%diag(c(u0)*tal[,j])%*%x

        Abetas[,j] <- solve(soma2)%*%soma1
        mu[,j]     <- x%*%Abetas[,j]
      }

      pii[g]       <- 1 - (sum(pii) - pii[g])
      zero.pos     <- NULL
      zero.pos     <- which(pii == 0)
      if(length(zero.pos) != 0)
      {
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }

      if(pii[1]< 0.5 & g==2)
      {
        mu         <- cbind(mu[,2],mu[,1])
        Abetas     <- cbind(Abetas[,2], Abetas[,1])
        pii        <- as.vector(c(pii[2], pii[1]))
        sigma2     <- as.vector(c(sigma2[2], sigma2[1]))
      }

      ft           <- function(nu)sum(log(d.mixedSL(cc, y, pii, mu, sigma2,nu)))
      nu           <- optimize(f=ft, interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,maximum=TRUE,tol=TOLERANCIA)$maximum

      auxlog       <- d.mixedSL(cc, y, pii, mu, sigma2,nu)
      if(length(which(auxlog == 0)) > 0) auxlog[which(auxlog == 0)] <- .Machine$double.xmin

      lk3          <- sum(log(auxlog))

      if(count<2){criterio <- abs(lk2 - lk3)/abs(lk3)
      }else {
        tmp        <- (lk3 - lk2)/(lk2 - lk1)
        tmp2       <- lk2 + (lk3 - lk2)/(1-tmp)
        criterio   <- abs(tmp2 - lk3)#; print(criterio)
      }

      lk2          <- lk3

    }
    lk             <- lk2

    #if (criteria == TRUE){
      cl <- apply(tal, 1, which.max)
      #icl <- 0
      #for (j in 1:g) icl<-icl+sum(log(pii[j]*dSL(cc, y, mu[,j], sigma2[j],nu)))
    #}
    end.time       <- Sys.time()
    time.taken     <- end.time - start.time

    EP              <- im.fmr.smn.cr(cc, y,x,Abetas,sigma2,pii,nu,family)$EP;print(EP)
    parameters      <- cbind(c(do.call(c, as.list(Abetas) ),pii[1:(g-1)],sigma2,nu))
    table           <- data.frame(parameters,c(EP,0))
    colnames(table) <- c("Estimate","Std. Error")

    k <- 0
    namesrowAbetas   <- c()
    for(j in 1:g)
      for(i in 1:length(Abetas[,j]))
      {
        k <- k + 1
        namesrowAbetas[k] <- paste("beta",i-1,j,sep="")
      }
    namesrowSigmas   <- c(); for(i in 1:g){namesrowSigmas[i] <- paste("sigma",i,sep="")}
    namesrowPii      <- c(); for(i in 1:(g-1)){namesrowPii[i]  <- paste("pii",i,sep="")}
    namesrowNu       <- paste("nu",sep="")
    rownames(table)  <- c(namesrowAbetas,namesrowPii,namesrowSigmas,namesrowNu)
  }

  ################################################################################
  ###                                     Normal-Contaminada
  ################################################################################

  if (family == "NormalC")
  {
    start.time     <- Sys.time()
    nu             <-  c(0.1,0.1)
    mu             <- matrix(0,n,g)

    for (k in 1:g){mu[,k]<- x%*%Abetas[,k]}

    criterio       <- 1
    count          <- 0

    lk = lk1 = lk2 <- sum(log(d.mixedCN(cc, y, pii, mu, sigma2,nu)))## log-likelihood

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
        sigma2[j]  <- sum(tal[,j]*(u2-2*u1*mu[,j]+u0*mu[,j]^2))/sum(tal[,j])
        soma1      <- t(x)%*%diag(tal[,j])%*%u1
        soma2      <- t(x)%*%diag(c(u0)*tal[,j])%*%x

        Abetas[,j] <- solve(soma2)%*%soma1
        mu[,j]     <- x%*%Abetas[,j]
      }

      pii[g]       <- 1 - (sum(pii) - pii[g])
      zero.pos     <- NULL
      zero.pos     <- which(pii == 0)
      if(length(zero.pos) != 0)
      {
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }


      if (pii[1]< 0.5 & g==2)
      {
        mu         <- cbind(mu[,2],mu[,1])
        Abetas     <- cbind(Abetas[,2], Abetas[,1])
        pii        <- as.vector(c(pii[2], pii[1]))
        sigma2     <- as.vector(c(sigma2[2], sigma2[1]))
      }

      ft <- function(nu)
      {
        a1         <- exp(nu[1])/(1+exp(nu[1]))
        a2         <- exp(nu[2])/(1+exp(nu[2]))
        resp       <- -sum(log(d.mixedCN(cc, y, pii, mu, sigma2,c(a1,a2))))
        return(resp)
      }

      #Art          <- optim(c(0.1,0.1),ft, method=c("BFGS"),control=list(maxit=20000))$par
      #nu3          <- min(round(exp(Art[1])/(1+exp(Art[1]))+0.05,1),0.9)
      #nu4          <- min(round(exp(Art[2])/(1+exp(Art[2]))+0.05,1),0.9)
      #nu           <- c(nu3,nu4)

      ft2          <- function(nu)sum(log(d.mixedCN(cc, y, pii, mu, sigma2,nu)))
      nu           <- optim(nu, ft2, control = list(fnscale = -1), method = "L-BFGS-B", lower = rep(0.01, 2), upper = rep(0.99,2))$par
      #print(nuu)
      auxlog       <- d.mixedCN(cc, y, pii, mu, sigma2,nu)
      if(length(which(auxlog == 0)) > 0) auxlog[which(auxlog == 0)] <- .Machine$double.xmin

      lk3          <- sum(log(auxlog))

      if(count<2){criterio <- abs(lk2 - lk3)/abs(lk3)
      }else {
        tmp        <- (lk3 - lk2)/(lk2 - lk1)
        tmp2       <- lk2 + (lk3 - lk2)/(1-tmp)
        criterio   <- abs(tmp2 - lk3)#; print(criterio)
      }

      lk2          <- lk3
    }
    lk             <- lk2
    #if (criteria == TRUE)
    #{
      cl     <- apply(tal, 1, which.max)
      #  icl <- 0
      #for (j in 1:g) icl<-icl+sum(log(pii[j]*dCN(cc, y, mu[,j], sigma2[j],nu)))
    #}
    #nu <- nuu
    end.time       <- Sys.time()
    time.taken     <- end.time - start.time

    EP              <- im.fmr.smn.cr(cc, y,x,Abetas,sigma2,pii,nu,family)$EP
    parameters      <- cbind(c(do.call(c, as.list(Abetas) ),pii[1:(g-1)],sigma2,nu))
    table           <- data.frame(parameters,c(EP,0,0))
    colnames(table) <- c("Estimate","Std. Error")

    k <- 0
    namesrowAbetas   <- c()
    for(j in 1:g)
      for(i in 1:length(Abetas[[j]]))
      {
        k <- k + 1
        namesrowAbetas[k] <- paste("beta",i-1,j,sep="")
      }
    namesrowSigmas   <- c(); for(i in 1:g){namesrowSigmas[i] <- paste("sigma",i,sep="")}
    namesrowPii      <- c(); for(i in 1:(g-1)){namesrowPii[i]  <- paste("pii",i,sep="")}
    namesrowNu       <- paste("nu",sep="")
    namesrowGama     <- paste("gamma",sep="")
    rownames(table)  <- c(namesrowAbetas,namesrowPii,namesrowSigmas,namesrowNu,namesrowGama)

  }

  if(family == "Normal")  d <- sum(p) + g + (g-1)    #Abeta + Sigma + pi
  if(family == "NormalC") d <- sum(p) + g + (g-1)+2  #Abeta + Sigma + pi+ (nu,gamma)
  if(family == "T")       d <- sum(p) + g + (g-1)+ 1 #Abeta + sigma +pi +nu
  if(family == "Slash")   d <- sum(p) + g + (g-1)+ 1 #Abeta + sigma +pi +nu
  aic <- -2*lk + 2*d
  bic <- -2*lk + log(n)*d
  edc <- -2*lk + 0.2*sqrt(n)*d
  #icl <- -2*icl + log(n)*d
  obj.out <- list(EP=EP,criterio=criterio,table=table,Abetas = Abetas, sigma2 = sigma2, pii = pii, nu=nu, lk=lk, aic = aic, bic = bic, edc = edc,  iter = count, n = length(y),time = time.taken,initial_values=initial_values,  convergence = criterio < error, group = cl)

  class(obj.out) <- family
  obj.out
}
#-----------------------------------------------------------------------------------------------------------------------#

