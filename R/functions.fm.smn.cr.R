##########################################################################################
#A funcao initial.Values permite obter os chutes iniciais para os parametros
#betas_j, sigma_j e nu (aqui estamos supongo que nu=nu1=...=nuG)
#A obtencao do chute inicial de nu e por medio de uma grade
#Aqui e considerando xi
##########################################################################################

initial.values.fm.smn.cr <- function(cc, y,x,g=2,algorithm="k-medoids",family="T",lower=1,upper=20,space=0.1,plotLog = TRUE,searchNU=TRUE,printNU=TRUE, saveFigure = FALSE)
{
  p <- ncol(x)
  n <- length(y)

  ######################################################################################
  #Segue os tres tipos de algoritmos de particao que podem ser utilizados para obter
  #os chutes iniciais para beta_j, sigma_j e pii_j

  if(algorithm == "trim-kmeans")
  {
    init           <- trimkmeans(y,k=g,trim=0.1,runs=3)
    cluster        <- init$classification[init$classification!=(g+1)]

    if(g==1)
    {
      pii          <- c(sum(cluster==1))/length(y)
      size         <- c(sum(cluster==1))
      sigma2       <- c()
      Abetas       <- matrix(0,p,g)

      for (j in 1:g)
      {
        dd         <- cluster
        Abetas[,j] <- solve(t(x[dd==j,])%*%x[dd==j,])%*%t(x[dd==j,])%*%y[dd==j]
        yr         <- y[dd==j]-x[dd==j,]%*%Abetas[,j]
        sigma2[j]  <- sum(yr^2)/size[j]
      }
    }

    if(g==2)
    {
      pii          <- c(sum(cluster==1),sum(cluster==2))/length(y)
      size         <- c(sum(cluster==1),sum(cluster==2))
      sigma2       <- c()
      Abetas       <- matrix(0,p,g)

      for (j in 1:g)
      {
        dd         <- cluster
        Abetas[,j] <- solve(t(x[dd==j,])%*%x[dd==j,])%*%t(x[dd==j,])%*%y[dd==j]
        yr         <- y[dd==j]-x[dd==j,]%*%Abetas[,j]
        sigma2[j]  <- sum(yr^2)/size[j]
      }
    }

    if(g==3)
    {
      pii          <- c(sum(cluster==1),sum(cluster==2),sum(cluster==3))/length(y)
      size         <- c(sum(cluster==1),sum(cluster==2),sum(cluster==3))
      sigma2       <- c()
      Abetas       <- matrix(0,p,g)

      for (j in 1:g)
      {
        dd         <- cluster
        Abetas[,j] <- solve(t(x[dd==j,])%*%x[dd==j,])%*%t(x[dd==j,])%*%y[dd==j]
        yr         <- y[dd==j]-x[dd==j,]%*%Abetas[,j]
        sigma2[j]  <- sum(yr^2)/size[j]
      }
    }

    if(g==4)
    {
      pii          <- c(sum(cluster==1),sum(cluster==2),sum(cluster==3),sum(cluster==4))/length(y)
      size         <- c(sum(cluster==1),sum(cluster==2),sum(cluster==3),sum(cluster==4))
      sigma2       <- c()
      Abetas       <- matrix(0,p,g)

      for (j in 1:g)
      {
        dd         <- cluster
        Abetas[,j] <- solve(t(x[dd==j,])%*%x[dd==j,])%*%t(x[dd==j,])%*%y[dd==j]
        yr         <- y[dd==j]-x[dd==j,]%*%Abetas[,j]
        sigma2[j]  <- sum(yr^2)/size[j]
      }
    }

    initial_values <- c(Abetas,sigma2,pii)
  }

  if(algorithm == "MinMax_kmeans")
  {
    N            <- nrow(matrix(y))#;print(N)

    k            <- g    #number of clusters
    p_init       <- 0    #initial p
    p_max        <- 0.5  #maximum p
    p_step       <- 0.01 #p step
    t_max        <- 500  #maximum number of iterations
    beta         <- 0.3  #amount of memory for the weight updates
    Restarts     <- 10   #number of MinMax k-means restarts

    tmp          <- sample(1:N,replace=F)
    M            <- cbind(y[tmp[1:k],])

    init         <- MinMax_kmeans(matrix(y),M,k,p_init,p_max,p_step,t_max,beta)
    if(g==1) {size <- c(sum(init$Cluster_elem==1))}
    if(g==2) {size <- c(sum(init$Cluster_elem==1),sum(init$Cluster_elem==2))}
    if(g==3) {size <- c(sum(init$Cluster_elem==1),sum(init$Cluster_elem==2),sum(init$Cluster_elem==3))}
    if(g==4) {size <- c(sum(init$Cluster_elem==1),sum(init$Cluster_elem==2),sum(init$Cluster_elem==3),sum(init$Cluster_elem==4))}
    pii          <- size/N
    sigma2       <- c()
    Abetas       <- matrix(0,p,g)

    for (j in 1:g)
    {
      dd         <- init$Cluster_elem
      Abetas[,j] <- solve(t(x[dd==j,])%*%x[dd==j,])%*%t(x[dd==j,])%*%y[dd==j]
      yr         <- y[dd==j]-x[dd==j,]%*%Abetas[,j]
      sigma2[j]  <- sum(yr^2)/size[j]
    }

    initial_values <- c(Abetas,sigma2,pii)
  }

  if(algorithm == "k-means")
  {
    if(length(g) == 0) stop("g is not specified correctly.\n")

    if(g > 1)
    {
      init         <- kmeans(y,g,iter.max=50,nstart=1,algorithm="Hartigan-Wong")
      pii          <- init$size/length(y)
      sigma2       <- c()
      Abetas       <- solve(t(x)%*%x)%*%t(x)%*%y
      dd           <- init$cluster
      yr           <- y-x%*%Abetas
      sigma2       <- sum(yr^2)/init$size
      medj         <- as.vector(init$centers)
    }else{
      Abetas       <- solve(t(x)%*%x)%*%t(x)%*%y
      yr           <- y-x%*%Abetas
      sigma2       <- var(yr)
      medj         <- mean(yr)
      pii          <- 1
    }
    initial_values <- c(Abetas,medj,sigma2,pii)#;print(initial_values)
  }

  if (algorithm == "k-medoids")
  {
    if(length(g) == 0) {stop("g is not specified correctly.\n")}

    if(g > 1){
      Abetas       <- matrix(0,p)
      init         <- Cluster_Medoids(as.matrix(y),clusters=g, distance_metric="euclidean", swap_phase = TRUE, fuzzy = TRUE,seed=sample(1:100000,1))
      pii          <- init$medoid_indices/length(y)
      Abetas       <- solve(t(x)%*%x)%*%t(x)%*%y
      yr           <- y-x%*%Abetas

      sigma2       <- c()
      medj         <- c()
      for (j in 1:g)
      {
        dd         <- init$clusters
        sigma2[j]  <- sum(yr^2)/init$medoid_indices[j]
        medj[j]    <- mean(yr[dd=j])
      }
    } else{
      Abetas       <- solve(t(x)%*%x)%*%t(x)%*%y
      yr           <- y-x%*%Abetas
      sigma2       <- var(yr)
      medj         <- mean(yr)
      pii          <- 1
    }

    initial_values  <- c(Abetas,medj,sigma2,pii)#;print(initial_values)
  }
  ########################################################################################################################
  #-----------------------------------------------------------------------------------------------------------------------#
  #    As linhas de codigo de abaixo sao para poder obter o chute inicial para os graus de liberdade para a T, SL e CN
  #-----------------------------------------------------------------------------------------------------------------------#
  ########################################################################################################################

  if(searchNU ==TRUE)
  {
    mu             <- matrix(0,n,g)
    for(k in 1:g){mu[,k]<- medj[k] + x%*%Abetas}

    if(family == "T")
    {
      nu           <- seq(lower,upper,space)
      vero         <- c()
      for(i in 1:length(nu))
      {
        vero[i]    <- sum(log(d.mixedT(cc, y, pii, mu, sigma2,nu[i])))## log-likelihood
        if(printNU==TRUE)
          cat(i,"to",length(nu),"\n")
      }
      m           <- cbind(nu,vero)
      maxi        <- max(vero)
      for(i in 1:length(vero))
      {
        if(maxi==m[i,2])
          nu0     <- m[i,1]
      }

      if(plotLog == TRUE && saveFigure == TRUE)
      {
        postscript("plotTnu.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
        plot(nu,vero,type="l",xlab=expression(paste(nu)),ylab="Log-likelihood")
        abline(v=nu0,lty=2)
        eq <- bquote(bold(nu[max] == .(nu0)))
        mtext(eq,side=3,cex=1.5)
        dev.off() #Fechando o dispositivo potscript
      }

      if(plotLog == TRUE  && saveFigure == FALSE)
      {
        plot(nu,vero,type="l",xlab=expression(paste(nu)),ylab="Log-likelihood")
        abline(v=nu0,lty=2)
        eq <- bquote(bold(nu[max] == .(nu0)))
        mtext(eq,side=3,cex=1.5)
      }

    }

    if(family == "Slash")
    {
      nu           <- seq(lower,upper,space)
      vero         <- c()
      for(i in 1:length(nu))
      {
        vero[i]    <- sum(log(d.mixedSL(cc, y, pii, mu, sigma2,nu[i])))## log-likelihood
        if(printNU==TRUE)
          cat(i,"to",length(nu),"\n")
      }
      m            <- cbind(nu,vero)
      maxi         <- max(vero)
      for(i in 1:length(vero))
      {
        if(maxi==m[i,2])
          nu0      <- m[i,1]
      }

      if(plotLog == TRUE && saveFigure == TRUE)
      {
        postscript("plotSlashnu.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
        plot(nu,vero,type="l",xlab=expression(paste(nu)),ylab="Log-likelihood")
        abline(v=nu0,lty=2)
        eq <- bquote(bold(nu[max] == .(nu0)))
        mtext(eq,side=3,cex=1.5)
        dev.off() #Fechando o dispositivo potscript
      }

      if(plotLog == TRUE  && saveFigure == FALSE)
      {
        plot(nu,vero,type="l",xlab=expression(paste(nu)),ylab="Log-likelihood")
        abline(v=nu0,lty=2)
        eq <- bquote(bold(nu[max] == .(nu0)))
        mtext(eq,side=3,cex=1.5)
      }
    }

    if(family == "NormalC")
    {
      nuu          <- seq(lower,upper,space)
      comp         <- expand.grid(nuu,nuu)
      nu           <- cbind(comp[,1],comp[,2])
      vero         <- c()
      for(i in 1:nrow(nu))
      {
        vero[i]    <- sum(log(d.mixedCN(cc, y, pii, mu, sigma2,nu[i,])))
        if(printNU==TRUE)
          cat(i,"to",nrow(nu),"\n")

      }
      m       <- cbind(nu,vero)
      maxi    <- max(vero)
      nu1     <- 0
      nu2     <- 0
      for(i in 1:length(vero))
      {
        if(maxi==m[i,3])
        {
          nu1 <- m[i,1]
          nu2 <- m[i,2]
        }
      }
      nu0     <- c(nu1,nu2)

      if(plotLog == TRUE && saveFigure == TRUE)
      {
        z <- matrix(0,nrow=length(nuu),ncol=length(nuu))
        for(i in 1:length(nuu))
          for(j in 1:length(nuu))
            z[i,j] <- sum(log(d.mixedCN(cc, y, pii, mu, sigma2,c(nuu[i],nuu[j]))))
          postscript("persp.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
          res1 <- persp(x=nuu, y=nuu, z,theta = 95, phi = 15,xlab="nu",ylab="gamma")
          dev.off() #Fechando o dispositivo potscript

          postscript("NuGamma.eps", width=5.75, height=5.75, horizontal=FALSE, onefile=TRUE)
          plot(nu[,1],nu[,2],type="n",ylab=expression(paste(gamma)),xlab=expression(paste(nu)),ylim=c(0,1.2))
          abline(v=nu1,lty=2)
          abline(h=nu2,lty=2)
          legend('topleft',legend= c(as.expression(bquote(nu[max] == .(nu1))), as.expression(bquote(gamma[max] == .(nu2)))),bty="n",cex=1.5)
          dev.off() #Fechando o dispositivo potscript
      }

      if(plotLog == TRUE && saveFigure == FALSE)
      {
        z <- matrix(0,nrow=length(nuu),ncol=length(nuu))
        for(i in 1:length(nuu))
          for(j in 1:length(nuu))
            z[i,j] <- sum(log(d.mixedCN(cc, y, pii, mu, sigma2,c(nuu[i],nuu[j]))))
          par(mfrow=c(1,2))
          res1 <- persp(x=nuu, y=nuu, z,theta = 95, phi = 15,xlab="nu",ylab="gamma")
          plot(nu[,1],nu[,2],type="n",ylab=expression(paste(gamma)),xlab=expression(paste(nu)),ylim=c(0,1.2))
          abline(v=nu1,lty=2)
          abline(h=nu2,lty=2)
          legend('topright',legend= c(as.expression(bquote(nu[max] == .(nu1))), as.expression(bquote(gamma[max] == .(nu2)))),bty="n",cex=1.5)
          par(mfrow=c(1,1))
      }
    }

    ###############################################################################################################
    #Objetos de saida
    ###############################################################################################################
    if(family == "Normal"){
      nu0     <- 0
      obj.out <- list(Abetas=Abetas,medj=medj,sigma2=sigma2,pii=pii,nu=nu0)
    }else{
      obj.out <- list(vero=vero, Abetas=Abetas,medj=medj,sigma2=sigma2,pii=pii, nu=nu0)
    }
  }else{
    if(family == "Normal"){
      nu0     <- 0
      obj.out <- list(Abetas=Abetas,medj=medj, sigma2=sigma2,pii=pii,nu=nu0)
    }else{
      obj.out <- list(Abetas=Abetas,medj=medj, sigma2=sigma2,pii=pii)
    }

  }
}

############################################################################################

## Densidade das SMN com locacao-escala #######

dSlash    <- function(y,mu,sigma2,nu)
{
  resp     <- z <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y))
  {
    z[i]    <- (y[i]-mu)/sqrt(sigma2)
    f2      <- function(u) nu*u^(nu-0.5)*dnorm(z[i]*sqrt(u))/sqrt(sigma2)
    resp[i] <- integrate(f2,0,1)$value
  }
  return(resp)
}

dNormalC  <- function(y,mu,sigma2,nu)
{
  Acum     <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y))
  {
    eta     <- nu[1]
    gama    <- nu[2]
    Acum[i] <- eta*dnorm(y[i],mu,sqrt(sigma2/gama)) + (1-eta)*dnorm(y[i],mu,sqrt(sigma2))
  }
  return(Acum)
}

################################################################################


MomN <- function(mu,sigma2,a)
{
  z   <- (a-mu)/sqrt(sigma2)
  aux <-(1-pnorm(z))
  if(length(which(aux == 0)) > 0) aux[which(aux == 0)] <- .Machine$double.xmin
  Ey  <- mu+sqrt(sigma2)*dnorm(z)/aux
  Ey2 <- mu^2+sigma2+(dnorm(z)/aux)*(mu+a)*sqrt(sigma2)
  return(list(Ey=Ey,Ey2=Ey2))
}

MomT <- function(mu,sigma2,nu,a)
{
  a1  <- (a-mu)/sqrt(sigma2)
  aux <- 1-pt(a1,nu)
  if(length(which(aux == 0)) > 0){ aux[which(aux == 0)] <- .Machine$double.xmin }
  G1  <- 0.5*(gamma((nu-1)/2)*nu^(nu/2))/(aux*gamma(nu/2)*gamma(1/2))
  Ex  <- G1*(nu+a1^2)^(-(nu-1)/2)
  Ey  <- mu+sqrt(sigma2)*Ex
  Ex2 <- nu/(nu-2)+G1*(a1*(nu+a1^2)^(-(nu-1)/2))
  Ey2 <- mu^2+sigma2*Ex2+2*mu*sqrt(sigma2)*Ex
  return(list(Ey=Ey,Ey2=Ey2))
}

cdfNI  <- function(x,mu,sigma2,nu,type)
{
  resp <- matrix(0,length(x),1)

  if(type=="Normal"){ resp<-pnorm(x,mu,sqrt(sigma2)) }

  if(type=="T")
  {
    z    <- (x-mu)/sqrt(sigma2)
    resp <- pt(z,df=nu)
  }
  return(resp)
}

MomNIT <- function(mu,sigma2,nu,a,type)
{
  if(type=="Normal")
  {
    z   <- (a-mu)/sqrt(sigma2)
    aux <- (1-pnorm(z))
    if(length(which(aux == 0)) > 0) aux[which(aux == 0)] <- .Machine$double.xmin
    Ey  <- mu+sqrt(sigma2)*dnorm(z)/aux
    Ey2 <- mu^2+sigma2+(dnorm(z)/aux)*(mu+a)*sqrt(sigma2)
  }

  if(type=="T")
  {
    a1   <- (a-mu)/sqrt(sigma2)
    aux1 <- (1-cdfNI(a1,0,1,nu,type))
    if(length(which(aux1 == 0)) > 0) {  aux1[which(aux1 == 0)] <- .Machine$double.xmin }
    Aux  <- (1-cdfNI(a1,0,nu/(nu-2),nu-2,type))/aux1
    G1   <- 0.5*(gamma((nu-1)/2)*nu^(nu/2))/(aux1*gamma(nu/2)*gamma(1/2))
    Ex   <- G1*(nu+a1^2)^(-(nu-1)/2)
    Ey   <- mu+sqrt(sigma2)*Ex
    Ex2  <- Aux*(nu/(nu-2))+G1*(a1*(nu+a1^2)^(-(nu-1)/2))
    Ey2  <- mu^2+sigma2*Ex2+2*mu*sqrt(sigma2)*Ex
  }

  return(list(Ey=Ey,Ey2=Ey2))
}


CalMom <- function(mu,sigma2,nu,a,type)
{
  n<-length(mu)

  if(type=="Normal")
  {
    Mom <- MomNIT(mu=mu,sigma2=sigma2,a=a,type=type)
    Ey0 <- rep(1,n)
    Ey1 <- Mom$Ey
    Ey2 <- Mom$Ey2
  }

  if(type=="T")
  {
    sigma2a <- sigma2*nu/(nu+2)
    aux0    <- (1-cdfNI(a,mu,sigma2,nu,type))
    if(length(which(aux0 == 0)) > 0) aux0[which(aux0 == 0)] <- .Machine$double.xmin
    aux1    <- (1-cdfNI(a,mu,sigma2a,nu+2,type))/aux0
    aux2    <- gamma((nu+1)/2)*gamma((nu+2)/2)*(nu+1)/(nu*gamma(nu/2)*gamma((nu+3)/2))
    Ey0     <- aux1*aux2*rep(1,n)
    Mom     <- MomNIT(mu,sigma2a,nu+2,a,type)
    Ey1     <- aux1*aux2*Mom$Ey
    Ey2     <- aux1*aux2*Mom$Ey2
  }

  return(list(Ey0=Ey0,Ey1=Ey1,Ey2=Ey2))
}


#################################################################################################################
imm.fm.smn.cr = function(cc, y,x1,model)
{
  #if((class(model) != "T") && (class(model) != "Normal")) stop(paste("Family",class(model),"not recognized.",sep=" "))
  if(class(model)=="Normal")
  {
    n           <- length(y)

    Lim1        <- y #Esta sendo considerado o censura tipo 2
    Lim2        <- rep(Inf,n)

    p           <- ncol(x1)-1
    g           <- length(model$res$pii)

    Abetas      <- model$res$Abetas
    medj        <- model$res$medj
    sigma2      <- model$res$sigma2
    pii         <- model$res$pii

    beta0       <- Abetas[1]
    betas       <- as.matrix(Abetas[2:(p+1)])   # parameters of regression dimension "p"
    varphi      <- rep(0,g) # beta0 + mu_j
    mu = mu1    <- matrix(0,n,g)
    x           <- as.matrix(x1[,2:(p+1)])

    #u0          <- model$u0
    #u1          <- model$u1
    #u2          <- model$u2

    for (k in 1:g)
    {
      varphi[k] <- beta0+medj[k]
      mu1[,k]   <- x%*%betas
      mu[,k]    <- mu1[,k]+varphi[k]
    }

    tal            <- matrix(0,n,g)
    Sibeta         <- matrix(0,p,n)
    Sivarphi       <- matrix(0,g,n)
    Sisigma        <- matrix(0,g,n)
    Sip            <- matrix(0,g,n)

    for(j in 1:g)
    {
      d            <- ((y-mu[,j])^2)/sigma2[j]

      NCensEUY  <- NCensurEsperUY(y,mu[,j],sigma2[j],nu=NULL,0,type="Normal")
      u0        <- NCensEUY$EUY0
      u1        <- NCensEUY$EUY1
      u2        <- NCensEUY$EUY2
      #aux1<-MomN(mu[,j],sigma2[j],y)

      CensEUY   <- CensEsperUY1(mu[cc==1,j],sigma2=sigma2[j],nu=0,delta=0,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type="Normal", cens="2")
      u0[cc==1] <- CensEUY$EUY0
      u1[cc==1] <- CensEUY$EUY1
      u2[cc==1] <- CensEUY$EUY2

      d1           <- dNormal(cc, y, mu[,j], sigma2[j])
      d2           <- d.mixedN(cc, y, pii, mu, sigma2)
      tal[,j]      <- d1*pii[j]/d2

      Sibeta       <- Sibeta + t(x)%*%diag(tal[,j]*c(u1)/sigma2[j]) - t(x)%*%diag(tal[,j]*c(u0)*(mu1[,j] + varphi[j])/sigma2[j])
      Sip[j,]      <- (1/pii[j])*tal[,j] - (1/pii[g])*tal[,g]
      Sivarphi[j,] <- (1/sigma2[j])*(tal[,j]*c(u1) - tal[,j]*c(u0)*(mu1[,j] + varphi[j]))
      Sisigma[j,]  <- -0.5*(1/sigma2[j]^2)*tal[,j]*(sigma2[j] - c(u2) + 2*c(u1)*(mu1[,j] + varphi[j]) - c(u0)*(mu1[,j] + varphi[j])^2)
    }

    soma      <- matrix(0, length(Abetas)  + 2*g, length(Abetas)  + 2*g)
    si        <- matrix(0, n, (length(Abetas) ) + 2*g)

    for(i in 1:n)
    {
      si[i,]  <- c(Sibeta[,i],Sip[1:(g-1),i],Sivarphi[,i],Sisigma[,i])
      soma    <- soma + cbind(si[i,])%*%si[i,]
    }
    #print(Sinu[,i])

    #print(sqrt(diag(soma)))
    #Jacobian

    p      <- length(betas)
    j11    <- diag(rep(1,p))
    j12    <- matrix(0, nrow=p, ncol=g-1)
    j13    <- matrix(0, nrow=p, ncol=g)
    j14    <- matrix(0, nrow=p, ncol=g)
    J1     <- cbind(j11,j12,j13,j14)

    j21    <- matrix(0, nrow=g-1, ncol=p)
    j22    <- diag(rep(1,g-1))
    j23    <- matrix(0, nrow=g-1, ncol=g)
    j24    <- matrix(0, nrow=g-1, ncol=g)
    J2     <- cbind(j21,j22,j23,j24)

    j31    <- matrix(0, nrow=1, ncol=p)
    j32    <- matrix(0, nrow=1, ncol=g-1)
    j33    <- matrix(1, nrow=1, ncol=g)
    j34    <- matrix(0, nrow=1, ncol=g)
    J3     <- cbind(j31,j32,j33,j34)

    j41    <- matrix(0, nrow=g, ncol=p)
    j42    <- as.matrix(pii^(-1))%*%((medj[g] - medj[1:(g-1)]))
    j43    <- diag(rep(1,g))
    j44    <- matrix(0, nrow=g, ncol=g)
    J4     <- cbind(j41,j42,j43,j44)

    j51    <- matrix(0, nrow=g, ncol=p)
    j52    <- matrix(0, nrow=g, ncol=g-1)
    j53    <- matrix(0, nrow=g, ncol=g)
    j54    <- diag(rep(1,g))
    J5     <- cbind(j51,j52,j53,j54)

    Jacobian       <- rbind(J1,J2,J3,J4,J5)
    IM             <- Jacobian%*%solve(soma)%*%t(Jacobian)
    EP             <- as.matrix(sqrt(diag(Jacobian%*%solve(soma)%*%t(Jacobian))))
  }

  if(class(model)=="T")
  {
    n           <- length(y)
    p           <- ncol(x1)-1
    g           <- length(model$res$pii)

    Abetas      <- model$res$Abetas
    medj        <- model$res$medj
    sigma2      <- model$res$sigma2
    nu          <- model$res$nu
    pii         <- model$res$pii

    beta0       <- Abetas[1]
    betas       <- as.matrix(Abetas[2:(p+1)])   # parameters of regression dimension "p"
    varphi      <- rep(0,g) # beta0 + mu_j
    mu = mu1    <- matrix(0,n,g)
    x           <- as.matrix(x1[,2:(p+1)])

    u0          <- model$res$u0
    u1          <- model$res$u1
    u2          <- model$res$u2

    for (k in 1:g)
    {
      varphi[k] <- beta0+medj[k]
      mu1[,k]   <- x%*%betas
      mu[,k]    <- mu1[,k]+varphi[k]
    }


    tal            <- matrix(0,n,g)
    Sibeta         <- matrix(0,p,n)
    Sivarphi       <- matrix(0,g,n)
    Sisigma        <- matrix(0,g,n)
    Sip            <- matrix(0,g,n)

    for(j in 1:g)
    {
      d            <- ((y-mu[,j])^2)/sigma2[j]
      #u0           <- (nu[j]+1)/(nu[j]+d)
      #u1           <- y*(nu[j]+1)/(nu[j]+d)
      #u2           <- y^2*(nu[j]+1)/(nu[j]+d)

      d1           <- dT(cc, y, mu[,j], sigma2[j],nu)
      d2           <- d.mixedT(cc, y, pii, mu, sigma2,nu)
      tal[,j]      <- d1*pii[j]/d2

      Sibeta       <- Sibeta + t(x)%*%diag(tal[,j]*c(u1)/sigma2[j]) - t(x)%*%diag(tal[,j]*c(u0)*(mu1[,j] + varphi[j])/sigma2[j])
      Sip[j,]      <- (1/pii[j])*tal[,j] - (1/pii[g])*tal[,g]
      Sivarphi[j,] <- (1/sigma2[j])*(tal[,j]*c(u1) - tal[,j]*c(u0)*(mu1[,j] + varphi[j]))
      Sisigma[j,]  <- -0.5*(1/sigma2[j]^2)*tal[,j]*(sigma2[j] - c(u2) + 2*c(u1)*(mu1[,j] + varphi[j]) - c(u0)*(mu1[,j] + varphi[j])^2)
      #Sinu         <- (1/2)*(tal[,j]*(log(nu[j]/2) + 1 - digamma(nu[j]/2) + digamma((nu[j] + 1)/2) - log((nu[j] + d)/2)) - tal[,j]*c(u0))
      #ui           <- rgamma(10000, shape = nu[j]/2, rate = nu[j]/2)
      #montecarlo   <- mean(log(ui)*(1/gamma(nu[j]/2))*(ui/2)^(nu[j]/2)*ui^(nu[j]/2 - 1)*exp(-nu[j]*ui/2))
      #Sinu[j,]     <- tal[,j]*(-digamma(nu[j]/2) + 0.5*(log(nu[j]/2) + 1) + 0.5*(montecarlo - tal[,j]*c(u0)))
    }

    soma      <- matrix(0, length(Abetas)  + 2*g, length(Abetas)  + 2*g)
    si        <- matrix(0, n, (length(Abetas) ) + 2*g)

    for(i in 1:n)
    {
      si[i,]  <- c(Sibeta[,i],Sip[1:(g-1),i],Sivarphi[,i],Sisigma[,i])
      soma    <- soma + cbind(si[i,])%*%si[i,]
    }
    #print(Sinu[,i])

    #print(sqrt(diag(soma)))
    #Jacobian

    p      <- length(betas)
    j11    <- diag(rep(1,p))
    j12    <- matrix(0, nrow=p, ncol=g-1)
    j13    <- matrix(0, nrow=p, ncol=g)
    j14    <- matrix(0, nrow=p, ncol=g)
    #j15    <- matrix(0, nrow=p, ncol=g)
    J1     <- cbind(j11,j12,j13,j14)

    j21    <- matrix(0, nrow=g-1, ncol=p)
    j22    <- diag(rep(1,g-1))
    j23    <- matrix(0, nrow=g-1, ncol=g)
    j24    <- matrix(0, nrow=g-1, ncol=g)
    #j25    <- matrix(0, nrow=g-1, ncol=g)
    J2     <- cbind(j21,j22,j23,j24)

    j31    <- matrix(0, nrow=1, ncol=p)
    j32    <- matrix(0, nrow=1, ncol=g-1)
    j33    <- matrix(1, nrow=1, ncol=g)
    j34    <- matrix(0, nrow=1, ncol=g)
    #j35    <- matrix(0, nrow=1, ncol=g)
    J3     <- cbind(j31,j32,j33,j34)

    j41    <- matrix(0, nrow=g, ncol=p)
    j42    <- as.matrix(pii^(-1))%*%((medj[g] - medj[1:(g-1)]))
    j43    <- diag(rep(1,g))
    j44    <- matrix(0, nrow=g, ncol=g)
    #j45    <- matrix(0, nrow=g, ncol=g)
    J4     <- cbind(j41,j42,j43,j44)

    j51    <- matrix(0, nrow=g, ncol=p)
    j52    <- matrix(0, nrow=g, ncol=g-1)
    j53    <- matrix(0, nrow=g, ncol=g)
    j54    <- diag(rep(1,g))
    #j55    <- matrix(0, nrow=g, ncol=g)
    J5     <- cbind(j51,j52,j53,j54)

    Jacobian       <- rbind(J1,J2,J3,J4,J5)
    IM             <- Jacobian%*%solve(soma)%*%t(Jacobian)
    EP             <- as.matrix(sqrt(diag(Jacobian%*%solve(soma)%*%t(Jacobian))))
  }


  if(class(model)=="Slash")
  {
    n           <- length(y)
    p           <- ncol(x1)-1
    g           <- length(model$res$pii)

    Abetas      <- model$res$Abetas
    medj        <- model$res$medj
    sigma2      <- model$res$sigma2
    nu          <- c(3,3)
    pii         <- model$res$pii

    beta0       <- Abetas[1]
    betas       <- as.matrix(Abetas[2:(p+1)])   # parameters of regression dimension "p"
    varphi      <- rep(0,g) # beta0 + mu_j
    mu = mu1    <- matrix(0,n,g)
    x           <- as.matrix(x1[,2:(p+1)])

    u0          <- model$u0
    u1          <- model$u1
    u2          <- model$u2

    for (k in 1:g)
    {
      varphi[k] <- beta0+medj[k]
      mu1[,k]   <- x%*%betas
      mu[,k]    <- mu1[,k]+varphi[k]
    }


    tal            <- matrix(0,n,g)
    Sibeta         <- matrix(0,p,n)
    Sivarphi       <- matrix(0,g,n)
    Sisigma        <- matrix(0,g,n)
    Sip            <- matrix(0,g,n)
    Sinu           <- matrix(0,g,n)

    for(j in 1:g)
    {
      d            <- ((y-mu[,j])^2)/sigma2[j]
      #u0           <- (nu[j]+1)/(nu[j]+d)
      #u1           <- y*(nu[j]+1)/(nu[j]+d)
      #u2           <- y^2*(nu[j]+1)/(nu[j]+d)

      d1           <- dSL(cc, y, mu[,j], sigma2[j],nu[j])
      d2           <- d.mixedSL(cc, y, pii, mu, sigma2,nu)
      tal[,j]      <- d1*pii[j]/d2

      Sibeta       <- Sibeta + t(x)%*%diag(tal[,j]*c(u1)/sigma2[j]) - t(x)%*%diag(tal[,j]*c(u0)*(mu1[,j] + varphi[j])/sigma2[j])
      Sip[j,]      <- (1/pii[j])*tal[,j] - (1/pii[g])*tal[,g]
      Sivarphi[j,] <- (1/sigma2[j])*(tal[,j]*c(u1) - tal[,j]*c(u0)*(mu1[,j] + varphi[j]))
      Sisigma[j,]  <- -0.5*(1/sigma2[j]^2)*tal[,j]*(sigma2[j] - c(u2) + 2*c(u1)*(mu1[,j] + varphi[j]) - c(u0)*(mu1[,j] + varphi[j])^2)
      #Sinu[j,]      <- (1/2)*(tal[,j]*(log(nu[j]/2) + 1 - digamma(nu[j]/2) + digamma((nu[j] + 1)/2) - log((nu[j] + d)/2)) - tal[,j]*c(u0))
      ui            <- rbeta(10000, nu[j]/2, nu[j]/2)
      montecarlo   <- mean(log(ui)*(1/gamma(nu[j]/2))*(ui/2)^(nu[j]/2)*ui^(nu[j]/2 - 1)*exp(-nu[j]*ui/2))
      Sinu[j,]     <- tal[,j]*(-digamma(nu[j]/2) + 0.5*(log(nu[j]/2) + 1) + 0.5*(montecarlo - tal[,j]*c(u0)))
    }

    soma      <- matrix(0, length(Abetas)  + 3*g, length(Abetas)  + 3*g)
    si        <- matrix(0, n, (length(Abetas) ) + 3*g)

    for(i in 1:n)
    {
      si[i,]  <- c(Sibeta[,i],Sip[1:(g-1),i],Sivarphi[,i],Sisigma[,i],Sinu[,i])
      soma    <- soma + cbind(si[i,])%*%si[i,]
    }
    #print(Sinu[,i])

    #print(sqrt(diag(soma)))
    #Jacobian

    p      <- length(betas)
    j11    <- diag(rep(1,p))
    j12    <- matrix(0, nrow=p, ncol=g-1)
    j13    <- matrix(0, nrow=p, ncol=g)
    j14    <- matrix(0, nrow=p, ncol=g)
    j15    <- matrix(0, nrow=p, ncol=g)
    J1     <- cbind(j11,j12,j13,j14,j15)

    j21    <- matrix(0, nrow=g-1, ncol=p)
    j22    <- diag(rep(1,g-1))
    j23    <- matrix(0, nrow=g-1, ncol=g)
    j24    <- matrix(0, nrow=g-1, ncol=g)
    j25    <- matrix(0, nrow=g-1, ncol=g)
    J2     <- cbind(j21,j22,j23,j24,j25)

    j31    <- matrix(0, nrow=1, ncol=p)
    j32    <- matrix(0, nrow=1, ncol=g-1)
    j33    <- matrix(1, nrow=1, ncol=g)
    j34    <- matrix(0, nrow=1, ncol=g)
    j35    <- matrix(0, nrow=1, ncol=g)
    J3     <- cbind(j31,j32,j33,j34,j35)

    j41    <- matrix(0, nrow=g, ncol=p)
    j42    <- as.matrix(pii^(-1))%*%((medj[g] - medj[1:(g-1)]))
    j43    <- diag(rep(1,g))
    j44    <- matrix(0, nrow=g, ncol=g)
    j45    <- matrix(0, nrow=g, ncol=g)
    J4     <- cbind(j41,j42,j43,j44,j45)

    j51    <- matrix(0, nrow=g, ncol=p)
    j52    <- matrix(0, nrow=g, ncol=g-1)
    j53    <- matrix(0, nrow=g, ncol=g)
    j54    <- diag(rep(1,g))
    j55    <- matrix(0, nrow=g, ncol=g)
    J5     <- cbind(j51,j52,j53,j54,j55)

    j61    <- matrix(0, nrow=g, ncol=p)
    j62    <- matrix(0, nrow=g, ncol=g-1)
    j63    <- matrix(0, nrow=g, ncol=g)
    j64    <- matrix(0, nrow=g, ncol=g)
    j65    <- diag(rep(1,g))
    J6     <- cbind(j61,j62,j63,j64,j65)

    Jacobian       <- rbind(J1,J2,J3,J4,J5,J6)
    IM             <- Jacobian%*%solve(soma)%*%t(Jacobian)
    EP             <- as.matrix(sqrt(diag(Jacobian%*%solve(soma)%*%t(Jacobian))))
  }

  return(list(IM=IM,class=class(model),EP=EP))
}
#model=fitT_g22
#imm.fm.smn.cr(cc, y,x1,model)


#To obtain the observed information matrix
imm.fm.smn.cr3 = function(cc, y,x1,model)
{
  #if((class(model) != "T") && (class(model) != "Normal")) stop(paste("Family",class(model),"not recognized.",sep=" "))
  n           <- length(y)
  p           <- ncol(x1)-1
  g           <- length(model$res$pii)

  Abetas      <- model$res$Abetas
  medj        <- model$res$medj
  sigma2      <- model$res$sigma2
  nu          <- model$res$nu
  pii         <- model$res$pii

  beta0       <- Abetas[1]
  betas       <- as.matrix(Abetas[2:(p+1)])   # parameters of regression dimension "p"
  varphi      <- rep(0,g) # beta0 + mu_j
  mu = mu1    <- matrix(0,n,g)
  x           <- as.matrix(x1[,2:(p+1)])

  for (k in 1:g)
  {
    varphi[k] <- beta0+medj[k]
    mu1[,k]   <- x%*%betas
    mu[,k]    <- mu1[,k]+varphi[k]
  }

  tal         <- matrix(0,n,g)
  Sibeta      <- matrix(0,p,n)
  Sivarphi    <- matrix(0,2,n)
  Sisigma     <- matrix(0,2,n)
  Sip         <- matrix(0,g,n)

  for(j in 1:g)
  {
    d            <- ((y-mu[,j])^2)/sigma2[j]
    u0           <- (nu+1)/(nu+d)
    u1           <- y*(nu+1)/(nu+d)
    u2           <- y^2*(nu+1)/(nu+d)

    d1           <- dT(cc, y, mu[,j], sigma2[j],nu)
    d2           <- d.mixedT(cc, y, pii, mu, sigma2,nu)
    tal[,j]      <- d1*pii[j]/d2

    Sibeta       <- Sibeta + t(x)%*%diag(tal[,j]*c(u1)/sigma2[j]) - t(x)%*%diag(tal[,j]*c(u0)*(mu1[,j] + varphi[j])/sigma2[j])
    Sip[j,]      <- (1/pii[j])*tal[,j] - (1/pii[g])*tal[,g]
    Sivarphi[j,] <- (1/sigma2[j])*(tal[,j]*c(u1) - tal[,j]*c(u0)*(mu1[,j] + varphi[j]))
    Sisigma[j,]  <- -0.5*(1/sigma2[j]^2)*tal[,j]*(sigma2[j] - c(u2) + 2*c(u1)*(mu1[,j] + varphi[j]) - c(u0)*(mu1[,j] + varphi[j])^2)
  }

  soma      <- matrix(0, (length(Abetas) -1) + 3*g -1, (length(Abetas) -1) + 3*g -1)
  si        <- matrix(0, n, (length(Abetas) -1) + 3*g -1)

  for(i in 1:n)
  {
    si[i,]  <- c(Sibeta[,i],Sip[1:(g-1),i],Sivarphi[,i],Sisigma[,i])
    soma    <- soma + cbind(si[i,])%*%si[i,]
  }

  #Jacobian

  p      <- length(betas)
  j11    <- diag(rep(1,p))
  j12    <- matrix(0, nrow=p, ncol=g-1)
  j13    <- matrix(0, nrow=p, ncol=1)
  j14    <- matrix(0, nrow=p, ncol=g)
  j15    <- matrix(0, nrow=p, ncol=g)
  J1     <- cbind(j11,j12,j13,j14,j15)

  j21    <- matrix(0, nrow=g-1, ncol=p)
  j22    <- diag(rep(1,g-1))
  j23    <- matrix(0, nrow=g-1, ncol=1)
  j24    <- -as.matrix((medj[1:(g-1)] - medj[g])^(-1))%*%t(as.matrix(pii))
  j24    <- matrix(j24,ncol=g,nrow=g-1)
  j25    <- matrix(0, nrow=g-1, ncol=g)
  J2     <- cbind(j21,j22,j23,j24,j25)

  j31    <- matrix(0, nrow=g, ncol=p)
  j32    <- matrix(0, nrow=g, ncol=g-1)
  j33    <- matrix(1, nrow=g, ncol=1)
  j34    <- diag(rep(1,g))
  j35    <- matrix(0, nrow=g, ncol=g)
  J3     <- cbind(j31,j32,j33,j34,j35)

  j41    <- matrix(0, nrow=g, ncol=p)
  j42    <- matrix(0, nrow=g, ncol=g-1)
  j43    <- matrix(0, nrow=g, ncol=1)
  j44    <- matrix(0, nrow=g, ncol=g)
  j45    <- diag(rep(1,g))
  J4     <- cbind(j41,j42,j43,j44,j45)

  Jacobian       <- rbind(J1,J2,J3,J4)
  IM             <- t(Jacobian)%*%solve(soma)%*%Jacobian
  EP             <- as.matrix(sqrt(diag(t(Jacobian)%*%solve(soma)%*%Jacobian)))

  return(list(IM=IM,class=class(model),EP=EP, Jacobian=Jacobian))
}

#imm.fm.smn.cr(y,x1,cc,fit)


#########################################################################################
#Programa esta sendo trasladado de Matlab para R
#Por Luis Benites
MinMax_kmeans <- function(X,M,k,p_init,p_max,p_step,t_max,beta)
{
  N        <- nrow(X)
  # This function implements the MinMax k-means algorithm as described in
  # G.Tzortzis and A.Likas, "The MinMax k-means Clustering Algorithm",
  # Pattern Recognition, 2014.

  #Function Inputs
  #===============
  #
  # X is an Nxd data matrix, where each row corresponds to an instance.
  # M is a kxd matrix of the initial cluster centers. Each row corresponds to a center.
  # k is the number of clusters.
  # p_init is the initial p value (0<=p_init<1).
  # p_max is the maximum admissible p value (0<=p_max<1).
  # p_step is the step for increasing p (p_step>=0).
  # t_max is the maximum number of iterations (necessary as convergence cannot
  #                                        be guaranteed for MinMax k-means).
  # beta controls the amount of memory for the weight updates (0<=beta<=1).


  #Function Outputs
  #================
  #
  # Cluster_elem is an N-dimensional row vector containing the final cluster assignments.
  # Clusters are indexed 1,2,...,k.
  # M is a kxd matrix of the final cluster centers. Each row corresponds to a center.
  # Var is a k-dimensional row vector containing the final variance of each cluster.
  # Courtesy of G. Tzortzis (Matla code)

  # Luis Benites <lbenitesanchez@gmail.com> translated the code into R

  if(p_init <0 || p_init >= 1)
    print('p_init must take a value in [0,1)')

  if(p_max <0 || p_max >= 1)
    print('p_max must take a value in [0,1)')

  if(p_max < p_init)
    print('p_max must be greater or equal to p_init')

  if(p_step < 0)
    print('p_step must be a non-negative number')

  if(beta<0 || beta>1)
    print('beta must take a value in [0,1]')




  if(p_init==p_max){
    if(p_step!=0)
      cat('p_step reset to zero, since p_max equals p_init \n\n')

    p_flag <- 0
    p_step <- 0

  }else if(p_step==0){

    if(p_init!=p_max)
      cat('p_s=max reset to equal p_init, since p_step=0\n\n')

    p_flag <- 0
    p_max  <- p_init

  }else{
    p_flag <- 1 #p_flag indicates whether p will be increased during the iterations.
  }

  #--------------------------------------------------------------------------
  #Weights are uniformly initialized.
  W         <- rep(1,k)/k

  #Other initializations
  p         <- p_init #Initial p value
  p_prev    <- p - 10^{-8} #Dummy value
  empty     <- 0 #Count the number of iterations for which an empty or singleton cluster is detected
  Iter      <- 1 #Numbers of iterations
  E_w_old   <- 10^{20} #Previous iteration objective (used to check convergence)
  Var       <- rep(0,k)
  Cluster_elem_history <- list()
  W_history <- list()

  #--------------------------------------------------------------------------

  #cat("\n Start of MinMax k-means iterations \n")
  #cat("------------------------------------- \n \n")

  #Calculate the squared Euclidean distances between the instances and the
  #initial cluster centers.
  Dist <- sqdist(t(M),t(X))

  #The MinMax k-means iterative procedure.
  kk   <- 0
  while(Iter < t_max)
  {

    #Calculate the weighted distances.
    #Each cluster is multipied by its weight.
    WDist <- Dist
    for(i in 1:k)
      WDist[i,] <- W[i]^p*Dist[i,]

    #Update the cluster assignments.
    min_WDist    <- apply(WDist,2,min)
    Cluster_elem <- c()
    for(i in 1:N) {Cluster_elem[i] <- which.min(WDist[,i])}

    #Calculate the MinMax k-means objective.
    E_w <- sum(min_WDist)

    #If empty or singleton clusters are detected after the update.
    for(i in 1:k)
    {
      I <- which(Cluster_elem==i)
      if(vector.is.empty(I) || length(I)==1)
      {
        #cat('Empty or singleton clusters detected for p=%g.\n',p)
        #cat('Reverting to previous p value.\n\n')

        E_w    <- NaN  #Current objective undefined.
        empty  <- empty + 1

        #Reduce p when empty or singleton clusters are detected.
        if (empty>1){
          p    <- p-p_step
          #The last p increase may not correspond to a complete p_step,
          #if the difference p_max-p_init is not an exact multiple of p_step.
        }else{
          p    <- p_prev
        }

        p_flag <-0  #Never increase p again.

        #p is not allowed to drop out of the given range.
        if(p<p_init || p_step==0)
        {
          #cat('\n+++++++++++++++++++++++++++++++++++++++++\n\n')
          #cat('p cannot be decreased further.\n')
          #cat('Either p_step=0 or p_init already reached.\n')
          #cat('Aborting Execution\n')
          #cat('\n+++++++++++++++++++++++++++++++++++++++++\n\n')

          #Return NaN to indicate that no solution is produced
          Cluster_elem <- rep(NaN,nrow(X))
          M            <- matrix( rep(do.call("rbind", rep(list(NA), k)),ncol(X)), ncol=ncol(X), byrow=TRUE)
          Var          <- rep(NaN,k)
        }
        #Continue from the assignments and the weights corresponding
        #to the decreased p value.
        Cluster_elem <- Cluster_elem_history[[empty]]
        W            <- W_history[[empty]]
        break
      }
    }

    if(!is.nan(E_w))
    {
      #cat('p =',p,'\n')
      #cat('The MinMax k-means objective is E_w=',E_w,'\n\n')
    }

    #Check for convergence. Never converge if in the current (or previous)
    #iteration empty or singleton clusters were detected.
    if (!is.nan(E_w) && !is.nan(E_w_old) && (abs(1-E_w/E_w_old)<1e-6 || Iter>=t_max))
    {
      #Calculate the cluster variances.
      Var      <- c()
      for(i in 1:k)
      {
        I      <- Cluster_elem==i
        Var[i] <- sum(Dist[i,I])
      }

      #cat('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n')
      #cat('Converging for p =', p,'after',Iter,'iterations.','\n')
      #cat('The final MinMax k-means objective is E_w =',E_w,'\n')
      #cat('The maximum cluster variance is E_max =',max(Var),'\n')
      #cat('The sum of the cluster variances is E_sum =',sum(Var),'\n')
      #cat('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n')

      break
    }

    E_w_old <- E_w

    #Update the cluster centers
    for(i in 1:k)
    {
      I     <- Cluster_elem==i
      M[i,] <- mean(X[I,]) #apply(X[I,],2,mean)
    }

    #Increase the p value
    if(p_flag==1)
    {
      #Keep the assignments-weights corresponding to the current p.
      #These are needed when empty or singleton clusters are found in
      #subsequent iterations.
      kk                         <- kk + 1
      Cluster_elem_history[[kk]] <- Cluster_elem
      W_history[[kk]]            <- W

      p_prev   <- p
      p        <- p+p_step

      if(p>=p_max)
      {
        p      <- p_max
        p_flag <- 0
        #cat('p_max reached\n\n');
      }
    }

    #Recalculate the distances between the instances and the cluster centers.
    Dist = sqdist(t(M),t(X))

    #Calculate the cluster variances
    Var <- c()
    for(i in 1:k)
    {
      I      <- Cluster_elem==i
      Var[i] <- sum(Dist[i,I])
    }

    W_old    <- W

    #Update the weights
    for(i in 1:k)
    {
      W[i]   <- 0
      for(j in 1:k)
        W[i] <- W[i]+(Var[i]/Var[j])^(1/(p-1))
      W[i]   <- 1/W[i]
    }

    #Memory effect.
    W     <- (1-beta)*W+beta*W_old
    Iter  <- Iter+1
  }
  object <- list(Cluster_elem=Cluster_elem,W=W,M=M,Var=Var)
  object
}


sqdist <- function(a,b)
{
  #sqdist - computes pairwise squared Euclidean distances between points
  #original version by Roland Bunschoten, 1999

  aa <- apply(a^2,2,sum)
  bb <- apply(b^2,2,sum)
  ab <- t(a)%*%b
  d  <- abs(repmat(cbind(aa),1,length(bb)) + t(repmat(cbind(bb),1,length(aa))) - 2*t(a)%*%b)
  return(d)
}

vector.is.empty <- function(x) return(length(x) ==0 )


repmat = function(X,m,n)
{
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}
