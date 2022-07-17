##########################################################################################################################
##########################################################################################################################
#
#              A funcao initial.Values permite obter os chutes iniciais para os parametros
#              betas_j, sigma_j e nu (aqui estamos supongo que nu=nu1=...=nuG)
#              A obtencao do chute inicial de nu e por medio de uma grade
#
##########################################################################################################################

#A funcao initial.Values permite obter os chutes iniciais para os parametros
#betas_j, sigma_j e nu (aqui estamos supongo que nu=nu1=...=nuG)
#A obtencao do chute inicial de nu e por medio de uma grade
#Aqui e considerando xij
initial.values.fmr.smn.cr <- function(cc, y,x,g=2,algorithm="k-medoids",family="T",lower=1,upper=20,space=0.1,plotLog = TRUE,searchNU=TRUE,printNU=TRUE, saveFigure = FALSE)
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

    k.iter.max     <- 50
    n.start        <- 1
    algorithm      <- "Hartigan-Wong"

    if(g > 1)
    {
      init         <- kmeans(y,g,k.iter.max,n.start,algorithm)
      pii          <- init$size/length(y)
      sigma2       <- c()
      Abetas       <- matrix(0,p,g)

      for (j in 1:g)
      {
        dd         <- init$cluster
        Abetas[,j] <- solve(t(x[dd==j,])%*%x[dd==j,])%*%t(x[dd==j,])%*%y[dd==j]
        yr         <- y[dd==j]-x[dd==j,]%*%Abetas[,j]
        sigma2[j]  <- sum(yr^2)/init$size[j]
      }
    }else{
      Abetas       <- solve(t(x)%*%x)%*%t(x)%*%y
      yr           <- y-x%*%Abetas
      sigma2       <- var(yr)
      pii          <- 1
    }
    initial_values <- c(Abetas,sigma2,pii)#;print(initial_values)
  }

  if (algorithm == "k-medoids")
  {
    if(length(g) == 0) {stop("g is not specified correctly.\n")}

    if(g > 1){
      Abetas       <- matrix(0,p,g)
      init         <- Cluster_Medoids(as.matrix(y),clusters=g, distance_metric="euclidean", swap_phase = TRUE, fuzzy = TRUE,seed=sample(1:100000,1))
      pii          <- init$medoid_indices/length(y)
      sigma2       <- c()
      for (j in 1:g)
      {
        dd         <- init$clusters
        Abetas[,j] <- solve(t(x[dd==j,])%*%x[dd==j,])%*%t(x[dd==j,])%*%y[dd==j]
        yr         <- y[dd==j]-x[dd==j,]%*%Abetas[,j]
        sigma2[j]  <- sum(yr^2)/init$medoid_indices[j]
      }
    } else{
      Abetas       <- solve(t(x)%*%x)%*%t(x)%*%y
      yr           <- y-x%*%Abetas
      sigma2       <- var(yr)
      pii          <- 1
    }

    initial_values  <- c(Abetas,sigma2,pii)#;print(initial_values)
  }
  ########################################################################################################################
  #-----------------------------------------------------------------------------------------------------------------------#
  #    As linhas de codigo de abaixo sao para poder obter o chute inicial para os graus de liberdade para a T, SL e CN
  #-----------------------------------------------------------------------------------------------------------------------#
  ########################################################################################################################

  if(searchNU ==TRUE)
  {
    mu             <- matrix(0,n,g)
    for(k in 1:g){mu[,k]<- x[[k]]%*%Abetas[[k]]}

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
          #eq <- bquote(bold(nu[max] == .(nu0)),bold(nu[max] == .(nu0)))
          #mtext(eq,side=3,cex=1.5)
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
      obj.out <- list(Abetas=Abetas,sigma2=sigma2,pii=pii,nu=nu0)
    }else{
      obj.out <- list(vero=vero, Abetas=Abetas,sigma2=sigma2,pii=pii, nu=nu0)
    }
  }else{
    if(family == "Normal"){
      nu0     <- 0
      obj.out <- list(Abetas=Abetas,sigma2=sigma2,pii=pii,nu=nu0)
    }else{
      obj.out <- list(Abetas=Abetas,sigma2=sigma2,pii=pii)
    }

  }
}

##########################################################################################################################

genmixsmsn = function(n, family, mu,Sigma,nu,pii)
{
  y       <- c()
  g       <- length(pii)

  if(family=="Skew.normal")
    for (i in 1:n)
    {
      arg1  <- list(mu[1], Sigma=Sigma[1], shape=c(0))
      arg2  <- list(mu[2], Sigma=Sigma[2], shape=c(0))
      y[i]  <- rmix(1, pii, family, list(arg1,arg2))
    }

  if(family=="Skew.t" || family=="Skew.slash")
  {
    arg1  <- list(mu[1], Sigma[1], c(0), nu=nu)
    arg2  <- list(mu[2], Sigma[2], c(0), nu=nu)
    for (i in 1:n){y[i]  <- rmix(1, pii, family, list(arg1,arg2))}
  }


  if(family=="Skew.cn")
    for (i in 1:n)
    {
      y[i]  <- rmix(1, pii, family, list(c(mu=mu[1], sigma2=Sigma[1], shape=0, nu=nu),c(mu=mu[2], sigma2=Sigma[2], shape=0, nu=nu)))
    }

  return(y)
}


##########################################################################################################################
##########################################################################################################################
## Densidade das SMN com locacao-escala #######
dNormal <- function(cc, y, mu, sigma2 = 1)
{
  densN        <- vector(mode = "numeric", length = length(y))
  densN[cc==0] <- dnorm(y[cc==0], mu[cc==0], sqrt(sigma2))
  densN[cc==1] <- 1-pnorm((y[cc==1] - mu[cc==1])/sqrt(sigma2))
  if(length(which(densN == 0)) > 0) densN[which(densN == 0)] <- .Machine$double.xmin
  return(densN)
}

####
dT <- function(cc, y, mu, sigma2 = 1,nu=4)
{
  densN        <- vector(mode = "numeric", length = length(y))
  aux          <- (y-mu)/sqrt(sigma2)
  densN[cc==0] <- dt(aux[cc==0],nu)/sqrt(sigma2)
  densN[cc==1] <- 1-pt(aux[cc==1],nu)
  if(length(which(densN == 0)) > 0) densN[which(densN == 0)] <- .Machine$double.xmin
  return(densN)
}

####

dSL <- function(cc, y, mu, sigma2 = 1,nu=2)
{
  densN        <- vector(mode = "numeric", length = length(y))
  aux          <-(y-mu)/sqrt(sigma2)
  densN[cc==0] <- dSlash (aux[cc==0],0,1,nu)/sqrt(sigma2)
  densN[cc==1] <- 1-AcumSlash(aux[cc==1],0,1,nu)
  if(length(which(densN == 0)) > 0) densN[which(densN == 0)] <- .Machine$double.xmin
  return(densN)
}

#####

dCN <- function(cc, y, mu, sigma2 = 1,nu=c(0.1,0.1))
{
  densN        <- vector(mode = "numeric", length = length(y))
  aux          <- (y-mu)/sqrt(sigma2)
  densN[cc==0] <- dNormalC(aux[cc==0],0,1,nu)/sqrt(sigma2)
  densN[cc==1] <- 1-AcumNormalC(aux[cc==1],0,1,nu)
  if(length(which(densN == 0)) > 0) densN[which(densN == 0)] <- .Machine$double.xmin
  return(densN)
}

####

###########    Densidades das Misturas de SMN ##################

d.mixedN <- function(cc, x, pi1, mu, sigma2)
{
  # x: o vetor de dados
  ## mu[,] uma matrix
  # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
  g    <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*dNormal(cc, x, mu[, j], sigma2[j])
  return(dens)
}


d.mixedT <- function(cc, x, pi1, mu, sigma2, nu)
{
  # x: o vetor de dados
  ## mu[,] uma matrix
  # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
  g    <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*dT(cc, x, mu[, j], sigma2[j],nu)
  return(dens)
}

d.mixedT2 <- function(cc, x, pi1, mu, sigma2, nu)
{
  # x: o vetor de dados
  ## mu[,] uma matrix
  # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
  g    <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*dT(cc, x, mu[, j], sigma2[j],nu[j])
  return(dens)
}


d.mixedSL <- function(cc, x, pi1, mu, sigma2, nu)
{
  # x: o vetor de dados
  ## mu[,] uma matrix
  # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
  g    <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*dSL(cc, x, mu[, j], sigma2[j],nu)
  return(dens)
}

d.mixedSL2 <- function(cc, x, pi1, mu, sigma2, nu)
{
  # x: o vetor de dados
  # mu[,] uma matrix
  # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
  g    <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*dSL(cc, x, mu[, j], sigma2[j],nu[j])
  return(dens)
}


d.mixedCN <- function(cc, x, pi1, mu, sigma2, nu)
{
  # x: o vetor de dados
  #  mu[,] uma matrix
  # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
  g    <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*dCN(cc, x, mu[, j], sigma2[j],nu)
  return(dens)
}


##########     FIM   das Densidades da Normal e T Censuras  ############
################################################################################


AcumPearsonVII <- function(y,mu,sigma2,nu,delta)
{
  Acum <- z <- vector(mode = "numeric", length = length(y))
  sigma2a   <- sigma2*(delta/nu)
  for (i in 1:length(y))
  {
    z[i]    <- (y[i]-mu)/sqrt(sigma2a)
    Acum[i] <- pt(z[i],df=nu)
  }
  return(Acum)
}

P <- function(y,mu,sigma2,nu,delta)
{
  A       <- z <- vector(mode = "numeric", length = length(y))
  sigma2a <- sigma2*(delta/nu)
  n       <- length(y)
  i       <- 0
  while (i<n)
  {
    i     <- i +1
    z[i]  <- (y[i]-mu)/sqrt(sigma2a)
    A[i]  <- pt(z[i],df=nu)
  }
  return(A)
}


AcumSlash <- function(y,mu,sigma2,nu)
{
  Acum   <- z <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y))
  {
    z[i] <- (y[i]-mu)/sqrt(sigma2)
    f1   <- function(u) nu*u^(nu-1)*pnorm(z[i]*sqrt(u))
    Acum[i]<- integrate(f1,0,1)$value
  }
  return(Acum)
}

AcumNormalC <- function(y,mu,sigma2,nu)
{
  Acum   <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y))
  {
    eta  <- nu[1]
    gama <- nu[2]
    Acum[i] <- eta*pnorm(y[i],mu,sqrt(sigma2/gama)) + (1-eta)*pnorm(y[i],mu,sqrt(sigma2))
  }
  return(Acum)
}

dPearsonVII<- function(y,mu,sigma2,nu,delta)
{
  f <- z <- vector(mode = "numeric", length = length(y))
  sigma2a <- sigma2*(delta/nu)
  for (i in 1:length(y))
  {
    z[i] <- (y[i]-mu)/sqrt(sigma2a)
    f[i] <- dt(z[i],df=nu)/sqrt(sigma2a)
  }
  return(f)
}




NCensurEsperUY <-  function(y,mu,sigma2,nu,delta,type)
{
  EUY0   <- EUY1 <- EUY2 <- c()
  d      <- (y-mu)^2/sigma2
  n      <- length(y)
  if(type=="T")
  {
    EUY0 <- (nu+1)/(nu+d)
    EUY1 <- y*(nu+1)/(nu+d)
    EUY2 <- (y^2)*(nu+1)/(nu+d)
  }
  if(type=="Normal")
  {
    EUY0 <- rep(1,n)
    EUY1 <- y
    EUY2 <- (y^2)
  }
  if(type=="PearsonVII")
  {
    d    <- (y-mu)^2/sigma2
    EUY0 <- (nu+1)/(delta+d)
    EUY1 <- y*(nu+1)/(delta+d)
    EUY2 <- (y^2)*(nu+1)/(delta+d)
  }
  if(type=="Slash")
  {
    Num  <- GamaInc(nu+1.5,0.5*d)*(0.5*d)^(-nu-1.5)
    Den  <- GamaInc(nu+0.5,0.5*d)*(0.5*d)^(-nu-0.5)
    EUY0 <- Num/Den
    EUY1 <- y*Num/Den
    EUY2 <- (y^2)*Num/Den
  }
  if(type=="NormalC")
  {
    Num  <- 1-nu[1]+nu[1]*(nu[2])^(1.5)*exp(0.5*d*(1-nu[2]))
    Den  <- 1-nu[1]+nu[1]*(nu[2])^(0.5)*exp(0.5*d*(1-nu[2]))
    EUY0 <- Num/Den
    EUY1 <- y*Num/Den
    EUY2 <- (y^2)*Num/Den
  }
  return(list(EUY0=EUY0,EUY1=EUY1,EUY2=EUY2))
}

TN <- function(mu,sigma2,a,b)
{
  sigma2 <- as.vector(sigma2)
  lim1   <- (a-mu)/sqrt(sigma2)
  lim2   <- (b-mu)/sqrt(sigma2)
  EX     <-   (dnorm(lim2)- dnorm(lim1))/(pnorm(lim2)- pnorm(lim1))
  EX2    <- 1-(lim2*dnorm(lim2)- lim1*dnorm(lim1))/(pnorm(lim2)- pnorm(lim1))
  EX3    <- 2*EX - (lim2^2*dnorm(lim2)- lim1^2*dnorm(lim1))/(pnorm(lim2)- pnorm(lim1))
  EX4    <- 3*EX2 - (lim2^3*dnorm(lim2)- lim1^3*dnorm(lim1))/(pnorm(lim2)- pnorm(lim1))
  EY     <- mu + sqrt(sigma2)*EX
  EY2    <- mu^2 + 2*mu*sqrt(sigma2)*EX + sigma2*(EX2)
  EY3    <- mu^3 + 3*mu^2*sqrt(sigma2)*EX + 3*mu*sigma2*(EX2) + EX3
  EY4    <- mu^4 + 4*mu^3*sqrt(sigma2)*EX + 6*mu^2*sigma2*(EX2) + 4*mu*sigma2^1.5*EX3 + EX4*sigma2^2
  return(list(EY=EY, EY2=EY2, EY3=EY3, EY4=EY4))
}

GamaInc <- function(a,xx)
{
  res <- vector(mode = "numeric", length = length(xx))
  f <-function(t) exp(-t)*t^(a-1)
  for  (i in 1:length(xx))
  {
    res[i] <- integrate(f,0,xx[i])$value
  }
  return(res)
}


E_phi <- function(r,a,nu,delta,type=type,cens=cens)
{
  n       <- length(a)
  b       <- rep(Inf,n)
  b1      <- rep(-Inf,n)

  if(setequal(a,b)== TRUE | setequal(a,b1)== TRUE)
  {
    resp  <- rep(0,n)
  }
  else
  {
    if(type=="Normal")
    {
      resp <- dnorm(a)
    }
    if(type=="T")
    {
      Aux0 <- gamma(0.5*(nu+2*r))
      Aux1 <- gamma(nu/2)*sqrt(2*pi)
      Aux2 <- Aux0/Aux1
      Aux3 <- (0.5*nu)^(nu/2)
      Aux4 <- (0.5*(a^2+nu))^(-0.5*(nu+2*r))
      resp <- Aux2*Aux3*Aux4
    }
    if(type=="PearsonVII")
    {
      Aux0 <- gamma(0.5*(nu+2*r))
      Aux1 <- gamma(nu/2)*sqrt(2*pi)
      Aux2 <- Aux0/Aux1
      Aux3 <- (0.5*delta)^(nu/2)
      Aux4 <- (0.5*(a^2+delta))^(-0.5*(nu+2*r))
      resp <- Aux2*Aux3*Aux4
    }
    if(type=="Slash")
    {
      Aux0 <- nu/sqrt(2*pi)
      Aux1 <- (0.5*a^2)^(-(nu+r))
      Aux2 <- GamaInc(nu+r,0.5*a^2)
      resp <- Aux0*Aux1*Aux2
    }
    if(type=="NormalC")
    {
      Aux0 <- nu[1]*nu[2]^(r)*dnorm(a*sqrt(nu[2]))
      Aux1 <- (1-nu[1])*dnorm(a)
      resp <- Aux0 + Aux1
    }
  }
  return(resp)
}

E_Phi <- function(r,a,nu,delta,type=type)
{
  n      <- length(a)
  if(type=="Normal")
  {
    resp <- pnorm(a)
  }
  if(type=="T")
  {
    Aux0 <- gamma(0.5*(nu+(2*r)))
    Aux1 <- gamma(nu/2)
    Aux2 <- Aux0/Aux1
    Aux3 <- (0.5*nu)^(-r)
    Aux4 <- AcumPearsonVII(a,0,1,nu+(2*r),nu)
    resp <- Aux2*Aux3*Aux4
  }
  if(type=="PearsonVII")
  {
    Aux0 <- gamma(0.5*(nu+(2*r)))
    Aux1 <- gamma(nu/2)
    Aux2 <- Aux0/Aux1
    Aux3 <- (0.5*delta)^(-r)
    Aux4 <- AcumPearsonVII(a,0,1,nu+(2*r),delta)
    resp <- Aux2*Aux3*Aux4
  }
  if(type=="Slash")
  {
    Aux0 <- nu/(nu+r)
    Aux1 <- AcumSlash(a,0,1,nu+r)
    resp <- Aux0*Aux1
  }
  if(type=="NormalC")
  {
    Aux0 <- nu[2]^(r)*AcumNormalC(a,0,1,nu)
    Aux1 <- (1-nu[2]^(r))*(1-nu[1])*pnorm(a)
    resp <- Aux0 + Aux1
  }
  return(resp)
}

CensEsperUY1 <- function(mu,sigma2,nu,delta,Lim1,Lim2,type=type,cens=cens)
{
  Lim11   <- (Lim1-mu)/sqrt(sigma2)
  Lim21   <- (Lim2-mu)/sqrt(sigma2)
  n       <- length(Lim11)
  if(type=="Normal")
  {
    EU    <-  1
    FNIb  <-  pnorm(Lim21)
    FNIa  <-  pnorm(Lim11)
  }
  if(type=="T")
  {
    EU    <-  1
    FNIb  <-  pt(Lim21,nu)
    FNIa  <-  pt(Lim11,nu)
  }
  if(type=="PearsonVII")
  {
    FNIb  <-  AcumPearsonVII(Lim21,0,1,nu,delta)
    FNIa  <-  AcumPearsonVII(Lim11,0,1,nu,delta)
  }
  if(type=="Slash")
  {
    FNIb  <- AcumSlash(Lim21,0,1,nu)
    FNIa  <- AcumSlash(Lim11,0,1,nu)
  }
  if(type=="NormalC")
  {
    EU    <-  (nu[1]*nu[2]) + (1-nu[1])
    FNIb  <- AcumNormalC(Lim21,0,1,nu)
    FNIa  <- AcumNormalC(Lim11,0,1,nu)
  }
  if (cens=="1")
  {
    Aux11 <- rep(0,n)
  }else
  {
    Aux11 <- Lim11
  }
  if (cens=="2")
  {
    Aux22 <- rep(0,n)
  }else
  {
    Aux22 <- Lim21
  }
  Kaux    <- (FNIb-FNIa)

  if(length(which(Kaux == 0)) > 0) Kaux[which(Kaux == 0)] <- .Machine$double.xmin
  K     <- 1/Kaux
  EUX0  <- K*(E_Phi(1,Lim21, nu,delta,type)- E_Phi(1,Lim11, nu,delta,type))
  EUX1  <- K*(E_phi(0.5,Lim11,nu,delta,type,cens)- E_phi(0.5,Lim21,nu,delta,type,cens))
  EUX2  <- K*(E_Phi(0,Lim21, nu,delta,type)- E_Phi(0,Lim11, nu,delta,type) + Aux11*E_phi(0.5,Lim11,nu,delta,type,cens) - Aux22*E_phi(0.5,Lim21,nu,delta,type,cens))           # Neste c aso r =2
  EUX20 <- K*(E_Phi(2,Lim21, nu,delta,type)- E_Phi(2,Lim11, nu,delta,type))
  EUX21 <- K*(E_phi(1.5,Lim11,nu,delta,type,cens)- E_phi(1.5,Lim21,nu,delta,type,cens))
  EUY0  <- EUX0
  EUY1  <- mu*EUX0 + sqrt(sigma2)*EUX1
  EUY2  <- EUX0*mu^2 + 2*mu*sqrt(sigma2)*EUX1 + sigma2*EUX2
  EUY20 <- EUX20
  EUY21 <- mu*EUX20 + sqrt(sigma2)*EUX21
  return(list(EUY0=EUY0,EUY1=EUY1,EUY2=EUY2, EUY20=EUY20, EUY21=EUY21))
}

##########################################################################################################################
#
#                          As seguintes funcoes sao para calcular os EP utilizando a aproximacao de Basford
#
##########################################################################################################################

#Esta funcao para obter os EP e considerando nu1=nu2=....=nuG
#e tambem e considerando xij
im.fmr.smn.cr = function(cc, y,x,Abetas,sigma2,pii,nu,family)
{
  if((family != "T") && (family != "Normal")  && (family != "Slash") && (family != "NormalC")) stop(paste("Family",family,"not recognized.",sep=" "))

  if(family=="Normal")
  {
    n             <- length(y)

    Lim1          <- y #Esta sendo considerado o censura tipo 2
    Lim2          <- rep(Inf,n)

    p             <- ncol(x)
    g             <- length(pii)

    Abetas        <- Abetas
    sigma2        <- sigma2
    pii           <- pii

    mu            <- matrix(0,n,g)

    for (k in 1:g){mu[,k] <- x%*%Abetas[,k]}

    tal           <- matrix(0,n,g)
    Sibeta        <- rep(list(matrix(0,p,n)),2)
    Sisigma       <- matrix(0,g,n)
    Sip           <- matrix(0,g,n)

    for(j in 1:g)
    {
      NCensEUY    <- NCensurEsperUY(y,mu[,j],sigma2[j],nu=NULL,0,type="Normal")
      u0          <- NCensEUY$EUY0
      u1          <- NCensEUY$EUY1
      u2          <- NCensEUY$EUY2

      CensEUY     <- CensEsperUY1(mu[cc==1,j],sigma2=sigma2[j],nu=0,delta=0,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type="Normal", cens="2")
      u0[cc==1]   <- CensEUY$EUY0
      u1[cc==1]   <- CensEUY$EUY1
      u2[cc==1]   <- CensEUY$EUY2

      d1          <- dNormal(cc, y, mu[,j], sigma2[j])
      d2          <- d.mixedN(cc, y, pii, mu, sigma2)
      tal[,j]     <- d1*pii[j]/d2

      Sibeta[[j]] <- Sibeta[[j]] + t(x)%*%diag(tal[,j]*c(u1)/sigma2[j]) - t(x)%*%diag(tal[,j]*c(u0)*mu[,j]/sigma2[j])
      Sip[j,]     <- (1/pii[j])*tal[,j] - (1/pii[g])*tal[,g]
      Sisigma[j,] <- -0.5*(1/sigma2[j]^2)*tal[,j]*(sigma2[j] - c(u2) + 2*c(u1)*mu[,j] - c(u0)*mu[,j]^2)
    }

    soma          <- matrix(0, g*nrow(Abetas)  + 2*g - 1, g*nrow(Abetas)  + 2*g - 1)
    si            <- matrix(0, n, g*nrow(Abetas)  + 2*g - 1)

    if(g==1)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    if(g==2)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sibeta[[2]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    if(g==3)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sibeta[[2]][,i],Sibeta[[3]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    if(g==4)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sibeta[[2]][,i],Sibeta[[3]][,i],Sibeta[[4]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    IM            <- solve(soma)
    EP            <- as.matrix(sqrt(diag(IM)))
  }

  #%-----------------------------------------------------------------------------------------------------------%
  #%-----------------------------------------------------------------------------------------------------------%
  if(family=="T")
  {
    n             <- length(y)

    Lim1          <- y #Esta sendo considerado o censura tipo 2
    Lim2          <- rep(Inf,n)

    p             <- ncol(x)
    g             <- length(pii)

    Abetas        <- Abetas
    sigma2        <- sigma2
    pii           <- pii
    nu            <- nu

    mu            <- matrix(0,n,g)

    for (k in 1:g){mu[,k] <- x%*%Abetas[,k]}

    tal           <- matrix(0,n,g)
    Sibeta        <- rep(list(matrix(0,p,n)),2)
    Sisigma       <- matrix(0,g,n)
    Sip           <- matrix(0,g,n)

    for(j in 1:g)
    {
      NCensEUY    <- NCensurEsperUY(y,mu[,j],sigma2[j],nu=nu,0,type="T")
      u0          <- NCensEUY$EUY0
      u1          <- NCensEUY$EUY1
      u2          <- NCensEUY$EUY2

      CensEUY     <- CensEsperUY1(mu[cc==1,j],sigma2=sigma2[j],nu=nu,delta=0,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type="T", cens="2")
      u0[cc==1]   <- CensEUY$EUY0
      u1[cc==1]   <- CensEUY$EUY1
      u2[cc==1]   <- CensEUY$EUY2

      d1          <- dT(cc, y, mu[,j], sigma2[j],nu)
      d2          <- d.mixedT(cc, y, pii, mu, sigma2,nu)
      tal[,j]     <- d1*pii[j]/d2

      Sibeta[[j]] <- Sibeta[[j]] + t(x)%*%diag(tal[,j]*c(u1)/sigma2[j]) - t(x)%*%diag(tal[,j]*c(u0)*mu[,j]/sigma2[j])
      Sip[j,]     <- (1/pii[j])*tal[,j] - (1/pii[g])*tal[,g]
      Sisigma[j,] <- -0.5*(1/sigma2[j]^2)*tal[,j]*(sigma2[j] - c(u2) + 2*c(u1)*mu[,j] - c(u0)*mu[,j]^2)
    }

    soma          <- matrix(0, g*nrow(Abetas)  + 2*g - 1, g*nrow(Abetas)  + 2*g - 1)
    si            <- matrix(0, n, g*nrow(Abetas)  + 2*g - 1)

    if(g==1)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    if(g==2)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sibeta[[2]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    if(g==3)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sibeta[[2]][,i],Sibeta[[3]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    if(g==4)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sibeta[[2]][,i],Sibeta[[3]][,i],Sibeta[[4]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    IM            <- solve(soma)
    EP            <- as.matrix(sqrt(diag(IM)))
  }

  #%-----------------------------------------------------------------------------------------------------------%
  #%-----------------------------------------------------------------------------------------------------------%
  if(family=="Slash")
  {
    n             <- length(y)

    Lim1          <- y #Esta sendo considerado o censura tipo 2
    Lim2          <- rep(Inf,n)

    p             <- ncol(x)
    g             <- length(pii)

    Abetas        <- Abetas
    sigma2        <- sigma2
    pii           <- pii
    nu            <- nu

    mu            <- matrix(0,n,g)

    for (k in 1:g){mu[,k] <- x%*%Abetas[,k]}

    tal           <- matrix(0,n,g)
    Sibeta        <- rep(list(matrix(0,p,n)),2)
    Sisigma       <- matrix(0,g,n)
    Sip           <- matrix(0,g,n)

    for(j in 1:g)
    {
      NCensEUY    <- NCensurEsperUY(y,mu[,j],sigma2[j],nu=nu,0,type="Slash")
      u0          <- NCensEUY$EUY0
      u1          <- NCensEUY$EUY1
      u2          <- NCensEUY$EUY2

      CensEUY     <- CensEsperUY1(mu[cc==1,j],sigma2=sigma2[j],nu=nu,delta=0,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type="Slash", cens="2")
      u0[cc==1]   <- CensEUY$EUY0
      u1[cc==1]   <- CensEUY$EUY1
      u2[cc==1]   <- CensEUY$EUY2

      d1          <- dSL(cc, y, mu[,j], sigma2[j],nu)
      d2          <- d.mixedSL(cc, y, pii, mu, sigma2,nu)
      tal[,j]     <- d1*pii[j]/d2

      Sibeta[[j]] <- Sibeta[[j]] + t(x)%*%diag(tal[,j]*c(u1)/sigma2[j]) - t(x)%*%diag(tal[,j]*c(u0)*mu[,j]/sigma2[j])
      Sip[j,]     <- (1/pii[j])*tal[,j] - (1/pii[g])*tal[,g]
      Sisigma[j,] <- -0.5*(1/sigma2[j]^2)*tal[,j]*(sigma2[j] - c(u2) + 2*c(u1)*mu[,j] - c(u0)*mu[,j]^2)
    }

    soma          <- matrix(0, g*nrow(Abetas)  + 2*g - 1, g*nrow(Abetas)  + 2*g - 1)
    si            <- matrix(0, n, g*nrow(Abetas)  + 2*g - 1)

    if(g==1)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    if(g==2)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sibeta[[2]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    if(g==3)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sibeta[[2]][,i],Sibeta[[3]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    if(g==4)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sibeta[[2]][,i],Sibeta[[3]][,i],Sibeta[[4]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    IM            <- solve(soma)
    EP            <- as.matrix(sqrt(diag(IM)))
  }

  #%-----------------------------------------------------------------------------------------------------------%
  #%-----------------------------------------------------------------------------------------------------------%
  if(family=="NormalC")
  {
    n             <- length(y)

    Lim1          <- y #Esta sendo considerado o censura tipo 2
    Lim2          <- rep(Inf,n)

    p             <- ncol(x)
    g             <- length(pii)

    Abetas        <- Abetas
    sigma2        <- sigma2
    pii           <- pii
    nu            <- nu

    mu            <- matrix(0,n,g)

    for (k in 1:g){mu[,k] <- x%*%Abetas[,k]}

    tal           <- matrix(0,n,g)
    Sibeta        <- rep(list(matrix(0,p,n)),2)
    Sisigma       <- matrix(0,g,n)
    Sip           <- matrix(0,g,n)

    for(j in 1:g)
    {
      NCensEUY    <- NCensurEsperUY(y,mu[,j],sigma2[j],nu=nu,0,type="NormalC")
      u0          <- NCensEUY$EUY0
      u1          <- NCensEUY$EUY1
      u2          <- NCensEUY$EUY2

      CensEUY     <- CensEsperUY1(mu[cc==1,j],sigma2=sigma2[j],nu=nu,delta=0,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type="NormalC", cens="2")
      u0[cc==1]   <- CensEUY$EUY0
      u1[cc==1]   <- CensEUY$EUY1
      u2[cc==1]   <- CensEUY$EUY2

      d1          <- dCN(cc, y, mu[,j], sigma2[j],nu)
      d2          <- d.mixedCN(cc, y, pii, mu, sigma2,nu)
      tal[,j]     <- d1*pii[j]/d2

      Sibeta[[j]] <- Sibeta[[j]] + t(x)%*%diag(tal[,j]*c(u1)/sigma2[j]) - t(x)%*%diag(tal[,j]*c(u0)*mu[,j]/sigma2[j])
      Sip[j,]     <- (1/pii[j])*tal[,j] - (1/pii[g])*tal[,g]
      Sisigma[j,] <- -0.5*(1/sigma2[j]^2)*tal[,j]*(sigma2[j] - c(u2) + 2*c(u1)*mu[,j] - c(u0)*mu[,j]^2)
    }

    soma          <- matrix(0, g*nrow(Abetas)  + 2*g - 1, g*nrow(Abetas)  + 2*g - 1)
    si            <- matrix(0, n, g*nrow(Abetas)  + 2*g - 1)

    if(g==1)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    if(g==2)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sibeta[[2]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    if(g==3)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sibeta[[2]][,i],Sibeta[[3]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    if(g==4)
    {
      for(i in 1:n)
      {
        si[i,]    <- c(Sibeta[[1]][,i],Sibeta[[2]][,i],Sibeta[[3]][,i],Sibeta[[4]][,i],Sip[1:(g-1),i],Sisigma[,i])
        soma      <- soma + cbind(si[i,])%*%si[i,]
      }
    }

    IM            <- solve(soma)
    EP            <- as.matrix(sqrt(diag(IM)))
  }

  return(list(IM=IM,class=family,EP=EP))
}


