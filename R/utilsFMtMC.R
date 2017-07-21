
################################################################################
### Derivada Sigma respecto alpha_jk
################################################################################

deriv.sigma <- function(Sigma,k,p){
  k     <- as.numeric(k)
  Sigma <- as.matrix(Sigma)
  p     <- dim(Sigma)[2]
  if (dim(Sigma)[1] != dim(Sigma)[2]) stop("Sigma is not square matrix\n")
  if (k > (p+1)*p/2)  stop("k out of bounds\n")
  e    <- rep(0,(p+1)*p/2)
  e[k] <- 1
  dev  <- matrix(0,ncol = p, nrow = p)
  dev[lower.tri(dev, diag = TRUE)] <- e
  dev  <- dev + t(dev)
  if(sum(diag(dev)) > 1) dev <- dev/2
  return(dev)
}

################################################################################
### label matrix de Info
################################################################################

nombres <- function(g,p, order = "TRUE", uni.Sigma = "FALSE"){

  NAME <- piiN <- c()
  if (uni.Sigma == "FALSE"){
    if (order == "TRUE") {
      for (i in 1:g) {
        SigmaN <- muN <- shapeN <- c()

        piiN <- paste("pii_", i, sep = "")

        for (k in 1:p) {
          muN <- c(muN, paste("mu", i, "_", k, sep = ""))
        }

        l <- m <- 1
        for (k in 1:((p + 1) * p/2)) {
          Vis <- FALSE
          SigmaN <- c(SigmaN, paste("Sigma", i, "_", l,
                                    m, sep = ""))
          if (((l * m - p * floor((l * m)/p)) == 0) &&
                (l != m)) {
            l <- l + 1
            m <- l
            Vis <- TRUE
          }
          if (!Vis)
            m <- m + 1
        } # fim do seg for k

        NAME <- c(NAME, muN, shapeN, SigmaN, piiN)

      } # fim do seg for i
      return(NAME[-length(NAME)])
    }

    if (order == "FALSE"){
      for (i in 1:g) {
        SigmaN <- muN <- shapeN <- c()

        for (k in 1:p) {
          muN <- c(muN, paste("mu", i, "_", k, sep = ""))
        }

        l <- m <- 1
        for (k in 1:((p + 1) * p/2)) {
          Vis <- FALSE
          SigmaN <- c(SigmaN, paste("Sigma", i, "_", l,
                                    m, sep = ""))
          if (((l * m - p * floor((l * m)/p)) == 0) &&
                (l != m)) {
            l <- l + 1
            m <- l
            Vis <- TRUE
          }
          if (!Vis)
            m <- m + 1
        } # fim do seg for k



        NAME <- c(NAME, muN, shapeN, SigmaN, piiN)

      } # fim do seg for i

      for (w in 1:(g - 1)) piiN <- c(piiN, paste("pii_", w, sep = ""))
      return(c(NAME, piiN))
    }


  } #fim uni.Sigma = FALSE


  if (uni.Sigma == "TRUE"){
    if (order == "TRUE") {

      for (i in 1:g) {
        SigmaN <- muN <- shapeN <- c()
        piiN <- paste("pii_", i, sep = "")

        for (k in 1:p)  muN <- c(muN, paste("mu", i, "_", k, sep = ""))
        NAME <- c(NAME, muN, piiN)
      } # fim do seg for i
      NAME <- NAME[-length(NAME)]

      l <- m <- 1
      for (k in 1:((p + 1) * p/2)) {
        Vis <- FALSE
        SigmaN <- c(SigmaN, paste("Sigma_", l, m, sep = ""))
        if (((l * m - p * floor((l * m)/p)) == 0) && (l != m)) {
          l <- l + 1
          m <- l
          Vis <- TRUE
        }
        if (!Vis)
          m <- m + 1
      } # fim do seg for k
      NAME <- c(NAME, SigmaN)
      return(NAME)
    }

    if (order == "FALSE") {

      for (i in 1:g) {
        SigmaN <- muN <- shapeN <- c()
        piiN <- paste("pii_", i, sep = "")

        for (k in 1:p)  muN <- c(muN, paste("mu", i, "_", k, sep = ""))
        NAME <- c(NAME, muN, piiN)
      } # fim do seg for i
      NAME <- NAME[-length(NAME)]

      l <- m <- 1
      for (k in 1:((p + 1) * p/2)) {
        Vis <- FALSE
        SigmaN <- c(SigmaN, paste("Sigma_", l, m, sep = ""))
        if (((l * m - p * floor((l * m)/p)) == 0) && (l != m)) {
          l <- l + 1
          m <- l
          Vis <- TRUE
        }
        if (!Vis)
          m <- m + 1
      } # fim do seg for k
      NAME <- c(NAME, SigmaN)
      return(NAME)
    }


  } #fim uni.Sigma = TRUE

}

################################################################################
### Raiz Cuadrada de una Matriz
################################################################################

matrix.sqrt <- function(A) {
  sva <- svd(A)
  if (min(sva$d)>=0)
    Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
  else
    stop("Matrix square root is not defined")
  return(Asqrt)
}

################################################################################
### Dados Simulados
################################################################################

gen.SN.multi <- function(n, mu, Sigma, shape, v=NULL){
  p <- length(mu)
  delta <- shape / (sqrt(1 + t(shape)%*%shape))
  y <- matrix(0,n,p)
  for (i in 1:n) y[i,] <- mu + matrix.sqrt(Sigma)%*%(delta*abs(rnorm(1)) + matrix.sqrt(diag(p) - delta%*%t(delta))%*%as.vector(rmvnorm(1, mean = rep(0,p), sigma = diag(p))))
  return(y)
}

gen.ST.multi <- function(n, mu, Sigma, shape, v){
  v <- as.numeric(v)
  y  <- matrix(rep(mu, n), n, length(mu), byrow = T) + (rgamma(n, v/2, v/2))^(-1/2)*gen.SN.multi(n, rep(0, length(mu)), Sigma, shape)
  return(y)
}

gen.SCN.multi  <- function(n, mu, Sigma, shape, v){
  x1 <- matrix(0,n,length(mu))
  for (i in 1:n){
    u <- runif(1)
    if (u < v[1]) x1[i,] <- gen.SN.multi(1, mu, Sigma/v[1], shape)
    if (u > v[1]) x1[i,] <- gen.SN.multi(1, mu, Sigma, shape)
  }
  return(x1)
}

gen.SS.multi  <- function(n, mu, Sigma, shape, v){
  u1 <- runif(n)
  u2 <- u1^(1/(v))   # formula 10 do artigo e metodo da inversao
  ys <- mu + (u2)^(-1/2)*gen.SN.multi(n, c(0,0), Sigma, shape)
  return(ys)
}


rmmixcr <- function(n, pii, mu, Sigma, shape, nu, percCensu, family){

  if(sum(pii) != 1) stop("Sum of pii diferent of one")
  if(length(percCensu) != length(mu)) stop("Censure in a component")
  if(sum(percCensu) > length(mu)) stop("Sum censure greather that one")

  w <- rmultinom(n, size = 1, prob = pii)
  G <- rowSums(w)
  z <- h <- NULL
  p <- length(mu[[1]])

  for (r in 1:length(G)){

    if(family[[r]] == "Normal")      y <- gen.SN.multi(n = G[r], mu = mu[[r]], Sigma = Sigma[[r]], shape = shape[[r]])
    if(family[[r]] == "Skew.normal") y <- gen.SN.multi(n = G[r], mu = mu[[r]], Sigma = Sigma[[r]], shape = shape[[r]], v = nu[r])
    if(family[[r]] == "Skew.t")      y <- gen.ST.multi(n = G[r], mu = mu[[r]], Sigma = Sigma[[r]], shape = shape[[r]], v = nu[r])
    if(family[[r]] == "t")           y <- gen.ST.multi(n = G[r], mu = mu[[r]], Sigma = Sigma[[r]], shape = shape[[r]], v = nu[r])
    if(family[[r]] == "Skew.cn")     y <- gen.SCN.multi(n = G[r], mu = mu[[r]], Sigma = Sigma[[r]], shape = shape[[r]], v = nu[r])
    if(family[[r]] == "Skew.slash")  y <- gen.SS.multi(n = G[r], mu = mu[[r]], Sigma = Sigma[[r]], shape = shape[[r]], v = nu[r])

    m   <- dim(y)[1]
    y1  <- y
    aa  <- sort(y)
    aa1 <- matrix(t(aa),m*p,1)
    bb    <- aa1[1:(percCensu[r]*m*p)];
    cutof <- bb[percCensu[r]*m*p];
    cc    <- (y < cutof)+0
    y[cc==1] <- cutof
    y <- t(t(y))
    z <- rbind(z,y)
    h <- rbind(h,cc)
  }

  return(list(y = z, cc = h, G = G))
}

################################################################################
### Densidades Normal e Misturas
################################################################################

dmvNCens<-function(cc, y, mu, Sigma){

  m <- nrow(y)
  p <- ncol(y)

  ver <- matrix(0,m,1)
  mu  <- matrix(mu,p,1)

  for (j in 1:m ){
    cc1=matrix(cc[j,],p,1)
    y1=matrix(y[j,],p,1)

    if(sum(cc1)==0){
      ver[j]<-dmvnorm(as.vector(y1), as.vector(mu), Sigma)
    }
    if(sum(cc1)>0){
      if(sum(cc1)==p){
        ver[j]<- pmnorm(as.vector(y1),as.vector(mu),Sigma)
      }
      else{
        muc<- mu[cc1==1]+Sigma[cc1==1,cc1==0]%*%solve(Sigma[cc1==0,cc1==0])%*%(y1[cc1==0]-mu[cc1==0])
        Sc<- Sigma[cc1==1,cc1==1]-Sigma[cc1==1,cc1==0]%*%solve(Sigma[cc1==0,cc1==0])%*%Sigma[cc1==0,cc1==1]
        ver[j]<-dmnorm(as.vector(y1[cc1==0]),as.vector(mu[cc1==0]),Sigma[cc1==0,cc1==0])*(pmnorm(as.vector(y1[cc1==1]),as.vector(muc),Sc))
      }
    }
  }

  return(ver)

}


### Mixturas

d.mixedNCens <- function(cc, y, pi1, mu, Sigma){
  #y: Matrix de dados m x p
  #pi1: deve ser do tipo vetor de dimensao g
  #mu: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
  #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p

  g    <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*dmvNCens(cc, y , mu[[j]], Sigma[[j]])
  return(dens)
}

################################################################################
### Densidades T Student e Misturas
################################################################################


dmvTCens<-function(cc, y, mu, Sigma, nu){

  # GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)

  m <- nrow(y)
  p <- ncol(y)

  ver<-matrix(0,m,1)

  mu<-matrix(mu,p,1)

  #    nu=ceiling(nu)             ###???

  for(j in 1:m ){

    cc1=matrix(cc[j,],p,1)
    y1=matrix(y[j,],p,1)

    if(sum(cc1)==0){
      ver[j]<- dmt(as.vector(y1),as.vector(mu),Sigma,df=nu)
    }
    if(sum(cc1)>0){

      if(sum(cc1)==p){

        ver[j]<-pmt(as.vector(y1),as.vector(mu),Sigma,df=nu)
      }

      else {
        nu1<-(nu+length(cc1[cc1==0]))
        muc<-mu[cc1==1]+Sigma[cc1==1,cc1==0]%*%solve(Sigma[cc1==0,cc1==0])%*%(y1[cc1==0]-mu[cc1==0])
        Sc <-Sigma[cc1==1,cc1==1]-Sigma[cc1==1,cc1==0]%*%solve(Sigma[cc1==0,cc1==0])%*%Sigma[cc1==0,cc1==1]
        Qy1<-t(y1[cc1==0]-mu[cc1==0])%*%solve(Sigma[cc1==0,cc1==0])%*%(y1[cc1==0]-mu[cc1==0])
        auxcte<-as.numeric((nu+Qy1)/(nu+length(cc1[cc1==0])))
        Sc22<-auxcte*Sc
        muUi<-muc
        SigmaUi<-Sc22
        SigmaUi<-(SigmaUi+t(SigmaUi))/2
        Sigma[cc1==0,cc1==0]<-(Sigma[cc1==0,cc1==0]+t(Sigma[cc1==0,cc1==0]))/2
        ver[j,]<-dmt(as.vector(y1[cc1==0]),as.vector(mu[cc1==0]),as.matrix(Sigma[cc1==0,cc1==0]),df=nu)*pmt(as.vector(y1[cc1==1]),as.vector(muUi),SigmaUi,df=nu1)


      }
    }
  }

  return(ver)
}


### Mixturas

d.mixedTCens <- function(cc, y, pi1, mu, Sigma,nu){
  #y: Matrix de dados m x p
  #pi1: deve ser do tipo vetor de dimensao g
  #mu: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
  #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p

  g <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*dmvTCens(cc, y , mu[[j]], Sigma[[j]], nu)
  return(dens)
}


################################################################################
## Momentos da Normal Truncada ##
################################################################################

MomemNT<-function(u=c(0,0),S=diag(2),qc=c(1,2)) {

  nic=length(u)

  if (nic==1) {

    qq <- (1/sqrt(S))*(-qc+u)
    R<-1
    alpha <- pnorm(-qq)

    dd <- dnorm(-qq)
    H <- qq*dd
    EX <- (1/alpha)*dd   # a vector with a length of nic
    EXX <- 1+1/alpha*H
    varX <- EXX-EX^2
    Eycens <- -sqrt(S)*EX+u
    varyic<- varX*S
    E2yy<-varyic+Eycens^2

  }

  else {

    qq <- diag(1/sqrt(diag(S)))%*%(-qc+u)
    R <-  diag(1/sqrt(diag(S)))%*%S%*%diag(1/sqrt(diag(S)))
    alpha <- pmvnorm(upper=as.vector(-qq), corr=R)
    #print(qq)
    dd <- rep(0, nic)   #derivative vector

    for (j in 1:nic){
      V <- R[-j, -j, drop=F]-R[-j,j, drop=F]%*%R[j,-j, drop=F]
      nu <- -qq[-j]+R[-j,j, drop=F]%*%qq[j]
      dd[j] <- dnorm(-qq[j])*pmvnorm(upper=as.vector(nu), sigma=V)
    }

    H <- matrix(rep(0, nic*nic), nrow=nic)
    RH <- matrix(rep(0, nic*nic), nrow=nic)

    if(nic==2)     {
      H[1,2] <- H[2,1] <- dmvnorm(-qq[c(1, 2)],sigma=matrix(c(1, R[1,2], R[2,1], 1), nrow=2))
      #sigma==R since qq is standardized
      RH[1,2] <- RH[2,1] <- R[1,2]*H[1,2]
    }

    else {
      for( s in 1:(nic-1)){
        for (t in (s+1):nic){
          invR <- solve(R[c(s,t), c(s,t), drop=F])
          nu <- -qq[-c(s,t)]+R[-c(s,t), c(s,t), drop=F]%*%invR%*%qq[c(s,t),,drop=F]
          V <-  R[-c(s,t), -c(s,t), drop=F]- R[-c(s,t), c(s,t), drop=F]%*%invR%*%R[c(s,t), -c(s,t), drop=F]
          H[s,t] <- H[t,s] <- pmvnorm(upper=as.vector(nu), sigma=V)*dmvnorm(-qq[c(s, t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))
          RH[s,t] <- RH[t,s] <- R[s,t]*H[s,t]
        }
      }
    }

    h <- qq*dd-apply(RH, 1, sum)
    diag(H) <- h
    EX <- (1/alpha)*R%*%dd   # a vector with a length of nic
    EXX <- R+1/alpha*R%*%H%*%R
    varX <- EXX-EX%*%t(EX)
    Eycens <- -diag(sqrt(diag(S)))%*%EX+u
    varyic <- diag(sqrt(diag(S)))%*%varX%*%diag(sqrt(diag(S)))
    E2yy <- varyic+Eycens%*%t(Eycens)

  }

  return(list(Ey=Eycens,Eyy=E2yy,Vary=varyic))

}


################################################################################
#################### Momentos da Distribui??o t Truncada #######################
################################################################################

####
cdfNI<-function(x,mu,sigma2,nu,type="Normal"){
  resp<-matrix(0,length(x),1)


  if(type=="Normal"){
    resp<-pnorm(x,mu,sqrt(sigma2))
  }

  if(type=="t"){
    z=(x-mu)/sqrt(sigma2)
    resp=pt(z,df=nu)
  }
  return(resp)
}

####

TT.moment = function(a,b,R,nu)
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)

  p = length(a)

  if(p==1){

    if(a== -Inf)  a <- -1e12
    if(b== Inf)   b <- 1e12

    G1<- 0.5*(gamma((nu-1)/2)*nu^(nu/2))/((cdfNI(b,0,1,nu,"t")-cdfNI(a,0,1,nu,"t"))*gamma(nu/2)*gamma(1/2))
    EX<- ((G1*((nu+a^2)^(-(nu-1)/2)-(nu+b^2)^(-(nu-1)/2))))
    EXX<- nu/(nu-2)+(G1*(a*(nu+a^2)^(-(nu-1)/2)-b*(nu+b^2)^(-(nu-1)/2)))

  }

  else{

    a = ifelse(a==-Inf,rep(-1e12,p),a)

    b = ifelse(b== Inf,rep( 1e12,p),b)

    al0 = pmvt(lower = a, upper = b, sigma = R, df = nu, algorithm = GB)[1]

    ### pdf & cdf

    la1 = (nu-2)/nu; la2 = (nu-4)/nu

    da = (nu-1)/(nu+a^2); db = (nu-1)/(nu+b^2)

    f1a = sqrt(la1)*dt(sqrt(la1)*a,df=nu-2)

    f1b = sqrt(la1)*dt(sqrt(la1)*b,df=nu-2)

    f2 = matrix(NA, p, p)

    G1a = G1b = rep(NA, p)

    G2 = matrix(NA, p, p)



    H = matrix(0,p,p)

    for(r in 1:(p-1))

    {

      temp = R[-r,r]

      S1 = R[-r,-r] - temp %*% t(R[r,-r])

      mua = temp * a[r]; low = a[-r]-mua; upp = b[-r]-mua

      G1a[r] = ifelse(p==2,pt(upp/sqrt(S1/da[r]),df=nu-1)-pt(low/sqrt(S1/da[r]),df=nu-1)

                      ,pmvt(lower = low, upper = upp, sigma = S1/da[r], df = nu-1, algorithm = GB)[1])

      mub = temp * b[r]; low = a[-r]-mub; upp = b[-r]-mub

      G1b[r] = ifelse(p==2,pt(upp/sqrt(S1/db[r]),df=nu-1)-pt(low/sqrt(S1/db[r]),df=nu-1)

                      ,pmvt(lower = low, upper = upp, sigma = S1/db[r], df = nu-1, algorithm = GB)[1])

      for(s in (r+1):p)

      {

        rs = c(r,s)

        pdf.aa = dmvt(c(a[r],a[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)

        pdf.ab = dmvt(c(a[r],b[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)

        pdf.ba = dmvt(c(b[r],a[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)

        pdf.bb = dmvt(c(b[r],b[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)

        if(p==2){cdf.aa=cdf.ab=cdf.ba=cdf.bb=1}

        if(p>2)

        {

          tmp = R[-rs,rs]%*%solve(R[rs,rs])

          mu.aa = c(tmp%*%c(a[r],a[s]))

          mu.ab = c(tmp%*%c(a[r],b[s]))

          mu.ba = c(tmp%*%c(b[r],a[s]))

          mu.bb = c(tmp%*%c(b[r],b[s]))

          daa = (nu-2)/(nu+(a[r]^2-2*R[r,s]*a[r]*a[s]+a[s]^2)/(1-R[r,s]^2))

          dab = (nu-2)/(nu+(a[r]^2-2*R[r,s]*a[r]*b[s]+b[s]^2)/(1-R[r,s]^2))

          dba = (nu-2)/(nu+(b[r]^2-2*R[r,s]*b[r]*a[s]+a[s]^2)/(1-R[r,s]^2))

          dbb = (nu-2)/(nu+(b[r]^2-2*R[r,s]*b[r]*b[s]+b[s]^2)/(1-R[r,s]^2))

          R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]

          cdf.aa = ifelse(p==3,pt((b[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)-pt((a[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)

                          ,pmvt(lower = a[-rs]-mu.aa, upper = b[-rs]-mu.aa, sigma = R21/daa, df=nu-2, algorithm = GB)[1])

          cdf.ab = ifelse(p==3,pt((b[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)-pt((a[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)

                          ,pmvt(lower = a[-rs]-mu.ab, upper = b[-rs]-mu.ab, sigma = R21/dab, df=nu-2, algorithm = GB)[1])

          cdf.ba = ifelse(p==3,pt((b[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)-pt((a[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)

                          ,pmvt(lower = a[-rs]-mu.ba, upper = b[-rs]-mu.ba, sigma = R21/dba, df=nu-2, algorithm = GB)[1])

          cdf.bb = ifelse(p==3,pt((b[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)-pt((a[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)

                          ,pmvt(lower = a[-rs]-mu.bb, upper = b[-rs]-mu.bb, sigma = R21/dbb, df=nu-2, algorithm = GB)[1])

        }

        H[r,s] = H[s,r] = pdf.aa*cdf.aa - pdf.ab*cdf.ab - pdf.ba*cdf.ba + pdf.bb*cdf.bb

      }

    }

    ##last part do loop
    r <- p
    temp = R[-r,r]

    S1 = R[-r,-r] - temp %*% t(R[r,-r])

    mua = temp * a[r]; low = a[-r]-mua; upp = b[-r]-mua

    G1a[r] = ifelse(p==2,pt(upp/sqrt(S1/da[r]),df=nu-1)-pt(low/sqrt(S1/da[r]),df=nu-1)

                    ,pmvt(lower = low, upper = upp, sigma = S1/da[r], df = nu-1, algorithm = GB)[1])

    mub = temp * b[r]; low = a[-r]-mub; upp = b[-r]-mub

    G1b[r] = ifelse(p==2,pt(upp/sqrt(S1/db[r]),df=nu-1)-pt(low/sqrt(S1/db[r]),df=nu-1)

                    ,pmvt(lower = low, upper = upp, sigma = S1/db[r], df = nu-1, algorithm = GB)[1])

    qa = f1a*G1a; qb = f1b*G1b

    EX = c(R %*% (qa-qb)) / al0 / la1

    H = H / la2

    D = matrix(0,p,p)

    diag(D) = a * qa - b * qb - diag(R%*%H)

    al1 = pmvt(lower = a, upper = b, sigma = R/la1, df=nu-2, algorithm = GB)[1]

    EXX = (al1 * R + R %*% (H + D) %*% R) / al0 / la1

  }

  return(list(EX=EX,EXX=EXX))

}


##Calculate the first to moments when mu not 0 and Sigma not R

Mtmvt <- function(mu,Sigma,nu,lower,upper){

  p=length(lower)

  if(p==1){

    if(lower== -Inf)  lower <- -1e12
    if(upper== Inf)  upper <- 1e12

    a1<-(lower-mu)/sqrt(Sigma)
    b1<-(upper-mu)/sqrt(Sigma)
    M <- TT.moment(a1, b1, 1, nu)
    Ey<- mu+sqrt(Sigma)*M$EX
    Eyy<- mu^2+Sigma*M$EXX+2*mu*sqrt(Sigma)*M$EX
    Vary<- Eyy - Ey^2
  }

  else{
    Lambda <- diag(1/sqrt(diag(Sigma)))

    if(length(which(upper == Inf)) != 0)  upper[which(upper == Inf)] <- 1e12
    b <- as.vector(diag(1/sqrt(diag(Sigma))) %*% (upper - mu))

    if(length(which(lower == -Inf)) != 0) lower[which(lower == -Inf)] <- -1e12
    a <- as.vector(diag(1/sqrt(diag(Sigma))) %*% (lower - mu))

    R <- Lambda %*% Sigma %*% Lambda

    M <- TT.moment(a, b, R, nu)

    Ey <- mu + solve(Lambda) %*% M$EX

    Eyy <- mu %*% t(mu) + solve(Lambda) %*% M$EX %*% t(mu) + mu %*% t(M$EX) %*% solve(Lambda) + solve(Lambda) %*% M$EXX %*% solve(Lambda)
    Vary<- Eyy- Ey%*%t(Ey)
  }

  return(list(Ey=Ey,Eyy=Eyy,Vary=Vary))
}

################################################################################
##
################################################################################
