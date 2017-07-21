Cens.MmixEM <- function(cc, y, nu=4, mu=NULL, Sigma = NULL, pii = NULL, g = NULL, get.init = TRUE,
                      criteria = TRUE, group = FALSE, family = "Normal", error = 0.0001,
                      iter.max = 300, uni.Sigma = FALSE, obs.prob= FALSE, kmeans.param = NULL)
{
#mu, Sigma, shape devem ser do tipo list(). O numero de entradas no list eh o numero g de componentes de misturas
#cada entrada do list deve ser de tamanho igual ao numero de colunas da matriz de dados y

  y           <- as.matrix(y)
  dimnames(y) <- NULL
  if(is.list(family)==TRUE) {
    if(family[[1]] != family[[2]]) stop("Diferent family\n")
    family = family[[1]]
  }

  if(!is.matrix(y)) stop("The response is not in a matrix format\n")
  #  if(ncol(y) <= 1) stop("For the univariate case use the smsn.mix function\n")
  if((family != "t") && (family != "Normal")) stop(paste("Family",family,"not recognized.",sep=" "))
  if((length(g) == 0) && ((length(mu)==0) || (length(Sigma)==0) || (length(pii)==0)))  stop("The model is not specified correctly.\n")

  if(get.init == FALSE){
    g <- length(mu)
    if((length(mu) != length(Sigma)) || (length(mu) != length(pii))) stop("The size of the initial values are not compatibles.\n")
    if(sum(pii) != 1) stop("probability of pii does not sum to 1.\n")
    for (j in 1:g){
      dimnames(Sigma[[j]]) <- NULL
      names(mu[[j]])       <- NULL
    }
  }

  if((length(g)!= 0) && (g < 1)) stop("g must be greater than 0.\n")

  n   <- nrow(y)         ## recebe uma matrix
  p   <- ncol(y)
  if(uni.Sigma == FALSE) r <- g*(((p+1)*p/2)+p+((g-1)/g))
  if(uni.Sigma == TRUE)  r <- g*(((p+1)*p/(2*g))+p+((g-1)/g))

  if (get.init == TRUE){

    if(length(g) == 0) stop("g is not specified correctly.\n")

    k.iter.max <- 50
    n.start    <- 1
    algorithm  <- "Hartigan-Wong"

    if(length(kmeans.param) > 0){
      if(length(kmeans.param$iter.max) > 0 )  k.iter.max <- kmeans.param$iter.max
      if(length(kmeans.param$n.start) > 0 )   n.start    <- kmeans.param$n.start
      if(length(kmeans.param$algorithm) > 0 ) algorithm  <- kmeans.param$algorithm
    }

    if(g > 1){
      init <- kmeans(y,g,k.iter.max,n.start,algorithm)
      pii  <- init$size/dim(y)[1] # sum(pii)==0.5 tem que ser 1
      mu   <- Sigma <- list()
      for (s in 1:g){
        mu[[s]]    <- as.vector(init$centers[s,])
        Sigma[[s]] <- var(y[init$cluster == s,])
        dimnames(Sigma[[s]]) <- NULL
        names(mu[[s]])       <- NULL
      }

      #Ordenar por tamanho das componentes de pii
      if(family == "t"){
        if( sqrt(sum((order(pii)-g:1)*(order(pii)-g:1))) != 0 ){
          musi  <- mu
          Sigsi <- Sigma
          for(l in 1:g){
            mu[[l]]    <- musi[[(order(pii, decreasing = TRUE)[l])]]
            Sigma[[l]] <- Sigsi[[(order(pii, decreasing = TRUE)[l])]]
          }
          pii <- pii[order(pii, decreasing = TRUE)]
        }
      } else {
        if( sqrt(sum((order(pii)-g:1)*(order(pii)-g:1))) != 0 ){
          musi  <- mu
          Sigsi <- Sigma
          for(l in 1:g){
            mu[[l]]    <- musi[[(order(pii, decreasing = TRUE)[l])]]
            Sigma[[l]] <- Sigsi[[(order(pii, decreasing = TRUE)[l])]]
          }
          pii <- pii[order(pii, decreasing = TRUE)]
        }
      }

    } else {
      mu  <- Sigma <- list()
      pii <- 1
      mu[[1]]    <- as.vector(colMeans(y))
      Sigma[[1]] <- var(y)
      dimnames(Sigma[[1]]) <- NULL
      names(mu[[1]])       <- NULL
    }
  }

  ################################################################################
  ####### Inicia Familia Normal
  ################################################################################

  if (family == "Normal"){
  start.time  <- Sys.time()
    if(uni.Sigma){
      Sigma.uni <- Sigma[[1]]
      if(g > 1) for(k in 2:g) Sigma.uni <- Sigma.uni + Sigma[[k]]
      Sigma.uni <- Sigma.uni / g
      for(k in 1:g) Sigma[[k]] <- Sigma.uni
    }


    criterio <- 1
    count    <- 0
    lkante   <- 1

    while((criterio > error) && (count <= iter.max)){

      MI      <- matrix(0, ncol = r, nrow = r)
      alpha.g <- matrix(0, ncol = (p+1)*p/2, nrow = n)
      si.pi   <- si.mu <- alpha <- alpha.2 <- si <- si.g <- sipi <- NULL

      count <- count + 1
      #print(paste("Count",count,sep=" "))

      tal   <- matrix(0, n, g)
      tuyi  <- matrix(0,n,p)
      tui   <- matrix(0,n,1)
      tuyyi <- list()

      for (j in 1:g){

        #Ordenar por tamanho das componentes de pii
        if( sqrt(sum((order(pii)-g:1)*(order(pii)-g:1))) != 0 ){
          musi  <- mu
          Sigsi <- Sigma
          for(l in 1:g){
            mu[[l]]    <- musi[[(order(pii, decreasing = TRUE)[l])]]
            Sigma[[l]] <- Sigsi[[(order(pii, decreasing = TRUE)[l])]]
          }
          pii <- pii[order(pii, decreasing = TRUE)]
        }

        soma1 <- matrix(0,p,1)
        soma2 <- matrix(0,p,p)

        d1 <- dmvNCens(cc, y, mu[[j]], Sigma[[j]])
        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin

        d2 <-  d.mixedNCens(cc,y, pii, mu, Sigma)

        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

        tal[,j] <- d1*pii[j] / d2

        ### M-step: atualizar mu, Delta, Gama, sigma2 ###

        pii[j] <- (1/n)*sum(tal[,j])
        #        pii[g] <- 1 - (sum(pii) - pii[g])
        #######

        mus    <- mu[[j]]
        Sigmas <- Sigma[[j]]

        for (i in 1:n ){

          cc1 <- matrix(cc[i,],p,1)
          y1  <- matrix(y[i,],p,1)

          if(sum(cc1)==0){

            uy  <- matrix(y1,p,1)
            uyy <- y1%*%t(y1)

          }

          if(sum(cc1)>=1){

            if(sum(cc1)==p){
              muc<-mus
              Sc<-Sigmas
              aux<- MomemNT(muc,Sc,y1)
              uy<-aux$Ey
              uyy<- aux$Eyy
            }

            else {
              muc = mus[cc1==1]+Sigmas[cc1==1,cc1==0]%*%solve(Sigmas[cc1==0,cc1==0])%*%(y1[cc1==0]-mus[cc1==0])
              Sc <-Sigmas[cc1==1,cc1==1]-Sigmas[cc1==1,cc1==0]%*%solve(Sigmas[cc1==0,cc1==0])%*%Sigmas[cc1==0,cc1==1]
              Sc<-(Sc+t(Sc))/2
              aux <- MomemNT(muc,Sc,y1[cc1==1])
              uy <- matrix(y1,p,1)
              uy[cc1==1]<- aux$Ey
              uyy<-matrix(0,p,p)
              uyy[cc1==1,cc1==1]<-aux$Vary
              uyy<- uyy+uy%*%t(uy)

            }

          }

          soma1<- soma1 +  (tal[i,j])*uy
          soma2<- soma2 +  (tal[i,j])*(uyy-mus%*%t(uy)-(uy)%*%t(mus)+mus%*%t(mus))

          #uyi[j,]<-t(uy)

          #### Matrix de Info

          ################################################
          tuyi[i,]   <- t((tal[i,j])*uy)
          tui[i,]    <- tal[i,j]
          tuyyi[[i]] <- (tal[i,j])*(uyy)
          #################################################

          # Sigma[[j]] <- matrix.sqrt(Sigma[[j]]) %*% matrix.sqrt(Sigma[[j]])
          # Sigma[[j]] <- (Sigma[[j]]+t(Sigma[[j]]))/2


          if(uni.Sigma == FALSE){
            si.mu <- solve(Sigma[[j]]) %*% (tuyi[i,] - (tui[i]*mu[[j]]))
            si.pi <- (tal[i,j]/as.numeric(pii[j])) - (tal[i,g]/as.numeric(pii[g]))
            for(t in 1:((p+1)*p/2)) {
              alpha[t] <- -1/2*sum(diag( ( (tal[i,j]*solve(Sigma[[j]])) %*% deriv.sigma(Sigma[[j]],t,p) ) -
                                           ((tuyyi[[i]] - mu[[j]]%*%t(tuyi[i,]) - tuyi[i,]%*%t(mu[[j]]) + ((tui[i]*mu[[j]])%*%t(mu[[j]]))) %*%
                                              solve(Sigma[[j]]) %*% deriv.sigma(Sigma[[j]],t,p) %*% solve(Sigma[[j]])) ))

            }
            si <- rbind(si,c(si.mu, alpha, si.pi))
          } else {
            si.mu <- solve(Sigma[[j]]) %*% (tuyi[i,] - (tui[i]*mu[[j]]))
            si.pi <- (tal[i,j]/as.numeric(pii[j])) - (tal[i,g]/as.numeric(pii[g]))
            for(t in 1:((p+1)*p/2)) {
              alpha[t] <- -1/2*sum(diag( ( (tal[i,j]*solve(Sigma[[j]])) %*% deriv.sigma(Sigma[[j]],t,p) ) -
                                           ((tuyyi[[i]] - mu[[j]]%*%t(tuyi[i,]) - tuyi[i,]%*%t(mu[[j]]) + ((tui[i]*mu[[j]])%*%t(mu[[j]]))) %*%
                                              solve(Sigma[[j]]) %*% deriv.sigma(Sigma[[j]],t,p) %*% solve(Sigma[[j]])) ))
            }
            alpha.2 <- rbind(alpha.2,alpha)
            si <- rbind(si,c(si.mu,si.pi))
          }

        } ## end for de i

        if(uni.Sigma == TRUE) alpha.g <- alpha.g + alpha.2
        si.g    <- cbind(si.g, si)
        si      <- NULL
        alpha.2 <- NULL

        mu[[j]]    <- soma1 / sum(tal[,j])
        Sigma[[j]] <- soma2 / sum(tal[,j])
        Sigma[[j]] <- (Sigma[[j]]+t(Sigma[[j]]) )/2


        if(uni.Sigma == TRUE){
          GS <- 0
          for (w in 1:g) GS <- GS+kronecker(tal[,w],Sigma[[w]])
          Sigma.uni <- t(rowSums(array(t(GS),dim=c(p,p,n)),dims=2))/n
          for (w in 1:g) Sigma[[w]] <- Sigma.uni
        }


        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }


      }   ##FIM DO FOR  de g


      if(uni.Sigma == FALSE){
        si.g <- si.g[,-dim(si.g)[2]]
        if(g==3) si.g <- cbind(si.g[,-c((p+(p+1)*p/2 + 1),(p+(p+1)*p/2 + 1)*2)],si.g[,c((p+(p+1)*p/2 + 1),(p+(p+1)*p/2 + 1)*2)])
        if(g==2) si.g <- cbind(si.g[,-(p+(p+1)*p/2 + 1)],si.g[,(p+(p+1)*p/2 + 1)])
        for(v in 1:n) MI = MI + si.g[v,]%*%t(si.g[v,])
        colnames(MI) <- rownames(MI) <- nombres(g, p, order = "FALSE")
        #round(MI); print(paste("det(MI) =",det(MI),sep=" "))
        imm.sdev <- sqrt(diag(solve(MI)))
      } else {
        si.g <- si.g[,-dim(si.g)[2]]
        si.g <- cbind(si.g,alpha.g)
        for(v in 1:n) MI = MI + si.g[v,]%*%t(si.g[v,])
        colnames(MI) <- rownames(MI) <- nombres(g, p, uni.Sigma = "TRUE")
        #round(MI); print(paste("det(MI) =",det(MI),sep=" "))
        imm.sdev <- sqrt(diag(solve(MI)))
      }

      lk <- sum(log(d.mixedNCens(cc,y, pii, mu, Sigma) ))
      criterio <- abs((lk/lkante-1))

      lkante <- lk
      #print(paste("Criterio",criterio, sep=" "))

      if (criteria == TRUE){
        cl <- apply(tal, 1, which.max)
        #icl <- 0
        #for (j in 1:g) icl <- icl+sum(log(pii[j]*dmvSN(y[cl==j,], media[[j]], Sigma[[j]])))
        #icl <- -2*icl+(3*g-1)*log(n)

      }

    } # fim do while
    end.time        <- Sys.time()
    time.taken      <- end.time - start.time
  } # fim family Normal

  ################################################################################
  ####### Inicia Familia t
  ################################################################################

  if (family == "t"){
    start.time  <- Sys.time()
    if(uni.Sigma == TRUE){
      Sigma.uni <- Sigma[[1]]
      if(g > 1) for(k in 2:g) Sigma.uni <- Sigma.uni + Sigma[[k]]
      Sigma.uni <- Sigma.uni / g
      for(k in 1:g) Sigma[[k]] <- Sigma.uni
    }

    criterio <- 1
    count    <- 0
    lkante   <- 1

    while((criterio > error) && (count <= iter.max)){

      MI      <- matrix(0, ncol = r, nrow = r)
      alpha.g <- matrix(0, ncol = (p+1)*p/2, nrow = n)
      si.pi   <- si.mu <- alpha <- alpha.2 <- si <- si.g <- sipi <- NULL

      count <- count + 1
#      print(paste("Count",count,sep=" "))

      tal   <- matrix(0, n, g)
      tuyi  <- matrix(0,n,p)
      tui   <- matrix(0,n,1)
      tuyyi <- list()


      for (j in 1:g){

        #Ordenar por tamanho das componentes de pii
        if( sqrt(sum((order(pii)-g:1)*(order(pii)-g:1))) != 0 ){
          musi  <- mu
          Sigsi <- Sigma
          for(l in 1:g){
            mu[[l]]    <- musi[[(order(pii, decreasing = TRUE)[l])]]
            Sigma[[l]] <- Sigsi[[(order(pii, decreasing = TRUE)[l])]]
          }
          pii <- pii[order(pii, decreasing = TRUE)]
        }

        soma1 <- matrix(0,p,1)
        soma2 <- 0
        soma3 <- matrix(0,p,p)

        # print(j)
        d1 <- dmvTCens(cc, y, mu[[j]], Sigma[[j]],nu)

        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin # de donde pertenece

        d2 <-  d.mixedTCens(cc, y, pii, mu, Sigma, nu)

        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin # de donde pertenece

        tal[,j] <- d1*pii[j] / d2

        ### M-step: atualizar mu, Delta, Gama, sigma2 ###

        pii[j] <- (1/n)*sum(tal[,j])
        #       pii[g] <- 1 - (sum(pii) - pii[g])
        #######
        mus    <- mu[[j]]
        Sigmas <- Sigma[[j]]
        #######

        for (i in 1:n){
          #print(i)

          cc1 <- matrix(cc[i,],p,1)
          y1  <- matrix(y[i,],p,1)

          dm  <- t(y1-mus)%*%solve(Sigmas)%*%(y1-mus)
          cdm <- as.numeric((nu+p)/(nu+dm))

          if(sum(cc1)==0){

            tuy  <- (matrix(y1,p,1))*cdm
            tuyy <- (y1%*%t(y1))*cdm
            tu   <- cdm

            #ver[j,]<- dmvt(as.vector(y1),as.vector(mu),as.matrix(Sigma),df=nu,log=FALSE)
          }

          if(sum(cc1)>0){

            if(sum(cc1) == p){

              muUi     <- mus
              SigmaUi  <- Sigmas
              SigmaUiA <- SigmaUi*nu/(nu+2)
              auxupper <- y1-muUi
              # auxU1<-pmvt(upper = c(auxupper), sigma = SigmaUiA, df = nu+2,algorithm = GB)[1]
              auxU1    <- pmt(as.vector(y1),as.vector(muUi),SigmaUiA,df=nu+2)
              # auxU2<-pmvt(upper = c(auxupper), sigma = SigmaUi, df = nu,algorithm = GB) [1]
              auxU2    <- pmt(as.vector(y1),as.vector(muUi), SigmaUi,df=nu)
              MoMT     <- Mtmvt(muUi,SigmaUiA,nu+2,rep(-Inf,p),y1)
              U0       <- as.numeric(auxU1/auxU2)
              U1       <- auxU1/auxU2*MoMT$Ey
              U2       <- auxU1/auxU2*MoMT$Eyy

              tuy  <- U1
              tuyy <- U2
              tu   <- U0

              #   ver[j,]<-pmvt(upper = c(auxupper), sigma = SigmaUi, df = nu,algorithm = GB) [1]


            }

            else {

              PsiA <- Sigmas*nu/(nu+2)
              nu1  <- (nu+length(cc1[cc1==0]))

              muc  <- mus[cc1==1]+Sigmas[cc1==1,cc1==0]%*%solve(Sigmas[cc1==0,cc1==0])%*%(y1[cc1==0]-mus[cc1==0])
              Sc   <- Sigmas[cc1==1,cc1==1]-Sigmas[cc1==1,cc1==0]%*%solve(Sigmas[cc1==0,cc1==0])%*%Sigmas[cc1==0,cc1==1]
              Sc   <- (Sc+t(Sc))/2
              ScA  <- nu/(nu+2)*Sc

              Qy1  <- t(y1[cc1==0]-mus[cc1==0])%*%solve(Sigmas[cc1==0,cc1==0])%*%(y1[cc1==0]-mus[cc1==0])
              Qy2  <- t(y1[cc1==0]-mus[cc1==0])%*%solve(PsiA[cc1==0,cc1==0])%*%(y1[cc1==0]-mus[cc1==0])

              auxcte  <- as.numeric((nu+Qy1)/(nu+length(cc1[cc1==0])))
              auxcte1 <- as.numeric((nu+2+Qy2)/(nu+2+length(cc1[cc1==0])))

              Sc22 <- auxcte*Sc

              muUi    <- muc
              SigmaUi <- Sc22

              SigmaUiA <- auxcte1*ScA
              auxupper <- y1[cc1==1]-muUi


              #auxU1<-pmvt(upper = c(auxupper), sigma = SigmaUiA, df = nu1+2, algorithm = GB)[1]
              auxU1 <- pmt(as.vector(y1[cc1==1]), as.vector(muUi), SigmaUiA,df=nu1+2)

              #auxU2<-pmvt(upper = c(auxupper), sigma = SigmaUi, df = nu1, algorithm = GB)[1]
              auxU2 <- pmt(as.vector(y1[cc1==1]), as.vector(muUi), SigmaUi, df=nu1)

              MoMT <- Mtmvt(muUi,SigmaUiA,nu1+2,rep(-Inf,length(cc1[cc1==1])),y1[cc1==1])

              U0 <- as.numeric(auxU1/auxU2)/auxcte
              U1 <- (U0)*(MoMT$Ey)
              U2 <- (U0)*(MoMT$Eyy)

              Auxtuy <- (matrix(y1,p,1))

              tuy <- Auxtuy*U0
              tuy[cc1==1]<- U1

              tuyy <- (Auxtuy%*%t(Auxtuy))

              AAx <- tuyy[cc1==0,cc1==0]*U0
              ABx <- Auxtuy[cc1==0]%*%t(U1)
              BAx <- t(ABx)
              BBx <- U2

              tuyy[cc1==0,cc1==0] <- AAx
              tuyy[cc1==0,cc1==1] <- ABx
              tuyy[cc1==1,cc1==0] <- BAx
              tuyy[cc1==1,cc1==1] <- BBx


              tu <- U0
            }

          }

          soma1<- soma1 + (tal[i,j])*tuy
          soma2<- soma2 + (tal[i,j])*tu
          soma3<- soma3 + (tal[i,j])*(tuyy-(tuy)%*%t(mus)-(mus)%*%t(tuy)+ tu*mus%*%t(mus))

          #   tuyi[j,]<-t(tuy)
          #   tui[j]<-tu

          #### Matrix de Info

          ################################################
          tuyi[i,]   <- t((tal[i,j])*tuy)
          tui[i]     <- (tal[i,j])*tu
          tuyyi[[i]] <- (tal[i,j])*(tuyy)
          #################################################

          # Sigma[[j]] <- matrix.sqrt(Sigma[[j]]) %*% matrix.sqrt(Sigma[[j]])
          # Sigma[[j]] <- (Sigma[[j]]+t(Sigma[[j]]))/2


          if(uni.Sigma == FALSE){
            si.mu <- solve(Sigma[[j]]) %*% (tuyi[i,] - (tui[i]*mu[[j]]))
            si.pi <- (tal[i,j]/as.numeric(pii[j])) - (tal[i,g]/as.numeric(pii[g]))
            for(t in 1:((p+1)*p/2)) {
              alpha[t] <- -1/2*sum(diag( ( (tal[i,j]*solve(Sigma[[j]])) %*% deriv.sigma(Sigma[[j]],t,p) ) -
                                           ((tuyyi[[i]] - mu[[j]]%*%t(tuyi[i,]) - tuyi[i,]%*%t(mu[[j]]) + ((tui[i]*mu[[j]])%*%t(mu[[j]]))) %*%
                                              solve(Sigma[[j]]) %*% deriv.sigma(Sigma[[j]],t,p) %*% solve(Sigma[[j]])) ))

            }
            si <- rbind(si,c(si.mu, alpha, si.pi))
          } else {
            si.mu <- solve(Sigma[[j]]) %*% (tuyi[i,] - (tui[i]*mu[[j]]))
            si.pi <- (tal[i,j]/as.numeric(pii[j])) - (tal[i,g]/as.numeric(pii[g]))
            for(t in 1:((p+1)*p/2)) {
              alpha[t] <- -1/2*sum(diag( ( (tal[i,j]*solve(Sigma[[j]])) %*% deriv.sigma(Sigma[[j]],t,p) ) -
                                           ((tuyyi[[i]] - mu[[j]]%*%t(tuyi[i,]) - tuyi[i,]%*%t(mu[[j]]) + ((tui[i]*mu[[j]])%*%t(mu[[j]]))) %*%
                                              solve(Sigma[[j]]) %*% deriv.sigma(Sigma[[j]],t,p) %*% solve(Sigma[[j]])) ))
            }
            alpha.2 <- rbind(alpha.2,alpha)
            si <- rbind(si,c(si.mu,si.pi))
          }

        } ## end for de i

        if(uni.Sigma == TRUE) alpha.g <- alpha.g + alpha.2
        si.g    <- cbind(si.g, si)
        si      <- NULL
        alpha.2 <- NULL

        mu[[j]]    <- soma1 / soma2
        Sigma[[j]] <- soma3 / sum(tal[,j])
        Sigma[[j]] <- (Sigma[[j]]+t(Sigma[[j]]))/2


        if(uni.Sigma == TRUE){
          GS <- 0
          for (w in 1:g) GS <- GS + kronecker(tal[,w],Sigma[[w]])
          Sigma.uni <- t(rowSums(array(t(GS),dim=c(p,p,n)),dims=2))/n
          for (w in 1:g)  Sigma[[w]] <- Sigma.uni
        }


        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }



      }   ##FIM DO FOR  de g

      if(uni.Sigma == FALSE){
        si.g <- si.g[,-dim(si.g)[2]]
        if(g==3) si.g <- cbind(si.g[,-c((p+(p+1)*p/2 + 1),(p+(p+1)*p/2 + 1)*2)],si.g[,c((p+(p+1)*p/2 + 1),(p+(p+1)*p/2 + 1)*2)])
        if(g==2) si.g <- cbind(si.g[,-(p+(p+1)*p/2 + 1)],si.g[,(p+(p+1)*p/2 + 1)])
        for(v in 1:n) MI = MI + si.g[v,]%*%t(si.g[v,])
        colnames(MI) <- rownames(MI) <- nombres(g, p, order = "FALSE")
        #round(MI); print(paste("det(MI) =",det(MI),sep=" "))
        imm.sdev <- sqrt(diag(solve(MI)))
      } else {
        si.g <- si.g[,-dim(si.g)[2]]
        si.g <- cbind(si.g,alpha.g)
        for(v in 1:n) MI = MI + si.g[v,]%*%t(si.g[v,])
        colnames(MI) <- rownames(MI) <- nombres(g, p, uni.Sigma = "TRUE")
        #round(MI); print(paste("det(MI) =",det(MI),sep=" "))
        imm.sdev <- sqrt(diag(solve(MI)))
      }

      lk <- sum(log(d.mixedTCens(cc,y, pii, mu, Sigma, nu) ))
      criterio <- abs((lk/lkante-1))
      lkante <- lk
      #print(paste("Criterio", criterio, sep=" "))

      if (criteria == TRUE){
        cl <- apply(tal, 1, which.max)
        #icl <- 0
        #for (j in 1:g) icl <- icl+sum(log(pii[j]*dmvSN(y[cl==j,], media[[j]], Sigma[[j]])))
        #icl <- -2*icl+(3*g-1)*log(n)

      }

    } # fim do while
    end.time        <- Sys.time()
    time.taken      <- end.time - start.time
  } # fim family T

#---------------------------------------------------------------------------------
  ## Shape zero para contour
  shape <- list()
  for(b in 1:g) shape[[b]] <- rep(0,p)
  if(criteria == TRUE){

    if(uni.Sigma){
      if((family == "t") | (family == "Normal")) d <- g*p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)])  + (g-1) #mu + Sigma + pi
      else d <- g*2*p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)]) + (g-1) #mu + shape + Sigma + pi
    } else {
      if((family == "t") | (family == "Normal")) d <- g*(p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)]) ) + (g-1) #mu + shape + Sigma + pi
      else d <- g*(2*p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)]) ) + (g-1) #mu + shape + Sigma + pi
    }
    aic <- -2*lk + 2*d
    bic <- -2*lk + log(n)*d
    edc <- -2*lk + 0.2*sqrt(n)*d
    obj.out <- list(mu = mu, Sigma = Sigma, shape = shape, pii = pii, nu=nu, logLik = lk, aic = aic, bic = bic, edc = edc, iter = count, n = nrow(y), MI = MI, imm.sdev = imm.sdev, group = cl,time = time.taken, convergence = criterio<error)
  }




  if(criteria == FALSE) obj.out <- list(mu = mu, Sigma = Sigma, shape = shape, pii = pii, nu=nu, iter = count, n = nrow(y), MI = MI, imm.sdev = imm.sdev, group = apply(tal, 1, which.max),time = time.taken, convergence = criterio<error, iter=count)

  if (group == FALSE) obj.out <- obj.out[-length(obj.out)]
  obj.out$uni.Sigma <- uni.Sigma

  if (obs.prob == TRUE){
    nam <- c()
    for (i in 1:ncol(tal)) nam <- c(nam,paste("Group ",i,sep=""))
    dimnames(tal)[[2]] <- nam
    obj.out$obs.prob   <- tal
    if((ncol(tal) - 1) > 1) obj.out$obs.prob[,ncol(tal)] <- 1 - rowSums(obj.out$obs.prob[,1:(ncol(tal)-1)])
    else obj.out$obs.prob[,ncol(tal)] <- 1 - obj.out$obs.prob[,1]
    obj.out$obs.prob[which(obj.out$obs.prob[,ncol(tal)] <0),ncol(tal)] <- 0.0000000000
    obj.out$obs.prob <- round(obj.out$obs.prob,10)
  }

  class(obj.out) <- family

  obj.out
} #aqui  termina


