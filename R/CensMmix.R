CensMmix <-  function(cc, y, nu=3, mu=NULL, Sigma = NULL, pii = NULL, g = NULL, get.init = TRUE,
                        criteria = TRUE, group = FALSE, family = "Normal", error = 0.0001,
                        iter.max = 300, uni.Sigma = FALSE, obs.prob= FALSE, kmeans.param = NULL)
{
  #Running the algorithm
  out <- Cens.MmixEM(cc,y,nu,mu, Sigma, pii, g, get.init, criteria, group, family, error,
                     iter.max, uni.Sigma, obs.prob, kmeans.param)

  #show result
  cat('\n')
  cat('---------------------------------------------------------\n')
  cat('        Censored multivariate finite mixture model       \n')
  cat('---------------------------------------------------------\n')
  cat('\n')
  cat('Family =',class(out))
  cat('\n')
  cat('Components =',g)
  cat('\n')
  cat('Row =',nrow(y),",","Columns=",ncol(y))
  cat('\n')
  cat('-----------\n')
  cat('Estimates\n')
  cat('-----------\n')
  cat('\n')
  #Mu vector
   namesrowMedj <- matrix("",ncol(y),g)
   ttable       <- data.frame(out$mu)
   for(i in 1:g){for(j in 1:ncol(y)){namesrowMedj[j,i]   <- paste("mu",i,j,sep="")}}
   muM <- matrix(0,ncol(y),1)
   for(j in 1:g)
   {
     muM <- as.matrix(ttable[,j])
     rownames(muM) <- c(namesrowMedj[,j])
     colnames(muM) <- paste("mu",j,sep="")
     print(muM)
     cat('\n')
   }
   cat('\n')
  #Sigma matirx
   for(k in 1:g)
   {
    a <- out$Sigma[[k]]
    a <- a[lower.tri(a, diag=TRUE)]
    b <- matrix(0, ncol(y), ncol(y))
    b[lower.tri(b, diag=TRUE)] <- a
    b <- t(b)
    b <- round(b,digits = 4)

    namesrowSigmas   <- c()
    for(j in 1:ncol(y)){namesrowSigmas[j] <- paste("sigma",j,sep="")}
    rownames(b) <- colnames(b) <- namesrowSigmas
    cat("Sigma",k)
    cat('\n')
    print(b)
    cat('\n')
   }

  if(family=="t")
  {
   cat("nu=",out$nu)
   cat('\n')
  }
  cat('------------------------\n')
  cat('Model selection criteria\n')
  cat('------------------------\n')
  cat('\n')
  critFin           <- c(out$logLik, out$aic, out$bic, out$edc)
  critFin           <- round(t(as.matrix(critFin)),digits=3)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","EDC"))
  print(critFin)
  cat('-------\n')
  cat('Details\n')
  cat('-------\n')
  cat('\n')
  cat('Iterations =',out$iter)
  cat('\n')
  cat("Processing time =",out$time,units(out$time))
  cat('\n')
  cat("Convergence =",out$convergence)
  cat('\n')
  res            <- list(SE = out$imm.sdev, iter = out$iter, mu=out$mu, Sigma=out$Sigma, pii=out$pii,  nu=out$nu, loglik=out$logLik, aic=out$aic, bic=out$bic, edc=out$edc, time = out$time, convergence = out$convergence)
  obj.out        <- list(res = res)
  class(obj.out) <-  "FMtMC"
  return(obj.out)
}


#CensMmix(cc, y, nu=3, mu=NULL, Sigma = NULL, pii = NULL, g = NULL, get.init = TRUE,criteria = TRUE, group = FALSE, family = "Normal", error = 0.0001,iter.max = 300, uni.Sigma = FALSE, obs.prob= FALSE, kmeans.param = NULL)
