sp <- function(n, p, ini=NA,
                 dist.str=rep("uniform",p), dist.param=vector("list",p),
                 dist.samp=NA, scale.ind=T, num_subsamp=10000,
                 num_iter=200, tol_out=1e-10*sqrt(p), bd=NA){
  
  #Standard distributions
  if (is.na(dist.samp)){
    dist.samp <- matrix()
    #Encoding distribution string
    dist.vec <- c("uniform","normal","exponential","gamma","lognormal","student-t","weibull","cauchy","beta")
    dist.ind <- rep(NA,p)
    for (i in 1:p){
      dist.ind[i] <- which(dist.vec==dist.str[i])
      if (!any(dist.vec==dist.str[i])){
        stop("Please provide a valid distribution!")
      }
    }
    
    #Setting bounds and distribution parameters
    ini.flg <- TRUE
    if (is.na(ini)){
      if (p==1){
        ini <- matrix(randtoolbox::sobol(n,p,scrambling=T),ncol=1)
      }else{
        ini <- randtoolbox::sobol(n,p,scrambling=T)
      }
      
      ini.flg <- FALSE
    }
    for (i in 1:p){
      if (is.null(dist.param[[i]])){
        switch(dist.ind[i],
               "1" = {dist.param[[i]] <- c(0,1)}, #uniform
               "2" = {dist.param[[i]] <- c(0,1)}, #normal
               "3" = {dist.param[[i]] <- c(1)}, #exponential
               "4" = {dist.param[[i]] <- c(1,1)}, #gamma
               "5" = {dist.param[[i]] <- c(0,1)}, #lognormal
               "6" = {dist.param[[i]] <- c(1)}, #student-t
               "7" = {dist.param[[i]] <- c(1,1)}, #weibull
               "8" = {dist.param[[i]] <- c(0,1)}, #cauchy
               "9" = {dist.param[[i]] <- c(2,4)} #beta
        )
      }
      if (!ini.flg){
        switch(dist.ind[i],
               "1" = {ini[,i] <- stats::qunif(ini[,i], dist.param[[i]][1],dist.param[[i]][2])},
               "2" = {ini[,i] <- stats::qnorm(ini[,i], dist.param[[i]][1],dist.param[[i]][2])},
               "3" = {ini[,i] <- stats::qexp(ini[,i], dist.param[[i]][1])},
               "4" = {ini[,i] <- stats::qgamma(ini[,i], shape=dist.param[[i]][1], scale=dist.param[[i]][2])},
               "5" = {ini[,i] <- stats::qlnorm(ini[,i], dist.param[[i]][1],dist.param[[i]][2])},
               "6" = {ini[,i] <- stats::qt(ini[,i], df=dist.param[[i]][1])},
               "7" = {ini[,i] <- stats::qweibull(ini[,i], dist.param[[i]][1],dist.param[[i]][2])},
               "8" = {ini[,i] <- stats::qcauchy(ini[,i], dist.param[[i]][1],dist.param[[i]][2])},
               "9" = {ini[,i] <- stats::qbeta(ini[,i], dist.param[[i]][1],dist.param[[i]][2])}
        )
      }
      if (is.na(bd)){
        bd <- matrix(NA,nrow=p,ncol=2,byrow=T)
        switch(dist.ind[i],
               "1" = {bd[i,] <- c(0,1);},
               "2" = {bd[i,] <- c(-1e8,1e8);},
               "3" = {bd[i,] <- c(0,1e8);},
               "4" = {bd[i,] <- c(0,1e8);},
               "5" = {bd[i,] <- c(0,1e8);},
               "6" = {bd[i,] <- c(-1e8,1e8);},
               "7" = {bd[i,] <- c(0,1e8);},
               "8" = {bd[i,] <- c(-1e8,1e8);},
               "9" = {bd[i,] <- c(0,1);}
        )
      }
    }
    
    des <- sp_cpp(n,p,ini,dist.ind,dist.param,dist.samp,FALSE,
                  bd,num_subsamp,num_iter,tol_out,parallel::detectCores())
    
  }else{
    
    #Standardize
    if (scale.ind==T){
      sdpts <- sqrt(apply(dist.samp,2,stats::var))
      mmpts <- apply(dist.samp,2,mean)
      dist.samp <- sweep(sweep(dist.samp,2,mmpts,"-"),2,sdpts,"/")
    }
    
    dist.ind <- c(NA)
    dist.param <- list(NA)
    #Setting bounds and initial point set
    ini.flg <- TRUE
    if (is.na(ini)){
      nn <- min(nrow(dist.samp),1e6)
      ini <- stats::kmeans(dist.samp[sample(1:nrow(dist.samp),nn,F),],centers=n)$centers
      ini.flg <- FALSE
    }
    if (is.na(bd)){
      bd <- matrix(NA,nrow=p,ncol=2,byrow=T)
      for (i in 1:p){
        bd[i,] <- range(dist.samp[,i])
      }
    }
    
    des <- thincpp(dist.samp,ini,num_subsamp,bd,num_iter,1e4,1e-2*sqrt(ncol(dist.samp)),tol_out,1e-4,0,1)
    
    #Scale back
    if (scale.ind==T){
      ini <- sweep(sweep(ini,2,sdpts,"*"),2,mmpts,"+")
      des <- sweep(sweep(des,2,sdpts,"*"),2,mmpts,"+")
    }
  }
  
  #Compute support points
  # if (asymp){
  # des <- std_largep_cpp(n,p,num_iter,num_inn_iter,tol_inn,tol_out,eps,dist.ind,dist.param,0,1)
  # }
  # else{
  
  return(list(sp=des,ini=ini))
}