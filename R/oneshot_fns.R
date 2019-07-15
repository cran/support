sp <- function(n, p, ini=NA,
               dist.str=NA, dist.param=vector("list",p),
               dist.samp=NA, scale.flg=T, wts=NA, bd=NA,
               num.subsamp=ifelse(any(is.na(dist.samp)),
               max(10000,10*n),min(10000,nrow(dist.samp))),
               iter.max=max(250,iter.min), iter.min=50,
               tol=1e-10, par.flg=TRUE, n0=n*p){
  
  #Set std.flg (standard distn or data reduction)
  if (!any(is.na(dist.samp))&&any(is.na(dist.str))){
    # data reduction
    std.flg = F;
  }
  else if(any(is.na(dist.samp))&&!any(is.na(dist.str))){
    # standard dist'n
    std.flg = T;
  }
  else{
    stop("Exactly one of 'dist.samp' or 'dist.str' should be NA.")
  }
  
  #Set weights
  if (any(is.na(wts))){
    if (std.flg){
      wts <- rep(1.0,num.subsamp)
    }else{
      wts <- rep(1.0,nrow(dist.samp))
    }
  }else{
    if (std.flg){
      stop("wts must be NA for standard distributions.")
    }else{
      wts <- nrow(dist.samp)*wts;
    }
  }
  
  # Set cores
  if (par.flg){
    num.cores <- parallel::detectCores()
  }else{
    num.cores <- 1
  }
  # print(num.cores)
  
  #Standard distributions
  if (std.flg){
    dist.samp <- matrix(-1,nrow=2,ncol=2)
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
    if (any(is.na(ini))){
      ini.flg <- FALSE
      if (p==1){
        ini <- matrix(randtoolbox::sobol(n,p,scrambling=T,seed=sample(1e6,1)),ncol=1)
      }else{
        ini <- randtoolbox::sobol(n,p,scrambling=T,seed=sample(1e6,1))
      }
    }
    
    if (any(is.na(bd))){
      bd <- matrix(NA,nrow=p,ncol=2,byrow=T)
      bd.flg <- FALSE
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
      if (!bd.flg){
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
                  bd,num.subsamp,iter.max,iter.min,tol,num.cores,n0,wts)
    
  }else{
    
    #Standardize
    if (scale.flg==T){
      sdpts <- sqrt(apply(dist.samp,2,stats::var))
      mmpts <- apply(dist.samp,2,mean)
      dist.samp <- sweep(sweep(dist.samp,2,mmpts,"-"),2,sdpts,"/")
    }
    
    #Jitter if repeats
    if (any(duplicated(dist.samp))){
      dist.samp <- jitter(dist.samp)
    }
    
    dist.ind <- c(NA)
    dist.param <- list(NA)
    #Setting bounds and initial point set
    if (any(is.na(ini))){
      ini <- jitter(dist.samp[sample(1:nrow(dist.samp),n,F),])
    }else{
      ini <- sweep(sweep(ini,2,mmpts,"-"),2,sdpts,"/")
    }
    if (p==1){
      ini <- matrix(ini,ncol=1)
    }
    
    if (any(is.na(bd))){
      bd <- matrix(NA,nrow=p,ncol=2,byrow=T)
      for (i in 1:p){
        bd[i,] <- range(dist.samp[,i])
      }
    }
    
    des <- sp_cpp(n,p,ini,dist.ind,dist.param,dist.samp,TRUE,
                  bd,num.subsamp,iter.max,iter.min,tol,num.cores,n0,wts)
    
    #Scale back
    if (scale.flg==T){
      ini <- sweep(sweep(ini,2,sdpts,"*"),2,mmpts,"+")
      des <- sweep(sweep(des,2,sdpts,"*"),2,mmpts,"+")
    }
    
    
  }

  return(list(sp=des,ini=ini))
}
