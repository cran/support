sp <- function(n, p, ini=NA, std.flg=T,
               dist.str=rep("normal",p), dist.param=vector("list",p),
               dist.samp=NA, scale.ind=T, bd=NA,
               num_subsamp=max(10000,50*n), max_iter=200, min_iter=50,
               tol_out=min(1e-2/n,1e-4)*p, warm.cl=F){
  
  #Standard distributions
  if (std.flg){
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
      ini.flg <- FALSE
      if (p==1){
        ini <- matrix(randtoolbox::sobol(n,p,scrambling=T,seed=sample(1e6,1)),ncol=1)
      }else{
        ini <- randtoolbox::sobol(n,p,scrambling=T,seed=sample(1e6,1))
      }
    }
    
    if (is.na(bd)){
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
    
    num_subsamp <- max(num_subsamp, 25*n)
    des <- sp_cpp(n,p,ini,dist.ind,dist.param,dist.samp,FALSE,
                  bd,num_subsamp,max_iter,min_iter,tol_out,parallel::detectCores())
    
  }else{
    
    #Set subsample size
    num_subsamp <- min(num_subsamp, nrow(dist.samp))
    
    #Standardize
    if (scale.ind==T){
      sdpts <- sqrt(apply(dist.samp,2,stats::var))
      mmpts <- apply(dist.samp,2,mean)
      dist.samp <- sweep(sweep(dist.samp,2,mmpts,"-"),2,sdpts,"/")
    }
    
    dist.ind <- c(NA)
    dist.param <- list(NA)
    #Setting bounds and initial point set
    if (is.na(ini)){
      if (warm.cl){
        nn <- min(nrow(dist.samp),1e6)
        ini <- stats::kmeans(dist.samp[sample(1:nrow(dist.samp),nn,F),],centers=n)$centers
      }else{
        ini <- jitter(dist.samp[sample(1:nrow(dist.samp),n,F),])
      }
    }else{
      ini <- sweep(sweep(ini,2,mmpts,"-"),2,sdpts,"/")
    }
    if (p==1){
      ini <- matrix(ini,ncol=1)
    }
    
    if (is.na(bd)){
      bd <- matrix(NA,nrow=p,ncol=2,byrow=T)
      for (i in 1:p){
        bd[i,] <- range(dist.samp[,i])
      }
    }
    
    num_subsamp <- max(num_subsamp, 25*n)
    des <- sp_cpp(n,p,ini,dist.ind,dist.param,dist.samp,TRUE,
                  bd,num_subsamp,max_iter,min_iter,tol_out,parallel::detectCores())
    
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
