e_dist <- function(D, dist.str=NA, dist.param=vector("list",ncol(D)),
                   nsamp=1e6, dist.samp=NA){
  
  p <- ncol(D)
  
  #Set std.flg (standard dist'n or data reduction)
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
    
    #Sample big sample
    bigsamp <- matrix(NA,nrow=nsamp,ncol=p)
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
      switch(dist.ind[i],
             "1" = {bigsamp[,i] <- stats::runif(nsamp, dist.param[[i]][1],dist.param[[i]][2])},
             "2" = {bigsamp[,i] <- stats::rnorm(nsamp, dist.param[[i]][1],dist.param[[i]][2])},
             "3" = {bigsamp[,i] <- stats::rexp(nsamp, dist.param[[i]][1])},
             "4" = {bigsamp[,i] <- stats::rgamma(nsamp, shape=dist.param[[i]][1], scale=dist.param[[i]][2])},
             "5" = {bigsamp[,i] <- stats::rlnorm(nsamp, dist.param[[i]][1],dist.param[[i]][2])},
             "6" = {bigsamp[,i] <- stats::rt(nsamp, df=dist.param[[i]][1])},
             "7" = {bigsamp[,i] <- stats::rweibull(nsamp, dist.param[[i]][1],dist.param[[i]][2])},
             "8" = {bigsamp[,i] <- stats::rcauchy(nsamp, dist.param[[i]][1],dist.param[[i]][2])},
             "9" = {bigsamp[,i] <- stats::rbeta(nsamp, dist.param[[i]][1],dist.param[[i]][2])}
      )
    }
    
    ret <- energycrit(bigsamp,D)
    
  }else{
    # dist.samp
    
    ret <- energycrit(dist.samp,D)
    
  }
  
  return(ret)
  
}