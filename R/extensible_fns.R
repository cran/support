sp_seq <- function(D, nseq, ini=NA, num.rep=1,
                   dist.str=NA, dist.param=vector("list",p),
                   dist.samp=NA, scale.flg=T, bd=NA, 
                   num.subsamp=ifelse(any(is.na(dist.samp)),
                   max(10000,10*(nseq+nrow(D))),
                   min(10000,nrow(dist.samp))),
                   iter.max=max(200,iter.min), iter.min=50,
                   tol=1e-10, par.flg=TRUE){

  cur <- D
  p <- ncol(cur)
  dist.samp.o <- dist.samp
  
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
  
  # Set cores
  if (par.flg){
    num.cores <- parallel::detectCores()
  }else{
    num.cores <- 1
  }
  
  seq.lst <- vector("list",num.rep) # random restarts
  for (jj in 1:num.rep){
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
      if (any(is.na(ini))){
        if (p==1){
          ini <- matrix(randtoolbox::sobol(nseq,p,scrambling=T,seed=sample(1e6,1)),ncol=1)
        }else{
          # print(nseq)
          ini <- randtoolbox::sobol(nseq,p,scrambling=T,seed=sample(1e6,1))
        }
        # ini <- matrix(runif(nseq*p),ncol=p)
        
        ini.flg <- FALSE
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
      
      des <- sp_seq_cpp(cur,nseq,ini,dist.ind,dist.param,dist.samp,FALSE,
                        bd,num.subsamp,iter.max,iter.min,tol,num.cores)
      
    }else{
      
      # #Set subsample size
      # num.subsamp <- min(num.subsamp, nrow(dist.samp))
      
      #Standardize
      if (scale.flg==T){
        sdpts <- sqrt(apply(dist.samp,2,stats::var))
        mmpts <- apply(dist.samp,2,mean)
        dist.samp <- sweep(sweep(dist.samp,2,mmpts,"-"),2,sdpts,"/")
        cur <- sweep(sweep(cur,2,mmpts,"-"),2,sdpts,"/")
      }
      
      #Jitter if repeats
      dist.ind <- c(NA)
      dist.param <- list(NA)
      if (any(duplicated(dist.samp))){
        dist.samp <- jitter(dist.samp)
      }
      
      #Set initial points
      ini <- jitter(dist.samp[sample(1:nrow(dist.samp),nseq,F),])
      if (p==1){
        ini <- matrix(ini,ncol=1)
      }
      if (any(is.na(bd))){
        bd <- matrix(NA,nrow=p,ncol=2,byrow=T)
        for (i in 1:p){
          bd[i,] <- range(dist.samp[,i])
        }
      }
      
      # print(rbind(cur,ini))
      # print(dist.samp)
      
      # num.subsamp <- max(num.subsamp, 25*nseq)
      des <- sp_seq_cpp(cur,nseq,ini,dist.ind,dist.param,dist.samp,TRUE,
                        bd,num.subsamp,iter.max,iter.min,tol,num.cores)    

      #Scale back
      if (scale.flg==T){
        des <- sweep(sweep(des,2,sdpts,"*"),2,mmpts,"+")
      }
      # print(des)
    }
    
    # print(des)
    
    # seq <- des[(nrow(cur)+1):(nrow(cur)+nseq),]
    seq.lst[[jj]] <- des
    print('')
  }
  
  # Pick the best sequential points
  if (num.rep > 1){
    seq.crit <- rep(NA,num.rep)
    for (jj in 1:num.rep){
      seq.crit[jj] <- e_dist(seq.lst[[jj]],dist.str=dist.str,dist.param=dist.param,dist.samp=dist.samp.o)
    }
    des <- seq.lst[[which.min(seq.crit)]]
    seq <- des[(nrow(D)+1):(nrow(D)+nseq),]
  }else{
    des <- seq.lst[[1]]
    seq <- des[(nrow(D)+1):(nrow(D)+nseq),]
  }

  return(list(D=D,seq=seq))
  
}
