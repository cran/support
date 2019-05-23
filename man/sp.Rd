\name{sp}
\alias{sp}
\docType{package}
\title{
  Computing support points using difference-of-convex programming
}
\description{
  \code{sp} is the main function for computing the support points in Mak and Joseph (2018). Current options include support points on standard distributions (specified via \code{dist.str}) or support points for reducing big data (specified via \code{dist.samp}). For big data reduction, weights on each data point can be specified via \code{wts}.
}
\usage{
  sp(n, p, ini=NA,
    dist.str=NA, dist.param=vector("list",p),
    dist.samp=NA, scale.flg=T, wts=NA, bd=NA,
    num.subsamp=max(10000,50*n), iter.max=max(250,iter.min), iter.min=50,
    tol.out=1e-6*sqrt(p), warm.cl=F, n0=n*p)
}
\arguments{
  \item{n}{Number of support points.}
  \item{p}{Dimension of sample space.}
  \item{ini}{An \eqn{n} x \eqn{p} matrix for initialization.}
  \item{dist.str}{A \eqn{p}-length string vector for desired distribution (assuming independence). Current options include uniform, normal, exponential, gamma, lognormal, student-t, weibull, cauchy and beta. Exactly one of \code{dist.str} or \code{dist.samp} should be \code{NA}.}
  \item{dist.param}{A \eqn{p}-length list for desired distribution parameters in \code{dist.str}. The following parameters are supported:
  \itemize{
  \item \bold{Uniform}: Minimum, maximum;
  \item \bold{Normal}: Mean, standard deviation;
  \item \bold{Exponential}: Rate parameter;
  \item \bold{Gamma}: Shape parameter, scale parameter;
  \item \bold{Lognormal}: Log-mean, Log-standard deviation;
  \item \bold{Student-t}: Degree-of-freedom;
  \item \bold{Weibull}: Shape parameter, scale parameter;
  \item \bold{Cauchy}: Location parameter, scale parameter;
  \item \bold{Beta}: Shape parameter 1, shape parameter 2.
  }
  }
  \item{dist.samp}{An \eqn{N} x \eqn{p} matrix for the big dataset (e.g., MCMC chain) to be reduced. Exactly one of \code{dist.str} or \code{dist.samp} should be \code{NA}.}
  \item{scale.flg}{Should the big data \code{dist.samp} be normalized to unit variance before processing?}
  \item{wts}{Weights on each data point in \code{dist.samp}. Uniform weights are assigned if \code{NA}.}
  \item{bd}{A \eqn{p} x \eqn{2} matrix for the lower and upper bounds of each variable.}
  \item{num.subsamp}{Batch size for resampling.}
  \item{iter.max}{Maximum iterations for optimization.}
  \item{iter.min}{Minimum iterations for optimization.}
  \item{tol.out}{Error tolerance for optimization.}
  \item{warm.cl}{Should the initial solution be pre-processed via clustering?}
  \item{n0}{Momentum parameter for optimization.}
}
\value{
  \item{sp}{An \eqn{n} x \eqn{p} matrix for support points.}
  \item{ini}{An \eqn{n} x \eqn{p} matrix for initial points.}
}
\references{
Mak, S. and Joseph, V. R. (2018). Support points. \emph{Annals of Statistics}, 46(6A):2562-2592.
}
\examples{
    ## Support points on distributions
    
    #Generate 25 SPs for the 2-d i.i.d. N(0,1) distribution
    n <- 25 #number of points
    p <- 2 #dimension
    D <- sp(n,p,dist.str=rep("normal",p))
    
    x1 <- seq(-3.5,3.5,length.out=100) #Plot contours
    x2 <- seq(-3.5,3.5,length.out=100)
    z <- exp(-outer(x1^2,x2^2,FUN="+")/2)
    contour.default(x=x1,y=x2,z=z,drawlabels=FALSE,nlevels=10)
    points(D$sp,pch=16,cex=1.25,col="red")
    
  \dontrun{
    #Generate 50 SPs for the 2-d i.i.d. Beta(2,4) distribution
    n <- 50
    p <- 2
    dist.param <- vector("list",p)
    for (l in 1:p){
      dist.param[[l]] <- c(2,4)
    }
    D <- sp(n,p,dist.str=rep("beta",p),dist.param=dist.param)
    
    x1 <- seq(0,1,length.out=100) #Plot contours
    x2 <- seq(0,1,length.out=100)
    z <- matrix(NA,nrow=100,ncol=100)
    for (i in 1:100){
      for (j in 1:100){
        z[i,j] <- dbeta(x1[i],2,4) * dbeta(x2[j],2,4)
      }
    }
    contour.default(x=x1,y=x2,z=z,drawlabels=FALSE,nlevels=10 )
    points(D$sp,pch=16,cex=1.25,col="red")
    
    #Generate 100 SPs for the 3-d i.i.d. Exp(1) distribution
    n <- 100
    p <- 3
    D <- sp(n,p,dist.str=rep("exponential",p))
    pairs(D$sp,xlim=c(0,5),ylim=c(0,5),pch=16)
    
    ## Support points for big data reduction: Franke's function
    library(MHadaptive)
    
    #Use modified Franke's function as posterior
    franke2d <- function(xx){
      if ((xx[1]>1)||(xx[1]<0)||(xx[2]>1)||(xx[2]<0)){
        return(-Inf)
      }
      else{
        x1 <- xx[1]
        x2 <- xx[2]
        
        term1 <- 0.75 * exp(-(9*x1-2)^2/4 - (9*x2-2)^2/4)
        term2 <- 0.75 * exp(-(9*x1+1)^2/49 - (9*x2+1)/10)
        term3 <- 0.5 * exp(-(9*x1-7)^2/4 - (9*x2-3)^2/4)
        term4 <- -0.2 * exp(-(9*x1-4)^2 - (9*x2-7)^2)
        
        y <- term1 + term2 + term3 + term4
        return(2*log(y))
      }
    }
    
    #Generate MCMC samples
    li_func <- franke2d #Desired log-posterior
    ini <- c(0.5,0.5) #Initial point for MCMc
    NN <- 1e5 #Number of MCMC samples desired
    burnin <- NN/2 #Number of burn-in runs
    mcmc_r <- Metro_Hastings(li_func, pars=ini, prop_sigma=0.05*diag(2),
                   iterations=NN, burn_in=burnin)
    
    #Compute n SPs
    n <- 100
    D <- sp(n,2,dist.samp=mcmc_r$trace)
    
    #Plot SPs
    par(mfrow=c(1,2))
    x1 <- seq(0,1,length.out=100) #contours
    x2 <- seq(0,1,length.out=100)
    z <- matrix(NA,nrow=100,ncol=100)
    for (i in 1:100){
      for (j in 1:100){
        z[i,j] <- franke2d(c(x1[i],x2[j]))
      }
    }
    plot(mcmc_r$trace,pch=4,col="gray",cex=0.75,
      xlab="",ylab="",xlim=c(0,1),ylim=c(0,1)) #big data
    points(D$sp,pch=16,cex=1.25,col="red")
    contour.default(x=x1,y=x2,z=z,drawlabels=TRUE,nlevels=10) #contour
    points(D$sp,pch=16,cex=1.25,col="red")

    ## Support points for big data: Rosenbrock distribution
    
    #Use Rosenbrock function as posterior
    rosen2d <- function(x) {
      B <- 0.03
      -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
    }
    
    #Generate MCMC samples
    li_func <- rosen2d #Desired log-posterior
    ini <- c(0,1) #Initial point for MCMc
    NN <- 1e5 #Number of MCMC samples desired
    burnin <- NN/2 #Number of burn-in runs
    mcmc_r <- Metro_Hastings(li_func, pars=ini, prop_sigma=0.25*diag(2),
                   iterations=NN, burn_in=burnin)
    
    #Compute n SPs
    n <- 100
    D <- sp(n,2,dist.samp=mcmc_r$trace,warm.cl=T)
    
    #Plot SPs
    x1 <- seq(-25,25,length.out=100) #contours
    x2 <- seq(-15,6,length.out=100)
    z <- matrix(NA,nrow=100,ncol=100)
    for (i in 1:100){
      for (j in 1:100){
        z[i,j] <- rosen2d(c(x1[i],x2[j]))
      }
    }
    plot(mcmc_r$trace,pch=4,col="gray",cex=0.75,
      xlab="",ylab="",xlim=c(-25,25),ylim=c(-15,6)) #big data
    points(D$sp,pch=16,cex=1.25,col="red")
    contour.default(x=x1,y=x2,z=z,drawlabels=TRUE,nlevels=10) #contour
    points(D$sp,pch=16,cex=1.25,col="red")
  }
}
