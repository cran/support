\name{sp}
\alias{sp}
\docType{package}
\title{
  Generating support points using difference-of-convex programming
}
\description{
  \code{sp} is the main function for generating the support points proposed in Mak and Joseph (2017) <\url{https://arxiv.org/abs/1609.01811}>. This implementation uses difference-of-convex optimization tools to efficiently generate support points with \eqn{n <= 10000} in up to \eqn{p <= 500} dimensions. To compute support points on standard distributions, the parameters \code{dist.str} and \code{dist.param} should be specified for the desired distribution and their parameters. To compute support points for big data reduction (e.g., Markov-chain Monte Carlo chains or large datasets), the data should be provided in the parameter \code{dist.samp}.
}
\usage{
  sp(n, p, ini=NA,
      dist.str=rep("uniform",p), dist.param=vector("list",p),
      dist.samp=NA, scale.ind=T, num_subsamp=10000,
      num_iter=200, tol_out=1e-10*sqrt(p), bd=NA)
}
\arguments{
  \item{n}{Number of support points.}
  \item{p}{Dimension of sample space.}
  \item{ini}{An \eqn{n} x \eqn{p} matrix for initializing the optimization algorithm.}
  \item{dist.str}{A \eqn{p}-length vector of strings for the desired marginal distributions (independence is assumed). Current choices for distributions include uniform, normal, exponential, gamma, lognormal, student-t, weibull, cauchy and beta. Specify \emph{only} for standard distributions.}
  \item{dist.param}{A \eqn{p}-length list for the desired distribution parameters. The parameter settings for each distribution are as follows:
  \itemize{
  \item \bold{Uniform}: Left bound, right bound;
  \item \bold{Normal}: Mean, standard deviation;
  \item \bold{Exponential}: Rate parameter;
  \item \bold{Gamma}: Shape parameter, scale parameter;
  \item \bold{Lognormal}: Log-mean, Log-standard deviation;
  \item \bold{Student-t}: Degree-of-freedom;
  \item \bold{Weibull}: Shape parameter, scale parameter;
  \item \bold{Cauchy}: Location parameter, scale parameter;
  \item \bold{Beta}: Shape parameter 1, shape parameter 2.
  }
Specify \emph{only} for standard distributions.
  }
  \item{dist.samp}{An \eqn{N} x \eqn{p} matrix for the MCMC chain (or large dataset) to reduce using support points. Specify \emph{only} for big data reduction.}
  \item{scale.ind}{Should the data be normalized to unit variance before processing? Specify \emph{only} for big data reduction.}
  \item{num_subsamp}{Batch sample size for resampling from desired distribution.}
  \item{num_iter}{Number of iterations for the optimization algorithm.}
  \item{tol_out}{Error tolerance for terminating the optimization algorithm.}
  \item{bd}{A \eqn{p} x \eqn{2} matrix for the lower and upper bounds of each distribution.}
}
\value{
  \item{sp}{An \eqn{n} x \eqn{p} matrix for \eqn{n} support points in \eqn{p} dimensions.}
  \item{ini}{An \eqn{n} x \eqn{p} matrix for the initialization point set used for optimization.}
}
\references{
Mak, S. and Joseph, V. R. (2017). Support points. \emph{Annals of Statistics}. To appear.
}
\examples{
  \dontshow{
    #Generate 25 support points for the 2-d i.i.d. N(0,1) distribution
    n <- 25 #number of points
    p <- 2 #dimension
    des <- sp(n,p,dist.str=rep("normal",p))
  }
  \donttest{
    ## Support points on standard distributions
    
    #Generate 25 support points for the 2-d i.i.d. N(0,1) distribution
    n <- 25 #number of points
    p <- 2 #dimension
    des <- sp(n,p,dist.str=rep("normal",p))
    
    x1 <- seq(-3.5,3.5,length.out=100) #Plot contours
    x2 <- seq(-3.5,3.5,length.out=100)
    z <- exp(-outer(x1^2,x2^2,FUN="+")/2)
    contour(x=x1,y=x2,z=z,drawlabels=F,nlevels=10)
    points(des$sp,pch=16,cex=1.25,col="red")
    
    #Generate 50 support points for the 2-d i.i.d. Beta(2,4) distribution
    n <- 50
    p <- 2
    dist.param <- vector("list",p)
    for (l in 1:p){
      dist.param[[l]] <- c(2,4)
    }
    des <- sp(n,p,dist.str=rep("beta",p),dist.param=dist.param)
    
    x1 <- seq(0,1,length.out=100) #Plot contours
    x2 <- seq(0,1,length.out=100)
    z <- matrix(NA,nrow=100,ncol=100)
    for (i in 1:100){
      for (j in 1:100){
        z[i,j] <- dbeta(x1[i],2,4) * dbeta(x2[j],2,4)
      }
    }
    contour(x=x1,y=x2,z=z,drawlabels=F,nlevels=10 )
    points(des$sp,pch=16,cex=1.25,col="red")
    
    #Generate 100 support points for the 5-d i.i.d. Exp(1) distribution
    n <- 100
    p <- 5
    des <- sp(n,p,dist.str=rep("exponential",p))
    pairs(des$sp,xlim=c(0,5),ylim=c(0,5),pch=16)
    
    ## Support points for big data reduction
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
    
    #Compute 100 support points
    n <- 100
    des <- sp(n,2,dist.samp=mcmc_r$trace,num_iter=1000)
    
    #Plot support points
    x1 <- seq(0,1,length.out=100) #contours
    x2 <- seq(0,1,length.out=100)
    z <- matrix(NA,nrow=100,ncol=100)
    for (i in 1:100){
      for (j in 1:100){
        z[i,j] <- franke2d(c(x1[i],x2[j]))
      }
    }
    contour(x=x1,y=x2,z=z,drawlabels=T,nlevels=10 )
    points(des$sp,pch=16,cex=1.25,col="red")

  }
}
