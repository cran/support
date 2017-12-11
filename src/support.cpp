// [[Rcpp::depends(BH)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <float.h>
#include <boost/any.hpp>
// #include <omp.h>
#include <time.h>
#include <limits.h>
#include <sstream>
#include <string>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace std;
using namespace arma;

double m_eps = pow(10.0,-10.0);

//-------------------------------------------------------------------------------
// Random number generators for different distributions
//-------------------------------------------------------------------------------
namespace sftrabbit {

template <typename RealType = double>
class beta_distribution
{
public:
  typedef RealType result_type;
  
  class param_type
  {
  public:
    typedef beta_distribution distribution_type;
    
    explicit param_type(RealType a = 2.0, RealType b = 2.0)
      : a_param(a), b_param(b) { }
    
    RealType a() const { return a_param; }
    RealType b() const { return b_param; }
    
    bool operator==(const param_type& other) const
    {
      return (a_param == other.a_param &&
              b_param == other.b_param);
    }
    
    bool operator!=(const param_type& other) const
    {
      return !(*this == other);
    }
    
  private:
    RealType a_param, b_param;
  };
  
  explicit beta_distribution(RealType a = 2.0, RealType b = 2.0)
    : a_gamma(a), b_gamma(b) { }
  explicit beta_distribution(const param_type& param)
    : a_gamma(param.a()), b_gamma(param.b()) { }
  
  void reset() { }
  
  param_type param() const
  {
    return param_type(a(), b());
  }
  
  void param(const param_type& param)
  {
    a_gamma = gamma_dist_type(param.a());
    b_gamma = gamma_dist_type(param.b());
  }
  
  template <typename URNG>
  result_type operator()(URNG& engine)
  {
    return generate(engine, a_gamma, b_gamma);
  }
  
  template <typename URNG>
  result_type operator()(URNG& engine, const param_type& param)
  {
    gamma_dist_type a_param_gamma(param.a()),
    b_param_gamma(param.b());
    return generate(engine, a_param_gamma, b_param_gamma); 
  }
  
  result_type min() const { return 0.0; }
  result_type max() const { return 1.0; }
  
  result_type a() const { return a_gamma.alpha(); }
  result_type b() const { return b_gamma.alpha(); }
  
  bool operator==(const beta_distribution<result_type>& other) const
  {
    return (param() == other.param() &&
            a_gamma == other.a_gamma &&
            b_gamma == other.b_gamma);
  }
  
  bool operator!=(const beta_distribution<result_type>& other) const
  {
    return !(*this == other);
  }
  
private:
  typedef std::gamma_distribution<result_type> gamma_dist_type;
  
  gamma_dist_type a_gamma, b_gamma;
  
  template <typename URNG>
  result_type generate(URNG& engine,
                       gamma_dist_type& x_gamma,
                       gamma_dist_type& y_gamma)
  {
    result_type x = x_gamma(engine);
    return x / (x + y_gamma(engine));
  }
};

template <typename CharT, typename RealType>
std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os,
                                      const beta_distribution<RealType>& beta)
{
  os << "~Beta(" << beta.a() << "," << beta.b() << ")";
  return os;
}

template <typename CharT, typename RealType>
std::basic_istream<CharT>& operator>>(std::basic_istream<CharT>& is,
                                      beta_distribution<RealType>& beta)
{
  std::string str;
  RealType a, b;
  if (std::getline(is, str, '(') && str == "~Beta" &&
      is >> a && is.get() == ',' && is >> b && is.get() == ')') {
    beta = beta_distribution<RealType>(a, b);
  } else {
    is.setstate(std::ios::failbit);
  }
  return is;
}

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector csample_num( NumericVector x,
                           int size,
                           bool replace) {
  //Sampling from vector
  NumericVector prob(x.size());
  for (int i=0;i<x.size();i++){
    prob[i] = 1.0/x.size();
  }
  NumericVector ret = Rcpp::RcppArmadillo::sample(x, size, replace, prob);
  return (ret);
}

//-------------------------------------------------------------------------------
// Functions for computing energy criterion
//-------------------------------------------------------------------------------

double sgn(double val) {
  int ret = (val > 0.0) - (val < 0.0);
  return( (double) ret);
}

// Computing the energy criterion using a Monte-Carlo sample

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
NumericVector energycrit(NumericMatrix& Rcpp_point, NumericMatrix& Rcpp_des, NumericMatrix& cn, int num_proc) {
  //Rcpp_point - Approximating points in projected subspace
  //Rcpp_des   - Design in full design space
  //cn         - Combination matrix
  //num_proc   - Number of processors available

  int num_combn = cn.nrow();
  NumericVector ret(num_combn);
  // omp_set_num_threads(num_proc);

  //Computes the energy criterion under L_pw norm

  int dim_num = Rcpp_point.ncol(); //dimension of projected subspace
  int point_num = Rcpp_point.nrow(); //number of approximating points
  int des_num = Rcpp_des.nrow(); //number of design points

  // #pragma omp parallel for
  for (int l=0; l<num_combn; l++){
    double runcrit = 0.0;
    for (int i=0; i<des_num; i++){
      for (int j=0; j<point_num; j++){
        double runcrittmp = 0.0;
        for (int k=0; k<dim_num; k++){
          runcrittmp += pow(abs(Rcpp_point(j,k) - Rcpp_des(i,cn(l,k))), 2.0);
        }
        runcrit += sqrt(runcrittmp);
      }
    }
    runcrit = runcrit * (2.0/(des_num*point_num));

    double runcrit2 = 0.0;
    for (int i=0; i<des_num; i++){
      for (int j=0; j<des_num; j++){
        double runcrittmp2 = 0.0;
        for (int k=0; k<dim_num; k++){
          runcrittmp2 += pow(abs(Rcpp_des(i,cn(l,k)) - Rcpp_des(j,cn(l,k))), 2.0);
        }
        runcrit2 += sqrt(runcrittmp2);
      }
    }
    ret(l) = runcrit2 / ((double)(des_num*des_num));
  }

  return(ret);

}

// A faster C++ implementation of L2-discrepancy

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
NumericVector starL2cpp(NumericMatrix& D, NumericMatrix& cn, int num_proc) {

  //int dimen = D.ncol();
  int n = D.nrow();
  int dimen = cn.ncol();
  NumericVector ret(cn.nrow());
  // omp_set_num_threads(num_proc);

// #pragma omp parallel for
  for (int l=0; l<cn.nrow(); l++){
    double term1 = pow(3.0,-dimen);
    double runsum = 0.0;
    for (int i=0; i<n; i++){
      double runprod = 1.0;
      for (int j=0; j<dimen; j++){
        runprod = runprod * (1 - pow(D(i,cn(l,j)-1),2) );
      }
      runsum += runprod;
    }
    double term2 = (pow(2.0,1-dimen)/n) * runsum;

    runsum = 0.0;
    for (int i=0; i<n; i++){
      for (int j=0; j<n; j++){
        double runprod = 1.0;
        for (int k=0; k<dimen; k++){
          runprod = runprod * (1-max(D(i,cn(l,k)-1),D(j,cn(l,k)-1)));
        }
        runsum += runprod;
      }
    }
    double term3 = (1/pow(n,2.0))*runsum;

    //             cout << term1 << ", " << term2 << ", " << term3 << endl;

    ret(l) = sqrt(term1-term2+term3);
  }

  return (ret);

}

// A faster C++ implementation of the energy criterion for the standard normal distribution

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
NumericVector energy_norm_cpp(NumericMatrix& yMat, NumericMatrix& cn, int num_proc)
{
  /*
  Computes the approximate eigenvalues of $W_infty$ using Monte-Carlo
  Modified from the R package 'energy' - Rizzo (2014)
  */

  //Vectorize yMat
  // int d = yMat.ncol(); //number of data points
  int n = yMat.nrow(); //number of clusters
  int d = cn.ncol();
  NumericVector retvec(cn.nrow());
  // omp_set_num_threads(num_proc);

// #pragma omp parallel for
  for (int l=0; l<cn.nrow(); l++){
    std::vector<double> y(d*n);
    for (int i=0; i<n; i++){
      for (int j=0; j<d; j++){
        y[i*d+j] = yMat(i,cn(l,j)-1);
      }
    }
    //    int    d=(*dim), n=(*nobs);
    int    i, j, k, p, maxterms=2000;
    double D=(double)d;
    double meanyy, meanzz;
    double delta, eps=1.0e-7;
    double normy, yy, dif, sum, sum0, term;
    double lg0, lg1,logak, loggk;

    //Second mean
    lg0 = lgamma(D/2.0);
    lg1 = lgamma((D+1.0)/2.0);
    meanzz = 2.0 * exp(lg1 - lg0);

    //Compting the vector of first mean as series
    arma::vec meanyz(n);
    for (i=0; i<n; i++) {
      yy = 0.0;
      for (j=0; j<d; j++) {
        dif = y[i*d+j] * y[i*d+j];
        yy += dif;
      }
      normy = sqrt(yy);
      delta = 1.0;
      sum = 0.0;
      k = 0;
      while (delta > eps && k < maxterms) {
        sum0 = sum;
        logak = (k+1)*log(yy) - lgamma(k+1) - k*M_LN2 -
          log(2*k+1) - log(2*k+2);
        loggk = lg1 + lgamma(k+1.5) - lgamma(k+D/2+1);
        term = exp(logak + loggk);
        if (k % 2 == 0)
          sum += term;
        else
          sum -= term;
        delta = fabs(sum - sum0);
        k++;
      }
      if (delta < eps)
        meanyz(i) = meanzz/M_SQRT2 + M_SQRT_2dPI * sum;
      else {
        // std::cout << 'I got here!' << std::endl;
        meanyz(i) = normy;
        //Rf_warning('E|y-Z| did not converge, replaced by %f', normy);
      }
    }

    //Computing the entire kernel matrix
    double ret = 0.0;
    arma::vec normyvec(n);
    for (int i=0; i<n; i++){
      for (int j=0; j<n; j++){
        yy = 0.0;
        for (k=0; k<d; k++) {
          dif = pow( y[i*d+k] - y[j*d+k], 2.0);
          yy += dif;
        }
        normy = sqrt(yy);
        //      if (i == j){
        //        normyvec[i] = normy;
        //      }
        ret += meanyz(i) + meanyz(j) - meanzz - normy;
      }
    }
    retvec(l) = ret / pow((double)n,2.0);
  }

  return (retvec);
}

// [[Rcpp::export]]
double obj_qsp(arma::vec& des, arma::mat& distsamp, double q){
  //des - vectorized design
  int NN = distsamp.n_rows;
  int pp = distsamp.n_cols;
  int nn = des.n_elem/pp;
  double runsum = 0.0;
  double runsum2 = 0.0;
  double tmp = 0.0;
  double ret = 0.0;
  
  for (int i=0; i<nn; i++){
    for (int m=0; m<NN; m++){
      tmp = 0.0;
      for (int j=0; j<pp; j++){
        tmp += pow( des(i*pp+j)-distsamp(m,j), 2);
      }
      runsum += pow(tmp, q/2.0);
    }
  }
  runsum = 2.0/( (double) nn * NN ) * runsum;
  // cout << "runsum: " << runsum << endl;
  
  for (int i=0; i<nn; i++){
    for (int k=0; k<nn; k++){
      tmp = 0.0;
      for (int j=0; j<pp; j++){
        tmp += pow( des(i*pp+j)-des(k*pp+j), 2);
      }
      runsum2 += pow(tmp, q/2.0);
    }
  }
  runsum2 = runsum2/( (double) pow(nn,2.0));
  // cout << "runsum2: " << runsum2 << endl;
  
  ret = runsum - runsum2;
  return(ret);
}

// [[Rcpp::export]]
arma::vec grad_qsp(arma::vec& des, arma::mat& distsamp, double q){
  int NN = distsamp.n_rows;
  int pp = distsamp.n_cols;
  int nn = des.n_elem/pp;
  double tmp = 0.0;
  arma::vec ret(nn*pp); //point-by-point
  arma::vec ret2(nn*pp);
  for (int i=0; i<(nn*pp); i++){
    ret(i) = 0.0;
    ret2(i) = 0.0;
  }
  arma::vec tmpvec(pp); //for a single point
  
  for (int i=0; i<nn; i++){
    for (int m=0; m<NN; m++){
      
      //Reset tmpvec
      for (int j=0; j<pp; j++){
        tmpvec(j) = 0.0;
      }
      tmp = 0.0;
      
      for (int j=0; j<pp; j++){
        tmpvec(j) += des(i*pp+j)-distsamp(m,j);
        tmp += pow( tmpvec(j), 2);
      }
      tmpvec = tmpvec / pow(tmp, (2.0-q)/2.0);
      
      //Update ret
      for (int j=0; j<pp; j++){
        ret(i*pp+j) += tmpvec(j);
      }
      
    }
  }
  ret = (2.0 * q)/( (double) nn * NN ) * ret;
  
  for (int i=0; i<nn; i++){
    for (int k=0; k<nn; k++){
      if (k != i){
        //Reset tmpvec
        for (int j=0; j<pp; j++){
          tmpvec(j) = 0.0;
        }
        tmp = 0.0;
        
        for (int j=0; j<pp; j++){
          tmpvec(j) += des(i*pp+j)-des(k*pp+j);
          tmp += pow( tmpvec(j), 2 );
        }
        tmpvec = tmpvec / pow(tmp, (2.0-q)/2.0);
        
        //Update ret
        for (int j=0; j<pp; j++){
          ret2(i*pp+j) += tmpvec(j);
        }
      }
    }
  }
  ret2 = (2.0 * q) / ( (double) nn * nn ) * ret2;
  
  // cout << ret(0) << endl;
  // cout << ret2(0) << endl;
  // cout << tmpvec(0) << endl;
  return(ret-ret2);
}

//-------------------------------------------------------------------------------
// SP and PSP for thinning
//-------------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix thincpp(NumericMatrix& Rcpp_point, arma::mat& Rcpp_inides, int num_subsamp,
             NumericMatrix& bound, int it_max, int inn_it_max, double innertol, double outertol, double epsilon, double rate,
             int num_proc) 
  {
  // Description of Inputs:
  // Rcpp_point          - Sample data
  // Rcpp_inides         - Initial design
  // bound               - Upper and lower bound
  // it_max              - Maximum number of iterations
  // tol                 - Tolerance for stopping blockwise updates
  
  //Closed-form updates for p == 2
  
  int it_num = 0; //keeps track of iteration number
  int dim_num = Rcpp_point.ncol(); //dimension of data points
  int point_num = Rcpp_point.nrow(); //number of data points
  int des_num = Rcpp_inides.n_rows; //number of clusters
  bool cont = true; //stopping boolean for outside BCD loop
  bool cont2 = true; //stopping boolean for inside CCCP loop
  
  //Containers for optimization
  std::vector<double> prevdes(des_num*dim_num);
  std::vector<double> xprime(dim_num);
  std::vector<double> tmpvec(dim_num);
  // omp_set_num_threads(num_proc);
  
  //Vectorize design and points  
  std::vector<double> des(des_num*dim_num);
  for (int i=0; i<des_num; i++){
    for (int j=0; j<dim_num; j++){
      des[j+i*dim_num] = Rcpp_inides(i,j);
    }
  }
  std::vector<double> point(point_num*dim_num);
  for (int i=0; i<point_num; i++){
    for (int j=0; j<dim_num; j++){
      point[j+i*dim_num] = Rcpp_point(i,j);
    }
  }
  
  //Blockwise coordinate descent
  while (cont){
    
    Rcout << "SP: Iteration " << it_num << "/" << it_max << endl;
    time_t start = time(0);
    
    //Update innertol
    // innertol = innertol * exp(-(double)rate*it_num);
    
    //Update prevdes
    for (int i=0; i<des_num; i++){
      for (int j=0; j<dim_num; j++){
        prevdes[j+i*dim_num] = des[j+i*dim_num];
      }
    }
    
    //BCD for each point
    for (int m=0; m<des_num; m++){
      
      //        Rcout << "Point: " << m << endl;
      
      //Copy current design point as initial point
      for (int n=0; n<dim_num; n++){
        xprime[n] = des[n+m*dim_num];
      }
      
      //Reset cont2 flag
      cont2 = true;
      int it_num2 = 0;
      
      //Concave-convex updates until convergence
      while( (cont2) && (it_num2 <= inn_it_max) ){
        
        //          Rcout << "Inner iteration: " << it_num2 << endl;
        
        //Reset xprime
        for (int n=0; n<dim_num; n++){
          xprime[n] = 0.0;
        }
        
        //Update using closed-form formula
        //... for the design
        for (int o=0; o<des_num; o++){
          
          if (o != m){
            double tmptol = 0.0; //Running total for norm in denominator
            for (int n=0; n<dim_num; n++){
              tmpvec[n] = des[n+m*dim_num] - des[n+o*dim_num];
              tmptol += pow(tmpvec[n],2);
            }
            tmptol = sqrt(tmptol+epsilon);
            
            for (int n=0; n<dim_num; n++){
              xprime[n] += tmpvec[n]/tmptol;
            }
          }
          
        }
        
        for (int n=0; n<dim_num; n++){
          xprime[n] = xprime[n] * (num_subsamp/des_num);
        }
        
        //... for the samples
        
        // Subsample:
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, point_num-1);
        ivec samp(num_subsamp);
        for (int o=0; o<num_subsamp; o++){
          samp(o) = dis(gen);
        }
        
        // Compute sample-level 
        double tmpconst = 0.0; //For the inverse distance sum
        for (int o=0; o<num_subsamp; o++){
          
          double tmptol = 0.0; //Running total for norm in denominator
          for (int n=0; n<dim_num; n++){
            tmpvec[n] = pow(point[n+samp(o)*dim_num] - des[n+m*dim_num],2);
            tmptol += tmpvec[n];
          }
          //            Rcout << "tmptol: " << tmptol << endl;
          tmptol = sqrt(tmptol+epsilon);
          tmpconst += 1.0/tmptol;
          
          for (int n=0; n<dim_num; n++){
            xprime[n] += point[n+samp(o)*dim_num]/tmptol;
          }
          
        } 
        
        //Scale by inverse distances
        for (int n=0; n<dim_num; n++){
          xprime[n] = xprime[n]/tmpconst;
        }    
        
        //Update bounds
        for (int n=0; n<dim_num; n++){
          xprime[n] = min( max( xprime[n], bound(n,0)), bound(n,1));
        }
        
        //Compute running tolerance
        double runtol = 0.0;
        for (int n=0; n<dim_num; n++){
          runtol += pow(des[n+m*dim_num] - xprime[n],2);
        }
        //          Rcout << "runtol: " << runtol << endl;
        runtol = sqrt(runtol);
        
        //Update point
        for (int n=0; n<dim_num; n++){
          des[n+m*dim_num] = xprime[n];
        }
        
        //          Rcout << "runtol: " << runtol << endl;; //Temporary
        
        //Update cont2 flag
        if (runtol < innertol){
          cont2 = false;
        }
        //          cont2 = false; //////
        it_num2 ++;
      }
      
      //Update design point 
      for (int n=0; n<dim_num; n++){
        des[n+m*dim_num] = xprime[n];
      }
      
    }
    
    //Output time
    time_t end = time(0);
    double sec = difftime(end, start);
    Rcout << "Iteration time: " << sec * 1000.0 << endl;
    
    //Increment and update cont flag 
    it_num++;
    //Check maximum iterations
    if (it_num >= it_max){
      cont = false;
    }
    //Check convergence
    double rundiff = 0.0;
    for (int n=0; n<des_num; n++){
      for (int o=0; o<dim_num; o++){
        rundiff += abs(des[o+n*dim_num]-prevdes[o+n*dim_num]);
      }
    }
    if (rundiff < outertol){
      cont = false;
    }
    
    //      cont = false;//////
  }
  
  //Output the final NumericMatrix 
  NumericMatrix retdes(des_num, dim_num);
  for (int i=0; i<des_num; i++){
    for (int j=0; j<dim_num; j++){
      retdes(i,j) = des[j+i*dim_num];
    }
  }
  return(retdes);
  
}

//-------------------------------------------------------------------------------
// SP and PSP for standard distributions (including asymptotic approximation for large p)
//-------------------------------------------------------------------------------
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
NumericMatrix sp_cpp(int des_num, int dim_num, NumericMatrix& ini,
                             NumericVector& distind, List distparam, 
                             NumericMatrix& distsamp, bool thin, NumericMatrix& bd,
                             int point_num, int it_max, double tol, int num_proc){
  
  // Description of Inputs:
  // Rcpp_point          - Sample data
  // Rcpp_inides         - Initial design
  // bound               - Upper and lower bound
  // it_max              - Maximum number of iterations
  // tol                 - Tolerance for stopping blockwise updates
  
  int it_num = 0; //keeps track of iteration number
  bool cont = true; //stopping boolean for outside BCD loop
  bool cont2 = true; //stopping boolean for inside CCCP loop
  // Rcout << "Hi 1" << endl;
  
  //  Containers for optimization and sampling
  std::vector<double> prevdes(des_num*dim_num);
  // std::vector<double> tmpvec(dim_num);
  std::default_random_engine generator;
  int distint = 0;
  omp_set_num_threads(num_proc);
  
  //Vectorize design and points  
  std::vector<double> des(des_num*dim_num);
  for (int i=0; i<des_num; i++){
    for (int j=0; j<dim_num; j++){
      des[j+i*dim_num] = ini(i,j);
    }
  }
  
  //CCP
  while (cont){
    
    Rcout << "SP: Iteration " << it_num << "/" << it_max << endl;
    // time_t start = time(0);
  
    //Update prevdes
    for (int i=0; i<des_num; i++){
      for (int j=0; j<dim_num; j++){
        prevdes[j+i*dim_num] = des[j+i*dim_num];
      }
    }
    
    // Closed-form updates
    
    // Generate sample from F
    arma::mat rnd(point_num,dim_num);
    if (thin){
      std::default_random_engine generator;
      std::uniform_int_distribution<int> uddist(0,distsamp.nrow()-1);
      for (int i=0; i<point_num; i++){
        for (int j=0; j<dim_num; j++){
          rnd(i,j) = distsamp(uddist(generator),j);
        }
      }
    }else{
      for (int n=0; n<dim_num; n++){
        distint = distind(n);
        switch( distint ){
        case 1:
        {
          SEXP unill = distparam[n];
          NumericVector uniyy(unill);
          std::uniform_real_distribution<double> uni_dist (uniyy(0),uniyy(1));
          for (int o=0; o<point_num; o++){
            rnd(o,n) = uni_dist(generator);
          }
          break;
        }
        case 2:
        {
          SEXP normll = distparam[n]; 
          NumericVector normyy(normll);  
          std::normal_distribution<double> norm_dist (normyy(0),normyy(1));
          for (int o=0; o<point_num; o++){
            rnd(o,n) = norm_dist(generator);
          }
          break;
        }
        case 3:
        {
          SEXP expll = distparam[n]; 
          NumericVector expyy(expll);  
          std::exponential_distribution<double> exp_dist (expyy(0));
          for (int o=0; o<point_num; o++){
            rnd(o,n) = exp_dist(generator);
          }
          break;
        }
        case 4:
        {
          SEXP gamll = distparam[n]; 
          NumericVector gamyy(gamll);  
          std::gamma_distribution<double> gam_dist (gamyy(0),gamyy(1));
          for (int o=0; o<point_num; o++){
            rnd(o,n) = gam_dist(generator);
          }
          break;
        }
        case 5:
        {
          SEXP lnll = distparam[n]; 
          NumericVector lnyy(lnll);  
          std::lognormal_distribution<double> ln_dist (lnyy(0),lnyy(1));
          for (int o=0; o<point_num; o++){
            rnd(o,n) = ln_dist(generator)/ ( exp(pow(lnyy(1),2.0)-1.0) * exp(2*lnyy(0)+pow(lnyy(1),2.0)) );
          }
          break;
        }
        case 6:
        {
          SEXP tll = distparam[n]; 
          NumericVector tyy(tll);  
          std::student_t_distribution<double> t_dist (tyy(0));
          for (int o=0; o<point_num; o++){
            rnd(o,n) = t_dist(generator);
          }
          break;
        }
        case 7:
        {
          SEXP wbll = distparam[n]; 
          NumericVector wbyy(wbll);  
          std::weibull_distribution<double> wb_dist (wbyy(0),wbyy(1));
          for (int o=0; o<point_num; o++){
            rnd(o,n) = wb_dist(generator);
          }
          break;
        }
        case 8:
        {
          SEXP cauchll = distparam[n]; 
          NumericVector cauchyy(cauchll);  
          std::cauchy_distribution<double> cauch_dist (cauchyy(0),cauchyy(1));
          for (int o=0; o<point_num; o++){
            rnd(o,n) = cauch_dist(generator);
          }
          break;
        }
        case 9:
        {
          SEXP betall = distparam[n]; 
          NumericVector betayy(betall);  
          sftrabbit::beta_distribution<> beta_dist (betayy(0),betayy(1));
          for (int o=0; o<point_num; o++){
            rnd(o,n) = beta_dist(generator);
          }
          break;
        }
        }
      }
    }
    
    
    //Parallelize computation
    #pragma omp parallel for
    for (int m=0; m<des_num; m++){
      
      arma::vec xprime(dim_num);
      xprime.fill(0.0);
      arma::vec tmpvec(dim_num);
      tmpvec.fill(0.0);
      
      // if (m == 0){
      //   Rcout << xprime(0) << endl;
      // }
      
      //Summation for design part
      for (int o=0; o<des_num; o++){
        
        if (o != m){
          double tmptol = 0.0; //Running total for norm in denominator
          for (int n=0; n<dim_num; n++){
            tmpvec(n) = prevdes[n+m*dim_num] - prevdes[n+o*dim_num];
            tmptol += pow(tmpvec(n),2.0);
          }
          tmptol = sqrt(tmptol);
          
          for (int n=0; n<dim_num; n++){
            xprime(n) += tmpvec(n)/tmptol;
          }
        }
      }
      
      // if (m == 0){
      //   Rcout << xprime(0) << endl;
      // }
      
      for (int n=0; n<dim_num; n++){
        xprime(n) = xprime(n) * (point_num/des_num);
      }
      
      // if (m == 0){
      //   Rcout << xprime(0) << endl;
      // }
      
      //Summation for sample side
      double tmpconst = 0.0; //Running total for inverse distance sum
      for (int o=0; o<point_num; o++){
        double tmptol = 0.0; //Running total for norm in denominator
        for (int n=0; n<dim_num; n++){
          tmpvec(n) = pow(rnd(o,n) - prevdes[n+m*dim_num],2);
          tmptol += tmpvec(n);
        }
        tmptol = sqrt(tmptol);
        tmpconst += 1.0/tmptol;
        
        for (int n=0; n<dim_num; n++){
          xprime(n) += rnd(o,n)/tmptol;
        }
      }
      
      // if (m == 0){
      //   Rcout << xprime(0) << endl;
      // }
      
      //Scale by inverse distances
      for (int n=0; n<dim_num; n++){
        xprime(n) = xprime(n)/tmpconst;
      }    
      
      // if (m == 0){
      //   Rcout << xprime(0) << endl;
      // }
      
      //Update bounds
      for (int n=0; n<dim_num; n++){
        xprime(n) = min( max( xprime(n), bd(n,0)), bd(n,1));
      }

      //Update point
      for (int n=0; n<dim_num; n++){
        des[n+m*dim_num] = xprime(n);
      }
    }

    // //Output time
    // time_t end = time(0);
    // double sec = difftime(end, start);
    // Rcout << "Iteration time: " << sec * 1000.0 << endl;
    
    //Increment and update cont flag 
    it_num++;
    //Check maximum iterations
    if (it_num >= it_max){
      cont = false;
    }
    //Check convergence (point with maximum movement)
    double maxdiff = 0.0;
    double rundiff = 0.0;
    for (int n=0; n<des_num; n++){
      rundiff = 0.0;
      for (int o=0; o<dim_num; o++){
        rundiff += std::pow(des[o+n*dim_num]-prevdes[o+n*dim_num],2.0);
      }
      maxdiff = max(maxdiff,rundiff);
    }
    // Rcout << "maxdiff: " << maxdiff << endl;
    
    if (maxdiff < tol){
      cont = false;
    }
    
    //      cont = false;//////
  }
  
  //Output the final NumericMatrix 
  // Rcout << "One value: " << des[0] << endl;
  NumericMatrix retdes(des_num, dim_num);
  for (int j=0; j<dim_num; j++){
    // if (distind(j)==5){//Rescale lognormal
    //   SEXP lnll = distparam[j]; 
    //   NumericVector lnyy(lnll);  
    //   for (int i=0; i<des_num; i++){
    //     retdes(i,j) = des[j+i*dim_num] * ( exp(pow(lnyy(1),2.0)-1.0) * exp(2*lnyy(0)+pow(lnyy(1),2.0)) );
    //   }
    // }
    // else{
      for (int i=0; i<des_num; i++){
        retdes(i,j) = des[j+i*dim_num];
      }
    // }
  }
  
  return(retdes);
  
}

// // [[Rcpp::export]]
// NumericMatrix unif_largep(NumericMatrix& Rcpp_inides,
//                           int it_max, int inn_it_max, double innertol, double outertol, double epsilon, double rate,
//                           int num_proc) {
//   int des_num = Rcpp_inides.nrow();
//   int dim_num = Rcpp_inides.ncol();
//   arma::mat inides(des_num, dim_num);
//   for (int i=0; i<des_num; i++){
//     for (int j=0; j<dim_num; j++){
//       inides(i,j) = Rcpp_inides(i,j);
//     }
//   }
//   
//   //Do SBD
//   unif_largep_cpp(inides, it_max, inn_it_max, innertol, outertol, epsilon, rate, num_proc);
//   
//   //Return matrix
//   NumericMatrix retdes(des_num, dim_num);
//   for (int i=0; i<des_num; i++){
//     for (int j=0; j<dim_num; j++){
//       retdes(i,j) = inides(i,j);
//     }
//   }
//   return(retdes);
// }



//-------------------------------------------------------------------------------
// Old code
//-------------------------------------------------------------------------------

//// [[Rcpp::export]]
//List SBD_PSO(NumericMatrix& Rcpp_point, NumericMatrix& Rcpp_inipart, NumericVector& bound, int it_max, double innertol, double outertol, int pw) {
//  //PSO implementation of SBD_foo
//  //Rcpp_point - Samples
//  //Rcpp_inides - Cube of initial particles
//  //it_max - maximum PSO iterations
//  //innertol - tolerance for blockwise updates
//  //pw - p in Lp
//  
//  int it_num = 0; //keeps track of iteration number
//  int dim_num = Rcpp_point.ncol(); //dimension of data points
//  int des_num = Rcpp_inipart.nrow(); //number of clusters
//  int part_num = Rcpp_inipart.ncol()/dim_num; // # of PSo particles
//  
////  Rcout << "Here 0" << endl;
//  
//  //PSO objects
//  arma::mat des = zeros(des_num,dim_num); //Gbes design to return
//  arma::mat prevdes = zeros(des_num,dim_num); //Gbes design to return
//  arma::cube des_part = zeros(des_num,dim_num,part_num); // Design particles
//  arma::cube des_vel = zeros(des_num,dim_num,part_num);//current velocity
//  arma::cube des_lbes = zeros(des_num,dim_num,part_num);//lbes positions
//  arma::vec cur_obj(part_num); //current energy
//  arma::vec lbes_obj = DBL_MAX*ones(part_num); //local best energy
//  double gbes_obj = DBL_MAX; //global minimum energy over all particles
//  
////  Rcout << "Here 1" << endl;
//  
//  //Copy into cube object for slices
//  for (int i=0; i<des_num; i++){
//    for (int j=0; j<dim_num; j++){
//      for (int k=0; k<part_num; k++){
//        des_part(i,j,k) = Rcpp_inipart(i,j+(k)*dim_num);
//      }
//    }
//  }
//  
////  for (int j=0; j<dim_num; j++){
////    Rcout << "start_pt: " << Rcpp_inipart(0,j+(k-1)*dim_num) << endl;
////  }
//
//  
////  Rcout << "Here 2" << endl;
//  
//  //PSO constants
//  double w = 0.72;//inertia constant (from Merwe and Engelbrecht, 2003)
//  double c1 = 1.49;//acceleration constants
//  double c2 = 1.49;
//  bool cont = true; //flag for continuing
//  
//  while (cont){
//    
//    Rcout << "PSO iteration: " << it_num << "..." << endl;
//    
//    //Update prevdes
//    for (int i=0; i<des_num; i++){
//      for (int j=0; j<dim_num; j++){
//        prevdes(i,j) = des(i,j);
//      }
//    }
//    
////    Rcout << "Here 0" << endl;
//    
//    //Do one step of PSO
//    for (int i=0; i<part_num; i++){
//      if (pw==1){
//        SBD_L1(Rcpp_point, des_part.slice(i), 1, innertol); //1e-2 doesn't matter
//      }
//      else{
////        Rcout << "Here 0.0" << endl;
//        SBD_Lp(Rcpp_point, des_part.slice(i), bound, 1, innertol, outertol, pw); //1e-2 doesn't matter
////        Rcout << "Here 0.1" << endl;
//      }
//      
//      cur_obj(i) = energycrit(Rcpp_point,des_part.slice(i),pw);
////      Rcout << "Here 0.2" << endl;
//    }
//    
////    Rcout << "Here 1" << endl;     
//    
//    //Update lbes and gbes
//    for (int i=0;i<part_num;i++){
//      if (cur_obj(i)<lbes_obj(i)){
//        des_lbes.slice(i) = des_part.slice(i);
//        lbes_obj(i) = cur_obj(i);
//      }
//      if (cur_obj(i)<gbes_obj){
//        des = des_part.slice(i);
//        gbes_obj = cur_obj(i);
//        Rcout <<  "gbes_obj changed to:" << endl;
//        Rcout << gbes_obj << endl;
//      }
//    }
//    
////    Rcout << "Here 2" << endl;
//    
//    //Update velocities
//    for (int i=0; i<part_num; i++){
//      //update velocity  
//      des_vel.slice(i) = w*des_vel.slice(i) + c1*(randu(des_num,dim_num)%(des_lbes.slice(i)-des_part.slice(i))) + c2*(randu(des_num,dim_num)%(des-des_part.slice(i)));
//      //update particle
//      des_part.slice(i) += des_vel.slice(i);
//  //    enfBox(des_part.slice(i),lb,ub);
//    }
//    
////    Rcout << "Here 6" << endl;
//    
//    //Check convergence and update flag
//    it_num++;
//    if (it_num>=it_max){
//      cont = false;
//    }
//    //Check convergence
//    double rundiff = 0.0;
//    for (int n=0; n<des_num; n++){
//      for (int o=0; o<dim_num; o++){
//        rundiff += abs(des[o+n*dim_num]-prevdes[o+n*dim_num]);
//      }
//    }
//    Rcout << "PSO rundiff: " << rundiff << endl;
////    if (rundiff <= outertol){
////      cont = false;
////    }
//    
//  }
//  
////  Rcout << "Here 7" << endl;
//  
//  //Return globally optimal design
//  NumericMatrix retdes(des_num,dim_num);
//  for (int i=0; i<des_num; i++){
//    for (int j=0; j<dim_num; j++){
//      retdes(i,j) = des(i,j);
//    }
//  }
//  
////  Rcout << "Here 8" << endl;
//  
//  return List::create(Named("cur_energy")=cur_obj,Named("lbes_energy")=lbes_obj,Named("des")=retdes);
//
//}

//// [[Rcpp::export]]
//NumericMatrix SBD_L1(NumericMatrix& Rcpp_point, NumericMatrix& Rcpp_inides, int it_max, double tol) {
//  int des_num = Rcpp_inides.nrow();
//  int dim_num = Rcpp_inides.ncol();
//  arma::mat inides(des_num, dim_num);
//  for (int i=0; i<des_num; i++){
//    for (int j=0; j<dim_num; j++){
//      inides(i,j) = Rcpp_inides(i,j);
//    }
//  }
//  
//  //Do SBD
//  SBD_L1(Rcpp_point, inides, it_max, tol);
//  
//  //Return matrix
//  NumericMatrix retdes(des_num, dim_num);
//  for (int i=0; i<des_num; i++){
//    for (int j=0; j<dim_num; j++){
//      retdes(i,j) = inides(i,j);
//    }
//  }
//  return(retdes);
//}

//void SBD_Lp(NumericMatrix& Rcpp_point, arma::mat& Rcpp_inides, NumericVector& bound, int it_max, double innertol, double outertol, int pw) {
//  //Routing function for SBD_L1, SBD_L2 and SBD_Lp
//  
//  // Description of Inputs:
//  // Rcpp_point          - Sample data
//  // Rcpp_inides         - Initial design
//  // bound               - Upper and lower bound
//  // it_max              - Maximum number of iterations
//  // tol                 - Tolerance for stopping blockwise updates
//  // pw                  - p in Lp  
//  
//  if (pw == 1){
//    SBD_L1(Rcpp_point, Rcpp_inides, it_max, outertol);
//  }
//  else if (pw == 2){
//    SBD_L2(Rcpp_point, Rcpp_inides, bound, it_max, innertol, outertol);
//  }
//  else{
//    //Non-linear opt for p>2
//    
//    int it_num = 0; //keeps track of iteration number
//    int dim_num = Rcpp_point.ncol(); //dimension of data points
//    int point_num = Rcpp_point.nrow(); //number of data points
//    int des_num = Rcpp_inides.n_rows; //number of clusters
//    bool cont = true; //stopping boolean for outside BCD loop
//    bool cont2 = true; //stopping boolean for inside CCCP loop
////    double innertol = 1e-4;
//  
//    //Containers for optimization
//    std::vector<double> prevdes(des_num*dim_num);
////    std::vector<double> thindes(thin_num*dim_num);
//    column_vector start_pt(dim_num);
//    
//    //Vectorize design and points  
//    std::vector<double> des(des_num*dim_num);
//    for (int i=0; i<des_num; i++){
//      for (int j=0; j<dim_num; j++){
//        des[j+i*dim_num] = Rcpp_inides(i,j);
//      }
//    }
//    std::vector<double> point(point_num*dim_num);
//    for (int i=0; i<point_num; i++){
//      for (int j=0; j<dim_num; j++){
//        point[j+i*dim_num] = Rcpp_point(i,j);
//      }
//    }
//    
//    //Blockwise coordinate descent
//    while (cont){
//      
////      Rcout << "Iteration " << it_num << "..." << endl;
//      
//      //Update prevdes
//      for (int i=0; i<des_num; i++){
//        for (int j=0; j<dim_num; j++){
//          prevdes[j+i*dim_num] = des[j+i*dim_num];
//        }
//      }
//      
//      //BCD for each point
//      for (int m=0; m<des_num; m++){
//        
////        Rcout << "m = " << m << endl;
//      
//        //Copy onto start_pt as initial point
//        for (int n=0; n<dim_num; n++){
//          start_pt(n) = des[n+m*dim_num];
//        }
//        
//        
////        //Reset cont2 flag
////        cont2 = true;
//        
//        //Concave-convex updates until convergence
////        while(cont2){
//          //Optimization
//          objective obj(point, des, point_num, des_num, dim_num, m, pw, bound);
//          objective_grad obj_grad(point, des, point_num, des_num, dim_num, m, pw);
//          
////          clock_t start = clock();
////          Rcout << "Evaluation: " << obj(start_pt) << endl;
////          Rcout << double( clock() - start ) / (double)CLOCKS_PER_SEC<< " seconds." << endl;
////          
////          clock_t start1 = clock();
////          for (int pp=0; pp<dim_num; pp++){
////            Rcout << "start_pt: " << start_pt(pp) << endl;
////          }
////          Rcout << "Analytic derivative: " << obj_grad(start_pt) << endl;
////          Rcout << double( clock() - start1 ) / (double)CLOCKS_PER_SEC<< " seconds." << endl;
////          Rcout << "Numeric derivative: " << dlib::derivative(obj)(start_pt) << endl;
//          
//          find_min(dlib::lbfgs_search_strategy(3), dlib::objective_delta_stop_strategy(innertol), obj, obj_grad, start_pt, -1);
//
////          find_min_using_approximate_derivatives(dlib::lbfgs_search_strategy(3),
////                                               dlib::objective_delta_stop_strategy(tol),
////                                               obj, start_pt, -1); 
//                                               
////          //Update convergence flag
////          double runtol = 0.0;
////          for (int n=0; n<dim_num; n++){
////            runtol += abs(des[m+n*des_num] - start_pt(n));         
////          }
////          if (runtol < innertol){
////            cont2 = false;
////          }
////        }
//
//        //Update design point 
//        for (int n=0; n<dim_num; n++){
//          des[n+m*dim_num] = start_pt(n);
//        }
//        
//      }
// 
//      //Increment and update cont flag 
//      it_num++;
//      //Check maximum iterations
//      if (it_num >= it_max){
//        cont = false;
//      }
//      //Check convergence
//      double rundiff = 0.0;
//      for (int n=0; n<des_num; n++){
//        for (int o=0; o<dim_num; o++){
//          rundiff += abs(des[o+n*dim_num]-prevdes[o+n*dim_num]);
//        }
//      }
//      if (rundiff < outertol){
//        cont = false;
//      }
//    }
//    
//    //Output the final NumericMatrix 
//    for (int i=0; i<des_num; i++){
//      for (int j=0; j<dim_num; j++){
//        Rcpp_inides(i,j) = des[j+i*dim_num];
//      }
//    }
//    
//  }
//  
//}

//void SBD_L1(NumericMatrix& Rcpp_point, arma::mat& Rcpp_inides, int it_max, double tol) {
//  // Description of Inputs:
//  // Rcpp_point          - Sample data
//  // Rcpp_inides         - Initial design
//  // it_max              - Maximum number of iterations
//  // tol                 - Tolerance for stopping blockwise updates
//  // sorted              - true if each dimension in Rcpp_point is sorted
//  
//  int it_num = 0; //keeps track of iteration number
//  int dim_num = Rcpp_point.ncol(); //dimension of data points
//  int point_num = Rcpp_point.nrow(); //number of data points
//  int des_num = Rcpp_inides.n_rows; //number of clusters
//  bool cont = true; //stopping boolean for outside BCD loop
//  bool cont2 = true; //stopping boolean for inside CCCP loop
//  double innertol = 1e-4;
//  clock_t start, cdiff;
//  int msec;
//
//  //Containers to assign
//  vector<double> desvec(des_num);
//  vector<double> desmi(des_num-1); //sorted design without index i
//  vector<double> desmimean(des_num-1); //initial points to determine best starting pt
//  vector<double> desmival(des_num-1);
//  double prevsol;
//  double cursol;
//  double slp = 0.0; //slope for linearization
//  int num_sm = 0; //container for counts
//  int max_idx = 0; //maximum index for subgradients
//  std::vector<double> prevdes(des_num*dim_num);
//  
//  //Vectorize design and points  
//  std::vector<double> des(des_num*dim_num);
//  for (int i=0; i<des_num; i++){
//    for (int j=0; j<dim_num; j++){
//      des[i+j*des_num] = Rcpp_inides(i,j);
//    }
//  }
//  std::vector<double> point(point_num*dim_num);
//  for (int i=0; i<point_num; i++){
//    for (int j=0; j<dim_num; j++){
//      point[i+j*point_num] = Rcpp_point(i,j);
//    }
//  }
//  
////  //Sort each dimension of points if not sorted
////  if (!sorted){
////    sort(point.begin(), point.end());
////  }
//  
//  //Blockwise coordinate descent
//  while (cont){
//    
//    Rcout << "Iteration " << it_num << "..." << endl;
//    
//    //Update prevdes
//    for (int i=0; i<des_num; i++){
//      for (int j=0; j<dim_num; j++){
//        prevdes[i+j*des_num] = des[i+j*des_num];
//      }
//    }
//    
//    //BCD for each dimension and point
//    for (int l=0; l<dim_num; l++){
//    
//      //Copy current design and compute distance from points
//      copy(des.begin()+l*des_num, des.begin()+(l+1)*des_num, desvec.begin());
//      sort(desvec.begin(), desvec.end());
//      
//      double runsum = 0.0;
//      double runsum2 = 0.0;
//      
//      for (int n=0; n<(des_num-1); n++){
//        
//        runsum = 0.0;
//        runsum2 = 0.0;
//        desmimean[n] = 0.5*(desvec[n] + desvec[n+1]);
//        
//        for (int o=0; o<point_num; o++){
//          runsum += abs(point[o+l*point_num]-desmimean[n]);
//        }
//        runsum = (1.0/point_num)*runsum;
//        for (int o=0; o<des_num; o++){
//          runsum2 -= abs(desvec[o]-desmimean[n]);
//        }
//        runsum2 = (1.0/des_num)*runsum2;          
//        desmival[n] = runsum + runsum2;
//        
//      }
//      
//      for (int m=0; m<des_num; m++){
//        
//        //Add extra part
//        for (int n=0; n<(des_num-1); n++){
//          desmival[n] += (1.0/des_num)*abs(des[m+l*des_num]-desmimean[n]);
//        }
//        
//        //Reset CCCP flag
//        cont2 = true;
//        
//        //Find the minimum point as starting point
//        int min_idx = -1;
//        double cur_min = DBL_MAX;
//        for (int n=0; n<(des_num-1); n++){
//          if (n != m){
////            Rcout << "n: " << n << ", " << "desmival: " << desmival[n] << ", cur_min " << cur_min << endl;
//            if (desmival[n] <= cur_min){
//              min_idx = n;
//              cur_min = desmival[n];
//            }
//          }
//        }
////        vector<double>::iterator min_idx = min_element(desmival.begin(),desmival.end());
//        prevsol = desmimean[ min_idx ];
//        cursol = prevsol;
//        
//        //Majorize-minimize in CCCP
//        while(cont2){
//          //Compute slope for linearization
//          prevsol = cursol;
//          num_sm = 0;
//          for (int n=0; n<des_num; n++){
//            if (n != m){
//              if (des[n+l*des_num]<cursol){
//                num_sm++;
//              }
//            }
//          }
//          slp = -1.0/(des_num)*( (des_num-num_sm-1) - num_sm );
//          
//          //Solve KKT subgradient conditions
//          max_idx = (slp*point_num + point_num)/2-1;
////          cursol = 0.5*(point[max_idx+l*point_num] + point[max_idx+l*point_num+1]);
//          cursol = point[max_idx+l*point_num];
//          
//          //Check convergence for cont2 flag
//          if ( abs(cursol - prevsol) < innertol){
//            cont2 = false;
//          }
//          
//        }
//        
//        //Update point
////        Rcout << "desmival: ";
////        for (int o=0; o<(des_num-1); o++){
////          Rcout << desmival[o] << " ";
////        }
////        Rcout << endl;
////        Rcout << "desmival optimal: " << desmival[ min_idx ] << endl;
////        Rcout << "num_sm: " << num_sm << endl;
////        Rcout << "slp: " << slp << endl;
////        Rcout << "maxidx: " << max_idx << endl;
////        Rcout << "prevsol: " << prevsol << endl;
////        Rcout << "cursol: " << cursol << endl;
////        Rcout << "------------------------------" << endl;
//        des[m+l*des_num] = cursol;
//        
//        //Subtract extra part
//        for (int n=0; n<(des_num-1); n++){
//          desmival[n] -= (1.0/des_num)*abs(cursol-desmimean[n]);
//        }
//        
//      }
//
//    }
//    
//    //Increment and update cont flag 
//    it_num++;
//    //Check maximum iterations
//    if (it_num >= it_max){
//      cont = false;
//    }
//    //Check convergence
//    double rundiff = 0.0;
//    for (int n=0; n<des_num; n++){
//      for (int o=0; o<dim_num; o++){
//        rundiff += abs(des[n+o*des_num]-prevdes[n+o*des_num]);
//      }
//    }
//    if (rundiff < tol){
//      cont = false;
//    }
//  }
//  
//  //Output the final NumericMatrix 
////  NumericMatrix retdes(des_num,dim_num);
//  for (int i=0; i<des_num; i++){
//    for (int j=0; j<dim_num; j++){
//      Rcpp_inides(i,j) = des[i+j*des_num];
//    }
//  }
//  
////  return(retdes);
//}
