// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// herdobj_seq
double herdobj_seq(NumericVector& xx, NumericMatrix& des, NumericMatrix& distsamp, double sigma);
RcppExport SEXP _support_herdobj_seq(SEXP xxSEXP, SEXP desSEXP, SEXP distsampSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type des(desSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(herdobj_seq(xx, des, distsamp, sigma));
    return rcpp_result_gen;
END_RCPP
}
// herdgrad_seq
NumericVector herdgrad_seq(NumericVector& xx, NumericMatrix& des, NumericMatrix& distsamp, double sigma);
RcppExport SEXP _support_herdgrad_seq(SEXP xxSEXP, SEXP desSEXP, SEXP distsampSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type des(desSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(herdgrad_seq(xx, des, distsamp, sigma));
    return rcpp_result_gen;
END_RCPP
}
// herdobj_full
double herdobj_full(NumericVector& xx, NumericMatrix& des, int idx, NumericMatrix& distsamp, double sigma);
RcppExport SEXP _support_herdobj_full(SEXP xxSEXP, SEXP desSEXP, SEXP idxSEXP, SEXP distsampSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type des(desSEXP);
    Rcpp::traits::input_parameter< int >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(herdobj_full(xx, des, idx, distsamp, sigma));
    return rcpp_result_gen;
END_RCPP
}
// herdgrad_full
NumericVector herdgrad_full(NumericVector& xx, NumericMatrix& des, int idx, NumericMatrix& distsamp, double sigma);
RcppExport SEXP _support_herdgrad_full(SEXP xxSEXP, SEXP desSEXP, SEXP idxSEXP, SEXP distsampSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type des(desSEXP);
    Rcpp::traits::input_parameter< int >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(herdgrad_full(xx, des, idx, distsamp, sigma));
    return rcpp_result_gen;
END_RCPP
}
// pspobj_seq
double pspobj_seq(NumericVector& xx, NumericMatrix& des, NumericMatrix& distsamp, double lambda, double nu);
RcppExport SEXP _support_pspobj_seq(SEXP xxSEXP, SEXP desSEXP, SEXP distsampSEXP, SEXP lambdaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type des(desSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(pspobj_seq(xx, des, distsamp, lambda, nu));
    return rcpp_result_gen;
END_RCPP
}
// pspgrad_seq
NumericVector pspgrad_seq(NumericVector& xx, NumericMatrix& des, NumericMatrix& distsamp, double lambda, double nu);
RcppExport SEXP _support_pspgrad_seq(SEXP xxSEXP, SEXP desSEXP, SEXP distsampSEXP, SEXP lambdaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type des(desSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(pspgrad_seq(xx, des, distsamp, lambda, nu));
    return rcpp_result_gen;
END_RCPP
}
// pspobj_full
double pspobj_full(NumericVector& xx, NumericMatrix& des, int idx, NumericMatrix& distsamp, double lambda, double nu);
RcppExport SEXP _support_pspobj_full(SEXP xxSEXP, SEXP desSEXP, SEXP idxSEXP, SEXP distsampSEXP, SEXP lambdaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type des(desSEXP);
    Rcpp::traits::input_parameter< int >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(pspobj_full(xx, des, idx, distsamp, lambda, nu));
    return rcpp_result_gen;
END_RCPP
}
// pspgrad_full
NumericVector pspgrad_full(NumericVector& xx, NumericMatrix& des, int idx, NumericMatrix& distsamp, double lambda, double nu);
RcppExport SEXP _support_pspgrad_full(SEXP xxSEXP, SEXP desSEXP, SEXP idxSEXP, SEXP distsampSEXP, SEXP lambdaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type des(desSEXP);
    Rcpp::traits::input_parameter< int >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(pspgrad_full(xx, des, idx, distsamp, lambda, nu));
    return rcpp_result_gen;
END_RCPP
}
// pspobj_seq2
double pspobj_seq2(NumericVector& xx, NumericMatrix& des, NumericMatrix& distsamp, double lambda, double nu);
RcppExport SEXP _support_pspobj_seq2(SEXP xxSEXP, SEXP desSEXP, SEXP distsampSEXP, SEXP lambdaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type des(desSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(pspobj_seq2(xx, des, distsamp, lambda, nu));
    return rcpp_result_gen;
END_RCPP
}
// pspgrad_seq2
NumericVector pspgrad_seq2(NumericVector& xx, NumericMatrix& des, NumericMatrix& distsamp, double lambda, double nu);
RcppExport SEXP _support_pspgrad_seq2(SEXP xxSEXP, SEXP desSEXP, SEXP distsampSEXP, SEXP lambdaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type des(desSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(pspgrad_seq2(xx, des, distsamp, lambda, nu));
    return rcpp_result_gen;
END_RCPP
}
// pspobj_full2
double pspobj_full2(NumericVector& xx, NumericMatrix& des, int idx, NumericMatrix& distsamp, double lambda, double nu);
RcppExport SEXP _support_pspobj_full2(SEXP xxSEXP, SEXP desSEXP, SEXP idxSEXP, SEXP distsampSEXP, SEXP lambdaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type des(desSEXP);
    Rcpp::traits::input_parameter< int >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(pspobj_full2(xx, des, idx, distsamp, lambda, nu));
    return rcpp_result_gen;
END_RCPP
}
// pspgrad_full2
NumericVector pspgrad_full2(NumericVector& xx, NumericMatrix& des, int idx, NumericMatrix& distsamp, double lambda, double nu);
RcppExport SEXP _support_pspgrad_full2(SEXP xxSEXP, SEXP desSEXP, SEXP idxSEXP, SEXP distsampSEXP, SEXP lambdaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type des(desSEXP);
    Rcpp::traits::input_parameter< int >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(pspgrad_full2(xx, des, idx, distsamp, lambda, nu));
    return rcpp_result_gen;
END_RCPP
}
// printBar
void printBar(double prop);
RcppExport SEXP _support_printBar(SEXP propSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type prop(propSEXP);
    printBar(prop);
    return R_NilValue;
END_RCPP
}
// csample_num
NumericVector csample_num(NumericVector x, int size, bool replace);
RcppExport SEXP _support_csample_num(SEXP xSEXP, SEXP sizeSEXP, SEXP replaceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type replace(replaceSEXP);
    rcpp_result_gen = Rcpp::wrap(csample_num(x, size, replace));
    return rcpp_result_gen;
END_RCPP
}
// energycrit
double energycrit(NumericMatrix& Rcpp_point, NumericMatrix& Rcpp_des);
RcppExport SEXP _support_energycrit(SEXP Rcpp_pointSEXP, SEXP Rcpp_desSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type Rcpp_point(Rcpp_pointSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type Rcpp_des(Rcpp_desSEXP);
    rcpp_result_gen = Rcpp::wrap(energycrit(Rcpp_point, Rcpp_des));
    return rcpp_result_gen;
END_RCPP
}
// starL2cpp
NumericVector starL2cpp(NumericMatrix& D, NumericMatrix& cn, int num_proc);
RcppExport SEXP _support_starL2cpp(SEXP DSEXP, SEXP cnSEXP, SEXP num_procSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type D(DSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type cn(cnSEXP);
    Rcpp::traits::input_parameter< int >::type num_proc(num_procSEXP);
    rcpp_result_gen = Rcpp::wrap(starL2cpp(D, cn, num_proc));
    return rcpp_result_gen;
END_RCPP
}
// energy_norm_cpp
NumericVector energy_norm_cpp(NumericMatrix& yMat, NumericMatrix& cn, int num_proc);
RcppExport SEXP _support_energy_norm_cpp(SEXP yMatSEXP, SEXP cnSEXP, SEXP num_procSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type yMat(yMatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type cn(cnSEXP);
    Rcpp::traits::input_parameter< int >::type num_proc(num_procSEXP);
    rcpp_result_gen = Rcpp::wrap(energy_norm_cpp(yMat, cn, num_proc));
    return rcpp_result_gen;
END_RCPP
}
// obj_qsp
double obj_qsp(arma::vec& des, arma::mat& distsamp, double q);
RcppExport SEXP _support_obj_qsp(SEXP desSEXP, SEXP distsampSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type des(desSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(obj_qsp(des, distsamp, q));
    return rcpp_result_gen;
END_RCPP
}
// grad_qsp
arma::vec grad_qsp(arma::vec& des, arma::mat& distsamp, double q);
RcppExport SEXP _support_grad_qsp(SEXP desSEXP, SEXP distsampSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type des(desSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_qsp(des, distsamp, q));
    return rcpp_result_gen;
END_RCPP
}
// sp_cpp
NumericMatrix sp_cpp(int des_num, int dim_num, NumericMatrix& ini, NumericVector& distind, List distparam, NumericMatrix& distsamp, bool thin, NumericMatrix& bd, int point_num, int it_max, int it_min, double tol, int num_proc, double n0, NumericVector& wts, bool rnd_flg);
RcppExport SEXP _support_sp_cpp(SEXP des_numSEXP, SEXP dim_numSEXP, SEXP iniSEXP, SEXP distindSEXP, SEXP distparamSEXP, SEXP distsampSEXP, SEXP thinSEXP, SEXP bdSEXP, SEXP point_numSEXP, SEXP it_maxSEXP, SEXP it_minSEXP, SEXP tolSEXP, SEXP num_procSEXP, SEXP n0SEXP, SEXP wtsSEXP, SEXP rnd_flgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type des_num(des_numSEXP);
    Rcpp::traits::input_parameter< int >::type dim_num(dim_numSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type ini(iniSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type distind(distindSEXP);
    Rcpp::traits::input_parameter< List >::type distparam(distparamSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< bool >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type bd(bdSEXP);
    Rcpp::traits::input_parameter< int >::type point_num(point_numSEXP);
    Rcpp::traits::input_parameter< int >::type it_max(it_maxSEXP);
    Rcpp::traits::input_parameter< int >::type it_min(it_minSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type num_proc(num_procSEXP);
    Rcpp::traits::input_parameter< double >::type n0(n0SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type wts(wtsSEXP);
    Rcpp::traits::input_parameter< bool >::type rnd_flg(rnd_flgSEXP);
    rcpp_result_gen = Rcpp::wrap(sp_cpp(des_num, dim_num, ini, distind, distparam, distsamp, thin, bd, point_num, it_max, it_min, tol, num_proc, n0, wts, rnd_flg));
    return rcpp_result_gen;
END_RCPP
}
// sp_seq_cpp
NumericMatrix sp_seq_cpp(NumericMatrix& cur, int nseq, NumericMatrix& ini, NumericVector& distind, List distparam, NumericMatrix& distsamp, bool thin, NumericMatrix& bd, int point_num, int it_max, int it_min, double tol, int num_proc);
RcppExport SEXP _support_sp_seq_cpp(SEXP curSEXP, SEXP nseqSEXP, SEXP iniSEXP, SEXP distindSEXP, SEXP distparamSEXP, SEXP distsampSEXP, SEXP thinSEXP, SEXP bdSEXP, SEXP point_numSEXP, SEXP it_maxSEXP, SEXP it_minSEXP, SEXP tolSEXP, SEXP num_procSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type cur(curSEXP);
    Rcpp::traits::input_parameter< int >::type nseq(nseqSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type ini(iniSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type distind(distindSEXP);
    Rcpp::traits::input_parameter< List >::type distparam(distparamSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< bool >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type bd(bdSEXP);
    Rcpp::traits::input_parameter< int >::type point_num(point_numSEXP);
    Rcpp::traits::input_parameter< int >::type it_max(it_maxSEXP);
    Rcpp::traits::input_parameter< int >::type it_min(it_minSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type num_proc(num_procSEXP);
    rcpp_result_gen = Rcpp::wrap(sp_seq_cpp(cur, nseq, ini, distind, distparam, distsamp, thin, bd, point_num, it_max, it_min, tol, num_proc));
    return rcpp_result_gen;
END_RCPP
}
// gamma_eval
double gamma_eval(arma::vec& dd, arma::vec& theta);
RcppExport SEXP _support_gamma_eval(SEXP ddSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type dd(ddSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(gamma_eval(dd, theta));
    return rcpp_result_gen;
END_RCPP
}
// omega
arma::vec omega(arma::vec& theta_vec, arma::vec& gamma_vec, int max_ord);
RcppExport SEXP _support_omega(SEXP theta_vecSEXP, SEXP gamma_vecSEXP, SEXP max_ordSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type theta_vec(theta_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type gamma_vec(gamma_vecSEXP);
    Rcpp::traits::input_parameter< int >::type max_ord(max_ordSEXP);
    rcpp_result_gen = Rcpp::wrap(omega(theta_vec, gamma_vec, max_ord));
    return rcpp_result_gen;
END_RCPP
}
// opt_hess
arma::mat opt_hess(arma::vec zz, arma::vec omega_vec);
RcppExport SEXP _support_opt_hess(SEXP zzSEXP, SEXP omega_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type zz(zzSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega_vec(omega_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(opt_hess(zz, omega_vec));
    return rcpp_result_gen;
END_RCPP
}
// psp_mi
arma::vec psp_mi(arma::vec& xx, arma::mat& omega_mat, arma::mat& samp_mat, std::vector<double>& des_mat, int des_num, int ii);
RcppExport SEXP _support_psp_mi(SEXP xxSEXP, SEXP omega_matSEXP, SEXP samp_matSEXP, SEXP des_matSEXP, SEXP des_numSEXP, SEXP iiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type omega_mat(omega_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type samp_mat(samp_matSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type des_mat(des_matSEXP);
    Rcpp::traits::input_parameter< int >::type des_num(des_numSEXP);
    Rcpp::traits::input_parameter< int >::type ii(iiSEXP);
    rcpp_result_gen = Rcpp::wrap(psp_mi(xx, omega_mat, samp_mat, des_mat, des_num, ii));
    return rcpp_result_gen;
END_RCPP
}
// psp_cpp
NumericMatrix psp_cpp(NumericMatrix& Rcpp_inides, NumericVector& distind, List distparam, NumericVector& gam_param, arma::vec& gamma_vec, int max_ord, NumericMatrix& distsamp, bool thinind, NumericMatrix& gamsamp, bool gamind, NumericMatrix& bd, int point_num, int gam_point_num, int it_max, double tol, int num_proc);
RcppExport SEXP _support_psp_cpp(SEXP Rcpp_inidesSEXP, SEXP distindSEXP, SEXP distparamSEXP, SEXP gam_paramSEXP, SEXP gamma_vecSEXP, SEXP max_ordSEXP, SEXP distsampSEXP, SEXP thinindSEXP, SEXP gamsampSEXP, SEXP gamindSEXP, SEXP bdSEXP, SEXP point_numSEXP, SEXP gam_point_numSEXP, SEXP it_maxSEXP, SEXP tolSEXP, SEXP num_procSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type Rcpp_inides(Rcpp_inidesSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type distind(distindSEXP);
    Rcpp::traits::input_parameter< List >::type distparam(distparamSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type gam_param(gam_paramSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type gamma_vec(gamma_vecSEXP);
    Rcpp::traits::input_parameter< int >::type max_ord(max_ordSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< bool >::type thinind(thinindSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type gamsamp(gamsampSEXP);
    Rcpp::traits::input_parameter< bool >::type gamind(gamindSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type bd(bdSEXP);
    Rcpp::traits::input_parameter< int >::type point_num(point_numSEXP);
    Rcpp::traits::input_parameter< int >::type gam_point_num(gam_point_numSEXP);
    Rcpp::traits::input_parameter< int >::type it_max(it_maxSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type num_proc(num_procSEXP);
    rcpp_result_gen = Rcpp::wrap(psp_cpp(Rcpp_inides, distind, distparam, gam_param, gamma_vec, max_ord, distsamp, thinind, gamsamp, gamind, bd, point_num, gam_point_num, it_max, tol, num_proc));
    return rcpp_result_gen;
END_RCPP
}
// psp_seq_cpp
NumericMatrix psp_seq_cpp(NumericMatrix& cur, int nseq, NumericMatrix& ini, NumericVector& distind, List distparam, NumericVector& gam_param, arma::vec& gamma_vec, int max_ord, NumericMatrix& distsamp, bool thin, NumericMatrix& gamsamp, bool gamind, NumericMatrix& bd, int point_num, int gam_point_num, int it_max, double tol, int num_proc);
RcppExport SEXP _support_psp_seq_cpp(SEXP curSEXP, SEXP nseqSEXP, SEXP iniSEXP, SEXP distindSEXP, SEXP distparamSEXP, SEXP gam_paramSEXP, SEXP gamma_vecSEXP, SEXP max_ordSEXP, SEXP distsampSEXP, SEXP thinSEXP, SEXP gamsampSEXP, SEXP gamindSEXP, SEXP bdSEXP, SEXP point_numSEXP, SEXP gam_point_numSEXP, SEXP it_maxSEXP, SEXP tolSEXP, SEXP num_procSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type cur(curSEXP);
    Rcpp::traits::input_parameter< int >::type nseq(nseqSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type ini(iniSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type distind(distindSEXP);
    Rcpp::traits::input_parameter< List >::type distparam(distparamSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type gam_param(gam_paramSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type gamma_vec(gamma_vecSEXP);
    Rcpp::traits::input_parameter< int >::type max_ord(max_ordSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type distsamp(distsampSEXP);
    Rcpp::traits::input_parameter< bool >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type gamsamp(gamsampSEXP);
    Rcpp::traits::input_parameter< bool >::type gamind(gamindSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type bd(bdSEXP);
    Rcpp::traits::input_parameter< int >::type point_num(point_numSEXP);
    Rcpp::traits::input_parameter< int >::type gam_point_num(gam_point_numSEXP);
    Rcpp::traits::input_parameter< int >::type it_max(it_maxSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type num_proc(num_procSEXP);
    rcpp_result_gen = Rcpp::wrap(psp_seq_cpp(cur, nseq, ini, distind, distparam, gam_param, gamma_vec, max_ord, distsamp, thin, gamsamp, gamind, bd, point_num, gam_point_num, it_max, tol, num_proc));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_support_herdobj_seq", (DL_FUNC) &_support_herdobj_seq, 4},
    {"_support_herdgrad_seq", (DL_FUNC) &_support_herdgrad_seq, 4},
    {"_support_herdobj_full", (DL_FUNC) &_support_herdobj_full, 5},
    {"_support_herdgrad_full", (DL_FUNC) &_support_herdgrad_full, 5},
    {"_support_pspobj_seq", (DL_FUNC) &_support_pspobj_seq, 5},
    {"_support_pspgrad_seq", (DL_FUNC) &_support_pspgrad_seq, 5},
    {"_support_pspobj_full", (DL_FUNC) &_support_pspobj_full, 6},
    {"_support_pspgrad_full", (DL_FUNC) &_support_pspgrad_full, 6},
    {"_support_pspobj_seq2", (DL_FUNC) &_support_pspobj_seq2, 5},
    {"_support_pspgrad_seq2", (DL_FUNC) &_support_pspgrad_seq2, 5},
    {"_support_pspobj_full2", (DL_FUNC) &_support_pspobj_full2, 6},
    {"_support_pspgrad_full2", (DL_FUNC) &_support_pspgrad_full2, 6},
    {"_support_printBar", (DL_FUNC) &_support_printBar, 1},
    {"_support_csample_num", (DL_FUNC) &_support_csample_num, 3},
    {"_support_energycrit", (DL_FUNC) &_support_energycrit, 2},
    {"_support_starL2cpp", (DL_FUNC) &_support_starL2cpp, 3},
    {"_support_energy_norm_cpp", (DL_FUNC) &_support_energy_norm_cpp, 3},
    {"_support_obj_qsp", (DL_FUNC) &_support_obj_qsp, 3},
    {"_support_grad_qsp", (DL_FUNC) &_support_grad_qsp, 3},
    {"_support_sp_cpp", (DL_FUNC) &_support_sp_cpp, 16},
    {"_support_sp_seq_cpp", (DL_FUNC) &_support_sp_seq_cpp, 13},
    {"_support_gamma_eval", (DL_FUNC) &_support_gamma_eval, 2},
    {"_support_omega", (DL_FUNC) &_support_omega, 3},
    {"_support_opt_hess", (DL_FUNC) &_support_opt_hess, 2},
    {"_support_psp_mi", (DL_FUNC) &_support_psp_mi, 6},
    {"_support_psp_cpp", (DL_FUNC) &_support_psp_cpp, 16},
    {"_support_psp_seq_cpp", (DL_FUNC) &_support_psp_seq_cpp, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_support(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
