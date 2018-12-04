#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector log_ps_cpp (NumericVector ks, double alpha, double delta) {
  NumericVector vec = alpha * log(delta) + lgamma(ks+alpha) - lgamma(alpha)
    - lgamma(ks+1) - (ks+alpha)*log(delta+1);
  return vec;
}

// [[Rcpp::export]]
NumericVector marginal_negative_binomial_cpp (NumericVector k,
      double alpha, double delta) {
  return dnbinom(k, alpha, delta/(delta + 1));
}

// [[Rcpp::export]]
NumericVector log_marginal_negative_binomial_cpp (NumericVector k,
                                                             double alpha,
                                                             double delta) {
  double p = delta/(1+delta);
  NumericVector vec = lgamma(k+alpha) - lgamma(alpha) - lgamma(k+1) + alpha*log(p) +
    k*log(1-p);
  return vec;
}

// [[Rcpp::export]]
double loglike_cond_cpp (double alpha, double delta, int cc,
                                NumericVector ks, NumericVector fs) {
  NumericVector zero = NumericVector::create(0);
  double p0 = marginal_negative_binomial_cpp(zero,alpha,delta)[0];
  NumericVector ps = marginal_negative_binomial_cpp(ks, alpha, delta);
  NumericVector l_ps = log_ps_cpp(ks,alpha, delta);
  double result = sum(fs * (l_ps - log(1-p0)));
  return result;
}

// [[Rcpp::export]]
double log_factorial_cpp (int x) {
  if ((x % 1) > 0) {
    stop("Input must be an integer");
  }
  NumericVector OneToX(x);
  for (int i = 0; i < x; i++) {
    OneToX[i] = i+1;
  }
  return sum(log(OneToX));
}

// [[Rcpp::export]]
double loglike_unpenalized_cpp (NumericVector x, int ccc, int cc,
                                NumericVector ks, NumericVector fs) {
  double alpha = x[0];
  double delta = x[1];
  NumericVector zero = NumericVector::create(0);
  double p0 = marginal_negative_binomial_cpp(zero,alpha,delta)[0];
  double l_p0 = log_marginal_negative_binomial_cpp(zero,alpha,delta)[0];
  double lb = cc * log(1-p0) + (ccc-cc)*l_p0 + log_factorial_cpp(ccc)
    - log_factorial_cpp(ccc-cc);
  double lc = loglike_cond_cpp(alpha,delta,cc,ks,fs);
  return lb+lc;
}
