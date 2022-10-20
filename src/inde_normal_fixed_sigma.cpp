//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppGSL)]]
#include <cmath>
#include <iostream>
#include <RcppArmadillo.h>
//#include <RcppGSL.h>
//#include <gsl_rng.h>
//#include <gsl_randist.h>
//#include <gsl_cdf.h>
//#include <gsl_statistics_double.h>
//'@name Normal_fixed_sigma_marginal
//'@title Marginal likelihood of Normal prior when varinace is fixed
//‘@usage Normal_fixed_sigma_marginal(index_vec_input, Sigma, z, zSigmaz, p, tau, p_S, y_sigma)
//'@return marginal likelihood

double ind_normal_sigma_fixed_marginal_fun_indi(const double &zSigmaz_S, const double &tau,
                         const double & p_S, const double & det_S){
  double result=p_S/2.00*log(tau)-0.50*log(det_S)+(zSigmaz_S/(2.00)) ;
  return  result;
}

// [[Rcpp::export]]
double ind_Normal_fixed_sigma_marginal(const arma::uvec & index_vec_input, const arma::mat & Sigma,
                           const arma::vec & z, const double &tau,  const double & p_S, const double & y_sigma){
  
  arma::uvec index_vec=index_vec_input-1;
  
  arma::mat Sigma_S=Sigma.submat(index_vec,index_vec);
  arma::mat A=tau*arma::eye(p_S,p_S);
  double det_S=det(Sigma_S+A);
  arma::mat Sigma_S_inv=pinv(Sigma_S+A,0.00001,"std");
  arma::mat sub_z=z.rows(index_vec);
  arma::mat zSigmaz_S=sub_z.t()*Sigma_S_inv*sub_z;
    double b;
  
    b=ind_normal_sigma_fixed_marginal_fun_indi(zSigmaz_S(0,0),tau,p_S,det_S);
  
  double results=(b);
  
  return(results);
}
