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



// [[Rcpp::export]]
double outlier_Normal_fixed_sigma_marginal(const arma::uvec & index_vec_input, const arma::mat & Sigma,
                           const arma::vec & z, const double &tau,  const double & p_S, const double & y_sigma){
  
  arma::uvec index_vec=index_vec_input-1;
  
  arma::mat Sigma_S=Sigma.submat(index_vec,index_vec);
  arma::mat Sigma_S_inv=pinv(Sigma_S,0.00001,"std");
  arma::mat sub_z=z.rows(index_vec);
    double det_S=det(Sigma_S);
  arma::mat zSigmaz_S=sub_z.t()*Sigma_S_inv*sub_z;
    double results=-0.5*log(det_S)-(tau*zSigmaz_S(0,0)/(2*y_sigma*(1.00+tau)));
  
  return(results);
}
