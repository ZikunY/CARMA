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
//'@name Normal_marginal
//'@title Marginal likelihood of Normal prior
//‘@usage Normal_marginal(index_vec_input, Sigma, z, zSigmaz, p, tau, p_S)
//'@return marginal likelihood

double normal_marginal_fun_indi(const double &zSigmaz_S, const double &tau, const double &p, const double & zSigmaz,const double & p_S){
    double result=p_S/2.00*log(tau/(1.00+tau))-p/2.00*log(1-zSigmaz_S/((1.00+tau)*zSigmaz) ) ;
  return  result;
}

// [[Rcpp::export]]
double Normal_marginal(const arma::uvec & index_vec_input, const arma::mat & Sigma,
                           const arma::vec & z, const double & zSigmaz,  const double &tau, const double &p, const double & p_S){
  
  arma::uvec index_vec=index_vec_input-1;
  arma::mat Sigma_S=Sigma.submat(index_vec,index_vec);
  arma::mat Sigma_S_inv=pinv(Sigma_S);
  arma::mat sub_z=z.rows(index_vec);
  arma::mat zSigmaz_S=sub_z.t()*Sigma_S_inv*sub_z;
  double b;
  
  
  b=normal_marginal_fun_indi(zSigmaz_S(0,0),tau,p,zSigmaz,p_S);
  
  double results=b;
  
  return(results);
}
