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
//'@name Cauchy_fixed_sigma_marginal
//'@title Marginal likelihood of Cauchy prior when varinace is fixed
//‘@usage Cauchy_fixed_sigma_marginal(index_vec_input, Sigma, z, zSigmaz, p, tau, p_S, y_sigma)
//'@return marginal likelihood

double outlier_marginal_fun_indi(const double &zSigmaz_S, const double &tau,
                         const double & p_S, const double & y_sigma){
  double result=pow(tau/(1.00+tau),p_S/2.00)*exp(-tau*zSigmaz_S/(2*y_sigma*(1.00+tau)) ) ;
  return  result;
}

// [[Rcpp::export]]
double outlier_Cauchy_fixed_sigma_marginal(const arma::uvec & index_vec_input, const arma::mat & Sigma,
                           const arma::vec & z, const arma::vec &tau,  const double & p_S, const double & y_sigma){
  
  unsigned int S=tau.n_elem;
  arma::uvec index_vec=index_vec_input-1;
  
  arma::mat Sigma_S=Sigma.submat(index_vec,index_vec);
  arma::mat Sigma_S_inv=pinv(Sigma_S,0.00001,"std");
  arma::mat sub_z=z.rows(index_vec);
    
    double det_S=det(Sigma_S_inv);
    det_S=fabs(det_S);
    
  arma::mat zSigmaz_S=sub_z.t()*Sigma_S_inv*sub_z;
  //std::cout<<zSigmaz_S<<std::endl;
  arma::vec b(S);
  for(unsigned int j=0;j<S;j++){
    b[j]=outlier_marginal_fun_indi(zSigmaz_S(0,0),tau[j],p_S,y_sigma);
  }
  double results=log(mean(b))+0.5*log(det_S);
  
  return(results);
}
