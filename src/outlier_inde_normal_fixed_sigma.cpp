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

//double outlier_marginal_fun_indi(const double &zSigmaz_S, const double &tau,
  //                       const double & p_S, const double & det_S){
//  double result=p_S/2.00*log(tau)-0.50*log(det_S)+(zSigmaz_S/(2.00)) ;
//  return  result;
//}

// [[Rcpp::export]]
double outlier_ind_Normal_marginal(const arma::uvec & index_vec_input, const arma::mat & Sigma,
                           const arma::vec & z, const double &tau,  const double & p_S, const double & y_sigma){
  
  arma::uvec index_vec=index_vec_input-1;
  
  arma::mat Sigma_S=Sigma.submat(index_vec,index_vec);
  arma::mat A=tau*arma::eye(p_S,p_S);

  arma::mat Sigma_S_I_inv=pinv(Sigma_S+A,0.00001,"std");
  arma::mat Sigma_S_inv=pinv(Sigma_S,0.00001,"std");
  
    double det_S=det(Sigma_S_inv);
    double det_I_S=det(Sigma_S_I_inv);
    det_I_S=fabs(det_I_S);
    det_S=fabs(det_S);
//    std::cout<<det_I_S <<std::endl;
//    std::cout<<det_S <<std::endl;
 
  arma::mat sub_z=z.rows(index_vec);
  arma::mat zSigmaz_S=sub_z.t()*Sigma_S_inv*sub_z;
  arma::mat zSigmaz_I_S=sub_z.t()*Sigma_S_I_inv*sub_z;
  double b;
//  std::cout<<zSigmaz_S(0,0) <<std::endl;
//  std::cout<<zSigmaz_I_S(0,0) <<std::endl;
  
  b=+0.5*(log(det_S)+log(det_I_S))-0.5*(zSigmaz_S(0,0)-zSigmaz_I_S(0,0));
  double results=(b);
  
  return(results);
}
