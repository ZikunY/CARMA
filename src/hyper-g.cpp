//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppGSL)]]
#include <cmath>
#include <iostream>
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_sf_hyperg.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>

//'@name hyper_g_marginal
//'@title Marginal likelihood of hyper-g prior
//‘@usage hyper_g_marginal(index_vec_input, Sigma, z, zSigmaz, p, tau, p_S)
//'@return marginal likelihood

// [[Rcpp::export]]
double hyper_g_marginal(const arma::uvec & index_vec_input, const arma::mat & Sigma,
                           const arma::vec & z, const double & zSigmaz,const double &p,const arma::vec &tau,  const double & p_S){
    
    arma::uvec index_vec=index_vec_input-1;
    arma::mat Sigma_S=Sigma.submat(index_vec,index_vec);
    arma::mat Sigma_S_inv=pinv(Sigma_S);
    arma::mat sub_z=z.rows(index_vec);
    arma::mat zSigmaz_S=sub_z.t()*Sigma_S_inv*sub_z;
   
    double result=gsl_sf_hyperg_2F1_renorm(p/2.00,1,p_S/2.00+1.5,zSigmaz_S(0,0)/zSigmaz)/(p_S+1.00);
    return log(result);
}
