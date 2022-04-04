//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppGSL)]]
#include <cmath>
#include <iostream>
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_sf_hyperg.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>

//'@name hyper_g_fixed_sigma_marginal
//'@title Marginal likelihood of hyper-g prior when varinace is fixed
//‘@usage hyper_g_marginal(index_vec_input, Sigma, z, zSigmaz, p, tau, p_S, y_sigma)
//'@return marginal likelihood

// [[Rcpp::export]]
double hyper_g_fixed_sigma_marginal(const arma::uvec & index_vec_input, const arma::mat & Sigma,
                           const arma::vec & z,const arma::vec &tau, const double & p_S, const double & y_sigma){
    
    arma::uvec index_vec=index_vec_input-1;
    arma::mat Sigma_S=Sigma.submat(index_vec,index_vec);
    arma::mat Sigma_S_inv=pinv(Sigma_S);
    arma::mat sub_z=z.rows(index_vec);
    arma::mat zSigmaz_S=sub_z.t()*Sigma_S_inv*sub_z;
    
    double result=gsl_sf_hyperg_1F1(1.00,p_S/2.00+1.5,zSigmaz_S(0,0)/(2.00*y_sigma))/(p_S/2.00+.5);

    return log(result);
}
