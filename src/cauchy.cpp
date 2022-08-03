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
//'@name Cauchy_marginal
//'@title Marginal likelihood of Cauchy prior
//‘@usage Cauchy_marginal(index_vec_input, Sigma, z, zSigmaz, p, tau, p_S)
//'@return marginal likelihood 

double marginal_fun_indi(const double &zSigmaz_S, const double &tau, const double &p, const double & zSigmaz,const double & p_S){
  double result=pow(tau/(1.00+tau),p_S/2.00)*exp(p/2.00*log(zSigmaz)-p/2.00*log(zSigmaz-zSigmaz_S/(1.00+tau) ) );
  return  result;
}

// [[Rcpp::export]]
double Cauchy_marginal(const arma::uvec & index_vec_input, const arma::mat & Sigma,
                           const arma::vec & z, const double & zSigmaz, const arma::vec &tau, const double &p, const double & p_S){
  
  unsigned int S=tau.n_elem;
  arma::uvec index_vec=index_vec_input-1;
  
  arma::mat Sigma_S=Sigma.submat(index_vec,index_vec);
  arma::mat Sigma_S_inv=pinv(Sigma_S);
  arma::mat sub_z=z.rows(index_vec);
  arma::mat zSigmaz_S=sub_z.t()*Sigma_S_inv*sub_z;
  //std::cout<<zSigmaz_S<<std::endl;
  arma::vec b(S);
  for(unsigned int j=0;j<S;j++){
    b[j]=marginal_fun_indi(zSigmaz_S(0,0),tau[j],p,zSigmaz,p_S);
  }
  double results=log(mean(b));
  
  return(results);
}
