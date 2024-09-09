#include <iostream>
#include <RcppArmadillo.h>
//#include <RcppEnsmallen.h>  
#include <BBP_functions.h>
/*
// [[Rcpp::depends(RcppEnsmallen)]]
  
  class BBPgamea
  {
  public:
    
    BBPgamea(arma::mat& logaltruism, arma::vec& income) : income(income), logaltruism(logaltruism) {}
    
    // Return the objective function for model parameters beta.
    double EvaluateWithGradient(const arma::mat& theta, arma::mat& gradient)
    {
      //Rcpp::Rcout << income.n_elem; 
      gradient = gr_logpotential_fast_cpp(theta, logaltruism, income);
      //Rcpp::Rcout << 5; 
      double lp = logpotential_fast_cpp(theta, logaltruism, income); 
      //Rcpp::Rcout << 6; 
      return lp;
    }
    
  private:
    
    const arma::mat& logaltruism;
    const arma::vec& income;
  };
  
  // [[Rcpp::export]]
  arma::mat solve_BBP(arma::mat& altruism, arma::vec& income)
  {
    BBPgame btest(altruism, income);
    
    // create a Limited-memory BFGS optimizer object with default parameters
    ens::L_BFGS opt;
    opt.MaxIterations() = 10;
    
    // initial point (uniform random)
    arma::mat theta(income.n_elem,income.n_elem, arma::fill::zeros);
    Rcpp::Rcout << 1;
    opt.Optimize(btest, theta);
    Rcpp::Rcout << 2;
    
    return theta;
  }
  
    */
