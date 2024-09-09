#include <Eigen/Core>
#include <iostream>
#include <LBFGSB.h>  
#include <RcppArmadillo.h>
#include <BBP_functions.h>
//#include <RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace LBFGSpp;

Eigen::MatrixXd example_cast_eigen(arma::mat arma_A) {
  
  Eigen::MatrixXd eigen_B = Eigen::Map<Eigen::MatrixXd>(arma_A.memptr(),
                                                        arma_A.n_rows,
                                                        arma_A.n_cols);
  
  return eigen_B;
}

arma::mat matrixxd_to_armamat2(Eigen::MatrixXd& eigen_A) {
  arma::mat arma_B = arma::mat(eigen_A.data(), eigen_A.rows(), eigen_A.cols(),
                               false, false);
  return arma_B;
}
arma::mat matrixxd_to_armamat(Eigen::MatrixXd eigen_A) {
  arma::mat arma_B = arma::mat(eigen_A.data(), eigen_A.rows(), eigen_A.cols(),
                               true,   // changed from false to true.
                               false); 
  return arma_B;
}


class BBPgame
{
private:
  const arma::mat& logaltruism;
  const arma::vec& income;
public:
  BBPgame(arma::mat& logaltruism, arma::vec& income) : income(income), logaltruism(logaltruism) {}
  
  double operator()(const VectorXd& x, VectorXd& grad)
  {
    Rcpp::Rcout << "tic"; 
    double lp = logpotential_fast_cpp(matrixxd_to_armamat(x), logaltruism, income);
    arma::mat gr = gr_logpotential_fast_cpp(matrixxd_to_armamat(x), logaltruism, income);
    arma::mat gra = arma::reshape(gr, income.n_elem*income.n_elem,1);
    grad = example_cast_eigen(gra);
    return lp;
  }
  
  
};

// [[Rcpp::export]]
arma::mat solve_BBP2(arma::mat& altruism, arma::vec& income)
  {
  int n = income.n_elem*income.n_elem;
  // Set up parameters
  LBFGSBParam<double> param;  // New parameter class
  param.max_linesearch = 5000;
  param.max_iterations = 500;

  // Create solver and function object
  LBFGSBSolver<double> solver(param);  // New solver class
  BBPgame fun(altruism, income);
  
  // Bounds
  VectorXd lb = VectorXd::Constant(n, 0);
  VectorXd ub = VectorXd::Constant(n, 9999);
  
  // Initial guess
  VectorXd x = VectorXd::Constant(n, 0.1);
  
  // x will be overwritten to be the best point found
  double fx;
  int niter = solver.minimize(fun, x, fx, lb, ub);
  
  Rcpp::Rcout << niter << " iterations" << std::endl;
  Rcpp::Rcout << "x = \n" << x.transpose() << std::endl;
  Rcpp::Rcout << "f(x) = " << fx << std::endl;
  
  return 0;
}
