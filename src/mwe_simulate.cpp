#include <RcppArmadillo.h>
#include <iostream>
#include <RcppParallel.h>
#include <random>

using namespace arma;


// [[Rcpp::export]]
mat berlin_mwe(int n, double delta0,double delta1,double sigma, mat distance, mat kinship,  mat capacity, vec income,int reps,int seed,int rounds) {
  mat error = zeros(n,n);
  mat altruism = 1/(1+exp(-(1+error)));
  return(error);
}

// [[Rcpp::export]]
mat berlin_mwe2(int n) {
  mat error = zeros(n,n);
  mat altruism = 1/(1+exp(-(1+error)));
  return(error);
}
