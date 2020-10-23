// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <tictoc.h>
#include <BBP_functions.h>
#include <network_functions.h>
#include <RcppParallel.h>
#include <random>

using namespace arma;


// [[Rcpp::export]]
vec random_normal_seed(int n, double mean, double sd, int seed) {
  if (sd<0) {std::cout <<"negative sigma"<<endl;}
  std::mt19937 rng(seed);
  std::normal_distribution<double> distribution(mean,sd);
  vec ret(n);
  for (int i=0; i<n; ++i) {
    ret(i)=distribution(rng);
  }
  return(ret);
}
//draws a symmetric random matrix of width n, where the off-diagonal elements are iid normal with mean an sd.
// [[Rcpp::export]]
mat symmetrix_normal_error_matrix(int n, double mean, double sd, int seed) {
 uvec tri_coords = find(trimatu(ones(n,n),1));
 mat ret(n,n);
 ret.elem(tri_coords)= random_normal_seed(n*(n-1)/2,mean,sd,seed);
 return(symmatu(ret));
}
// [[Rcpp::export]]
mat normal_error_matrix(int n, double mean, double sd, int seed) {
  if (sd<0) {std::cout <<"negative sigma"<<endl;}
  std::mt19937 rng(seed);
  std::normal_distribution<double> distribution(mean,sd);
  mat ret(n,n);
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j) {
      if (i!=j) {
        ret(i,j)=distribution(rng);
      }
    }
  }
  return(ret);
}



//this is to ensure that consequtive indices dont result in consequctive seeds
// [[Rcpp::export]]
int seedfromindex(int index) {
  for (int i=0; i<10; ++i) {
  std::mt19937 mt_rand(index);
  return(mt_rand());
  }
}

// [[Rcpp::export]]
mat simulate_BBP_cpp(int n, double delta0,double delta1,double sigma, mat distance, mat kinship,  mat capacity, vec income,vec theta,int reps,int seed,int rounds) {
  if (seed==0) seed=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  
  mat finalMatrix = zeros(reps,12);
  vec converged(reps);
  vec roundsc(reps);
  
  for (int i=0; i<reps;i++) {
    //std::cout << "s:" <<i;
    mat error = normal_error_matrix(n,0,sigma,seedfromindex(seed)+i);
    mat altruism = 1/(1+exp(-(delta0+delta1*kinship+error)));
    altruism.diag().ones();
    int rounds_;
    bool converged_;
    mat eqtrans2=equilibrate_cpp_fast8_debug(altruism,income,capacity,true, rounds_,converged_,rounds);
    eqtrans2.elem( find(eqtrans2) ).ones();
    vec moments = compute_moments_cpp(eqtrans2,kinship,distance,income,theta);

    
    for(int mom=0;mom<12;mom++)
      finalMatrix(i,mom) = moments(mom);
    roundsc(i) = rounds_;
    converged(i) = converged_;
  }
  if (mean(converged)<0.9) {
    std::cout << "(" << floor(mean(converged)*100) <<"%c)" ;
    finalMatrix.zeros();
  }
  return(finalMatrix);
}

// [[Rcpp::export]]
mat simulate_BBP_symmetric_cpp(int n, double delta0,double delta1,double sigma, mat distance, mat kinship,  mat capacity, vec income,vec theta,int reps,int seed,int rounds) {
  if (seed==0) seed=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  
  mat finalMatrix = zeros(reps,12);
  vec converged(reps);
  vec roundsc(reps);
  
  for (int i=0; i<reps;i++) {
    //std::cout << "s:" <<i;
    mat error = symmetrix_normal_error_matrix(n,0,sigma,seedfromindex(seed)+i);
    mat altruism = 1/(1+exp(-(delta0+delta1*kinship+error)));
    altruism.diag().ones();
    int rounds_;
    bool converged_;
    mat eqtrans2=equilibrate_cpp_fast8_debug(altruism,income,capacity,true, rounds_,converged_,rounds);
    eqtrans2.elem( find(eqtrans2) ).ones();
    vec moments = compute_moments_cpp(eqtrans2,kinship,distance,income,theta);

    
    for(int mom=0;mom<12;mom++)
      finalMatrix(i,mom) = moments(mom);
    roundsc(i) = rounds_;
    converged(i) = converged_;
  }
  if (mean(converged)<0.9) {
    std::cout << "(" << floor(mean(converged)*100) <<"%c)" ;
    finalMatrix.zeros();
  }
  return(finalMatrix);
}




struct simulate_BBP_worker : public RcppParallel::Worker
{
  // source matrix
  int n,reps,seed,rounds;
  double delta0,delta1,sigma;
  const mat distance;
  const mat kinship;
  const mat capacity;
  const vec income;
  const vec theta;

  // destination matrix
  RcppParallel::RMatrix<double> output; //it seems that using a arma::mat here is problematic https://stackoverflow.com/questions/27523979/stack-imbalance-with-rcppparallel
  RcppParallel::RMatrix<double> convergance; //it seems that using a arma::mat here is problematic https://stackoverflow.com/questions/27523979/stack-imbalance-with-rcppparallel
  
  // initialize with source and destination
  simulate_BBP_worker(Rcpp::NumericMatrix output,Rcpp::NumericMatrix convergance,  int n, double delta0,double delta1,double sigma, const mat& distance, const mat& kinship,const  mat& capacity,const vec& income, const vec& theta,  int reps,int seed,int rounds) 
    : output(output), convergance(convergance), n(n), delta0(delta0), delta1(delta1), sigma(sigma), distance(distance), kinship(kinship), capacity(capacity), income(income), theta(theta), reps(reps), seed(seed), rounds(rounds)
  {}
  
  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i = begin; i < end; i++) {
      mat error = normal_error_matrix(n,0,sigma,seedfromindex(seed)+i);
      mat altruism = 1/(1+exp(-(delta0+delta1*kinship+error)));
      altruism.diag().ones();
      
      int rounds_;
      bool converged_;
      mat eqtrans2=equilibrate_cpp_fast8_debug(altruism,income,capacity,true, rounds_,converged_,rounds);
      eqtrans2.elem( find(eqtrans2) ).ones();
      
      vec moments=compute_moments_cpp(eqtrans2,kinship,distance,income,theta);
      for(int mom=0;mom<12;mom++)
        output(i,mom) = moments(mom);
      convergance(i,0) = converged_;
      
    }
    
  }
  
};
// [[Rcpp::export]]
Rcpp::NumericMatrix simulate_BBP_cpp_parallel(int n, double delta0,double delta1,double sigma, const mat& distance, const mat& kinship, const mat& capacity,const vec& income, const vec& theta,int reps,int seed,int rounds)  {
  if (seed==0) seed=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  
    
  
  // allocate the output matrix
  Rcpp::NumericMatrix output(reps,12);
  Rcpp::NumericMatrix convergance(reps,1);
  

    
  // SquareRoot functor (pass input and output matrixes)
  simulate_BBP_worker sim(output, convergance, n, delta0, delta1, sigma, distance, kinship, capacity, income, theta,reps, seed, rounds);
  

  // call parallelFor to do the work
  parallelFor(0, reps, sim);
  
  if (mean(convergance)<0.85) {
    Rcpp::Rcout  << "(" << floor(mean(convergance)*100) <<"%)";
  }
  if (mean(convergance)<0.20) {
    Rcpp::NumericMatrix z(reps,12);
    return z;
  }
  
  // return the output matrix
  return output;
}
