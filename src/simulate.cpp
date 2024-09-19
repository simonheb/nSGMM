#define ARMA_DONT_USE_OPENMP
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
  if (sd<0) {
    Rcpp::Rcout <<"negative sigma"<<endl;
  }
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
  if (sd<0) { Rcpp::Rcout <<"negative sigma"<<endl;}
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



//this is to ensure that consecutive indices don't result in consecutive seeds
// [[Rcpp::export]]
int seedfromindex(int index) {
  for (int i=0; i<10; ++i) {
  std::mt19937 mt_rand(index);
  return(mt_rand());
  }
}

// [[Rcpp::export]]
mat simulate_BBP_cpp(int n, double delta0,double delta1,double sigma, mat distance, mat kinship,  mat capacity, vec income,int reps,int seed,int rounds) {
  if (seed==0) seed=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  
  mat finalMatrix = zeros(reps,13);
  vec converged(reps);
  vec roundsc(reps);
  
  for (int i=0; i<reps;i++) {
    //std::cout << "s:" <<i;
    mat error = normal_error_matrix(n,0,sigma,seedfromindex(seed)+i);
    mat altruism = 1/(1+exp(-(delta0+delta1*kinship+error)));
    altruism.diag().ones();
    int rounds_;
    bool converged_;
    mat eqtrans2=equilibrate_cpp_fast8_debug(altruism,income,capacity,zeros(income.n_elem,income.n_elem),true, rounds_,converged_,rounds);
    eqtrans2.elem( find(eqtrans2) ).ones();
    vec moments = compute_moments_cpp(eqtrans2,kinship,distance,income);

    for(int mom=0;mom<13;mom++)
      finalMatrix(i,mom) = moments(mom);
    roundsc(i) = rounds_;
    converged(i) = converged_;
  }
  if (mean(converged)<0.9) {
    //Rcpp::Rcout << "" << floor(mean(converged)*100) <<"%" ;
    finalMatrix.zeros();
  }
  return(finalMatrix);
}

struct simulate_BBP_worker_smart : public RcppParallel::Worker
{
  // source matrix
  int n,reps,seed,rounds;
  double delta0,delta1,sigma;
  const mat distance;
  const mat kinship;
  const mat capacity;
  const vec income;

  // destination matrix
  RcppParallel::RMatrix<double> output; //it seems that using a arma::mat here is problematic https://stackoverflow.com/questions/27523979/stack-imbalance-with-rcppparallel
  RcppParallel::RMatrix<double> convergance; //it seems that using a arma::mat here is problematic https://stackoverflow.com/questions/27523979/stack-imbalance-with-rcppparallel
  
  // initialize with source and destination
  simulate_BBP_worker_smart(Rcpp::NumericMatrix output,Rcpp::NumericMatrix convergance,  int n, double delta0,double delta1,double sigma, const mat& distance, const mat& kinship,const  mat& capacity,const vec& income, int reps,int seed,int rounds) 
    : output(output), convergance(convergance), n(n), delta0(delta0), delta1(delta1), sigma(sigma), distance(distance), kinship(kinship), capacity(capacity), income(income), reps(reps), seed(seed), rounds(rounds)
  {}
  
  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
    size_t i = begin;
    while (i < end) {
      mat error = normal_error_matrix(n,0,sigma,seedfromindex(seed)+i);
      mat altruism = 1/(1+exp(-(delta0+delta1*kinship+error)));
      altruism.diag().ones();
      
      int rounds_;
      bool converged_;
      mat eqtrans2=equilibrate_cpp_fast8_debug(altruism,income,capacity,zeros(income.n_elem,income.n_elem),true,rounds_,converged_,rounds);
      eqtrans2.elem( find(eqtrans2) ).ones();
      
      vec moments=compute_moments_cpp(eqtrans2,kinship,distance,income);
      for(int mom=0;mom<13;mom++)
        output(i,mom) = moments(mom);
      convergance(i,0) = converged_;
      //as a "trick" to speed up computations in situations with low convergence, we can skip rounds, whenever there convergence ratio in the past was too low
      //We achieve this, by skippnig (using the result from some previous round) two iterations for every one that did not converge
      if (i>begin+10 & i < end-10 & !converged_) {
        i = i + 1;
        for(int mom=0;mom<13;mom++)
          output(i,mom) = output(begin+round((i-begin)*1/3),mom);
        convergance(i,0) = convergance(begin+round((i-begin)*1/3),0) ;
        i = i + 1;
        for(int mom=0;mom<13;mom++)
          output(i,mom) = output(begin+round((i-begin)*2/3),mom);
        convergance(i,0) = convergance(begin+round((i-begin)*2/3),0) ;
      }
      i = i+1;
    }
    
  }
  
};

// [[Rcpp::export]]
Rcpp::NumericMatrix simulate_BBP_cpp_parallel(int n, double delta0,double delta1,double sigma, const mat& distance, const mat& kinship, const mat& capacity,const vec& income,int reps,int seed,int rounds)  {
  if (seed==0) seed=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  
  
  // allocate the output matrix
  Rcpp::NumericMatrix output(reps,13);
  Rcpp::NumericMatrix convergance(reps,1);
  
  
  // functor (pass input and output matrixes)
  simulate_BBP_worker_smart sim(output, convergance, n, delta0, delta1, sigma, distance, kinship, capacity, income, reps, seed, rounds);
  
  
  // call parallelFor to do the work
  parallelFor(0, reps, sim);
  
  if (mean(convergance)<0.85) {
    //Rcpp::Rcout  << "" << floor(mean(convergance)*100) <<"%";
  }
  if (mean(convergance)<0.20) {
    Rcpp::NumericMatrix z(reps,13);
    return z;
  }
  
  // return the output matrix
  return output;
}







// [[Rcpp::export]]
mat simulate_BBP_cpp73(int n, double delta0,double delta1,double sigma, mat distance, mat kinship,  mat capacity, vec income,int reps,int seed,int rounds) {
  mat error = zeros(n,n);
  mat altruism = 1/(1+exp(-(1+error)));
  return(error);
}



// [[Rcpp::export]]
mat simulate_BBP_cpp1(int n, double delta0,double delta1,double sigma, mat distance, mat kinship,  mat capacity, vec income,int reps,int seed,int rounds) {
  if (seed==0) seed=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  
  mat finalMatrix = zeros(reps,13);
  vec converged(reps);
  vec roundsc(reps);
  
  for (int i=0; i<reps;i++) {
    //std::cout << "s:" <<i;
    mat error = normal_error_matrix(n,0,sigma,seedfromindex(seed)+i);
    mat altruism = 1/(1+exp(-(delta0+delta1*kinship+error)));
    altruism.diag().ones();
    int rounds_;
    bool converged_;
    mat eqtrans2=equilibrate_cpp_fast8_debug(altruism,income,capacity,zeros(income.n_elem,income.n_elem),true, rounds_,converged_,rounds);
    eqtrans2.elem( find(eqtrans2) ).ones();
    vec moments = compute_moments_cpp(eqtrans2,kinship,distance,income);
  }    
  return(finalMatrix);
}




// [[Rcpp::export]]
mat simulate_BBP_cpp2(int n, double delta0,double delta1,double sigma, mat distance, mat kinship,  mat capacity, vec income,int reps,int seed,int rounds) {
  if (seed==0) seed=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  
  mat finalMatrix = zeros(reps,13);
  vec converged(reps);
  vec roundsc(reps);
  
  for (int i=0; i<reps;i++) {
    //std::cout << "s:" <<i;
    mat error = normal_error_matrix(n,0,sigma,seedfromindex(seed)+i);
    mat altruism = 1/(1+exp(-(delta0+delta1*kinship+error)));
    altruism.diag().ones();
    int rounds_;
    bool converged_;
    mat eqtrans2=equilibrate_cpp_fast8_debug(altruism,income,capacity,zeros(income.n_elem,income.n_elem),true, rounds_,converged_,rounds);
    eqtrans2.elem( find(eqtrans2) ).ones();
  }    
  return(finalMatrix);
}

// [[Rcpp::export]]
mat simulate_BBP_cpp3(int n, double delta0,double delta1,double sigma, mat distance, mat kinship,  mat capacity, vec income,int reps,int seed,int rounds) {
  if (seed==0) seed=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  
  mat finalMatrix = zeros(reps,13);
  vec converged(reps);
  vec roundsc(reps);
  
  for (int i=0; i<reps;i++) {
    //std::cout << "s:" <<i;
    mat error = normal_error_matrix(n,0,sigma,seedfromindex(seed)+i);
    mat altruism = 1/(1+exp(-(delta0+delta1*kinship+error)));
    altruism.diag().ones();
  }    
  return(finalMatrix);
}
// [[Rcpp::export]]
mat simulate_BBP_cpp4(int n, double delta0,double delta1,double sigma, mat distance, mat kinship,  mat capacity, vec income,int reps,int seed,int rounds) {
  if (seed==0) seed=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  
  mat finalMatrix = zeros(reps,13);
  vec converged(reps);
  vec roundsc(reps);
  
  for (int i=0; i<reps;i++) {
    //std::cout << "s:" <<i;
    mat error = normal_error_matrix(n,0,sigma,seedfromindex(seed)+i);
  }    
  return(finalMatrix);
}
// [[Rcpp::export]]
mat simulate_BBP_cpp5(int n, double delta0,double delta1,double sigma, mat distance, mat kinship,  mat capacity, vec income,int reps,int seed,int rounds) {
  if (seed==0) seed=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  
  mat finalMatrix = zeros(reps,13);
  vec converged(reps);
  vec roundsc(reps);
  
  for (int i=0; i<reps;i++) {
    //std::cout << "s:" <<i;
    mat error = normal_error_matrix(n,0,sigma,seedfromindex(seed)+i);
  }    
  return(finalMatrix);
}
// [[Rcpp::export]]
mat simulate_BBP_cpp6(int n, double delta0,double delta1,double sigma, mat distance, mat kinship,  mat capacity, vec income,int reps,int seed,int rounds) {
  if (seed==0) seed=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  
  mat finalMatrix = zeros(reps,13);
  vec converged(reps);
  vec roundsc(reps);
  
  return(finalMatrix);
}