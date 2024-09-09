#ifndef BBP_functions_H
#define BBP_functions_H

#include <RcppArmadillo.h>
#include <iostream>
#include <tictoc.h>

using namespace arma;

//mat equilibrate_cpp_fast7_smarter(const mat& altruism, const vec& income, const mat& capacity,int modmode=5);
mat equilibrate_cpp_fast8_smarter(const mat& altruism, const vec& income, const mat& capacity, int maxrounds=500);
mat equilibrate_cpp_fast8_debug(const mat& altruism, const vec& income, const mat& capacity, mat transfers,bool verbose, int& r, bool& NE, int maxrounds=500);

mat gr_logpotential_fast_cpp(const vec& transfers,const mat& logaltruism, const vec& income);
double logpotential_fast_cpp(const vec& transfers,const mat& logaltruism, const vec& income);
mat gr_logpotential_cpp(const mat& transfers,const mat& logaltruism, const vec& income);
double logpotential_cpp(const mat& transfers,const mat& logaltruism, const vec& income);


    
    

        
#endif

    
