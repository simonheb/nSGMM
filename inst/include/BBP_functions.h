#ifndef BBP_functions_H
#define BBP_functions_H

#include <RcppArmadillo.h>
#include <iostream>
#include <tictoc.h>

using namespace arma;

mat equilibrate_cpp_fast7_smarter(const mat& altruism, const vec& income, const mat& capacity,int modmode=5);
mat equilibrate_cpp_fast8_smarter(const mat& altruism, const vec& income, const mat& capacity,int modmode=5);

#endif
