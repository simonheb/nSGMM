#ifndef network_functions_H
#define network_functions_H

#include <RcppArmadillo.h>
#include <iostream>

using namespace arma;

double Ccomponents(const mat& transferstructure,uvec & t_components_csize,uvec & t_components_membership, mat & t_components_matrix, mat & t_conmat);
double component_counts(const mat& transferstructure);
vec compute_moments_cpp(const mat& btransfers,const mat& kinship,const mat& distance,const vec& income) ;

#endif