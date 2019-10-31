// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/nSGMM.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// utlity_cppvec
double utlity_cppvec(const vec& consumption, const rowvec& altruism);
RcppExport SEXP _nSGMM_utlity_cppvec(SEXP consumptionSEXP, SEXP altruismSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type consumption(consumptionSEXP);
    Rcpp::traits::input_parameter< const rowvec& >::type altruism(altruismSEXP);
    rcpp_result_gen = Rcpp::wrap(utlity_cppvec(consumption, altruism));
    return rcpp_result_gen;
END_RCPP
}
// mynegutility_cpp
double mynegutility_cpp(const vec& mytransfers, int i, mat transfers, const mat& altruism, const vec& income);
RcppExport SEXP _nSGMM_mynegutility_cpp(SEXP mytransfersSEXP, SEXP iSEXP, SEXP transfersSEXP, SEXP altruismSEXP, SEXP incomeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type mytransfers(mytransfersSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< mat >::type transfers(transfersSEXP);
    Rcpp::traits::input_parameter< const mat& >::type altruism(altruismSEXP);
    Rcpp::traits::input_parameter< const vec& >::type income(incomeSEXP);
    rcpp_result_gen = Rcpp::wrap(mynegutility_cpp(mytransfers, i, transfers, altruism, income));
    return rcpp_result_gen;
END_RCPP
}
// mynegutility_cpp_consumption
double mynegutility_cpp_consumption(const vec& mytransfers, int i, mat transfers, const mat& altruism, vec consumption);
RcppExport SEXP _nSGMM_mynegutility_cpp_consumption(SEXP mytransfersSEXP, SEXP iSEXP, SEXP transfersSEXP, SEXP altruismSEXP, SEXP consumptionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type mytransfers(mytransfersSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< mat >::type transfers(transfersSEXP);
    Rcpp::traits::input_parameter< const mat& >::type altruism(altruismSEXP);
    Rcpp::traits::input_parameter< vec >::type consumption(consumptionSEXP);
    rcpp_result_gen = Rcpp::wrap(mynegutility_cpp_consumption(mytransfers, i, transfers, altruism, consumption));
    return rcpp_result_gen;
END_RCPP
}
// get_BBP_BR_analytically_cpp_inline
mat get_BBP_BR_analytically_cpp_inline(uword i, mat transfers, const vec& income, const mat& altruism, const mat& capacities);
RcppExport SEXP _nSGMM_get_BBP_BR_analytically_cpp_inline(SEXP iSEXP, SEXP transfersSEXP, SEXP incomeSEXP, SEXP altruismSEXP, SEXP capacitiesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< uword >::type i(iSEXP);
    Rcpp::traits::input_parameter< mat >::type transfers(transfersSEXP);
    Rcpp::traits::input_parameter< const vec& >::type income(incomeSEXP);
    Rcpp::traits::input_parameter< const mat& >::type altruism(altruismSEXP);
    Rcpp::traits::input_parameter< const mat& >::type capacities(capacitiesSEXP);
    rcpp_result_gen = Rcpp::wrap(get_BBP_BR_analytically_cpp_inline(i, transfers, income, altruism, capacities));
    return rcpp_result_gen;
END_RCPP
}
// BBP_in_equilibrium_YaT_cpp
bool BBP_in_equilibrium_YaT_cpp(const mat& transfers, const vec& income, const mat& altruism, const mat& capacities, double tolerance);
RcppExport SEXP _nSGMM_BBP_in_equilibrium_YaT_cpp(SEXP transfersSEXP, SEXP incomeSEXP, SEXP altruismSEXP, SEXP capacitiesSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type transfers(transfersSEXP);
    Rcpp::traits::input_parameter< const vec& >::type income(incomeSEXP);
    Rcpp::traits::input_parameter< const mat& >::type altruism(altruismSEXP);
    Rcpp::traits::input_parameter< const mat& >::type capacities(capacitiesSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(BBP_in_equilibrium_YaT_cpp(transfers, income, altruism, capacities, tolerance));
    return rcpp_result_gen;
END_RCPP
}
// meannanrm
mat meannanrm(mat x, int dim);
RcppExport SEXP _nSGMM_meannanrm(SEXP xSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(meannanrm(x, dim));
    return rcpp_result_gen;
END_RCPP
}
// consumption_weights_cpp
mat consumption_weights_cpp(const mat& alphas, const mat transfers, const mat& t_conmat);
RcppExport SEXP _nSGMM_consumption_weights_cpp(SEXP alphasSEXP, SEXP transfersSEXP, SEXP t_conmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const mat >::type transfers(transfersSEXP);
    Rcpp::traits::input_parameter< const mat& >::type t_conmat(t_conmatSEXP);
    rcpp_result_gen = Rcpp::wrap(consumption_weights_cpp(alphas, transfers, t_conmat));
    return rcpp_result_gen;
END_RCPP
}
// BBP_c_from_atY_cpp
vec BBP_c_from_atY_cpp(const mat& alphas, const mat& transferstructure, const vec& incomes, const mat& t_components_matrix, const uvec& t_components_csize, const mat& t_conmat);
RcppExport SEXP _nSGMM_BBP_c_from_atY_cpp(SEXP alphasSEXP, SEXP transferstructureSEXP, SEXP incomesSEXP, SEXP t_components_matrixSEXP, SEXP t_components_csizeSEXP, SEXP t_conmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const mat& >::type transferstructure(transferstructureSEXP);
    Rcpp::traits::input_parameter< const vec& >::type incomes(incomesSEXP);
    Rcpp::traits::input_parameter< const mat& >::type t_components_matrix(t_components_matrixSEXP);
    Rcpp::traits::input_parameter< const uvec& >::type t_components_csize(t_components_csizeSEXP);
    Rcpp::traits::input_parameter< const mat& >::type t_conmat(t_conmatSEXP);
    rcpp_result_gen = Rcpp::wrap(BBP_c_from_atY_cpp(alphas, transferstructure, incomes, t_components_matrix, t_components_csize, t_conmat));
    return rcpp_result_gen;
END_RCPP
}
// BBP_T_from_tYc_cpp
mat BBP_T_from_tYc_cpp(mat transferstructure, const vec& incomes, const vec& consumptions, mat Tr, uword depth);
RcppExport SEXP _nSGMM_BBP_T_from_tYc_cpp(SEXP transferstructureSEXP, SEXP incomesSEXP, SEXP consumptionsSEXP, SEXP TrSEXP, SEXP depthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type transferstructure(transferstructureSEXP);
    Rcpp::traits::input_parameter< const vec& >::type incomes(incomesSEXP);
    Rcpp::traits::input_parameter< const vec& >::type consumptions(consumptionsSEXP);
    Rcpp::traits::input_parameter< mat >::type Tr(TrSEXP);
    Rcpp::traits::input_parameter< uword >::type depth(depthSEXP);
    rcpp_result_gen = Rcpp::wrap(BBP_T_from_tYc_cpp(transferstructure, incomes, consumptions, Tr, depth));
    return rcpp_result_gen;
END_RCPP
}
// BBP_T_from_atY_plain_cpp
mat BBP_T_from_atY_plain_cpp(const mat& alphas, const mat& transferstructure, const vec& incomes);
RcppExport SEXP _nSGMM_BBP_T_from_atY_plain_cpp(SEXP alphasSEXP, SEXP transferstructureSEXP, SEXP incomesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const mat& >::type transferstructure(transferstructureSEXP);
    Rcpp::traits::input_parameter< const vec& >::type incomes(incomesSEXP);
    rcpp_result_gen = Rcpp::wrap(BBP_T_from_atY_plain_cpp(alphas, transferstructure, incomes));
    return rcpp_result_gen;
END_RCPP
}
// equilibrate_cpp_fast8_debug
mat equilibrate_cpp_fast8_debug(const mat& altruism, const vec& income, const mat& capacity, bool verbose, int& r, bool& NE, int maxrounds);
RcppExport SEXP _nSGMM_equilibrate_cpp_fast8_debug(SEXP altruismSEXP, SEXP incomeSEXP, SEXP capacitySEXP, SEXP verboseSEXP, SEXP rSEXP, SEXP NESEXP, SEXP maxroundsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type altruism(altruismSEXP);
    Rcpp::traits::input_parameter< const vec& >::type income(incomeSEXP);
    Rcpp::traits::input_parameter< const mat& >::type capacity(capacitySEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int& >::type r(rSEXP);
    Rcpp::traits::input_parameter< bool& >::type NE(NESEXP);
    Rcpp::traits::input_parameter< int >::type maxrounds(maxroundsSEXP);
    rcpp_result_gen = Rcpp::wrap(equilibrate_cpp_fast8_debug(altruism, income, capacity, verbose, r, NE, maxrounds));
    return rcpp_result_gen;
END_RCPP
}
// equilibrate_cpp_fast8_smarter
mat equilibrate_cpp_fast8_smarter(const mat& altruism, const vec& income, const mat& capacity, int maxrounds);
RcppExport SEXP _nSGMM_equilibrate_cpp_fast8_smarter(SEXP altruismSEXP, SEXP incomeSEXP, SEXP capacitySEXP, SEXP maxroundsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type altruism(altruismSEXP);
    Rcpp::traits::input_parameter< const vec& >::type income(incomeSEXP);
    Rcpp::traits::input_parameter< const mat& >::type capacity(capacitySEXP);
    Rcpp::traits::input_parameter< int >::type maxrounds(maxroundsSEXP);
    rcpp_result_gen = Rcpp::wrap(equilibrate_cpp_fast8_smarter(altruism, income, capacity, maxrounds));
    return rcpp_result_gen;
END_RCPP
}
// recip_cpp
double recip_cpp(const mat& adj, const mat& radj);
RcppExport SEXP _nSGMM_recip_cpp(SEXP adjSEXP, SEXP radjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type adj(adjSEXP);
    Rcpp::traits::input_parameter< const mat& >::type radj(radjSEXP);
    rcpp_result_gen = Rcpp::wrap(recip_cpp(adj, radj));
    return rcpp_result_gen;
END_RCPP
}
// component_counts
double component_counts(const mat& transferstructure);
RcppExport SEXP _nSGMM_component_counts(SEXP transferstructureSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type transferstructure(transferstructureSEXP);
    rcpp_result_gen = Rcpp::wrap(component_counts(transferstructure));
    return rcpp_result_gen;
END_RCPP
}
// Ccomponents
double Ccomponents(const mat& transferstructure, uvec& t_components_csize, uvec& t_components_membership, mat& t_components_matrix, mat& t_conmat);
RcppExport SEXP _nSGMM_Ccomponents(SEXP transferstructureSEXP, SEXP t_components_csizeSEXP, SEXP t_components_membershipSEXP, SEXP t_components_matrixSEXP, SEXP t_conmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type transferstructure(transferstructureSEXP);
    Rcpp::traits::input_parameter< uvec& >::type t_components_csize(t_components_csizeSEXP);
    Rcpp::traits::input_parameter< uvec& >::type t_components_membership(t_components_membershipSEXP);
    Rcpp::traits::input_parameter< mat& >::type t_components_matrix(t_components_matrixSEXP);
    Rcpp::traits::input_parameter< mat& >::type t_conmat(t_conmatSEXP);
    rcpp_result_gen = Rcpp::wrap(Ccomponents(transferstructure, t_components_csize, t_components_membership, t_components_matrix, t_conmat));
    return rcpp_result_gen;
END_RCPP
}
// zeroOneBFS
mat zeroOneBFS(const mat& m, int src);
RcppExport SEXP _nSGMM_zeroOneBFS(SEXP mSEXP, SEXP srcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type src(srcSEXP);
    rcpp_result_gen = Rcpp::wrap(zeroOneBFS(m, src));
    return rcpp_result_gen;
END_RCPP
}
// BFS_dist_all
mat BFS_dist_all(mat graph);
RcppExport SEXP _nSGMM_BFS_dist_all(SEXP graphSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type graph(graphSEXP);
    rcpp_result_gen = Rcpp::wrap(BFS_dist_all(graph));
    return rcpp_result_gen;
END_RCPP
}
// dijkstra
vec dijkstra(mat graph, int src);
RcppExport SEXP _nSGMM_dijkstra(SEXP graphSEXP, SEXP srcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type graph(graphSEXP);
    Rcpp::traits::input_parameter< int >::type src(srcSEXP);
    rcpp_result_gen = Rcpp::wrap(dijkstra(graph, src));
    return rcpp_result_gen;
END_RCPP
}
// dijkstra_all
mat dijkstra_all(mat graph);
RcppExport SEXP _nSGMM_dijkstra_all(SEXP graphSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type graph(graphSEXP);
    rcpp_result_gen = Rcpp::wrap(dijkstra_all(graph));
    return rcpp_result_gen;
END_RCPP
}
// compute_moments_cpp
vec compute_moments_cpp(const mat& btransfers, const mat& kinship, const mat& distance);
RcppExport SEXP _nSGMM_compute_moments_cpp(SEXP btransfersSEXP, SEXP kinshipSEXP, SEXP distanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type btransfers(btransfersSEXP);
    Rcpp::traits::input_parameter< const mat& >::type kinship(kinshipSEXP);
    Rcpp::traits::input_parameter< const mat& >::type distance(distanceSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_moments_cpp(btransfers, kinship, distance));
    return rcpp_result_gen;
END_RCPP
}
// random_normal_seed
vec random_normal_seed(int n, double mean, double sd, int seed);
RcppExport SEXP _nSGMM_random_normal_seed(SEXP nSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(random_normal_seed(n, mean, sd, seed));
    return rcpp_result_gen;
END_RCPP
}
// symmetrix_normal_error_matrix
mat symmetrix_normal_error_matrix(int n, double mean, double sd, int seed);
RcppExport SEXP _nSGMM_symmetrix_normal_error_matrix(SEXP nSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(symmetrix_normal_error_matrix(n, mean, sd, seed));
    return rcpp_result_gen;
END_RCPP
}
// seedfromindex
int seedfromindex(int index);
RcppExport SEXP _nSGMM_seedfromindex(SEXP indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type index(indexSEXP);
    rcpp_result_gen = Rcpp::wrap(seedfromindex(index));
    return rcpp_result_gen;
END_RCPP
}
// simulate_BBP_cpp
mat simulate_BBP_cpp(int n, double delta0, double delta1, double delta2, double sigma, mat distance, mat kinship, mat capacity, vec income, int reps, int seed, int rounds);
RcppExport SEXP _nSGMM_simulate_BBP_cpp(SEXP nSEXP, SEXP delta0SEXP, SEXP delta1SEXP, SEXP delta2SEXP, SEXP sigmaSEXP, SEXP distanceSEXP, SEXP kinshipSEXP, SEXP capacitySEXP, SEXP incomeSEXP, SEXP repsSEXP, SEXP seedSEXP, SEXP roundsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type delta0(delta0SEXP);
    Rcpp::traits::input_parameter< double >::type delta1(delta1SEXP);
    Rcpp::traits::input_parameter< double >::type delta2(delta2SEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< mat >::type distance(distanceSEXP);
    Rcpp::traits::input_parameter< mat >::type kinship(kinshipSEXP);
    Rcpp::traits::input_parameter< mat >::type capacity(capacitySEXP);
    Rcpp::traits::input_parameter< vec >::type income(incomeSEXP);
    Rcpp::traits::input_parameter< int >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type rounds(roundsSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_BBP_cpp(n, delta0, delta1, delta2, sigma, distance, kinship, capacity, income, reps, seed, rounds));
    return rcpp_result_gen;
END_RCPP
}
// simulate_BBP_cpp_parallel
Rcpp::NumericMatrix simulate_BBP_cpp_parallel(int n, double delta0, double delta1, double delta2, double sigma, const mat& distance, const mat& kinship, const mat& capacity, const vec& income, int reps, int seed, int rounds);
RcppExport SEXP _nSGMM_simulate_BBP_cpp_parallel(SEXP nSEXP, SEXP delta0SEXP, SEXP delta1SEXP, SEXP delta2SEXP, SEXP sigmaSEXP, SEXP distanceSEXP, SEXP kinshipSEXP, SEXP capacitySEXP, SEXP incomeSEXP, SEXP repsSEXP, SEXP seedSEXP, SEXP roundsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type delta0(delta0SEXP);
    Rcpp::traits::input_parameter< double >::type delta1(delta1SEXP);
    Rcpp::traits::input_parameter< double >::type delta2(delta2SEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type distance(distanceSEXP);
    Rcpp::traits::input_parameter< const mat& >::type kinship(kinshipSEXP);
    Rcpp::traits::input_parameter< const mat& >::type capacity(capacitySEXP);
    Rcpp::traits::input_parameter< const vec& >::type income(incomeSEXP);
    Rcpp::traits::input_parameter< int >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type rounds(roundsSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_BBP_cpp_parallel(n, delta0, delta1, delta2, sigma, distance, kinship, capacity, income, reps, seed, rounds));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nSGMM_utlity_cppvec", (DL_FUNC) &_nSGMM_utlity_cppvec, 2},
    {"_nSGMM_mynegutility_cpp", (DL_FUNC) &_nSGMM_mynegutility_cpp, 5},
    {"_nSGMM_mynegutility_cpp_consumption", (DL_FUNC) &_nSGMM_mynegutility_cpp_consumption, 5},
    {"_nSGMM_get_BBP_BR_analytically_cpp_inline", (DL_FUNC) &_nSGMM_get_BBP_BR_analytically_cpp_inline, 5},
    {"_nSGMM_BBP_in_equilibrium_YaT_cpp", (DL_FUNC) &_nSGMM_BBP_in_equilibrium_YaT_cpp, 5},
    {"_nSGMM_meannanrm", (DL_FUNC) &_nSGMM_meannanrm, 2},
    {"_nSGMM_consumption_weights_cpp", (DL_FUNC) &_nSGMM_consumption_weights_cpp, 3},
    {"_nSGMM_BBP_c_from_atY_cpp", (DL_FUNC) &_nSGMM_BBP_c_from_atY_cpp, 6},
    {"_nSGMM_BBP_T_from_tYc_cpp", (DL_FUNC) &_nSGMM_BBP_T_from_tYc_cpp, 5},
    {"_nSGMM_BBP_T_from_atY_plain_cpp", (DL_FUNC) &_nSGMM_BBP_T_from_atY_plain_cpp, 3},
    {"_nSGMM_equilibrate_cpp_fast8_debug", (DL_FUNC) &_nSGMM_equilibrate_cpp_fast8_debug, 7},
    {"_nSGMM_equilibrate_cpp_fast8_smarter", (DL_FUNC) &_nSGMM_equilibrate_cpp_fast8_smarter, 4},
    {"_nSGMM_recip_cpp", (DL_FUNC) &_nSGMM_recip_cpp, 2},
    {"_nSGMM_component_counts", (DL_FUNC) &_nSGMM_component_counts, 1},
    {"_nSGMM_Ccomponents", (DL_FUNC) &_nSGMM_Ccomponents, 5},
    {"_nSGMM_zeroOneBFS", (DL_FUNC) &_nSGMM_zeroOneBFS, 2},
    {"_nSGMM_BFS_dist_all", (DL_FUNC) &_nSGMM_BFS_dist_all, 1},
    {"_nSGMM_dijkstra", (DL_FUNC) &_nSGMM_dijkstra, 2},
    {"_nSGMM_dijkstra_all", (DL_FUNC) &_nSGMM_dijkstra_all, 1},
    {"_nSGMM_compute_moments_cpp", (DL_FUNC) &_nSGMM_compute_moments_cpp, 3},
    {"_nSGMM_random_normal_seed", (DL_FUNC) &_nSGMM_random_normal_seed, 4},
    {"_nSGMM_symmetrix_normal_error_matrix", (DL_FUNC) &_nSGMM_symmetrix_normal_error_matrix, 4},
    {"_nSGMM_seedfromindex", (DL_FUNC) &_nSGMM_seedfromindex, 1},
    {"_nSGMM_simulate_BBP_cpp", (DL_FUNC) &_nSGMM_simulate_BBP_cpp, 12},
    {"_nSGMM_simulate_BBP_cpp_parallel", (DL_FUNC) &_nSGMM_simulate_BBP_cpp_parallel, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_nSGMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
