#define ARMA_DONT_PRINT_FAST_MATH_WARNING
#include <armadillo> // Include Armadillo linear algebra library
#include <cmath>
#include <iostream>
#include <string>

arma::mat generateProbabilityMatrix(std::string proba_str, size_t N, float K, float* p) ;

arma::mat generateMVN(const arma::vec& mu, const arma::mat& Sigma, int num_samples) ;

arma::mat generate_ksi(int Nb) ;

arma::mat get_mat_LR(arma::mat ksi, int rank, int Nb) ;
