#include <armadillo> // Include Armadillo linear algebra library
#include <cmath>
#include <iostream>
#include <string>
// g++ -I ~/mambaforge/include -L ~/mambaforge/lib multivariate.cpp -llapack -lblas -larmadillo

arma::mat generateProbabilityMatrix(std::string proba_str, size_t N, float K, float* p) {

    arma::vec theta(N), phi(N);
    for(size_t i = 0; i < N; ++i) {
        theta(i) = 2.0 * M_PI * i / N;
        phi(i) = 2.0 * M_PI * i / N;
    }

    arma::mat P(N,N);

    if(proba_str == "cosine")
      for(size_t i = 0; i < N; ++i) {
        for(size_t j = 0; j < N; ++j) {
          P(i,j) = K / N * (1.0 + p[j] * std::cos(theta(i) - phi(j)));
        }
      }
    else
      for(size_t i = 0; i < N; ++i) {
        for(size_t j = 0; j < N; ++j) {
          P(i,j) = K / N ;
        }
      }
    
    return P
}

arma::mat generateRandomMatrix(const arma::mat& P) {

    arma::mat C = arma::randu<arma::mat>(P.n_rows, P.n_cols); // Uniform random matrix
    C.transform( [&P] (double val, arma::uword i, arma::uword j) { return ( val <= P(i,j) ) ? 1.0 : 0.0; } ); // Transform based on threshold

    return C;
}

arma::mat generateMVN(const arma::vec& mu, const arma::mat& Sigma, int num_samples) {
    int num_vars = mu.n_elem; // Number of variables (Dimension of the Gaussian)
    arma::mat Y = arma::randn(num_samples, num_vars); // Independent standard normals
    arma::mat C; // Cholesky Decomposition
    arma::chol(C, Sigma); // Perform Cholesky Decomposition
    return arma::repmat(mu, 1, num_samples).t() + Y * C;
}

int main() {
    arma::vec mu = {0, 0, 0, 0}; // mean vector
    arma::mat Sigma = {{1, 0.5, 0.3, 0.1}, {0.5, 1, 0.2, 0.4}, {0.3, 0.2, 1, 0.7}, {0.1, 0.4, 0.7, 1}}; // covariance matrix
    int num_samples = 10000;

    arma::mat X = generateMVN(mu, Sigma, num_samples);
    // X.print();

    arma::mat sample_cov = arma::cov(X); // The samples are in X columns
    sample_cov.print("Sample Covariance:");

    double cov00 = arma::as_scalar(arma::cov(X.col(0), X.col(0)));
    double cov01 = arma::as_scalar(arma::cov(X.col(0), X.col(1)));
    double cov02 = arma::as_scalar(arma::cov(X.col(0), X.col(2)));
    double cov03 = arma::as_scalar(arma::cov(X.col(0), X.col(3)));

    std::cout << "Covariance of first and first columns: " << cov00 << std::endl;
    std::cout << "Covariance of first and second columns: " << cov01 << std::endl;
    std::cout << "Covariance of first and third columns: " << cov02 << std::endl;
    std::cout << "Covariance of first and fourth columns: " << cov03 << std::endl;

    // arma::rowvec x1 = X.row(0);
    // arma::rowvec x2 = X.row(1);
    // arma::rowvec x3 = X.row(2);
    // arma::rowvec x4 = X.row(3);

    return 0;
}
