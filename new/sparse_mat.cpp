#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include "globals.hpp"
#include "sparse_mat.hpp"

float genConProb(int i, int j) {

  int pres_pop = which_pop[j];
  int post_pop = which_pop[i];

  float theta_i = (2.0 * M_PI * i) / (float) Na[post_pop];
  float theta_j = (2.0 * M_PI * j) / (float) Na[pres_pop];

  float proba = K / Na[pres_pop];

  if (PROBA == "cos")
    proba *= (1.0 + KAPPA[pres_pop + N_POP * post_pop] * cos(theta_i - theta_j));
  else if (PROBA == "spec")
    proba *= (1.0 + KAPPA[pres_pop + N_POP * post_pop] / sqrt(K) * cos(theta_i - theta_j));

  return proba;
}


void genSparseMatCSC(unsigned long*& colptr, int*& indices) {
    std::mt19937 rng;
    rng.seed(42);
    std::uniform_real_distribution<float> dist(0.0, 1.0);

    unsigned long nnz = 0;

    std::cout << " " << PROBA << std::endl;

    colptr[0] = 0;
    for (int j = 0; j < N; j++) { // presynaptic
      for (int i = 0; i < N; i++) { // postsynaptic

        if (dist(rng) < genConProb(i,j)) {
          indices[nnz] = i;
          nnz++;
        }
      }
      colptr[j+1] = nnz;
    }
}

void cscToDense(unsigned long* colptr, int* indices, int** dense) {
    // Initialize dense matrix with zeros
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            dense[i][j]=0;
        }
    }

    // Fill dense matrix with non-zero values from CSC format
    for (int j = 0; j < N; j++) {
        for (int i = colptr[j]; i < colptr[j+1]; i++) {
            dense[indices[i]][j] = 1 ;
        }
    }
}
