#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>

#include "globals.hpp"
#include "utils.hpp"
#include "sparse_mat.hpp"

float genConProb(int i, int j) {

  int pres_pop = which_pop[j];
  int post_pop = which_pop[i];

  float theta_i = (2.0 * M_PI * (i - cNa[post_pop])) / (float) Na[post_pop];
  float theta_j = (2.0 * M_PI * (j - cNa[pres_pop])) / (float) Na[pres_pop];

  float proba = Ka[pres_pop] / Na[pres_pop];

  if (PROBA == "cos")
    proba *= (1.0 + KAPPA[pres_pop + N_POP * post_pop] * cos(theta_i - theta_j));
  else if (PROBA == "spec")
    proba *= (1.0 + KAPPA[pres_pop + N_POP * post_pop] / sqrt(Ka[pres_pop]) * cos(theta_i - theta_j));

  return proba;
}

void genSparseMatCSC(size_t*& colptr, int*& indices) {

  std::cout << "Generating Sparse Matrix" ;

  std::mt19937 rng;
  rng.seed(42);
  std::uniform_real_distribution<float> unif(0.0, 1.0);

  size_t nnz = 0;

  std::cout << " " << PROBA << std::endl;

  colptr[0] = 0;
  for (int j = 0; j < N; j++) { // presynaptic
    for (int i = 0; i < N; i++) { // postsynaptic

      if (unif(rng) < genConProb(i,j)) {
        indices[nnz] = i;
        nnz++;
      }
    }
    colptr[j+1] = nnz;
  }

  if (IF_SAVE_MAT)
    saveSparseMatCSC(colptr, indices);

  std::cout << " Done" << std::endl;
}

void getSparseMatCSC(size_t *&colptr, int *&indices) {

  std::cout << "Loading Sparse Matrix from:" << MAT_PATH;

  std::ifstream colptrFile(MAT_PATH + "/colptr.txt");
  loadArrayFromFile(colptrFile, colptr, (size_t) N+1);
  colptrFile.close();

  std::ifstream indicesFile(MAT_PATH + "/indices.txt");
  loadArrayFromFile(indicesFile, indices, colptr[N]);
  indicesFile.close();

  std::cout << " Done" << std::endl;
}

void saveSparseMatCSC(size_t* colptr, int* indices){

  std::cout << "Saving Sparse Matrix to: " << MAT_PATH ;
  ensureDirExists(MAT_PATH);

  std::ofstream colptrFile(MAT_PATH + "/colptr.txt");
  saveArrayToFile(colptrFile, colptr, (size_t) N+1);
  colptrFile.close();

  std::ofstream indicesFile(MAT_PATH + "/indices.txt");
  saveArrayToFile(indicesFile, indices, colptr[N]);
  indicesFile.close();

  std::cout << " Done" << std::endl;
}

void cscToDense(size_t* colptr, int* indices, int** dense) {
    // Initialize dense matrix with zeros
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            dense[i][j]=0;
        }
    }

    // Fill dense matrix with non-zero values from CSC format
    for (int j = 0; j < N; j++) {
        for (size_t i = colptr[j]; i < colptr[j+1]; i++) {
            dense[indices[i]][j] = 1 ;
        }
    }
}
