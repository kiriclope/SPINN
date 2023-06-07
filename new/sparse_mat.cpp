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

  float theta_i = (2.0 * M_PI * i) / (float) Na[post_pop];
  float theta_j = (2.0 * M_PI * j) / (float) Na[pres_pop];

  float proba = Ka[pres_pop] / Na[pres_pop];

  if (PROBA == "cos")
    proba *= (1.0 + KAPPA[pres_pop + N_POP * post_pop] * cos(theta_i - theta_j));
  else if (PROBA == "spec")
    proba *= (1.0 + KAPPA[pres_pop + N_POP * post_pop] / sqrt(Ka[pres_pop]) * cos(theta_i - theta_j));

  return proba;
}

void getSparseMatCSC(size_t *&colptr, int *&indices) {

  std::cout << "Loading Sparse Matrix" ;

  std::ifstream inFile("./matrix/colptr.txt");
  loadArrayFromFile(inFile, colptr, N+1);
  inFile.close();

  std::ifstream inFile2("./matrix/indices.txt");
  loadArrayFromFile(inFile2, indices, colptr[N+1]);
  inFile2.close();

  std::cout << " Done" << std::endl;

}

void saveSparseMatCSC(size_t* colptr, int* indices){

  std::cout << "Saving Sparse Matrix" ;

  std::ofstream colptrFile("./matrix/colptr.txt");
  saveArrayToFile(colptrFile, colptr, (size_t) N+1);
  colptrFile.close();

  std::ofstream indicesFile("./matrix/indices.txt");
  saveArrayToFile(indicesFile, indices, colptr[N+1]);
  indicesFile.close();

  std::cout << " Done" << std::endl;

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
