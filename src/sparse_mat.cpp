#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <vector>
#include <utility>
#include <array>
#include <numeric> // for std::inner_product

#include "globals.hpp"
#include "utils.hpp"
#include "sparse_mat.hpp"
// #include "multivariate.hpp"

std::random_device rd_mat;
std::mt19937 rng(rd_mat());

// LOW RANK
std::random_device rd_lr;
std::mt19937 rng_lr(rd_lr());

std::array<float, 3> generateTrivariateGaussian(std::vector<float>means,
                                                std::vector<float> stds,
                                                std::vector<float> rhos) {
  
  // Generate three independent standard normal variables
  float Z1 = white(rng_lr);
  float Z2 = white(rng_lr);
  float Z3 = white(rng_lr);
  
  // Extract correlation coefficients for readability
  float rho_21 = rhos[0];
  float rho_31 = rhos[1];
  float rho_32 = rhos[2];
  
  // float deter = (1.0 - rho_yz*rho_yz) - rho_xy * (rho_xy - rho_xz*rho_yz) + rho_xz * (rho_xy*rho_yz - rho_xz);
  // std::cout << deter << std::endl;
    
  // Cholesky decomposition of the covariance matrix
  float L11 = 1.0;
  float L21 = rho_21;
  float L22 = sqrt(1.0 - L21 * L21);
  // float L31 = rho_xz / L11;
  float L31 = rho_31;
  // float L32 = (rho_32 - L31 * L21) / L22 ;
  float L32 = rho_32;
  float L33 = sqrt(1.0 - L31 * L31 - L32 * L32);
  
  // Generating correlated variables
  float X1 = L11 * Z1 + means[0];
  float X2 = L21 * Z1 + L22 * Z2 + means[1];
  float X3 = L31 * Z1 + L32 * Z2 + L33 * Z3 + means[2];
  
  return {X3, X1, X2};
}

std::array<float, 3> generateBivariateGaussian(std::vector<float> means, std::vector<float> stds, std::vector<float> correlation) {  
 
  float Z1 = white(rng_lr);
  float Z2 = white(rng_lr);
  
  float X1 = stds[0] * Z1 + means[0];
  float X2 = correlation[0] * stds[0] * Z1 + stds[1] * sqrt(1.0 - correlation[0] * correlation[0]) * Z2 + means[1];

  return {X1, X2, X2};
}

std::array<float, 3> generateGaussian(std::vector<float> means, std::vector<float> stds) {
    
  float Z1 = white(rng_lr);
  float X1 = stds[0] * Z1 + means[0];
    
  return {X1, X1, X1};
}

std::vector<std::vector<float>> outerProduct(std::vector<float>& vec) {
  size_t length = vec.size();
  std::vector<std::vector<float>> result(length, std::vector<float>(length));

  for (size_t i = 0; i < length; i++) {
    for (size_t j = 0; j < length; j++) {
      result[i][j] = vec[i] * vec[j];
    }
  }

  return result;
}

void generate_ksi(std::vector<float> &sample_0, std::vector<float> &sample_1, std::vector<float> &sample_2, int Nb) {

  std::cout << "Generating LR Vec ";
  
  if(LR_SEED!=0){
    std::cout << " LR_SEED " << LR_SEED ;
    rng_lr.seed(LR_SEED);
  }
  
  std::array<float, 3> samples;
  
  if(LR_RANK == 1)
    for (int i = 0; i < Nb; i++) {
      samples = generateGaussian(LR_MEAN, LR_STD);
      sample_0.push_back(samples[0]);
      sample_1.push_back(samples[0]);
    }
  
  if(LR_RANK == 2)
    for (int i = 0; i < Nb; i++) {
      samples = generateBivariateGaussian(LR_MEAN, LR_STD, LR_RHO);
      sample_1.push_back(samples[0]); // ksi_1
      sample_2.push_back(samples[1]); // ksi_2
    }
  
  if(LR_RANK==3) {    
    for (int i = 0; i < Nb; i++) {
      samples = generateTrivariateGaussian(LR_MEAN, LR_STD, LR_RHO);
      sample_0.push_back(samples[0]); // Input
      sample_1.push_back(samples[1]); // ksi_1
      sample_2.push_back(samples[2]); // ksi_2
    }
  }
  std::cout <<" ksi ";
  for (int i=0; i<5; i++)
    std::cout << sample_1[i] << " ";
  std::cout << std::endl;

  std::cout << "Saving LR Vec to:" << MAT_PATH;
  
  std::ofstream ksi0File(MAT_PATH + "/ksi_0.txt");
  saveVectorToFile(ksi0File, sample_0);
  ksi0File.close();
  
  std::ofstream ksi1File(MAT_PATH + "/ksi_1.txt");
  saveVectorToFile(ksi1File, sample_1);
  ksi1File.close();
  
  std::ofstream ksi2File(MAT_PATH + "/ksi_2.txt");
  saveVectorToFile(ksi2File, sample_2);
  ksi2File.close();

  std::cout << " Done" << std::endl;  
}


std::vector<std::vector<float>> gen_mat_LR(std::vector<float> sample_1, std::vector<float> sample_2) {
  std::vector<std::vector<float>> outer = outerProduct(sample_1);  
  return outer;
}

float genConProb(int i, int j) {

  int pres_pop = which_pop[j];
  int post_pop = which_pop[i];
  
  float theta_i = (2.0 * M_PI * (i - cNa[post_pop])) / (float) Na[post_pop];
  float theta_j = (2.0 * M_PI * (j - cNa[pres_pop])) / (float) Na[pres_pop];
  
  float proba = Ka[pres_pop] / Na[pres_pop];

  if (PROBA[pres_pop + N_POP * post_pop] == "cos")
    proba *= (1.0 + KAPPA[pres_pop + N_POP * post_pop] * cos(theta_i - theta_j));

  else if (PROBA[pres_pop + N_POP * post_pop] == "spec")
    proba *= (1.0 + KAPPA[pres_pop + N_POP * post_pop] / sqrt(Ka[pres_pop]) * cos(theta_i - theta_j));

  else if (PROBA[pres_pop + N_POP * post_pop] == "lr") {
    proba *= (1.0 + KAPPA[pres_pop + N_POP * post_pop]
              * (ksi_1[i - cNa[post_pop]] * ksi_1[j - cNa[pres_pop]]
                 + ksi_2[i - cNa[post_pop]] * ksi_2[j - cNa[pres_pop]])
              / sqrt(Ka[pres_pop]));
  }
  
  if (proba>1.0)
    proba = 1;
  if (proba<0.0)
    proba = 0;
  
  return proba;
}

void genSparseMatCSC(size_t*& colptr, int*& indices) {
  
  std::cout << "Generating Sparse Matrix" ;
  
  // arma::mat ksi_mat;
  // std::vector<float> ksi_1;
  // std::vector<float> ksi_2;
  // std::vector<std::vector<float>> outer;
  
  // if (PROBA[0] == "lr") {    
  //   // outer = gen_mat_LR(ksi_1, ksi_2);
  //   // ksi = generate_ksi(Na[0]);
  //   // ksi_mat = get_mat_LR(ksi, 2, Na[0]);
  // }
  
  size_t nnz = 0;
  
  for(int i=0;i<N_POP;i++)
    for(int j=0;j<N_POP;j++)
      std::cout << " " << PROBA[i + N_POP * j] ;
  std::cout << std::endl;
  
  colptr[0] = 0;
  for (int j = 0; j < N; j++) { // presynaptic
    for (int i = 0; i < N; i++) { // postsynaptic
      // std::cout << i << " " << j << std::endl;
      if (unif(rng) < genConProb(i, j)) {
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

void load_ksi(std::vector<float> &sample_0, std::vector<float> &sample_1, std::vector<float> &sample_2) {
  std::cout << "Loading LR vec from:" << MAT_PATH;
  
  std::ifstream ksi0File(MAT_PATH + "/ksi_0.txt");
  loadVectorFromFile(ksi0File, sample_0);
  ksi0File.close();
  
  std::ifstream ksi1File(MAT_PATH + "/ksi_1.txt");
  loadVectorFromFile(ksi1File, sample_1);
  ksi1File.close();
  
  std::ifstream ksi2File(MAT_PATH + "/ksi_2.txt");
  loadVectorFromFile(ksi2File, sample_2);
  ksi2File.close();

  std::cout << std::endl;
  std::cout <<" ksi ";  
  for(int i=0; i<5; i++)
    std::cout << sample_1[i] << " ";
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
