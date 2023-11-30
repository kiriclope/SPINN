#ifndef SPARSE_HPP
#define SPARSE_HPP

#include <vector>
#include <array>

float genConProb(int i, int j, std::vector<float> ksi_1, std::vector<float> ksi_2);

void getSparseMatCSC(size_t *&colptr, int *&indices);
void genSparseMatCSC(size_t*& colptr, int*& indices);
void saveSparseMatCSC(size_t* colptr, int* indices);
void cscToDense(size_t* colptr, int* indices, int** dense);

void generate_ksi(std::vector<float> &sample_0, std::vector<float> &sample_1, std::vector<float> &sample_2, int Nb);
void load_ksi(std::vector<float> &sample_0, std::vector<float> &sample_1, std::vector<float> &sample_2, int Nb);

std::vector<std::vector<float>> gen_mat_LR(std::vector<float> sample_1, std::vector<float> sample_2);
std::vector<std::vector<float>> outerProduct(std::vector<float>& vec);

std::array<float, 3> generateGaussian(std::vector<float> means, std::vector<float> stds, std::vector<float> correlation);
std::array<float, 3> generateBivariateGaussian(std::vector<float> means, std::vector<float> stds, std::vector<float> correlation);
std::array<float, 3> generate_trivariate_gaussian(std::vector<float> means, std::vector<float> stds, std::vector<float> rhos);

#endif
