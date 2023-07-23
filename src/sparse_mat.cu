#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cmath>
#include <string>
#include "globals.hpp"

// nvcc -arch=sm_75 -o output_file source_file.cu

__device__ float genConProb(int i, int j) {

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

__global__ void genSparseMatCSC(unsigned long* colptr, int* indices) {

  int j = blockIdx.x * blockDim.x + threadIdx.x;

  // Each CUDA thread generates connections for multiple pre-synaptic neurons
  for (int k = threadIdx.y; k < N; k += blockDim.y) {
    int i = k;
    float P = genConProb(i, j, which_pop, Na, Ka, N_POP, KAPPA, PROBA);
    int nnz = 0;
    if (P > 0) {
      nnz = atomicAdd(&colptr[j+1], 1);
      indices[colptr[j] + nnz] = i;
    }

    // Generate connections for the rest of the post-synaptic neurons
    for (i = k+blockDim.y; i < N; i += blockDim.y) {
      P = genConProb(i, j, which_pop, Na, Ka, N_POP, KAPPA, PROBA);
      if (P > 0) {
        nnz = atomicAdd(&colptr[j+1], 1);
        indices[colptr[j] + nnz] = i;
      }
    }
  }

  // Write NNZ to colptr[0] once all the threads have completed
  __syncthreads();

  if (threadIdx.y == 0 && threadIdx.x == 0) {
    unsigned long total_nnz = colptr[blockDim.x * blockDim.y];
    colptr[0] = 0;
    for (int i = 1; i <= N; i++) {
      unsigned long tmp = colptr[i];
      colptr[i] = colptr[i-1] + total_nnz;
      total_nnz = tmp;
    }
  }
}

int main() {
  
  std::string configname = argv[1] ;
  loadConfig(configname);
  
  dim3 blockSize(32, 8);  // 256 threads
  dim3 gridSize((N + blockSize.x - 1) / blockSize.x, 1);  // 1D grid
  
  genSparseMatCSC<<<gridSize, blockSize>>>(colptr, indices);
  
  cudaDeviceSynchronize();  // Wait for kernel to finish
  
}
