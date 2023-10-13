#include <iostream>
#include <fstream>

#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>

#include "globals.hpp"
#include "utils.hpp"

__global__ void initCurand(curandState* states, unsigned long seed)
{
  size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  curand_init(seed, idx, 0, &states[idx]);
}

__device__ float genConProb(int i, int j, int Na, int Nb, float Kb, float kappa) {
  
  float theta_i = (2.0 * M_PI * i) / (float) Na;
  float theta_j = (2.0 * M_PI * j) / (float) Nb;
  
  float proba = Kb / (float) Nb;
  proba *= (1.0 + kappa * cos(theta_i - theta_j));
  
  return proba;
}

__global__ void genSparseMatKernel(unsigned long long* nnz, size_t* colptr, int* indices, curandState* states, int Na, int Nb, float Kb, float kappa) {
  
  int j = threadIdx.x + blockIdx.x * blockDim.x;
  
  if (j < Nb) {
    float rand_num;
    float proba;
    unsigned long long local_nnz = 0;
    
    colptr[j] = *nnz;
    for (int i = 0; i < Na; i++) { // postsynaptic
      rand_num = curand_uniform(&states[j]);
      
      proba = genConProb(i, j, Na, Nb, Kb, kappa);
      
      if (rand_num < proba) {
        indices[*nnz + local_nnz] = i;
        local_nnz++;
      }
    }
    colptr[j+1] = local_nnz + *nnz ;
    atomicAdd(nnz, local_nnz);
  }
}

void genSparseMatCSC(size_t*& colptr, int*& indices, int Na, int Nb, float Kb, float kappa, size_t offset) {
  
  std::cout << "Generating Sparse Matrix" ;
  
  // Allocate memory on the device
  unsigned long long* nnz=0, nnzz=0;
  cudaMallocManaged(&nnz, sizeof(*nnz));
  *nnz = 0; // important: initialize since atomicAdd is an incrementing operation
  
  size_t *d_colptr;
  int *d_indices;
  curandState *d_states;  
  
  for(int i_pop=0; i_pop<2; ++i_pop) {
    cudaMalloc(&d_colptr, (Nb+1) * sizeof(size_t));
    cudaMalloc(&d_indices, Na * Nb * sizeof(int));
    
    // Allocate space for curand states
    cudaMalloc(&d_states, Nb * sizeof(curandState));
  
    // Initialize curand
    initCurand<<<Nb, 1>>>(d_states, i_pop*4);

    nnzz = *nnz;
    // Invoke the kernel    
    genSparseMatKernel<<<Nb, 1>>>(nnz, d_colptr, d_indices, d_states, Na, Nb, Kb, kappa);
    cudaDeviceSynchronize(); // make sure all device operations finish before accessing data on the host
    
    std::cout << " nnz " << *nnz << std::endl; 
    // cudaFree(nnz);
    
    // Copy result back to host
    cudaMemcpy(&colptr[offset], d_colptr, (Nb+1) * sizeof(size_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&indices[nnzz], d_indices, Na * Nb * sizeof(int), cudaMemcpyDeviceToHost);
    
    // Cleanup
    cudaFree(d_colptr);
    cudaFree(d_indices);
    cudaFree(d_states);
  }
  
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

int main(int argc , char** argv) {
  
  loadConfig("../conf/config_EI.yml");
  std::cout << " Done" << std::endl;
  
  size_t* colptr = new size_t[N+1]();
  int* indices = new int[(size_t) (N * N)]();

  for(int i_pop=0; i_pop<N_POP;++i_pop)
    for(int j_pop=0; j_pop<N_POP;++j_pop) {
      genSparseMatCSC(colptr, indices, Na[i_pop], Na[j_pop], Ka[j_pop], KAPPA[j_pop + N_POP * i_pop], cNa[j_pop]);
      
      std::cout << colptr[0] << " " << colptr[1] << std::endl;
    }
  
  if(IF_SAVE_MAT)
    saveSparseMatCSC(colptr, indices);
  
  // Remember to free dynamically allocated memory after you're done
  delete[] colptr;
  delete[] indices;
  
  return 0;
}
