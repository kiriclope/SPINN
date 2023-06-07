#ifndef SPARSE_HPP
#define SPARSE_HPP

void getSparseMatCSC(size_t *&colptr, int *&indices);
void genSparseMatCSC(size_t*& colptr, int*& indices);
void saveSparseMatCSC(size_t* colptr, int* indices, size_t len);
void cscToDense(size_t* colptr, int* indices, int** dense);

#endif
