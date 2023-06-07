#ifndef SPARSE_HPP
#define SPARSE_HPP

void getSparseMatCSC(unsigned long *&colptr, int *&indices);
void genSparseMatCSC(unsigned long*& colptr, int*& indices);
void saveSparseMatCSC(unsigned long* colptr, int* indices);
void cscToDense(unsigned long* colptr, int* indices, int** dense);

#endif
