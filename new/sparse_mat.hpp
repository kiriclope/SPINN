#ifndef SPARSE_HPP
#define SPARSE_HPP

void genSparseMatCSC(unsigned long*& colptr, int*& indices);
void cscToDense(unsigned long* colptr, int* indices, int** dense);

#endif
