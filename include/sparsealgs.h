#ifndef SPARSEALGS_H
#define SPARSEALGS_H

#include "sparsematrix.h"

// multiply CSR sparse matrix with dense vector
// and store result in input vector
template<int ROWS, int COLS, int NNZ, typename T> 
void spMV(const SparseCSR<ROWS, COLS, NNZ, T> &csr, const T inVector[COLS], T outVector[COLS]) {
	for (int i = 0; i < ROWS; i++) {
		T dot = 0;
		for (int j = csr.rowptr[i]; j < csr.rowptr[i + 1]; j++) {
			dot += csr.data[j] * inVector[csr.col[j]];
		}
		outVector[i] = dot;
	}
}

// multiply CSR sparse matrix with CSC sparse matrix
// based on "The Algorithms for FPGA Implementation of Sparse Matrices Multiplication"
template<int M, int K, int N, int NNZ1, int NNZ2, typename T>
void innerProductSpMM(const SparseCSR<M, K, NNZ1, T> &csr, const SparseCSC<K, N, NNZ2, T> &csc, T outMatrix[M][N]) {
	
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			
			T dot = 0;
			for (int ja = csr.rowptr[i], jb = csc.colptr[j]; ja < csr.rowptr[i + 1] && jb < csc.colptr[j + 1]; ) {
				if (csr.col[ja] < csc.row[jb]) {
					ja++;
				} else if (csr.col[ja] == csc.row[jb]) {
					dot += csr.data[ja] * csc.data[jb];
					ja++;
					jb++;
				} else {
					jb++;
				}
			}
			outMatrix[i][j] = dot;

		}
	}
}

#endif // SPARSEALGS_H
