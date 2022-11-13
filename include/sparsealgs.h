#ifndef SPARSEALGS_H
#define SPARSEALGS_H

#include "sparsematrix.h"

// multiply CSR sparse matrix with dense vector and store results in outVector
template<int ROWS, int COLS, int NNZ, typename T> 
void spMV(const SparseCSR<ROWS, COLS, NNZ, T> &csr, const T inVector[COLS], T outVector[ROWS]) {
	for (int i = 0; i < ROWS; i++) {
		T dot = 0;
		for (int j = csr.rowptr[i]; j < csr.rowptr[i + 1]; j++) {
			dot += csr.data[j] * inVector[csr.col[j]];
		}
		outVector[i] = dot;
	}
}

// multiply ELL sparse matrix with dense vector and store results in outVector
template<int ROWS, int COLS, int MAXNNZCOLS, typename T>
void spMV(const SparseELL<ROWS, COLS, MAXNNZCOLS, T> &ell, const T inVector[COLS], T outVector[ROWS]) {
	for (int i = 0; i < ROWS; i++) {
		T dot = 0;
		for (int j = 0; j < MAXNNZCOLS && ell.data[i * MAXNNZCOLS + j] != 0; j++) {
			dot += ell.data[i * MAXNNZCOLS + j] * inVector[ell.col[i * MAXNNZCOLS + j]];
		}
		outVector[i] = dot;
	}	
}

// multiply TJDS sparse matrix with dense vector and store results in outVector
// NOTE: assumes that inVector has been reordered as part of sparse matrix encoding
template<int ROWS, int COLS, int NNZ, int TJ_TILES, typename T>
void spMV(const SparseTJDS<ROWS, COLS, NNZ, TJ_TILES, T> &tjds, const T inVector[COLS], T outVector[ROWS]) {
	
	for (int i = 0; i < ROWS; i++) {
		outVector[i] = 0;
	}

	for (int i = 0; i < TJ_TILES; i++) {
		int vecIdx = 0;
		for (int j = tjds.start[i]; j < tjds.start[i + 1]; j++) {
			outVector[tjds.row_index[j]] += tjds.val[j] * inVector[vecIdx];
			vecIdx += 1;
		}
	}	
}

// multiply SSS sparse matrix with dense vector
// based on "Improving the Performance of the Symmetric Sparse Matrix-Vector Multiplication in Multicore (Alg. 2)"
template<int N, int LOWERNNZ, typename T>
void spMV(const SparseSSS<N, LOWERNNZ, T> &sss, const T inVector[N], T outVector[N]) {
	for (int r = 0; r < N; r++) {
		outVector[r] = sss.dvalues[r] * inVector[r];
		for (int j = sss.rowptr[r]; j < sss.rowptr[r + 1]; j++) {
			int c = sss.col[j];
			outVector[r] += sss.values[j] * inVector[c];
			outVector[c] += sss.values[j] * inVector[r];
		}
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

// outer product dataflow for SpMM, using CSC and CSR sparse matrices
// outer product dataflow eliminates index-matching
template<int M, int K, int N, int NNZ1, int NNZ2, typename T>
void outerProductSpMM(const SparseCSC<M, K, NNZ1, T> &csc, const SparseCSR<K, N, NNZ2, T> &csr, T outMatrix[M][N]) {
	for (int k = 0; k < K; k++) {
		
		for (int i = csc.colptr[k]; i < csc.colptr[k + 1]; i++) {
			for (int j = csr.rowptr[k]; j < csr.rowptr[k + 1]; j++) {
				outMatrix[csc.row[i]][csr.col[j]] += csc.data[i] * csr.data[j];
			}
		}
		
	}
}

// gustavson dataflow for SpMM, using csr inputs
template<int M, int K, int N, int NNZ1, int NNZ2, typename T>
void gustavsonProductSpMM(const SparseCSR<M, K, NNZ1, T> &a, const SparseCSR<K, N, NNZ2, T> &b, T outMatrix[M][N]) {
	
	for (int i = 0; i < M; i++) {
		for (int k = a.rowptr[i]; k < a.rowptr[i + 1]; k++) {
			for (int j = b.rowptr[a.col[k]]; j < b.rowptr[a.col[k] + 1]; j++) {
				outMatrix[i][b.col[j]] += a.data[k] * b.data[j];
			}
		}
	}
	
}

// column-wise dataflow for SpMM, using csc inputs
template<int M, int K, int N, int NNZ1, int NNZ2, typename T> 
void columnWiseProductSpMM(const SparseCSC<M, K, NNZ1, T> &a, const SparseCSC<K, N, NNZ2, T> &b, T outMatrix[M][N]) {
	
	for (int j = 0; j < N; j++) {
		for (int k = b.colptr[j]; k < b.colptr[j + 1]; k++) {
			for (int i = a.colptr[b.row[k]]; i < a.colptr[b.row[k] + 1]; i++) {
				outMatrix[a.row[i]][j] += a.data[i] * b.data[k];
			}
		}
	}
	
}

#endif // SPARSEALGS_H
