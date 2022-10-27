// sparsematrix.h
#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <iostream>

// csr sparse matrix format
// ROWS: number of rows of original dense matrix
// COLS: number of columns of original dense matrix
// NNZ: number of non-zero elements
template<int ROWS, int COLS, int NNZ, typename T>
class SparseCSR {
public:
	// constructor sparsifies dense matrix
	SparseCSR(T denseMatrix[ROWS][COLS]) {
		// scan input matrix row-wise for csr
		rowptr[0] = 0;
		int countNNZ = 0;
		for (int i = 0; i < ROWS; i++) {
			for (int j = 0; j < COLS; j++) {
				if (denseMatrix[i][j] != 0.0) {
					data[countNNZ] = denseMatrix[i][j];
					col[countNNZ] = j;
					countNNZ += 1;
				}
			}
			rowptr[i + 1] = countNNZ;
		}
	}

	T data[NNZ];
	int col[NNZ];
	int rowptr[ROWS + 1];
};

// print csr matrix  using << operator
template<int ROWS, int COLS, int NNZ, typename T>
std::ostream &operator<<(std::ostream &os, const SparseCSR<ROWS, COLS, NNZ, T> &csr) {
	std::cout << "data = [ ";
	for (int i = 0; i < NNZ; i++) {
		std::cout << csr.data[i] << " ";
	}
	std::cout << "]\ncol = [ ";
	for (int i = 0; i < NNZ; i++) {
		std::cout << csr.col[i] << " ";
	}
	std::cout << "]\nrowptr = [ ";
	for (int i = 0; i < ROWS + 1; i++) {
		std::cout << csr.rowptr[i] << " ";
	}
	std::cout << "]\n";

	return os;
}

// CSC sparse matrix format
// ROWS: number of rows of original dense matrix
// COLS: number of columns of original dense matrix
// NNZ: number of non-zero elements
template<int ROWS, int COLS, int NNZ, typename T> 
class SparseCSC {
public:
	// constructor sparsifies dense matrix
	SparseCSC(T denseMatrix[ROWS][COLS]) {
		// scan input matrix column-wise for csc
		colptr[0] = 0;
		int countNNZ = 0;
		for (int j = 0; j < COLS; j++) {
			for (int i = 0; i < ROWS; i++) {
				if (denseMatrix[i][j] != 0.0) {
					data[countNNZ] = denseMatrix[i][j];
					row[countNNZ] = i;
					countNNZ += 1;
				}
			}
			colptr[j + 1] = countNNZ;
		}

	}

	T data[NNZ];
	int row[NNZ];
	int colptr[COLS + 1];
};

// print csc matrix  using << operator
template<int ROWS, int COLS, int NNZ, typename T>
std::ostream &operator<<(std::ostream &os, const SparseCSC<ROWS, COLS, NNZ, T> &csc) {
	std::cout << "data = [ ";
	for (int i = 0; i < NNZ; i++) {
		std::cout << csc.data[i] << " ";
	}
	std::cout << "]\nrow = [ ";
	for (int i = 0; i < NNZ; i++) {
		std::cout << csc.row[i] << " ";
	}
	std::cout << "]\ncolptr = [ ";
	for (int i = 0; i < ROWS + 1; i++) {
		std::cout << csc.colptr[i] << " ";
	}
	std::cout << "]\n";

	return os;
}

#endif //SPARSEMATRIX_H
