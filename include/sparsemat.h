/// sparse_matrix_formats.h
#ifndef SPARSEMAT_H
#define SPARSEMAT_H

#include <iostream>
#include <vector>

#include "matrix.h"

/// a basic implementation of the csr sparse matrix format
class SparseCSR {
public:
	// constructor creates sparse csr matrix given a regular matrix m
	SparseCSR(const Matrix &m) {
		nnz = 0;
		matrixRows = m.getRows();
		rowptr = new int[m.getRows() + 1];
	
		// sparsify matrix m
		// using vectors avoids looping through m twice
        // rowptr array size is known so no vector needed
		std::vector<double> dataVector;
		std::vector<int> colVector;

		rowptr[0] = 0;
		for (int i = 0; i < m.getRows(); i++) {
			for (int j = 0; j < m.getCols(); j++) {
				if (m(i, j) != 0.0) {
					dataVector.push_back(m(i, j));
					colVector.push_back(j);
					nnz += 1;
				}
			}
			rowptr[i+1] = nnz;
		}
		
		// copy elements from vectors to corresponding arrays (?)
		data = new double[nnz];
		std::copy(dataVector.begin(), dataVector.end(), data);

		col = new int[nnz];
		std::copy(colVector.begin(), colVector.end(), col);
	}
		
	
	~SparseCSR() {
		delete[] data;
		delete[] col;
		delete[] rowptr;
	}
	
	// operators
	// for printing csr matrix using cout
	friend std::ostream &operator<<(std::ostream &os, const SparseCSR &csr);

	// getters
	int getNnz() const {
		return nnz;
	}

	int getMatrixRows() const {
		return matrixRows;
	}

private:
	double *data;
	int *col;
	int *rowptr;

	int nnz;
	int matrixRows;
};

std::ostream &operator<<(std::ostream &os, const SparseCSR &csr) {
	std::cout << "data = [ ";
	for (int i = 0; i < csr.nnz; i++) {
		std::cout << csr.data[i] << " ";
	}
	std::cout << "]\ncols = [ ";
	for (int i = 0; i < csr.nnz; i++) {
		std::cout << csr.col[i] << " ";
	}
	std::cout << "]\nrowptr = [ ";
	for (int i = 0; i < csr.matrixRows + 1; i++) {
		std::cout << csr.rowptr[i] << " ";
	}
	std::cout << "]\n";

	return os;
}

#endif //SPARSEMAT_H
