#include <iostream>

#include "sparsematrix.h"
#include "sparsealgs.h"

// a function that prints the contents of a 2D array
// with M rows and N columns
template<int M, int N, typename T>
void print2Darray(T arr[M][N]) {
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			std::cout << arr[i][j] << ' ';
		}
		std::cout << '\n';
	}
}

template<int M, typename T>
void print1Darray(T arr[M]) {
	for (int i = 0; i < M; i++) {
		std::cout << arr[i] << ' ';
	}
	std::cout << '\n';
}

int main() {
	// ---Sparse Matrix Vector Multiplication---
	// some dense matrix
	double denseData[5][5] = {
		{ 1.0, 0, 0, 2.0, 0 }, 
		{ 3.0, 4.0, 0, 5.0, 0 }, 
		{ 6.0, 0, 7.0, 8.0, 9.0 }, 
		{0, 0, 10.0, 11.0, 0 }, 
		{ 0, 0, 0, 0, 12.0 }
	};

	std::cout << "---Dense Matrix---\n";
	print2Darray<5, 5, double>(denseData);
	
	// create csr matrix
	SparseCSR<5, 5, 12, double> csrMatrix(denseData);
	std::cout << "\n---CSR Matrix---\n";
	std::cout << csrMatrix;
	
	// spMV matrix multiplication with csr sparse matrix
	double inVector[5] = { 1, 2, 3, 4, 5 };
	double outVector[5];
	
	spMV(csrMatrix, inVector, outVector);
	std::cout << "\n---Sparse Matrix Vector Multiplication (CSR)---\n";
	print1Darray<5, double>(outVector);

	// ---Sparse Matrix Matrix Multiplication (inner product, csr-csc)---
	double dense1[3][3] = {
		{ 1.0, 0.0, 1.2 },
		{ 0.0, 0.0, 2.0 }, 
		{ 0.0, 3.0, 0.0 }
	};

	double dense2[3][3] = {
		{ 0.0, 2.0, 3.2 }, 
		{ 1.0, 0.0, 0.0 }, 
		{ 0.0, 2.2, 0.0 }
	};

	// sparsify the above dense matrices
	SparseCSR<3, 3, 4, double> csr(dense1);
	SparseCSC<3, 3, 4, double> csc(dense2);
	
	double denseout[3][3];
	innerProductSpMM(csr, csc, denseout);
	std::cout << "\n ---inner product SpMM---\n";
	std::cout << "result: \n";
	print2Darray<3, 3, double>(denseout);

	return 0;
}
