#include <iostream>
#include <cstdlib>

#include "sparsematrix.h"
#include "sparsealgs.h"
#include "randommatrix.h"

// prints the contents of a 2D array with M rows and N columns
template<int M, int N, typename T>
void print2Darray(T arr[M][N]) {
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			std::cout << arr[i][j] << ' ';
		}
		std::cout << '\n';
	}
}

// prints the contents of an array of length M
template<int M, typename T>
void print1Darray(T arr[M]) {
	for (int i = 0; i < M; i++) {
		std::cout << arr[i] << ' ';
	}
	std::cout << '\n';
}

int main() {
	// ---Sparse Matrix Formats---
	/*
	// known example matrix
	double denseData[6][6] = {
		{ 1.0, 0.0, 0.0, 2.0, 0.0, 0.0 }, 
		{ 3.0, 4.0, 0.0, 5.0, 0.0, 0.0 }, 
		{ 6.0, 0.0, 7.0, 8.0, 9.0, 0.0 }, 
		{ 0.0, 0.0, 10.0, 11.0, 0.0, 0.0 }, 
		{ 0.0, 0.0, 0.0, 0.0, 12.0, 0.0 }, 
		{ 0.0, 0.0, 0.0, 0.0, 13.0, 14.0 }
	};
	
	std::cout << "---Dense Matrix---\n";
	print2Darray<6, 6, double>(denseData);

	// create and print all 6 sparse formats
	SparseCOO<6, 6, 14, double> cooMatrix(denseData);
	std::cout << "---COO Matrix---\n" << cooMatrix;

	SparseCSR<6, 6, 14, double> csrMatrix(denseData);
	std::cout << "---CSR Matrix---\n" << csrMatrix;
	
	SparseCSC<6, 6, 14, double> cscMatrix(denseData);
	std::cout << "---CSC Matrix---\n" << cscMatrix;

	SparseBSR<6, 6, 2, 6, double> bsr(denseData);
	std::cout << "---BSR Matrix---\n" << bsr;

	SparseELL<6, 6, 4, double> ell(denseData);
	std::cout << "---ELL Matrix---\n" << ell;
	
	SparseTJDS<6, 6, 14, 4, double> tjds(denseData);
	std::cout << "---TJDS Matrix---\n" << tjds;
	*/	

	// 1) coo, csr, csc formats using createSparseMat()
	srand(time(NULL));
	int randomDense[10][7] = {};
	createSparseMat<10, 7, 1, 9, int>(11, randomDense);

	std::cout << "\n---Random Dense Matrix---\n";
	print2Darray<10, 7, int>(randomDense);
	
	SparseCOO<10, 7, 11, int> cooRandom(randomDense);
	std::cout << "---COO Matrix---\n" << cooRandom;

	SparseCSR<10, 7, 11, int> csrRandom(randomDense);
	std::cout << "---CSR Matrix---\n" << csrRandom;
	
	SparseCSC<10, 7, 11, int> cscRandom(randomDense);
	std::cout << "---CSC Matrix---\n" << cscRandom;
	
	// 2) bsr format using createSparseMatBlock()
	int randomBlockDense[6][18] = {};
	createSparseMatBlock<6, 18, 2, 1, 9, int>(7, randomBlockDense);
	
	std::cout << "\n---Random Block Dense Matrix---\n";
	print2Darray<6, 18, int>(randomBlockDense);

	SparseBSR<6, 18, 2, 7, int> bsrRandom(randomBlockDense);
	std::cout << "---BSR Matrix---\n" << bsrRandom;
	
	// 3) ell format using createSparseMatEll()
	int randomEllDense[5][11] = {};
	createSparseMatEll<5, 11, 1, 9, int>(7, randomEllDense);
	
	std::cout << "\n---Random Dense Matrix (for ELL)---\n";
	print2Darray<5, 11>(randomEllDense);

	SparseELL<5, 11, 7, int> ellRandom(randomEllDense);
	std::cout << "---ELL Matrix---\n" << ellRandom;
	
	// 4) tjds format using createSparseMatTjds()
	int randomTjdsDense[5][4] = {};
	createSparseMatTjds<5, 4, 1, 9, int>(10, 4, randomTjdsDense);
	
	std::cout << "\n---Random Dense Matrix (for TJDS)---\n";
	print2Darray<5, 4>(randomTjdsDense);

	SparseTJDS<5, 4, 10, 4, int> tjdsRandom(randomTjdsDense);
	std::cout << "---TJDS Matrix---\n" << tjdsRandom;


	// ---Algorithms---
	/*
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
	*/
	
	return 0;
}
