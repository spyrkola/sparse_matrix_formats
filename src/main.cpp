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
	srand(time(NULL));
	
	// ---Sparse Matrix Vector Multiplication Algorithms
	std::cout << "======Matrix Vector Multiplication======\n";
	int vecIn[9] = {}; int vecOut[6] = {};
	for (int i = 0; i < 9; i++) {
		vecIn[i] = 1 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (8)));
	}
	std::cout << "Dense Vector: \n";
	print1Darray<9, int>(vecIn);
	
	// 0) coo SpMV
	std::cout << "---COO SpMV---\n";
	int denseCooMat[6][9] = {};
	createSparseMat<6, 9, 1, 9, int>(10, denseCooMat);
	std::cout << "matrix:\n";
	print2Darray<6, 9, int>(denseCooMat);

	SparseCOO<6, 9, 10, int> cooMat(denseCooMat);
	spMV(cooMat, vecIn, vecOut);
	std::cout << "result:\n";
	print1Darray<6, int>(vecOut);
	
	// 1) csr SpMV
	std::cout << "---CSR SpMV---\n";
	int denseMat[6][9] = {};
	createSparseMat<6, 9, 1, 9, int>(7, denseMat);
	std::cout << "matrix:\n";
	print2Darray<6, 9, int>(denseMat);
	
	SparseCSR<6, 9, 7, int> csrMat(denseMat);
	spMV(csrMat, vecIn, vecOut);
	std::cout << "result:\n";
	print1Darray<6, int>(vecOut);
	
	// 2) bsr SpMV
	std::cout << "---BSR SpMV---\n";
	int denseBsrMat[6][9] = {};
	createSparseMatBlock<6, 9, 3, 1, 9, int>(4, denseBsrMat);
	std::cout << "matrix:\n";
	print2Darray<6, 9, int>(denseBsrMat);
	
	SparseBSR<6, 9, 3, 4, int> bsrMat(denseBsrMat);
	spMV(bsrMat, vecIn, vecOut);
	std::cout << "result: \n";
	print1Darray<6, int>(vecOut);
	
	// 3) ell SpMV
	std::cout << "---ELL SpMV---\n";
	int denseEllMat[6][9] = {};
	createSparseMatEll<6, 9, 1, 9, int>(5, denseEllMat);
	std::cout << "matrix:\n";
	print2Darray<6, 9, int>(denseEllMat);

	SparseELL<6, 9, 5, int> ellMat(denseEllMat);
	spMV(ellMat, vecIn, vecOut);
	std::cout << "result:\n";
	print1Darray<6, int>(vecOut);

	// 4) tjds SpMV
	// create copy of original input vector for reordering
	int vecCpy[9];
	for (int i = 0; i < 9; i++) {
		vecCpy[i] = vecIn[i];
	}

	std::cout << "---TJDS SpMV---\n";
	int denseTjdsMat[6][9] = {};
	createSparseMatTjds<6, 9, 1, 9, int>(11, 3, denseTjdsMat);
	std::cout << "matrix:\n";
	print2Darray<6, 9, int>(denseTjdsMat);


	SparseTJDS<6, 9, 11, 3, int> tjdsMat(denseTjdsMat, vecCpy);
	spMV(tjdsMat, vecCpy, vecOut);
	std::cout << "results:\n";
	print1Darray<6, int>(vecOut);

	// 5) sss SpMV
	int vecOutSym[9];
	std::cout << "---SSS SpMV---\n";
	int denseSssMat[9][9] = {};
	createSparseMatSss<9, 1, 9, int>(7, denseSssMat);
	std::cout << "matrix:\n";
	print2Darray<9, 9, int>(denseSssMat);

	SparseSSS<9, 7, int> sssMat(denseSssMat);
	spMV(sssMat, vecIn, vecOutSym);
	std::cout << "result:\n";
	print1Darray<9, int>(vecOutSym);
	
	// ---Sparse Matrix Matrix Multiplication Algorithms---	
	std::cout << "\n======Matrix Matrix Multiplication======\n";
	double dense1[5][8] = {};
	createSparseMat<5, 8, 1, 9, double>(15, dense1);

	double dense2[8][4] = {};
	createSparseMat<8, 4, 1, 9, double>(7, dense2);

	std::cout << "matrix 1:\n";
	print2Darray<5, 8, double>(dense1);
	std::cout << "matrix 2:\n";
	print2Darray<8, 4, double>(dense2);
	
	// sparsify dense matrices
	SparseCSR<5, 8, 15, double> aCSR(dense1);
	SparseCSC<5, 8, 15, double> aCSC(dense1);

	SparseCSR<8, 4, 7, double> bCSR(dense2);
	SparseCSC<8, 4, 7, double> bCSC(dense2);
	
	// 1) inner product dataflow
	double denseOutInner[5][4] = { 0.0 };
	innerProductSpMM(aCSR, bCSC, denseOutInner);
	std::cout << "\n ---inner product SpMM---\n";
	print2Darray<5, 4, double>(denseOutInner);
	
	// 2) outer product dataflow
	double denseOutOuter[5][4] = { 0.0 };
	outerProductSpMM(aCSC, bCSR, denseOutOuter);
	std::cout << "---outer product SpMM---\n";
	print2Darray<5, 4, double>(denseOutOuter);
	
	// 3) gustavson product dataflow
	double denseOutGustavson[5][4] = { 0.0 };
	gustavsonProductSpMM(aCSR, bCSR, denseOutGustavson);
	std::cout << "---gustavson product SpMM---\n";
	print2Darray<5, 4, double>(denseOutGustavson);
	
	// 4) column-wise product dataflow
	double denseOutColumnWise[5][4] = { 0.0 };
	columnWiseProductSpMM(aCSC, bCSC, denseOutColumnWise);
	std::cout << "---column-wise product SpMM---\n";
	print2Darray<5, 4, double>(denseOutColumnWise);
	
	return 0;
}
