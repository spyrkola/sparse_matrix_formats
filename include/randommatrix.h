// randommatrix.h
#ifndef RANDOMMATRIX_H
#define RANDOMMATRIX_H

#include <cstdlib>
#include "sparsematrix.h"

// generate a random 2D array for coo, csr, csc formats
template<int K, int M, int LO, int HI, typename T>
void createSparseMat(int nnz, T mat[K][M]) {
	for (int j = 0; j < nnz;) {
		int index1 = (int) (K * ((double) rand() / (RAND_MAX + 1.0)));
		int index2 = (int) (M * ((double) rand() / (RAND_MAX + 1.0)));
		if (mat[index1][index2]) {  /* something already at this index */
			continue;  /*skip incrementing and try again*/
		}
		mat[index1][index2] = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
		++j;
	}
}

// generate a random 2D array for bsr format
template<int M, int N, int BLOCKSIZE, int LO, int HI, typename T>
void createSparseMatBlock(int nnzb, T mat[M][N]) {
	// randomly select the blocks that contain non-zero elements
	// then, randomly set some of the elements of these blocks
	bool isNonZero[M / BLOCKSIZE][N / BLOCKSIZE] = { false };
	for (int k = 0; k < nnzb;) {
		int index1 = (int) (M / BLOCKSIZE * ((double) rand() / (RAND_MAX + 1.0)));
		int index2 = (int) (N / BLOCKSIZE * ((double) rand() / (RAND_MAX + 1.0)));
		
		if (isNonZero[index1][index2]) {
			continue;
		}
		
		bool hasNonZero = false;
		do {
			for (int i = index1 * BLOCKSIZE; i < index1 * BLOCKSIZE + BLOCKSIZE; i++) {
				for (int j = index2 * BLOCKSIZE; j < index2 * BLOCKSIZE + BLOCKSIZE; j++) {
					if (rand() % 2 > 0) {
						mat[i][j] = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
						hasNonZero = true;
					}
				}
			}
		} while (!hasNonZero);
				
		isNonZero[index1][index2] = true;
		k += 1;
	}
}

// generate a random 2D array for ell format
template<int M, int N, int LO, int HI, typename T>
void createSparseMatEll(int maxnnzcols, T mat[M][N]) {
	// guarantee one row with the maximum number of non-zero elements
	// the rest of the rows can have any valid number of non-zero elements
	int maxedRow = (int) (M * ((double) rand() / (RAND_MAX + 1.0)));
	
	for (int i = 0; i < M; i++) {
		int numOfRowElements = (i == maxedRow) ? maxnnzcols : (rand() % (maxnnzcols + 1));
		for (int j = 0; j < numOfRowElements;) {
			int index = (int) (N * ((double) rand() / (RAND_MAX + 1.0)));
			if (mat[i][index] != 0) {
				continue;
			}
			
			mat[i][index] = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
			++j;
		}
	}
}

// generate a random 2D array for tjds format
template<int M, int N, int LO, int HI, typename T>
void createSparseMatTjds(int nnz, int tj_tiles, T mat[M][N]) {
	// guarantee one col with tj_tiles non-zero elements
	// the rest of the rows can have any valid num. of non-zero elements
	int maxedCol = (int) (N * ((double) rand() / (RAND_MAX + 1.0)));
	int remainingNNZ = nnz - tj_tiles;

	// first, deal with maxed column
	// then, create the rest of the matrix
	for (int i = 0; i < tj_tiles;) {
		int indexMaxed = (int) (M * ((double) rand() / (RAND_MAX + 1.0)));
		if (mat[indexMaxed][maxedCol] != 0) {
			continue;
		}
		
		mat[indexMaxed][maxedCol] = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
		++i;
	}
	
	int elementsPerCol[N] = { 0 };
	elementsPerCol[maxedCol] = tj_tiles;
	
	for (int i = 0; i < remainingNNZ;) {
		int index1 = (int) (M * ((double) rand() / (RAND_MAX + 1.0)));
		int index2 = (int) (N * ((double) rand() / (RAND_MAX + 1.0)));
		
		if (elementsPerCol[index2] == tj_tiles || mat[index1][index2] != 0) {
			continue;
		}

		mat[index1][index2] = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
		elementsPerCol[index2] += 1;
		i += 1;
	}
}

#endif // RANDOMMATRIX_H
