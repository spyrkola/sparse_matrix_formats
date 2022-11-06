// sparsematrix.h
#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <iostream>

// COO sparse matrix format
// ROWS: number of rows of original dense matrix
// COLS: number of columns of original dense matrix
// NNZ: number of non-zero elements
template<int ROWS, int COLS, int NNZ, typename T>
class SparseCOO {
public:
	// constructor sparsifies dense matrix by
	// scanning input matrix and storing non-zero elements
	// along with their corresponding rows and cols
	SparseCOO(T denseMatrix[ROWS][COLS]) {
		int countNNZ = 0;
		for (int i = 0; i < ROWS; i++) {
			for (int j = 0; j < COLS; j++) {
				if (denseMatrix[i][j] != 0.0) {
					data[countNNZ] = denseMatrix[i][j];
					row[countNNZ] = i;
					col[countNNZ] = j;
					countNNZ += 1;
				}
			}
		}
	}

	T data[NNZ];
	int row[NNZ];
	int col[NNZ];
};

// CSR sparse matrix format
// ROWS: number of rows of original dense matrix
// COLS: number of columns of original dense matrix
// NNZ: number of non-zero elements
template<int ROWS, int COLS, int NNZ, typename T>
class SparseCSR {
public:
	// constructor sparsifies dense matrix
	// scan input matrix row-wise for csr
	SparseCSR(T denseMatrix[ROWS][COLS]) {
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

// CSC sparse matrix format
// ROWS: number of rows of original dense matrix
// COLS: number of columns of original dense matrix
// NNZ: number of non-zero elements
template<int ROWS, int COLS, int NNZ, typename T> 
class SparseCSC {
public:
	// constructor sparsifies dense matrix
	// scan input matrix column-wise for csc
	SparseCSC(T denseMatrix[ROWS][COLS]) {
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

// BSR sparse matrix format
// IMPORTANT NOTE: assumes that matrix dimensions are multiples of BLOCKSIZE
// ROWS: num. of rows of original dense matrix
// COLS: num. of cols of original dense matrix
// BLOCKSIZE: dimension of blocks
// NNZBLOCKS: number of blocks that contain nonzero elements
template<int ROWS, int COLS, int BLOCKSIZE, int NNZBLOCKS, typename T>
class SparseBSR {
public: 
	SparseBSR(T denseMatrix[ROWS][COLS]) {
		blockRowptr[0] = 0;
		int countNNZblocks = 0;
		// move along dense matrix block by block
		// for each block, find if nonzero elements exist
		// if they do, store block vals. in data (in row-major ordrer), block column in blockCol and increment counter
		for (int i = 0; i < ROWS; i += BLOCKSIZE) {
			for (int j = 0; j < COLS; j += BLOCKSIZE) {

				bool foundNonZero = false;
				for (int block_i = i; (block_i < i + BLOCKSIZE) && (!foundNonZero); block_i++)  {
					for (int block_j = j; (block_j < j + BLOCKSIZE) && (!foundNonZero); block_j++) {
						if (denseMatrix[block_i][block_j] != 0.0) {
							foundNonZero = true;
						}
					}
				}
				
				if (foundNonZero) {
					blockCol[countNNZblocks] = j / BLOCKSIZE;
					for (int block_i = i; block_i < i + BLOCKSIZE; block_i++) {
						for (int block_j = j; block_j < j + BLOCKSIZE; block_j++) {
							int idx = countNNZblocks * BLOCKSIZE * BLOCKSIZE + (block_i % BLOCKSIZE) * BLOCKSIZE + block_j % BLOCKSIZE;
							data[idx] = denseMatrix[block_i][block_j];
						}
					}
					countNNZblocks += 1;
				}				
			}

			blockRowptr[i / BLOCKSIZE + 1] = countNNZblocks;
		}
	}
	
	int blockRowptr[ROWS  / BLOCKSIZE + 1];
	int blockCol[NNZBLOCKS];
	T data[NNZBLOCKS * BLOCKSIZE * BLOCKSIZE];
};

// ELL sparse matrix format
// ROWS: num. of rows of original dense matrix
// COLS: num. of cols of original dense matrix
// MAXNNZCOLS: the maximum number of non-zero elements in one row
template<int ROWS, int COLS, int MAXNNZCOLS, typename T>
class SparseELL {
public:
	SparseELL(T denseMatrix[ROWS][COLS]) {
		// scan dense matrix row-wise
		// and add non-zero elements to the correct spot in data, col arrays
		// if necessary, pad data, col arrays with 0s and -1s respectively
		for (int i = 0; i < ROWS; i++) {
			int currentCol = 0;
			for (int j = 0; j < COLS; j++) {
				if (denseMatrix[i][j] != 0.0) {
					data[i * MAXNNZCOLS + currentCol] = denseMatrix[i][j];
					col[i * MAXNNZCOLS + currentCol] = j;
					currentCol += 1;
				}
			}

			while (currentCol < MAXNNZCOLS) {
				data[i * MAXNNZCOLS + currentCol] = 0;
				col[i * MAXNNZCOLS + currentCol] = -1;
				currentCol += 1;
			}
		}	
	}
	
	T data[ROWS * MAXNNZCOLS];
	int col[ROWS * MAXNNZCOLS];
};

// TJDS sparse matrix format
// based on "Systolic Sparse Matrix Vector Multiply in the Age of TPUs and Accelerators"
// ROWS: num. of rows of original dense matrix
// COLS: num. of cols of original dense matrix
// NNZ: num. of non-zero elements
// TJ_TILES: equal to the number of coefficients in the most populated column
template<int ROWS, int COLS, int NNZ, int TJ_TILES, typename T>
class SparseTJDS {
public:
	SparseTJDS(T denseMatrix[ROWS][COLS]) {
		// create matrix of row indices (use -1 instead of 0 because C++ is 0-based)
		int Ids[ROWS][COLS];
		for (int i = 0; i < ROWS; i++) {
			for (int j = 0; j < COLS; j++) {
				Ids[i][j] = (denseMatrix[i][j] != 0) ? i : -1;
			}
		}
		
		// compress denseMatrix and Ids arrays column-wise
		for (int j = 0; j < COLS; j++) {
			int currentRow = 0;
			for (int i = 0; i < ROWS; i++) {
				if (denseMatrix[i][j] != 0.0 && currentRow == i) {
					currentRow += 1;
				} else if (denseMatrix[i][j] != 0.0 && currentRow < i) {
					denseMatrix[currentRow][j] = denseMatrix[i][j];
					denseMatrix[i][j] = 0.0;
					
					Ids[currentRow][j] = Ids[i][j];
					Ids[i][j] = -1;

					currentRow += 1;
				}
			}
		}
				
		// reorder columns based on number of non-zero elements per column
		// do this by scanning from bottom to top and from left to right
		// if you encounter nonzero element in a column, place that column
		// left and increment sortedCols counter
		T denseTemp;
		int IdsTemp;
		int sortedCols = 0;
		for (int i = ROWS - 1; i >= 0; i--) {
			for (int j = sortedCols; j < COLS && sortedCols < COLS - 1; j++) {
				if (denseMatrix[i][j] != 0.0) {
					for (int k = 0; k < ROWS; k++) {
						denseTemp = denseMatrix[k][j];
						denseMatrix[k][j] = denseMatrix[k][sortedCols];
						denseMatrix[k][sortedCols] = denseTemp;

						IdsTemp = Ids[k][j];
						Ids[k][j] = Ids[k][sortedCols];
						Ids[k][sortedCols] = IdsTemp;	
					}
					sortedCols += 1;
				}
			}
		}
		
		// finally, store the above arrays as rows
		start[0] = 0;
		int nonZeroCount = 0;
		for (int i = 0; i < TJ_TILES; i++) {
			for (int j = 0; j < COLS; j++) {
				if (denseMatrix[i][j] != 0.0) {
					val[nonZeroCount] = denseMatrix[i][j];
					row_index[nonZeroCount] =  Ids[i][j];
					nonZeroCount += 1;
				}
			}
			start[i + 1] = nonZeroCount;
		}

	}
	
	T val[NNZ];
	int row_index[NNZ];
	int start[TJ_TILES + 1];
};

// some << operator overloads
// print coo matrix  using << operator
template<int ROWS, int COLS, int NNZ, typename T>
std::ostream &operator<<(std::ostream &os, const SparseCOO<ROWS, COLS, NNZ, T> &coo) {
	std::cout << "data = [ ";
	for (int i = 0; i < NNZ; i++) {
		std::cout << coo.data[i] << " ";
	}
	std::cout << "]\nrow = [ ";
	for (int i = 0; i < NNZ; i++) {
		std::cout << coo.row[i] << " ";
	}
	std::cout << "]\ncol = [ ";
	for (int i = 0; i < NNZ; i++) {
		std::cout << coo.col[i] << " ";
	}
	std::cout << "]\n";

	return os;
}

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

// print bsr matrix  using << operator
template<int ROWS, int COLS, int BLOCKSIZE, int NNZBLOCKS, typename T>
std::ostream &operator<<(std::ostream &os, const SparseBSR<ROWS, COLS, BLOCKSIZE, NNZBLOCKS, T> &bsr) {
	std::cout << "data = [ ";
	for (int i = 0; i < NNZBLOCKS * BLOCKSIZE * BLOCKSIZE; i++) {
		std::cout << bsr.data[i] << " ";
	}
	std::cout << "]\nblockRowptr = [ ";
	for (int i = 0; i < ROWS / BLOCKSIZE + 1; i++) {
		std::cout << bsr.blockRowptr[i] << " ";
	}
	std::cout << "]\nblockCol = [ ";
	for (int i = 0; i < NNZBLOCKS; i++) {
		std::cout << bsr.blockCol[i] << " ";
	}
	std::cout << "]\n";

	return os;
}

// print ell matrix  using << operator
template<int ROWS, int COLS, int MAXNNZCOLS, typename T>
std::ostream &operator<<(std::ostream &os, const SparseELL<ROWS, COLS, MAXNNZCOLS, T> &ell) {
	std::cout << "data = [ ";
	for (int i = 0; i < ROWS * MAXNNZCOLS; i++) {
		std::cout << ell.data[i] << " ";
	}
	std::cout << "]\ncol = [ ";
	for (int i = 0; i < ROWS * MAXNNZCOLS; i++) {
		std::cout << ell.col[i] << " ";
	}
	std::cout << "]\n";

	return os;
}

// print tjds matrix  using << operator
template<int ROWS, int COLS, int NNZ, int TJ_TILES, typename T>
std::ostream &operator<<(std::ostream &os, const SparseTJDS<ROWS, COLS, NNZ, TJ_TILES, T> &tjds) {
	std::cout << "val = [ ";
	for (int i = 0; i < NNZ; i++) {
		std::cout << tjds.val[i] << " ";
	}
	std::cout << "]\nrow_index = [ ";
	for (int i = 0; i < NNZ; i++) {
		std::cout << tjds.row_index[i] << " ";
	}
	std::cout << "]\nstart = [ ";
	for (int i = 0; i < TJ_TILES + 1; i++) {
		std::cout << tjds.start[i] << " ";
	}
	std::cout << "]\n";

	return os;
}

#endif //SPARSEMATRIX_H
