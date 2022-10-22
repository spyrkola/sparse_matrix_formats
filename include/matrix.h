/// matrix.h
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <stdexcept>

/// basic implementation of a matrix class.
class Matrix {
public: 
	// constructor copies values from array data[]
	// data[] contains 2D array elements stored in row-major order
	Matrix(int rows, int cols, double data[]) {
		if (rows <= 0 || cols <= 0) {
			throw std::invalid_argument("Matrix dimensions must be greater than zero.");
		}

		this->rows = rows;
		this->cols = cols;
		
		this->data = new double[rows * cols];		
		for (int i = 0; i < rows * cols; i++) {
			this->data[i] = data[i];	
		}
	}

	Matrix(const Matrix &mat) {
		rows = mat.getRows();
		cols = mat.getCols();

		data = new double[rows * cols];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				data[i * cols + j] = mat(i, j);
			}
		}
	}

	~Matrix() {
		delete[] data;
	}

	// operators
	// get matrix elements using (i, j) notation
	double operator()(int i, int j) const {
		if (i < 0 || i >= rows || j < 0 || j >= cols) {
			throw std::invalid_argument("Index out of bounds.");
		}		
		return data[i * cols + j];
	}
	
	// for printing matrix using cout
	friend std::ostream &operator<<(std::ostream &os, const Matrix &mat);

	// getters
	int getRows() const {
		return rows;
	}

	int getCols() const {
		return cols;
	}
	
private: 
	int rows;
	int cols;
	double *data;
};

std::ostream &operator<<(std::ostream &os, const Matrix &mat) {
	for (int i = 0; i < mat.rows; i++) {
		for (int j = 0; j < mat.cols; j++) {
			os << mat.data[i * mat.cols + j] << " ";
		}
		os << "\n";
	}
	return os;
}

#endif // MATRIX_H
