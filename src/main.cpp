#include <iostream>

#include "matrix.h"
#include "sparsemat.h"

int main() {
	
	double data[] = {
		1.0, 0, 0, 2.0, 0, 
		3.0, 4.0, 0, 5.0, 0, 
		6.0, 0, 7.0, 8.0, 9.0, 
		0, 0, 10.0, 11.0, 0, 
		0, 0, 0, 0, 12.0
	};

	Matrix mat(5, 5, data);
	std::cout << mat;	

	SparseCSR csrmat(mat);
	std::cout << csrmat;

	return 0;
}
