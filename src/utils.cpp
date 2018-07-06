#include "utils.h"

#include <stdexcept>

namespace TCERC {

std::vector<Vector> gaussian_elimination(const std::vector<Vector> &vertices)
{
	int n = vertices.size();
	int m;

	if (n <= 0)
		throw std::invalid_argument("zero dimension is not allowed");

	m = vertices[0].ambient_dimension();

	if (m <= 0)
		throw std::invalid_argument("zero dimension is not allowed");

	// Transpose when n < m
	std::vector<Vector> matrix;
	bool transposed = n < m;
	transposed = false;
	if (transposed) {
		for (int j = 0; j < m; j++) {
			std::vector<iRRAM::REAL> row = std::vector<iRRAM::REAL>();
			for (int i = 0; i < n; i++) {
				row.push_back(vertices[i][j]);
			}
			Vector v(row);
			matrix.push_back(v);
		}
		// Swap dimension
		std::swap(n, m);
	} else {
		for (int i = 0; i < n; i++) {
			matrix.push_back(vertices[i]);
		}
	}

	for (int j = 0; j < m; j++) {
		if (j == n) { // m > n, so stop here
			//printf("# of columns are larger than # of rows\n");
			break;
		}
		// int nonzero_row = choose( one of { (*(matrix[i].get_coordinates()))[j] != 0 } ) - 1;
		// TODO : use choose with many many arguments
		int nonzero_row = 1;
		std::vector<iRRAM::LAZY_BOOLEAN> nonzero_lists;
		//printf("Find nonzero row\n--\n");
		
		for (int i = j; i < n; i++) {
			//iRRAM::cout << "[" << matrix[i][j] << "]" << std::endl;
			iRRAM::LAZY_BOOLEAN nonzero = (matrix[i][j] > 0) || (matrix[i][j] < 0);
			nonzero_lists.push_back(nonzero);
		}
		//printf("--\n");

		if (j == n - 1)
			nonzero_row = n;
		else {
			nonzero_row = choose(nonzero_lists);
			nonzero_row += j;
		}

		//printf("Found : %d\n", nonzero_row);
		if (nonzero_row == 0) {
			throw std::invalid_argument("gaussian elimination error");
			// Zero column
			//continue;
			// Or return empty matrix here? for linear independency check
		//	printf("Nonzero 0?\n");
			//return std::vector<Vector>();
			// return std::vector<Vector>();
		}

		nonzero_row--; // choose function returns index starting from 1

		//printf("Swap rows[%dth, %dth]!\n", j, nonzero_row);
		Vector temp = matrix[nonzero_row];
		matrix[nonzero_row] = matrix[j];
		matrix[j] = temp;
		nonzero_row = j;

		for (int i = 0; i < n; i++) {
		//	printf("Vector algebra, ");
			
			if (i == nonzero_row) {
				temp = matrix[i] / matrix[i][j]; // Normalize
			} else {
				temp = matrix[i] - (matrix[nonzero_row] * (matrix[i][j] / matrix[nonzero_row][j]));
			}
		//	printf("Completed\n");
			matrix[i] = temp;
		}
	}

	// Transpose back
	std::vector<Vector> final_matrix;
	if (transposed) {
		for (int j = 0; j < m; j++) {
			std::vector<iRRAM::REAL> row;
			for (int i = 0; i < n; i++) {
				row.push_back(matrix[i][j]);
			}
			Vector v(row);
			final_matrix.push_back(v);
		}
		// Really need to swap back?
		std::swap(n, m);
	} else {
		final_matrix = matrix;
	}

	return final_matrix;
}

Vector linear_solver(const std::vector<Vector> &matrix, const Vector &coefficient)
{
	int d = matrix.size();
	int D = matrix[0].ambient_dimension();
        if (d != coefficient.ambient_dimension())
                throw std::invalid_argument("matrix row should be equal to vector's dimension");
	std::vector<Vector> joint_matrix;
	for (int i=0; i<d; i++){
		std::vector<iRRAM::REAL> row = std::vector<iRRAM::REAL>();
		for (int j=0; j<D; j++){
			row.push_back(matrix[i][j]);
		}
		row.push_back(coefficient[i]);
		Vector v(row);
		joint_matrix.push_back(v);
	}
	std::vector<Vector> answer_matrix = gaussian_elimination(joint_matrix);
	std::vector<iRRAM::REAL> root = std::vector<iRRAM::REAL>();
	for (int i=0; i<d; i++){
		root.push_back(answer_matrix[i][D]);
	}
	return Vector(root);
}

Vector project_vector(const std::vector<Vector> &matrix, const Vector &vec)
{
	std::vector<Vector> coeff_matrix;
	int size = matrix.size();
	for (int i = 0; i < size; i++) {
		std::vector<iRRAM::REAL> row;
		for (int j = 0; j < size; j++) {
			row.push_back(matrix[i] * matrix[j]);
		}
		coeff_matrix.push_back(Vector(row));
	}
	return linear_solver(coeff_matrix, vec);
}


/* TODO */
void check_linear_independence(const std::vector<std::vector<iRRAM::REAL>> &vectors)
{
	std::vector<Vector> matrix={};
	for (int i=0; i < vectors.size(); i++){
		Vector vec = Vector(vectors[i]);
		matrix.push_back(vec);
	}
	gaussian_elimination(matrix);
}

iRRAM::REAL det(const std::vector<std::vector<iRRAM::REAL>> &mat)
{
	iRRAM::REAL res = 1;
	int n = mat.size();
	int m;

	if (n <= 0)
		throw std::invalid_argument("zero dimension is not allowed");

	m = mat[0].size();

	if (m != n)
		throw std::invalid_argument("matrix should be square");

	// Transpose when n < m
	std::vector<Vector> matrix;
	bool transposed = n < m;
	transposed = false;
	if (transposed) {
		for (int j = 0; j < m; j++) {
			std::vector<iRRAM::REAL> row = std::vector<iRRAM::REAL>();
			for (int i = 0; i < n; i++) {
				row.push_back(mat[i][j]);
			}
			Vector v(row);
			matrix.push_back(v);
		}
		// Swap dimension
		std::swap(n, m);
	} else {
		for (int i = 0; i < n; i++) {
			matrix.push_back(mat[i]);
		}
	}

	for (int j = 0; j < m; j++) {
		if (j == n) { // m > n, so stop here
			//printf("# of columns are larger than # of rows\n");
			break;
		}
		// int nonzero_row = choose( one of { (*(matrix[i].get_coordinates()))[j] != 0 } ) - 1;
		// TODO : use choose with many many arguments
		int nonzero_row = 1;
		std::vector<iRRAM::LAZY_BOOLEAN> nonzero_lists;
		//printf("Find nonzero row\n--\n");
		
		for (int i = j; i < n; i++) {
			//iRRAM::cout << "[" << matrix[i][j] << "]" << std::endl;
			iRRAM::LAZY_BOOLEAN nonzero = (matrix[i][j] > 0) || (matrix[i][j] < 0);
			nonzero_lists.push_back(nonzero);
		}
		//printf("--\n");

		if (j == n - 1)
			nonzero_row = n;
		else {
			nonzero_row = choose(nonzero_lists);
			nonzero_row += j;
		}

		//printf("Found : %d\n", nonzero_row);
		if (nonzero_row == 0) {
			throw std::invalid_argument("gaussian elimination error");
			// Zero column
			//continue;
			// Or return empty matrix here? for linear independency check
		//	printf("Nonzero 0?\n");
			//return std::vector<Vector>();
			// return std::vector<Vector>();
		}

		nonzero_row--; // choose function returns index starting from 1

		//printf("Swap rows[%dth, %dth]!\n", j, nonzero_row);
		if (j != nonzero_row)
			res = -res;
		Vector temp = matrix[nonzero_row];
		matrix[nonzero_row] = matrix[j];
		matrix[j] = temp;
		nonzero_row = j;

		for (int i = 0; i < n; i++) {
		//	printf("Vector algebra, ");
			
			if (i == nonzero_row) {
				res *= matrix[i][j];
				temp = matrix[i] / matrix[i][j]; // Normalize
			} else {
				temp = matrix[i] - (matrix[nonzero_row] * (matrix[i][j] / matrix[nonzero_row][j]));
			}
		//	printf("Completed\n");
			matrix[i] = temp;
		}
	}

	// Transpose back
	std::vector<Vector> final_matrix;
	if (transposed) {
		for (int j = 0; j < m; j++) {
			std::vector<iRRAM::REAL> row;
			for (int i = 0; i < n; i++) {
				row.push_back(matrix[i][j]);
			}
			Vector v(row);
			final_matrix.push_back(v);
		}
		// Really need to swap back?
		std::swap(n, m);
	} else {
		final_matrix = matrix;
	}

	return res;
}


}
