/**
 * @file utils.h
 * @author Joonhyung Shin (tonyshin@kaist.ac.kr)
 * @author Namjo Ahn
 * @author Seungwoo Lee
 * @date May, 2018
 * @brief This is the header file for the utility functions.
 */

#ifndef UTILS_H
#define UTILS_H

#include <vector>

#include "euclidean.h"

namespace TCERC {

/**
 * @brief Calculate the reduced row echolon form of a matrix.
 * The matrix should be non-singular.
 *
 * @param vertices the vectors composing the non-singular matrix
 *
 * @return the reduced row echolon form of the matrix
 */
std::vector<Vector> gaussian_elimination(const std::vector<Vector> &vertices);

/**
 * @brief Solve the linear system.
 * The solution should be unique.
 *
 * @param matrix a matrix defining the linear system
 * @param coefficient the constant part of the linear system
 *
 * @return the solution to the linear system
 */
Vector linear_solver(const std::vector<Vector> &matrix, const Vector &coefficient);

Vector project_vector(const std::vector<Vector> &matrix, const Vector &vec);

/**
 * @brief Ensure that the given vectors are linearly independent.
 * If the vectors are not linearly independent, the test will freeze.
 *
 * @param vectors linearly independent vectors
 *
 * TODO
 */
void check_linear_independence(const std::vector<std::vector<iRRAM::REAL>> &vectors);

/**
 * @brief Compute the determinant of a matrix.
 *
 * @param matrix a square matrix
 *
 * @return the determinant of the matrix
 */
iRRAM::REAL det(const std::vector<std::vector<iRRAM::REAL>> &matrix);

}


#endif // UTILS_H
