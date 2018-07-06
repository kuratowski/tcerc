/**
 * @file tensor.h
 * @author Joonhyung Shin (tonyshin@kaist.ac.kr)
 * @author Namjo Ahn
 * @author Seungwoo Lee
 * @date May, 2018
 * @brief This is the header file for Tensor.
 */

#ifndef TENSOR_H
#define TENSOR_H

#include <vector>

#include <iRRAM/lib.h>
#include "subspace.h"

namespace TCERC {

/**
 * @brief An algebraic tensor lying in a subspace.
 */
class Tensor {

private:
	/* TODO: Change the implementation of tensor. */
	std::vector<bool> pq;
	int p, q;
	std::vector<iRRAM::REAL> data;

	std::vector<std::vector<iRRAM::REAL>> Amatrix;

	Subspace subspace;

public:
	/**
	 * @brief Construct a tensor given the coefficients.
	 *
	 * @param subspace the subspace that the tensor lies in
	 * @param p the number of contravariant parts
	 * @param q the number of covariant parts
	 * @param tensor the coefficients of the tensor
	 */
	Tensor(const Subspace &subspace, int p, int q, const std::vector<iRRAM::REAL> &tensor);

	/**
	 * @brief Construct a tensor given the coefficients.
	 *
	 * @param subspace the subspace that the tensor lies in
	 * @param type the type of the tensor
	 * @param tensor the coefficients of the tensor
	 */
	Tensor(const Subspace &subspace, std::vector<bool> type, const std::vector<iRRAM::REAL> &tensor);

	/**
	 * @brief Construct a tensor given a real matrix.
	 *
	 * @param subspace the subspace that the tensor lies in
	 * @param matrix the real matrix of the coefficients
	 */
	Tensor(const Subspace &subspace, const iRRAM::REALMATRIX &matrix);

	/**
	 * @brief Check whether the given index of the tensor is covariant.
	 *
	 * @param k the index of the tensor
	 *
	 * @return true if kth index is covariant
	 */
	bool is_covariant(int k) const;

	/**
	 * @brief Get the type of the tensor.
	 *
	 * @return the type of the tensor
	 */
	const std::vector<bool> &get_type() const;

	/* Basic arithmetics. */
	Tensor operator+(const Tensor &t) const;
	Tensor &operator+=(const Tensor &t);
	Tensor operator-(const Tensor &t) const;
	Tensor &operator-=(const Tensor &t);
	Tensor operator*(const iRRAM::REAL &c) const;
	Tensor &operator*=(const iRRAM::REAL &c);
	Tensor operator/(const iRRAM::REAL &c) const;
	Tensor &operator/=(const iRRAM::REAL &c);

	/**
	 * @brief Calculate the tensor product of two tensors.
	 */
	Tensor operator *(const Tensor &t) const;

	/**
	 * @brief A coefficient of the tensor.
	 *
	 * @param indices the indices of the tensor
	 *
	 * @return the corresponding coefficient
	 */
	iRRAM::REAL coeff(const std::vector<int> &indices) const;

	/* TODO: Contraction, Dualize, Transpose */

	/**
	 * @brief Contract the kth contravariant part and the lth covariant parts.
	 *
	 * @param k the contravariant index
	 * @param l the covariant index
	 *
	 * @return this tensor
	 */
	Tensor &contraction(int k, int l);

	/**
	 * @brief Set a dualization isomorphism.
	 *
	 * @param matrix the matrix representing the isomorphism
	 */
	void set_isomorphism(const std::vector<std::vector<iRRAM::REAL>> &matrix);

	/**
	 * @brief Dualize one contravariant or one covariant part.
	 *
	 * @param k the index of contravariant or covariant part
	 *
	 * @return this tensor
	 */
	Tensor &dualize(int k);

	/**
	 * @brief Transpose the tensor
	 *
	 * @param permutation the permutation used to permute the tensor
	 *
	 * @return this tensor
	 */
	Tensor &transpose(const std::vector<int> &permutation);
};


}

#endif // TENSOR_H
