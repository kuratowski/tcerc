/**
 * @file subspace.h
 * @author Joonhyung Shin (tonyshin@kaist.ac.kr)
 * @author Namjo Ahn
 * @author Seungwoo Lee
 * @date May, 2018
 * @brief This is the header file for Subspace.
 */

#ifndef SUBSPACE_H
#define SUBSPACE_H

#include <iRRAM/lib.h>
#include <vector>
#include <memory>

#include "simplex.h"
#include "euclidean.h"

namespace TCERC {

/**
 * @brief A subspace of a Euclidean space.
 */
class Subspace {
private:
	int D;
	int d;
	std::shared_ptr<std::vector<Vector>> basis;

public:
	/**
	 * @brief Construct a Euclidean space of the given dimension.
	 *
	 * @param d the dimension of the subspace
	 */
	Subspace(int d);

	/**
	 * @brief Construct a subspace with the given basis.
	 *
	 * @param vectors a basis of the subspace
	 */
	Subspace(const std::vector<Vector> &vectors);

	/**
	 * @brief Construct a subspace the simplex lies in.
	 *
	 * @param simplex a simplex lies in the subspace
	 */
	Subspace(const Simplex &s);

	/**
	 * @brief Determine whether two subspaces are pointing to the same basis.
	 * The result may be false even if two subspaces are mathematically equal.
	 */
	friend bool operator==(const Subspace &v, const Subspace &w) {
		return v.basis == w.basis;
	}

	/**
	 * @brief Determine whether two subspaces are pointing to different bases.
	 */
	friend bool operator!=(const Subspace &v, const Subspace &w) {
		return v.basis != w.basis;
	}

	/**
	 * @brief Get a basis of the subspace.
	 *
	 * @return a basis of the subspace
	 */
	const std::vector<Vector> &get_basis() const;

	/**
	 * @brief The dimension of the subspace.
	 *
	 * @return the dimension of the subspace
	 */
	int dimension() const;

	/**
	 * @brief The ambient dimension of the subspace.
	 *
	 * @return the ambient dimension of the subspace
	 */
	int ambient_dimension() const;

	/**
	 * @brief The absolute coordinate of a vector.
	 *
	 * @param vector the coordinate with respect to the basis of the subspace
	 *
	 * @return the translated coordinate of the vector
	 *
	 * TODO
	 */
	Vector embed(const std::vector<iRRAM::REAL> &vector) const;

	/**
	 * @brief The projection of a vector to the subspace.
	 *
	 * @param vector the vector we wish to project in absolute coordinate
	 *
	 * @return the projected vector in relative coordinate
	 *
	 * TODO
	 */
	std::vector<iRRAM::REAL> project(const Vector &vector) const;

	/**
	 * @brief Find a normal vector orthogonal to the given hyperplane.
	 *
	 * @param vectors the vectors that spans the hyperplane
	 *
	 * @return a normal vector orthonogal to the given hyperplane
	 *
	 * TODO
	 */
	std::vector<iRRAM::REAL> normal(const std::vector<std::vector<iRRAM::REAL>> &vectors) const;
};

}

#endif // VERTEX_H
