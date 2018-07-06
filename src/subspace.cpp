#include "subspace.h"
#include "utils.h"

namespace TCERC {

Subspace::Subspace(const std::vector<Vector> &vectors)
{
	if (vectors.empty()) {
		return;
	}
	this->d = vectors.size();
	this->D = vectors[0].ambient_dimension();
	for (int i = 1; i < vectors.size(); i++) {
		if (vectors[i].ambient_dimension() != this->D) {
			throw std::invalid_argument("invalid basis");
		}
	}

	std::vector<Vector> result = gaussian_elimination(vectors);
	if (result.empty()) {
		throw std::invalid_argument("basis must be linearly independent");
	}

	this->basis = std::make_shared<std::vector<Vector>>(vectors);
}

Subspace::Subspace(int d)
{
	this->d = d;
	this->D = d;

	std::vector<Vector> base;
	for (int i = 0; i < d; i++) {
		std::vector<iRRAM::REAL> ei;
		for (int j = 0; j < d; j++) {
			if (i == j)
				ei.push_back(1);
			else
				ei.push_back(0);
		}
		base.push_back(Vector(ei));
	}

	this->basis = std::make_shared<std::vector<Vector>>(base);
}


Subspace::Subspace(const Simplex &s)
{
	std::vector<Vertex> vertices = s.vertex_list();
	Vertex main_vertex = vertices[0];
	this->D = s.ambient_dimension();
	this->d = s.dimension();

	std::vector<Vector> vectors;
	for (int i = 1; i <= this->d; i++) {
		vectors.push_back(vertices[i] - main_vertex);
	}
	this->basis = std::make_shared<std::vector<Vector>>(vectors);
}

const std::vector<Vector> &Subspace::get_basis() const
{
	return *this->basis;
}

int Subspace::dimension() const
{
	return this->d;
}

int Subspace::ambient_dimension() const
{
	return this->D;
}

/* TODO */
Vector Subspace::embed(const std::vector<iRRAM::REAL> &vector) const
{
	int ambient_dim = this->D;
	int dim = this->d;
	Vector abs_vector = Vector(this->D);
	if (dim != vector.size())
		throw std::invalid_argument("Vector dimension different from basis dimension");

	for (int i = 0; i < dim; i++) {
		abs_vector += ((*this->basis)[i] * vector[i]);
	}
	return abs_vector;
}

/* TODO */
std::vector<iRRAM::REAL> Subspace::project(const Vector &vector) const
{
	int ambient_dim = this->D;
	int dim = this->d;
	std::vector<iRRAM::REAL> result = {};
	std::vector<Vector> basis_matrix = {};
	
	if (ambient_dim != vector.ambient_dimension())
		throw std::invalid_argument("Vector dimension different from ambient dimension");

	for (int i = 0; i < dim; i++) {
		basis_matrix.push_back((*this->basis)[i]);
	}
	for (int i = 0; i < dim; i++) {
		result.push_back((*this->basis)[i]*vector);
	}
	return project_vector(basis_matrix, Vector(result)).get_coordinate();
}

std::vector<iRRAM::REAL> Subspace::normal(const std::vector<std::vector<iRRAM::REAL>> &vectors) const
{
	/*
	 * First, let T1 = basis constructed from simplex's vertices, T1 = [b1 b2 b3 ... bd]
	 * considering the vertex we provided as an origin. T1 is a (D, d)-matrix
	 * Then (the facet's vertex - our vertex(origin)) is the basis itself,
	 * and the facet's edge is represented using T1. Let each edge as vi, i = 1 ... d,
	 * each vi is d-std::vector, and it's value at ambient space is (T1 * vi) which is D-std::vector
	 * We need exactly v1 ~ v(d - 1) to fully define the facet's subspace,
	 * which lies on d-dimensional subspace defined by T1.
	 * Then, let V = [v1; v2; v3; ...; vi] which is a (d - 1, d)-matrix. Calculate gaussian-elimination, denote Vg
	 * Vg = [Diag(d - 1); e_d] where Diag(d - 1) is (d - 1, d - 1) diagonal matrix from elimination, denote D,
	 * and e_d is d-std::vector. Denote (e1, e2, e3..., e(d - 1))
	 * Then the normal std::vector of this facet is the null space if V, which is defined as Vn = 0;
	 * n = (D(1, 1)/e1, D(2, 2)/e2, D(3, 3)/e3, ... D(d - 1, d - 1)/e(d - 1), -1)
	 *	Here, vi = (0, 0, 0, ..., -1, 1, 0, ..., 0) where only (i - 1)th element is -1 and ith element is 1
	 *		( vi is edge connecting ith basis of T1 and (i - 1)th basis of T1 )
	 *	So, n = (1, 1, 1, 1, ... 1) without need to calculate V and V
	 * Now convert this normal std::vector to ambient space.
	 * Let the new normal std::vector as (T2 * n), and it's normal to all (T1 * vi)s, so,
	 * (T2 * n) (_dot_) (T1*vi) = 0 = (T2 * n)' * (T1 * vi) = n' * (T2' * T1) * vi
	 * We already know that (n' * vi) = n (_dot_) vi = 0, since n is orthogonal to vi.
	 * We want (T2' * T1) is identity matrix, so T2' = left inverse of T1 (exists since T1 is basis with full rank)
	 * T2' = (T1' * T1)^(-1) * T1', so T2 = T1 * ((T1' * T1)^(-1))' = T1 * ((T1' * T1)')^(-1).
	 * T2 = T1 * (T1' * T1) ^ (-1)   // (T'T)' = T'T
	 *	(Here (T')^(-1) = (T^(-1))' is used too)
	 *
	 * T1' * T1 's inverse is calculated using guass-jordan elimination
	 */


	// Calculate n = norm_facet
	std::vector<iRRAM::REAL> norm_facet;
	for (int i = 0; i < this->d; i++) {
		norm_facet.push_back(1.0);
	}

	std::vector<Vector> basisT;
	for (int i = 0; i < vectors.size(); i++) {
		basisT.push_back(Vector(vectors[i]));
	}

	// Calculate T1'*T1 where (i, j)th element is a dot product of ith and jth basis
	//printf("Calculating T1'*T1\n");
	std::vector<Vector> basisT_basis;
	for (int i = 0; i < this->d; i++) {
		std::vector<iRRAM::REAL> row;
		for (int j = 0; j < this->d; j++) {
			row.push_back(basisT[i] * basisT[j]);
		}
		Vector v(row);
		basisT_basis.push_back(v);
	}


	// Expand each row of T1'*T1 to make a jordan form. Each row is 2d-std::vector
	//printf("Constructing jordan form\n");
	std::vector<Vector> basisT_basis_jordan;
	for (int i = 0; i < this->d; i++) {
		std::vector<iRRAM::REAL> row = basisT_basis[i].get_coordinate();
		for (int j = 0; j < this->d; j++) {
			row.push_back(i == j ? 1.0 : 0.0);
		}
		Vector v(row);
		basisT_basis_jordan.push_back(v);
	}

	// Compuate and take (d + 1 ~ 2d)th element of each row as an inverse
	//printf("Compute inverse\n");
	std::vector<Vector> basisT_basis_inverse_jordan = gaussian_elimination(basisT_basis_jordan);
	std::vector<Vector> basisT_basis_inverse;
	for (int i = 0; i < this->d; i++) {
		std::vector<iRRAM::REAL> row;
		for (int j = 0; j < this->d; j++) {
			row.push_back(basisT_basis_inverse_jordan[i][d + j]);
		}
		Vector v(row);
		basisT_basis_inverse.push_back(v);
	}

	// Calculate T1 as a transpose of T1'
	//printf("Compute transpose of T1' = T1\n");
	std::vector<Vector> basis;
	for (int j = 0; j < this->D; j++) {
		std::vector<iRRAM::REAL> row;
		for (int i = 0; i < this->d; i++) {
			row.push_back(basisT[i][j]);
		}
		Vector v(row);
		basis.push_back(v);
	}

	// Calculate T2 = T1 * (T1'*T1)^(-1)
	//printf("Compute T2 = T1 * (T1'*T1)^(-1)\n");
	std::vector<Vector> basis_norm;
	for (int i = 0; i < this->D; i++) {
		std::vector<iRRAM::REAL> row;
		for (int j = 0; j < this->d; j++) {
			row.push_back(basis[i] * basisT_basis_inverse[j]); // BasisT_Basis is same with it's transpose, so using row instead of column is totally fine
		}
		Vector v(row);
		basis_norm.push_back(v);
	}

	// Finally calculate T2 * n, which ith element is a dot product of T2's ith row and n itself
	std::vector<iRRAM::REAL> norm;
	Vector norm_facet_vertex(norm_facet);

	for (int i = 0; i < this->D; i++) {
		norm.push_back(basis_norm[i] * norm_facet_vertex);
	}

	iRRAM::REAL sum = 0;
	for (int i = 0; i < this->D; i++) {
		sum += norm[i] * norm[i];
	}
	sum = iRRAM::sqrt(sum);

	for (int i = 0; i < this->D; i++) {
		norm[i] /= sum;
	}

	return norm;

}

}
