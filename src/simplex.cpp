#include "simplex.h"
#include "utils.h"

#include <stdexcept>
#include <algorithm>

namespace TCERC {

Simplex::Simplex(int d)
{
	this->D = d + 1;
	this->d = d;

	for (int i = 0; i < d + 1; i++) {
		std::vector<iRRAM::REAL> ei(d + 1, 0);
		ei[i] = 1;
		vertex_set.insert(ei);
		vertices.push_back(ei);
	}
}

Simplex::Simplex(const std::vector<Vertex> &vertices)
{
	int d = vertices.size() - 1;

	if (d < 0)
		throw std::invalid_argument("no vertex provided");

	int D = vertices.front().ambient_dimension();

	if (D < d)
		throw std::invalid_argument("ambient dimension less than simplex dimension");

	for (auto &vertex : vertices) {
		if (vertex.ambient_dimension() != D)
			throw std::invalid_argument("dimension of input must be equal to the ambient dimension");

		this->vertex_set.insert(vertex);
	}

	this->d = d;
	this->D = D;
	this->vertices = vertices;

	// TODO : Check linear independency of user-input vertices
	/* if(!is_linear_independent(vertices)){
	 *	fprintf(stderr,"ERROR : Simplex initialization failed\n");
	 *	fprintf(stderr,"\tInput vertices lies on lower dimension, i.e. linearly dependent std::vectors\n");
	 *	exit(-1);
	 * }
	 */

	//this->construct_facets();
}
/*
Simplex::Simplex(int D, int d, std::vector< std::vector<iRRAM::REAL>* > &coordinates){

	this->vertices_container.reserve(d + 1);
	for (int i = 0; i < d + 1; i++) {
		Point v = Point(D, *coordinates[i]);
		this->vertices_container.push_back(v);
	}
	fo
		this->r (int i = 0; i < d + 1; i++) {vertices.push_back(&(this->vertices_container[i]));
	}
	this->D = D;
	this->d = d;
}
*/
Simplex::~Simplex()
{

}

bool Simplex::is_member_vertex(const Vertex &vertex) const
{
	//int target;
	/*for(target = 0; target < this->vertices.size(); i++){
		if(vertex == &(this->vertices)[target]){
			return true;
		}
	}*/
	for (std::set<Vertex, VertexCompare>::iterator it = this->vertex_set.begin(); it != this->vertex_set.end(); ++it) {
	//	iRRAM::cout << " @ " << (*it)[0] << " | " << (*it)[1] << " | " << (*it)[2] << std::endl;


	}
	return this->vertex_set.find(vertex) != this->vertex_set.end();

	/*
	for (Vertex &v : this->vertices) {
		if (vertex == v) {
			return true;
		}
	}
	return false;
	*/
}


int Simplex::ambient_dimension() const
{
	return this->D;
}

int Simplex::dimension() const
{
	return this->d;
}

static bool equal_vector(std::vector<Point *> v1, std::vector<Point *> v2)
{
	std::sort(v1.begin(), v1.end());
	std::sort(v2.begin(), v2.end());
	return v1 == v2;
}

const std::vector<Vertex> &Simplex::vertex_list() const
{
	return this->vertices;
}

/* TODO: Remove construct_facets and make these const */

Simplex Simplex::facet(const Vertex &vertex) const
{
	std::vector<Vertex> facet_v;
	for (auto &v : this->vertex_set) {
		if (v != vertex)
			facet_v.push_back(v);
	}
	return Simplex(facet_v);
}

Simplex Simplex::facet(int i) const
{
	return facet(this->vertices[i]);
}

Vector Simplex::normal(const Vertex &vertex) const
{
	if (!this->is_member_vertex(vertex)) {
		throw std::invalid_argument("input vertex is not a member of this simplex");
		//fprintf(stderr, "WARNING : Simplex::normal(Point *vertex)\n");
		//fprintf(stderr, "\tInput vertex is not a member of this Simplex!\n");
		//fprintf(stderr, "\t\tRETURN empty Vector\n");
		//return Vector();
	}


	//Simplex *facet = this->facet(vertex);




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

	// Calculate T1', each row is (each vertex - our given vertex) ignoring ours
	//printf("Constructing T1'\n");
	std::vector<Vector> basisT;
	for (int i = 0; i < this->d + 1; i++) {
		if (this->vertices[i] == vertex) {
			continue;
		}
		/*
		printf("%dth row : ", i);
		for (int j = 0; j < this->D; j++) {
			iRRAM::cout << "[" << (*(this->vertices[i]) - *vertex)[j] << "] ";
		}
		printf("\n");
		*/
		basisT.push_back(this->vertices[i] - vertex);
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
		/*
		printf("%dth row : ", i);
		for (int j = 0; j < this->d; j++) {
			iRRAM::cout << "[" << v[j] << "] ";
		}
		printf("\n");
		*/
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
		/*
		printf("%dth row : ", i);
		for (int j = 0; j < this->d*2; j++) {
			iRRAM::cout << "[" << v[j] << "] ";
		}
		printf("\n");
		*/
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
		/*
		printf("%dth row : ", i);
		for (int j = 0 ; j < this->d; j++) {
			iRRAM::cout << "[" << v[j] << "] ";
		}
		printf("\n");
		*/
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
		/*
		printf("%dth row :", j);
		for (int i = 0; i < this->d; i++) {
			iRRAM::cout << "[" << v[i] << "] ";
		}
		printf("\n");
		*/
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
		/*
		printf("%dth row : ", i);
		for (int j = 0; j < this->d; j++) {
			iRRAM::cout << "[" << v[j] << "] ";
		}
		printf("\n");
		*/
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

	return Vector(norm);
}

Vector Simplex::normal(int i) const
{
	return normal(this->vertices[i]);
}

iRRAM::REAL Simplex::content() const
{
	std::vector<std::vector<iRRAM::REAL>> cm_mat;
	const std::vector<Vertex> &v = this->vertices;

	int d = this->d;
	cm_mat.resize(d + 2);
	for (int i = 0; i < d + 1; i++) {
		for (int j = 0; j < d + 1; j++) {
			cm_mat[i].push_back((v[i] - v[j]) * (v[i] - v[j]));
		}
		cm_mat[i].push_back(1);
	}
	for (int j = 0; j < d + 1; j++) {
		cm_mat[d + 1].push_back(1);
	}
	cm_mat[d + 1].push_back(0);

	iRRAM::REAL cm_det = det(cm_mat);
	iRRAM::REAL d_fac = 1;
	iRRAM::REAL pow_d = 1;
	iRRAM::REAL sign = -1;
	for (int i = 1; i <= d; i++) {
		d_fac *= i;
		pow_d *= 2;
		sign = -sign;
	}

	return iRRAM::sqrt(sign * cm_det / pow_d) / d_fac;
}

}
