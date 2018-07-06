/**
 * @file simplex.h
 * @author Joonhyung Shin (tonyshin@kaist.ac.kr)
 * @author Namjo Ahn
 * @author Seungwoo Lee
 * @date May, 2018
 * @brief This is the header file for Vertex and Simplex.
 */

#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <iRRAM/lib.h>
#include "euclidean.h"
#include <set>
#include <map>
#include <memory>

namespace TCERC {

class Simplex;

/**
 * @brief A wrapper class for Point.
 */
class Vertex {
private:
	std::shared_ptr<Point> point;

public:
	/* Constructors. */

	/**
	 * @brief Construct a vertex located at the origin.
	 *
	 * @param D the dimension of the vertex
	 */
	Vertex(int D) : point(new Point(D)) {}
	Vertex(const std::vector<iRRAM::REAL> &coordinate) : point(new Point(coordinate)) {}
	Vertex(const Point &p) : point(new Point(p)) {}
	Vertex() : Vertex(0) {}

	/* Basic arithmetics. (Only operator- is provided until now) */
	Vector operator-(const Vertex &v) const { return *point - *v.point; }

	/* Destructor. */
	~Vertex() {}

	/* Comparisons. */

	/**
	 * @brief Determines whether two vertices point at the same point.
	 * This can be false even if the two points are actually equal.
	 */
	friend bool operator==(const Vertex &v, const Vertex &w) { return v.point == w.point; }

	/**
	 * @brief Determines whether two vertices point at the different point.
	 */
	friend bool operator!=(const Vertex &v, const Vertex &w) { return v.point != w.point; }

	/* Basic functionalities. */
	const Point &operator*() const { return point.operator*(); }
	const Point *operator->() const { return point.operator->(); }
	const iRRAM::REAL &operator[](int i) const { return (*point)[i]; }
	/**
	 * @brief The dimension of the vertex.
	 *
	 * @return the dimension of the vertex.
	 */
	int ambient_dimension() const { return point->ambient_dimension(); }

	/* This is to make a comparator of Vertex. */
	friend class Simplex;
};

class SComplex;

/**
 * @brief A geometric simplex lying in a Euclidean space.
 */
class Simplex {
private:
	int D;			/* ambient dimension */
	int d;			/* actual dimension */

	/**
	 * @brief Comparator for Vertices.
	 */
	struct VertexCompare {
		bool operator()(const Vertex &v, const Vertex &w) const {
			return v.point < w.point;
		}
	};
	std::set<Vertex, VertexCompare> vertex_set;
	std::vector<Vertex> vertices;
	// std::vector<Point> vertices_container;	//Deprecated Because Point is maintained by user
	std::map<Vertex, Simplex, VertexCompare> facets;

public:
	/* Constructors. */

	/**
	 * @brief Construct a standard d-simplex.
	 *
	 * @param d the dimension of the simplex.
	 */
	Simplex(int d);			/* standard d-simplex */

	Simplex(const std::vector<Vertex> &vertices);
	//Simplex(int D, int d, std::vector<std::vector<iRRAM::REAL> *> &coordinate);

	Simplex() : Simplex(0) {}

	/* Destructor. */
	~Simplex();

	/* Comparisons. */

	/**
	 * @brief Determine whether two simplices have the same set of vertices.
	 * The result can be false even if the two simplices are geometrically equal.
	 */
	friend bool operator==(const Simplex &x, const Simplex &y) {
		return x.vertex_set == y.vertex_set;
	}

	friend bool operator!=(const Simplex &x, const Simplex &y) {
		return x.vertex_set != y.vertex_set;
	}

	/* Basic functionalities. */

	/**
	 * @brief Check whether a vertex is one of the vertices of the simplex.
	 *
	 * @param vertex the vertex willing to check
	 *
	 * @return true if the vertex is one of the vertices of the simplex
	 */
	bool is_member_vertex(const Vertex &vertex) const;

	/**
	 * @brief The ambient dimension of the simplex.
	 *
	 * @return the ambient dimension of the simplex.
	 */
	int ambient_dimension() const;

	/**
	 * @brief The dimension of the simplex.
	 *
	 * @return the dimension of the simplex.
	 */
	int dimension() const;

	/**
	 * @brief A list of the vertices of the simplex.
	 *
	 * @return A list of vertices of the simplex.
	 */
	const std::vector<Vertex> &vertex_list() const;

	/**
	 * @brief A facet of the simplex.
	 *
	 * @param vertex the vertex not in the desired facet.
	 *
	 * @return the simplex consisting of the vertices other than the one privided.
	 */
	Simplex facet(const Vertex &vertex) const;

	Simplex facet(int i) const;

	/**
	 * @brief A normal vector of the simplex.
	 *
	 * @param vertex one of the vertices of the simplex.
	 *
	 * @return the normal vector orthogonal to the given vertex.
	 */
	Vector normal(const Vertex &vertex) const;

	Vector normal(int i) const;

	/**
	 * @brief Compute the content (volume) of the simplex.
	 *
	 * @return the content of the simplex
	 */
	iRRAM::REAL content() const;

	/* This is to make a comparator of Simplex. */
	friend class SComplex;
};

}

#endif /* simplex.h */
