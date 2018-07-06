#ifndef SCOMPLEX_H
#define SCOMPLEX_H

#include <iRRAM/lib.h>
#include <vector>
#include <utility>
#include <set>
#include <map>
#include <algorithm>

#include "simplex.h"

namespace TCERC {

/**
 * @brief A homogeneous simplicial complex homeomorphic to some manifold, possibly with boundary.
 */
class SComplex {
private:
#ifndef OLD_SCOMPLEX
	std::vector<Simplex> simplices;
        int D;
        int d;
        int n;
        // We are mapping simplices and their index together, for efficient connection management
        // Connection contains n vector<int>, each contains list of its connected simplices(their corresponding indexes)
        // ith vector contains list of indexes - it means ith simplex and jth simplex share facet thus connected
        // When new simplex is added, its new connection vector is added and connection status is updated
        // facet owner contains facet, and its parent simplex in their index

        struct SimplexCompare {
                bool operator()(const Simplex &x, const Simplex &y) const {
                        auto &x_set = x.vertex_set;
                        auto &y_set = y.vertex_set; 
                        if (x_set.size() != y_set.size())
                                return x_set.size() < y_set.size();
                        for (auto it_x = x_set.begin(), it_y = y_set.begin();
                             it_x != x_set.end();
                             it_x++, it_y++) { 
                                if (*it_x != *it_y)
                                        return Simplex::VertexCompare()(*it_x, *it_y);
                        }
                        return false;
                }
        };   
	// For SComplex comparasion
	std::set<Simplex, SimplexCompare> simplices_set;
        std::map<Simplex, int, SimplexCompare> idx_simplex;
        std::set<Vertex, Simplex::VertexCompare> points;
        std::map<Simplex, std::vector<int>, SimplexCompare> facet_owner;
        friend bool operator==(const SComplex &x, const SComplex &y) {
                return x.simplices_set == y.simplices_set;
        }

#else
	int D;
	int d;
	int n;
        // We are mapping simplices and their index together, for efficient connection management
        // Connection contains n vector<int>, each contains list of its connected simplices(their corresponding indexes)
        // ith vector contains list of indexes - it means ith simplex and jth simplex share facet thus connected
        // When new simplex is added, its new connection vector is added and connection status is updated
        // facet owner contains facet, and its parent simplex in their index

	struct SimplexCompare {
		bool operator()(const Simplex &x, const Simplex &y) const {
			auto &x_set = x.vertex_set;
			auto &y_set = y.vertex_set;
			if (x_set.size() != y_set.size())
				return x_set.size() < y_set.size();
			for (auto it_x = x_set.begin(), it_y = y_set.begin();
			     it_x != x_set.end();
			     it_x++, it_y++) {
				if (*it_x != *it_y)
					return Simplex::VertexCompare()(*it_x, *it_y);
			}
			return false;
		}
	};
	std::vector<Simplex *> simplices;
	std::map<Simplex, int, SimplexCompare> idx_simplex;
	std::set<Vertex, Simplex::VertexCompare> points;
	std::vector<std::vector<int>> connections; // First Simplex represents 
	//std::vector<Simplex *> boundaries;
	std::map<Simplex, std::vector<int>, SimplexCompare> facet_owner;
#endif

public:
#ifndef OLD_SCOMPLEX
	/* TODO: EVERYTHING */

	/**
	 * @brief Construct an empty simplicial complex.
	 *
	 * @param D the dimension of the ambient space
	 * @param d the dimension of the simplicial complex.
	 */
	SComplex(int D, int d);

	/**
	 * @brief Construct a simplicial complex with one simplex
	 *
	 * @param simplex the simplex composing the complex
	 */
	SComplex(const Simplex &simplex);

	~SComplex();

	/**
	 * @brief The ambient dimension of the simplicial complex.
	 *
	 * @return the ambient dimension of the simplicial complex
	 */
	int ambient_dimension() const;

	/**
	 * @brief The dimension of the simplicial complex.
	 *
	 * @return the dimension of the simplicial complex
	 */
	int dimension() const;

	/**
	 * @brief The number of simplices in the simplicial complex.
	 *
	 * @return the number of simplices in the simplicial complex
	 */
	int num_simplex() const;

	/**
	 * @brief The list of simplices in the simplicial complex.
	 *
	 * @return the list of simplices in the simplicial complex
	 */
	std::vector<Simplex> simplex_list() const;

	/**
	 * @brief Find a simplex adjacent to a particular simplex.
	 * If no such simplex exists, return itself. (TODO)
	 *
	 * @param simplex the benchmark simplex
	 * @param vertex the vertex not contained in the desired simplex
	 *
	 * @return the corresponding adjacent simplex
	 */
	const Simplex &adj_simplex(const Simplex &simplex, const Vertex &vertex) const;

	const Simplex &adj_simplex(const Simplex &simplex, int i) const;

	/**
	 * @brief Find all the adjacent simplices of a particular simplex.
	 *
	 * @param simplex the benchmark simplex
	 *
	 * @return the list of adjacent simplices
	 */
	std::vector<Simplex> adj_simplex_list(const Simplex &simplex) const;

	/**
	 * @brief Compute the boundary of the simplicial complex.
	 *
	 * @return the boundary simplicial complex
	 */
	SComplex boundary();

	/**
	 * @brief Find a simplex that the point lies in in some precision.
	 *
	 * @param point the point for the test
	 * @param n the precision paramter
	 *
	 * @return some simplex with distance between the point at most 2^-n
	 */
	const Simplex &point_locator(const Point &point, int n) const;

	/**
	 * @brief Add a new simplex to the simplicial complex.
	 *
	 * @param facet the boundary facet to attach the new simplex
	 * @param point the point in the new simplex other than the facet
	 */
	void add_simplex(const Simplex &facet, const Vertex &point);

	/**
	 * @brief Add a new simplex to the simplicial complex.
	 *
	 * @param points the points composing the new simplex.
	 */
	void add_simplex(const std::vector<Vertex> &points);

	/**
	 * @brief Add a new simplex to the simplicial complex.
	 *
	 * @param new simplex.
	 */
	void add_simplex(const Simplex &simpex);
#else
	SComplex(int D, int d, Simplex *simplex);
	~SComplex();

	int ambient_dimension();
	int dimension();
	int num_simplex();
	Simplex * operator[](const int i) const;
	Simplex *adj_simplex(Simplex *simplex, Vertex &vertex);
	Simplex *adj_simplex(Simplex *simplex, int i);
	std::vector<Simplex *> list_adj_simplex(Simplex *simplex);
//	Simplex *boundary_facet(int i);			/* arbitrary order */
//	Vector *boundary_normal(int i);
	SComplex boundary();

	Simplex *point_locator(Point &point, int n);

	std::vector<Simplex *> simplex_list();

	void add_simplex(Simplex *facet, Vertex &point);
	void add_simplex(std::vector<Vertex> points);
#endif
};

}

#endif /* scomplex.h */
