
#include "scomplex.h"
#include "utils.h"

namespace TCERC {

#ifndef OLD_SCOMPLEX
/* TODO: EVERYTHING */

SComplex::SComplex(int D, int d)
{
        if (d > D)
                throw std::invalid_argument("Dimension mismatch");

        this->D = D;
        this->d = d;
        this->n = 0;
        // We are mapping simplices and their index together, for efficient connection management
        // Connection contains n vector<int>, each contains list of its connected simplices(their corresponding indexes)
        // ith vector contains list of indexes - it means ith simplex and jth simplex share facet thus connected
        // When new simplex is added, its new connection vector is added and connection status is updated
        // facet owner contains facet, and its parent simplex in their index

//        simplices.push_back(simplex);
//        idx_simplex.insert(std::pair<Simplex, int>(*simplex, this->n));
//        connections.push_back(std::vector<int> {});
//        for (int i = 0; i < d + 1; i++) {
//                std::pair<Simplex, std::vector<int>> p(simplex->facet(i), std::vector<int> {this->n});
//                this->points.insert(simplex->vertex_list()[i]);
//                this->facet_owner.insert(p);
//        }

//        this->n++;
}

SComplex::SComplex(const Simplex &simplex)
{
	this->D = simplex.ambient_dimension();
	this->d = simplex.dimension();
	this->n = 0;
        this->simplices.push_back(simplex);
	this->simplices_set.insert(simplex);
        this->idx_simplex.insert(std::pair<Simplex, int>(simplex, this->n));
        for (int i = 0; i < d + 1; i++) {
                std::pair<Simplex, std::vector<int>> p(simplex.facet(i), std::vector<int> {this->n});
                this->points.insert(simplex.vertex_list()[i]);
		this->facet_owner.insert(p);
        }
        this->n++;
}

SComplex::~SComplex()
{

}

int SComplex::ambient_dimension() const
{
	return this->D;
}

int SComplex::dimension() const
{
	return this->d;
}

int SComplex::num_simplex() const
{
	return this->n;
}

std::vector<Simplex> SComplex::simplex_list() const
{
	return this->simplices;
}

const Simplex &SComplex::adj_simplex(const Simplex &simplex, const Vertex &vertex) const
{
	auto it = this->idx_simplex.find(simplex);
	int idx = it->second;
	//std::vector<int> list_adj_simplex = this->connections[idx];
	std::vector<Vertex> vertices = simplex.vertex_list();
	for (int i = 0; i < vertices.size(); i++){
                if (this->facet_owner.count(simplex.facet(vertices[i])) == 0 || !simplex.facet(vertices[i]).is_member_vertex(vertex))
			continue;

		auto it = this->facet_owner.find(simplex.facet(vertices[i]));
		if (it->second.size() == 2) {
			int ret_idx = (idx == it->second[0]) ? it->second[1] : it->second[0];
			return this->simplices[ret_idx];
		} 
	}
	return simplex;
}

const Simplex &SComplex::adj_simplex(const Simplex &simplex, int i) const
{
	Vertex vertex = simplex.vertices[i];
	return SComplex::adj_simplex(simplex, vertex);
//	std::map<Simplex, int>::iterator it = this->idx_simplex.find(simplex);
//	int idx = it->second;
//	std::vector<int> list_adj_simplex = this->connections[idx];
//	for (int i = 0; i < list_adj_simplex.size(); i++){
//		if (this->simplices[i].is_member_vertex(vertex))
//			return this->simplices[i]; 
//	}
//	return NULL;
}

std::vector<Simplex> SComplex::adj_simplex_list(const Simplex &simplex) const
{
        std::vector<Simplex> adj_list = {};
	auto it = this->idx_simplex.find(simplex);
	int idx = it->second;
	std::vector<Vertex> vertices = simplex.vertex_list();
	for (int i = 0; i < vertices.size(); i++){
                if (this->facet_owner.count(simplex.facet(vertices[i])) == 0)
			continue;

		auto it = this->facet_owner.find(simplex.facet(vertices[i]));
		if (it->second.size() == 2) {
			int ret_idx = (idx == it->second[0]) ? it->second[1] : it->second[0];
			adj_list.push_back(this->simplices[ret_idx]);
		} 
	}
        return adj_list;
}

SComplex SComplex::boundary()
{
        std::vector<Simplex> bound = {};
        for (auto const &x : facet_owner) {
                if (x.second.size() == 1) {
                        bound.push_back(x.first);
                }
        }
        SComplex boundary_complex = SComplex(bound[0]);
        for (int i = 1; i < bound.size(); i++) {
                boundary_complex.add_simplex(bound[i]);
        }
        return boundary_complex;
}

const Simplex &SComplex::point_locator(const Point &point, int n) const
{

	iRRAM::REAL precision = 1.0;

	for (int i = 0; i < n; i++) {
		precision /= 2.0;
	}

	std::vector<iRRAM::LAZY_BOOLEAN> candidates;

	for (int i = 0; i < this->simplices.size(); i++) {
		std::vector<Vertex> p = simplices[i].vertex_list();

		for (int j = 0; j < p.size(); j++) {
			candidates.push_back((point - *p[j]) * (point - *p[j]) < precision * precision);
		}

		std::vector<int> vertex_counter;
		for (int j = 0; j < p.size(); j++) {
			vertex_counter.push_back(0);
		}
		while (true) {
			int j;
			for (j = 0; j < p.size(); j++) {
				vertex_counter[j] ++;
				if (vertex_counter[j] == 2) {
					vertex_counter[j] = 0;
				} else {
					break;
				}
			}
			if (j == p.size()) break;

			std::vector<Vertex> ps;
			for (int j = 0; j < p.size(); j++) {
				if (vertex_counter[j] == 1) {
					ps.push_back(p[j]);
				}
			}

			if (ps.size() == 1)
				continue;
			
			std::vector<Vector> vs;
			for (int j = 0; j < ps.size() - 1; j++) {
				vs.push_back(*ps[j + 1] - *ps[0]);
			}
			
			std::vector<iRRAM::REAL> result;
			for (int n = 0; n < vs.size(); n++) {
				result.push_back((point - *ps[0]) * vs[n]);
			}
			Vector answer = project_vector(vs, Vector(result));

			iRRAM::REAL cs_sum = 0;
			for (int n = 0; n < vs.size(); n++) {
				cs_sum += answer[n];
			}

			Vector w = Vector(this->D);
			for (int n = 0; n < vs.size(); n++) {
				w = w + vs[n] * answer[n];
			}
			Vector x = (point - *ps[0]) - w;
			
			iRRAM::LAZY_BOOLEAN interior;
			interior = cs_sum < iRRAM::REAL(1.0);
			for (int n = 0; n < vs.size(); n++) {
				interior = interior && answer[n] > 0;
			}
			interior = interior && ((x * x) < precision * precision);
			candidates.push_back(interior);
		}
	}
	int interior_result = choose(candidates) - 1;

	//printf("Choose result = %d\n", interior_result);
	int size_of_simplex = (2 << (this->simplices[0]).dimension()) - 1;
	//printf("Interior result [%lf, %lf, %lf] -> %d, %d\n", point[0].as_double(), point[1].as_double(), point[2].as_double(), interior_result, interior_result / size_of_simplex);
	//iRRAM::cout << "interior result " << interior_result << " simplex " << interior_result / size_of_simplex << std::endl;
	if (interior_result == -1)
		throw std::invalid_argument("no simplex containing the point");
	return this->simplices[interior_result / size_of_simplex];

}

void SComplex::add_simplex(const Simplex &facet, const Vertex &point)
{
        if (this->facet_owner.count(facet) == 0)
                throw std::invalid_argument("Facet Not in Simplicial Complex");

        std::map<Simplex, std::vector<int>>::iterator iter;
        std::vector<Vertex> original_points = facet.vertex_list();
        original_points.push_back(point);
        Simplex new_simplex = Simplex(original_points);
        int d = this->dimension();

        this->idx_simplex.insert(std::pair<Simplex, int>(new_simplex, this->n));
        this->simplices.push_back(new_simplex);
	this->simplices_set.insert(new_simplex);

        for (int i = 0; i < d + 1; i++) {
                this->points.insert(new_simplex.vertex_list()[i]);
                iter = this->facet_owner.find(new_simplex.facet(i));
                // Facet that never used before - create new facet-simplex pair
                if (this->facet_owner.count(new_simplex.facet(i)) == 0)
                        this->facet_owner.insert(std::pair<Simplex, std::vector<int>>(new_simplex.facet(i), std::vector<int> {this->n}));
                // Facet that already used twice - Unavailable to add more
                else if (this->facet_owner.count(new_simplex.facet(i)) != 0 && iter->second.size() == 2)
                        throw std::invalid_argument("Facet Cannot be Used more than Twice");

                else {  
                        int original_number = iter->second[0];
                        iter->second.push_back(this->n);
                }
        }
        this->n++;

}

void SComplex::add_simplex(const std::vector<Vertex> &points)
{
        // Check whether d-1 size subset of &points is a boundary
        // The one point left cannot be validated (equality test)
        int d = this->dimension();
        std::map<Simplex, std::vector<int>>::iterator iter;
        if (points.size() != d+1)
                throw std::invalid_argument("Simplex in different dimension");
        Simplex new_simplex = Simplex(points);

        this->idx_simplex.insert(std::pair<Simplex, int>(new_simplex, this->n));
        this->simplices.push_back(new_simplex);
	this->simplices_set.insert(new_simplex);

        for (int i = 0; i < d + 1; i++) {
                this->points.insert(new_simplex.vertex_list()[i]);
                iter = this->facet_owner.find(new_simplex.facet(i));
                // Facet that never used before - create new facet-simplex pair
                if (this->facet_owner.count(new_simplex.facet(i)) == 0){
                        this->facet_owner.insert(std::pair<Simplex, std::vector<int>>(new_simplex.facet(i), std::vector<int> {this->n}));
                }// Facet that already used twice - Unavailable to add more
                else if (this->facet_owner.count(new_simplex.facet(i)) != 0 && iter->second.size() == 2)
                        throw std::invalid_argument("Facet Cannot be Used more than Twice");

                else {  

                        int original_number = iter->second[0];
                        iter->second.push_back(this->n);
                }
        }
        this->n++;
}

void SComplex::add_simplex(const Simplex &simplex)
{
	SComplex::add_simplex(simplex.vertex_list());
}

#else
SComplex::SComplex(int D, int d, Simplex *simplex) 
{

	if (d > D)
		throw std::invalid_argument("Dimension mismatch");

	this->D = D;
	this->d = d;
	this->n = 0;
	// We are mapping simplices and their index together, for efficient connection management
	// Connection contains n vector<int>, each contains list of its connected simplices(their corresponding indexes)
	// ith vector contains list of indexes - it means ith simplex and jth simplex share facet thus connected
	// When new simplex is added, its new connection vector is added and connection status is updated
	// facet owner contains facet, and its parent simplex in their index

	simplices.push_back(simplex);
	idx_simplex.insert(std::pair<Simplex, int>(*simplex, this->n));
	this->connection.push_back(std::vector<int> P{});	
	for (int i = 0; i < d + 1; i++) {
		std::pair<Simplex, std::vector<int>> p(simplex->facet(i), std::vector<int> {this->n});
		this->points.insert(simplex->vertex_list()[i]);
		this->facet_owner.insert(p);
	}

	this->n++;
}

SComplex::~SComplex()
{

}

int SComplex::ambient_dimension() 
{
	return this->D;
}

int SComplex::dimension()
{
	return this->d;
}

int SComplex::num_simplex()
{
	return this->n;
}

Simplex* SComplex::operator[](const int i) const
{
	if (i > this->simplices.size())
		throw std::invalid_argument("Index Out of Range");
        return this->simplices[i];
}

Simplex* SComplex::adj_simplex(Simplex *simplex, Vertex &vertex)
{
	std::map<Simplex, int>::iterator it = this->idx_simplex.find(*simplex);
	int idx = it->second;
	std::vector<int> list_adj_simplex = this->connections[idx];
	for (int i = 0; i < list_adj_simplex.size(); i++){
		if (this->simplices[i]->is_member_vertex(vertex))
			return this->simplices[i]; 
	}
	return NULL;
}

std::vector<Simplex *> SComplex::list_adj_simplex(Simplex *simplex)
{
	std::vector<Simplex *> adj_list = {};
	std::map<Simplex, int>::iterator it = this->idx_simplex.find(*simplex);
	int idx = it->second;
	std::vector<int> list_adj_simplex = this->connections[idx];
	for (int i = 0; i < list_adj_simplex.size(); i++){
		adj_list.push_back(this->simplices[i]); 
	}
	return adj_list;
}
/*
Vector* boundary_normal(int i)
{

}

Simplex* boundary_facet(int i)
{

}
*/

SComplex SComplex::boundary()
{
	std::vector<Simplex> bound = {};
	for (auto const &x : facet_owner) {
		if (x.second.size() == 1) {
			bound.push_back(x.first);
		}
	}
	SComplex boundary_complex = SComplex(this->D, this->d-1, &bound[0]);
	for (int i = 1; i < bound.size(); i++) {	
		boundary_complex.add_simplex(bound[i].vertex_list());
	}
	return boundary_complex;
}

/*
Simplex* SComplex::point_locator(Point &point, int n)
{

	iRRAM::REAL precision = 1.0;

	for (int i = 0; i < n; i++) {
		precision /= 2.0;
	}

	std::vector<iRRAM::LAZY_BOOLEAN> candidates;

	for (int i = 0; i < this->simplices.size(); i++) {
		//printf("%dth simplex\n", i);
		std::vector<Vertex> p = simplices[i]->vertex_list();

		for (int j = 0; j < p.size(); j++) {
			//printf("Push each vertex[%d] check\n", j); 
			candidates.push_back((point - *p[j]) * (point - *p[j]) < precision * precision);
		}

		std::vector<int> vertex_counter;
		for (int j = 0; j < p.size(); j++) {
			vertex_counter.push_back(0);
		}
		
		while (true) {
			int j;
			for (j = 0; j < p.size(); j++) {
				vertex_counter[j] ++;
				if (vertex_counter[j] == 2) {
					vertex_counter[j] = 0;
				} else {
					break;
				}
			}
			if (j == p.size()) break;

			std::vector<Vertex> ps;
			for (int j = 0; j < p.size(); j++) {
				if (vertex_counter[j] == 1) {
					ps.push_back(p[j]);
				}
			}

			if (ps.size() == 1)
				continue;
			
			std::vector<Vector> vs;
			for (int j = 0; j < ps.size() - 1; j++) {
				vs.push_back(*ps[j + 1] - *ps[0]);
			}
			
			std::vector<Vector> coeff_matrix;
			//iRRAM::REALMATRIX lmatrix(vs.size(), vs.size());
			for (int n = 0; n < vs.size(); n++) {
				std::vector<iRRAM::REAL> row;
				for (int m = 0; m < vs.size(); m++) {
					row.push_back(vs[n] * vs[m]);
					//lmatrix(n, m) = vs[n] * vs[m];
				}
				coeff_matrix.push_back(Vector(row));
			}
			std::vector<iRRAM::REAL> result;
			//iRRAM::REALMATRIX rmatrix(vs.size(), 1);
			for (int n = 0; n < vs.size(); n++) {
				result.push_back((point - *ps[0]) * vs[n]);
				//rmatrix(n, 0) = (*point - *ps[0]) * vs[n];
			}
			Vector answer = linear_solver(coeff_matrix, Vector(result));
			//iRRAM::REALMATRIX cs = solve(lmatrix, rmatrix, 1);
			iRRAM::REAL cs_sum = 0;
			//iRRAM::cout << "Printing Cs" << std::endl;
			for (int n = 0; n < vs.size(); n++) {
				cs_sum += answer[n];
				//cs_sum += cs(n, 0);
			//	iRRAM::cout << " | " << cs(n, 0);
			}
			//iRRAM::cout << std::endl;

			Vector w = Vector(this->D);
			for (int n = 0; n < vs.size(); n++) {
				w = w + vs[n] * answer[n];
				//w = w + vs[n] * cs(n, 0);
			}
			Vector x = (point - *ps[0]) - w;
			
			iRRAM::LAZY_BOOLEAN interior;
			interior = cs_sum < iRRAM::REAL(1.0);
			for (int n = 0; n < vs.size(); n++) {
				interior = interior && answer[n] > 0;
				//interior = interior && cs(n, 0) > 0;
			}
			interior = interior && ((x * x) < precision * precision);
			//printf("Sum < 1 ? %s, Interior ? %s\n", cs_sum < iRRAM::REAL(1.0) ? "Yes" : "No", interior ? "Yes" : "No");
			candidates.push_back(interior);
		}
	}
	
	int interior_result = choose(candidates) - 1;
	//printf("Choose result = %d\n", interior_result);
	int size_of_simplex = (2 << (this->simplices[0])->dimension()) - 1;
	if (interior_result == -1)
		return NULL;
	return this->simplices[interior_result / size_of_simplex];

}
*/

void SComplex::add_simplex(Simplex *facet, Vertex &point)
{
	if (this->facet_owner.count(*facet) == 0)
		throw std::invalid_argument("Facet Not in Simplicial Complex");

	std::map<Simplex, std::vector<int>>::iterator iter;
	std::vector<Vertex> original_points = facet->vertex_list();
	original_points.push_back(point);
	Simplex new_simplex = Simplex(original_points);
	int d = this->dimension();

        this->idx_simplex.insert(std::pair<Simplex, int>(new_simplex, this->n));
	this->simplices.push_back(&new_simplex);
        this->connections.push_back(std::vector<int> {});
        
	for (int i = 0; i < d + 1; i++) {
                this->points.insert(new_simplex.vertex_list()[i]);
		iter = this->facet_owner.find(new_simplex.facet(i));
		// Facet that never used before - create new facet-simplex pair
		if (this->facet_owner.count(new_simplex.facet(i)) == 0)
                	this->facet_owner.insert(std::pair<Simplex, std::vector<int>>(new_simplex.facet(i), std::vector<int> {this->n}));
		// Facet that already used twice - Unavailable to add more
		else if (this->facet_owner.count(new_simplex.facet(i)) != 0 && iter->second.size() == 2)
			throw std::invalid_argument("Facet Cannot be Used more than Twice");

		else {
			int original_number = iter->second[0];
			iter->second.push_back(this->n);
			// Updating the connection informaton - Simplex that contains current facet are now linked within connections
			this->connections[this->n].push_back(original_number);
			this->connections[original_number].push_back(this->n);
		}
        }
	this->n++;
}

void SComplex::add_simplex(std::vector<Vertex> points)
{

	// Check whether d-1 size subset of &points is a boundary
	// The one point left cannot be validated (equality test)
	int d = this->dimension();
	std::map<Simplex, std::vector<int>>::iterator iter;
	if (points.size() != d+1)
		throw std::invalid_argument("Simplex in different dimension");
	Simplex new_simplex = Simplex(points);
	
        this->idx_simplex.insert(std::pair<Simplex, int>(new_simplex, this->n));
	this->simplices.push_back(&new_simplex);
        this->connections.push_back(std::vector<int> {});
        
	for (int i = 0; i < d + 1; i++) {
                this->points.insert(new_simplex.vertex_list()[i]);
		iter = this->facet_owner.find(new_simplex.facet(i));
		// Facet that never used before - create new facet-simplex pair
		if (this->facet_owner.count(new_simplex.facet(i)) == 0){
                	this->facet_owner.insert(std::pair<Simplex, std::vector<int>>(new_simplex.facet(i), std::vector<int> {this->n}));
		}// Facet that already used twice - Unavailable to add more
		else if (this->facet_owner.count(new_simplex.facet(i)) != 0 && iter->second.size() == 2)
			throw std::invalid_argument("Facet Cannot be Used more than Twice");

		else {
			
			int original_number = iter->second[0];
			iter->second.push_back(this->n);
			// Updating the connection informaton - Simplex that contains current facet are now linked within connections
			this->connections[this->n].push_back(original_number);
			this->connections[original_number].push_back(this->n);
		}
        }
	this->n++;
}
#endif

}
