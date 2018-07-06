#include <iRRAM.h>
#include "scomplex.h"

using namespace iRRAM;
using namespace TCERC;
using namespace std;

void print_vertex(const Vertex &v, const vector<Vertex> &vs);

void compute()
{
	Vertex v1({1, 0, 0});
	Vertex v2({-1, 0, 0});
	Vertex v3({0, 1, 0});
	Vertex v4({0, -1, 0});
	Vertex v5({0, 0, 1});
	Vertex v6({0, 0, -1});

	vector<Vertex> vs = {v1, v2, v3, v4, v5, v6};

	/* Triangulation of 3-ball */
	Simplex s1({v1, v3, v5, v6});
	Simplex s2({v1, v4, v5, v6});
	Simplex s3({v2, v3, v5, v6});
	Simplex s4({v2, v4, v5, v6});

	/* Simplicial complex structure of 3-ball */
	SComplex ball(s1);
	ball.add_simplex(s2);
	ball.add_simplex(s3);
	ball.add_simplex(s4);

	/* Boundary complex gives a triangulation of 2-sphere */
	SComplex sphere = ball.boundary();

	/* Three points that lie in the triangulated 2-sphere */
	Point p1({REAL(1) / 3, -REAL(1) / 3, REAL(1) / 3});
	Point p2({0, REAL(1) / 2, -REAL(1) / 2});
	Point p3({-1, 0, 0});

	/* Report the point-locating simplex */
	Simplex sp1 = sphere.point_locator(p1, 30);
	iRRAM::cout << "P1 is contained in a simplex with vertices";
	for (auto &v : sp1.vertex_list()) {
		iRRAM::cout << " ";
		print_vertex(v, vs);
	}
	iRRAM::cout << endl;

	Simplex sp2 = sphere.point_locator(p2, 30);
	iRRAM::cout << "P2 is contained in a simplex with vertices";
	for (auto &v : sp2.vertex_list()) {
		iRRAM::cout << " ";
		print_vertex(v, vs);
	}
	iRRAM::cout << endl;

	Simplex sp3 = sphere.point_locator(p3, 30);
	iRRAM::cout << "P3 is contained in a simplex with vertices";
	for (auto &v : sp3.vertex_list()) {
		iRRAM::cout << " ";
		print_vertex(v, vs);
	}
	iRRAM::cout << endl;
}

void print_vertex(const Vertex &v, const vector<Vertex> &vs)
{
	for (int i = 1; i <= vs.size(); i++) {
		if (v == vs[i - 1])
			iRRAM::cout << "V" << i;
	}
}
