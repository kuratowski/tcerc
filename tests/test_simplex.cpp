#include "test_main.h"
#include "simplex.h"

#include <cstring>

using namespace TCERC;
using namespace iRRAM;
using namespace std;

TEST(TestSimplex, Create)
{
	Vertex v1({0, 0, 0});
	Vertex v2({1, 0, 0});
	Vertex v3({1, 1, 0});
	Vertex v4({0, 0, 1});
	Vertex v5({0, 0, 1});

	Simplex s1({v1, v2, v3, v4});
	Simplex s2({v4, v2, v3, v1});
	Simplex s3({v1, v2, v3, v5});

	vector<Vertex> s_vertices = s1.vertex_list();

	Simplex f1 = s1.facet(s_vertices[0]);
	Simplex f2 = s1.facet(s_vertices[0]);
	Simplex f3 = s1.facet(s_vertices[1]);

	EXPECT_TRUE(f1 == f2);
	EXPECT_FALSE(f1 == f3);
	EXPECT_TRUE(s1 == s2);
	EXPECT_FALSE(s1 == s3);

	EXPECT_TRUE(s1.is_member_vertex(v3));
	EXPECT_FALSE(s1.is_member_vertex(v5));
}

TEST(TestSimplex, Content)
{
	Vertex v0({0, 0, 1});
	Vertex v1({0, 1, 0});
	Vertex v2({1, 0, 0});
	Vertex v3({0, 0, 0});

	Simplex s2({v0, v1, v2});
	Simplex s3({v0, v1, v2, v3});

	REAL cont2 = s2.content();
	REAL cont3 = s3.content();

	REAL cont2_sq = cont2 * cont2;
	EXPECT_EQ(choose(cont2_sq * 4 > 3, cont2_sq * 4 < 3.001), 2);
	EXPECT_EQ(choose(cont2_sq * 4 < 3, cont2_sq * 4 > 2.999), 2);

	EXPECT_EQ(choose(cont3 * 6 > 1, cont3 * 6 < 1.001), 2);
	EXPECT_EQ(choose(cont3 * 6 < 1, cont3 * 6 > 0.999), 2);
}

/*
TEST(TestSimplex, Create)
{
	printf("\n");
	for (int i = 0; i < 40; i++)
		printf("-");
	printf("\nTEST_SIMPLEX_CREATE\n");

	// Create 3-dimensional Simplex in 3D ambient space
	vector<iRRAM::REAL> p1 = {REAL("0.0"), REAL("0.0"), REAL("0.0")};
	vector<iRRAM::REAL> p2 = {REAL(1.0), REAL(0.0), REAL(0.0)};
	vector<iRRAM::REAL> p3 = {REAL(1.0), REAL(1.0), REAL(0.0)};
	vector<iRRAM::REAL> p4 = {REAL(0.0), REAL(0.0), REAL(1.0)};
	vector<iRRAM::REAL> p5 = {REAL(0.0), REAL(0.0), REAL(1.0)};

	Point v1 = Point(3, p1);
	Point v2 = Point(3, p2);
	Point v3 = Point(3, p3);
	Point v4 = Point(3, p4);
	Point v5 = Point(3, p5);

	// Point Operation Test
	Point add = v2 + v3;
	Point sub = v2 - v3;
	REAL dot = v2 * v3;
	Point s_mul = v2 * REAL(3.0);
	Point s_div = v2 / REAL(3.0);

	printf("Point operation Didn't Failed\n");
	printf("\nAdd Test : [");
	iRRAM::cout << add[0] << " " << add[1] << "" << add[2];
	printf("]\n");

	printf("Sub Test : [");
	iRRAM::cout << sub[0] << " " << sub[1] << "" << sub[2];
	printf("]\n");

	printf("Dot Test : [");
	iRRAM::cout << dot;
	printf("]\n");

	printf("Scalar Mul Test : [");
	iRRAM::cout << s_mul[0] << " " << s_mul[1] << "" << s_mul[2];
	printf("]\n");

	printf("Scalar Div Test : [");
	iRRAM::cout << s_div[0] << " " << s_div[1] << "" << s_div[2];
	printf("]\n");

	// Standard basis and Standard Simplex
	vector<Point *> vertices1 = {&v1, &v2, &v3, &v4};
	//vector<vector<iRRAM::REAL> *> points1 = {&p1, &p2, &p3, &p4};
	// Same simplex, created in different Point order
	vector<Point *> vertices2 = {&v4, &v2, &v3, &v1};
	//vector<vector<iRRAM::REAL> *> points2 = {&p4, &p2, &p3, &p1};
	// Same simplex(coordinate), but different vertex object
	vector<Point *> vertices3 = {&v1, &v2, &v3, &v5};
	//vector<vector<iRRAM::REAL> *> points3 = {&p1, &p2, &p3, &p5};

	Simplex s1 = Simplex(vertices1);
	Simplex s2 = Simplex(vertices2);
	Simplex s3 = Simplex(vertices3);

	ASSERT_TRUE(s1.dimension() != -1);

	// Print Properties
	printf("\nStandard 3/3 simplex created!\n");
	printf("Dimension[%d] imbedded on Ambient Space[%d]\n", s1.dimension(), s1.ambient_dimension());
	printf("Must consists dimension+1 points[%d], actual size of vertices vector[%zu]\n", s1.dimension()+1, s1.vertex_list().size());

	printf("\nNormal vectors of first simplex\n");
	printf("Must be\n[ 1/sqrt(2)] [         0] [ 1/sqrt(2)]\n[-1/sqrt(2)] [ 1/sqrt(2)] [         0]\n[         0] [        -1] [         0]\n[         0] [         0] [        -1]\n");
	printf("Result is\n");
	vector<Point *> s_vertices = s1.vertex_list();
	for (int i = 0; i < s_vertices.size(); i++){
		Vector norm = s1.normal(s_vertices[i]);
		for (int j = 0; j < 3; j++) {
			iRRAM::cout << "[" << norm[j] << "] ";
		}
		printf("\n");
	}
	
	vector<vector<iRRAM::REAL> *> points4 = {&p2, &p3, &p4};
	vector<Point *> vertices4 = {&v2, &v3, &v4};
	Simplex s4 = Simplex(3, 2, vertices4);
	printf("\nSpecial 3/2 simplex for normal vector testing\n");
	printf("\nNormal vectors of 4th simplex\n");
	printf("Must be\n[-1/sqrt(6)] [ 2/sqrt(6)] [-1/sqrt(6)]\n[         0] [        -1] [         0]\n[ 1/sqrt(2)] [         0] [-1/sqrt(2)]\n");
	printf("Result is\n");
	vector<Point *> s4_vertices = s4.vertex_list();
	for (int i = 0; i < s4_vertices.size(); i++) {
		Vector norm = s4.normal(s4_vertices[i]);
		for (int j = 0; j < 3; j++) {
			iRRAM::cout << "[" << norm[j] << "] ";
		}
		printf("\n");
	}

	printf("\nPrint facets' vertices\n");
	for(int i = 0; i < s_vertices.size(); i++){
		Point *opposite = s_vertices[i];
		Simplex *facet = s1.facet(opposite);
		printf("%dth facet\n", i);
		vector<Point *> f_vertices = facet->vertex_list();
		for(int j = 0; j < f_vertices.size(); j++){
			Point *vertex = f_vertices[j];
			printf("\t%dth vertex : [", j);
			iRRAM::cout << (*vertex->get_coordinate())[0] << " " << (*vertex->get_coordinate())[1] << " " << (*vertex->get_coordinate())[2];
			printf("]\n");
		}
	}
	// Simplex == Operation Test
	// Expected True, False, True, False
	Simplex f1 = s1.facet(s_vertices[0]);
	Simplex f2 = s1.facet(s_vertices[0]);
	Simplex f3 = s1.facet(s_vertices[1]);

	EXPECT_TRUE(f1 == f2);
	EXPECT_FALSE(f1 == f3);
	EXPECT_TRUE(s1 == s2);
	EXPECT_FALSE(s1 == s3);

	printf("\nEnd of Simplex Creation test!\n");
	for (int i = 0; i < 40; i++) printf("-");
	printf("\n");
}
*/
