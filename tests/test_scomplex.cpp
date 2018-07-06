/* WARNING: This test contains some bug. We will fix it as soon as possible. */

#include "test_main.h"

#include "scomplex.h"

#include <stdlib.h>
#include <time.h>
using namespace TCERC;
using namespace iRRAM;
using namespace std;


TEST(TestSComplex, Boundary_and_Adj_Simplex){
        std::vector<REAL> origin = {0.0, 0.0, 0.0};
        std::vector<REAL> xplus = {1.0, 0.0, 0.0};
        std::vector<REAL> xminus = {-1.0, 0.0, 0.0};
        std::vector<REAL> yplus = {0.0, 1.0, 0.0};
        std::vector<REAL> yminus = {0.0, -1.0, 0.0};
        std::vector<REAL> zplus = {0.0, 0.0, 1.0};
        std::vector<REAL> zminus = {0.0, 0.0, -1.0};
        Vertex porigin = Vertex(origin);
        Vertex pxplus = Vertex(xplus);
        Vertex pxminus = Vertex(xminus);
        Vertex pyplus = Vertex(yplus);
        Vertex pyminus = Vertex(yminus);
        Vertex pzplus = Vertex(zplus);
        Vertex pzminus = Vertex(zminus);

        std::vector<Vertex> vertices_initial = {porigin, pxplus, pyplus, pzplus};
        Simplex initial = Simplex(vertices_initial);

        SComplex sc = SComplex(initial);
        std::vector<Vertex> vertices4 = {porigin, pxminus, pyplus, pzplus};
        Simplex s4 = Simplex(vertices4);
        sc.add_simplex(s4);

        std::vector<Vertex> vertices2 = {porigin, pxplus, pyminus, pzplus};
        Simplex s2 = Simplex(vertices2);
        sc.add_simplex(s2);

        std::vector<Vertex> vertices1 = {porigin, pxplus, pyplus, pzminus};
        Simplex s1 = Simplex(vertices1);
        sc.add_simplex(s1);

        std::vector<Vertex> vertices6 = {porigin, pxminus, pyminus, pzplus};
        Simplex s6 = Simplex(vertices6);
        sc.add_simplex(s6);

        std::vector<Vertex> vertices5 = {porigin, pxminus, pyplus, pzminus};
        Simplex s5 = Simplex(vertices5);
        sc.add_simplex(s5);

        std::vector<Vertex> vertices3 = {porigin, pxplus, pyminus, pzminus};
        Simplex s3 = Simplex(vertices3);
        sc.add_simplex(s3);

        std::vector<Vertex> vertices7 = {porigin, pxminus, pyminus, pzminus};
        Simplex s7 = Simplex(vertices7);
        sc.add_simplex(s7);

// Boundary Testing

        SComplex actual_bd = sc.boundary();
        SComplex test_bd = SComplex(initial.facet(porigin));

        std::vector<Vertex> v4 = {pxminus, pyplus, pzplus};
        test_bd.add_simplex(v4);

        std::vector<Vertex> v2 = {pxplus, pyminus, pzplus};
        test_bd.add_simplex(v2);

        std::vector<Vertex> v1 = {pxplus, pyplus, pzminus};
        test_bd.add_simplex(v1);

        std::vector<Vertex> v6 = {pxminus, pyminus, pzplus};
        test_bd.add_simplex(v6);

        std::vector<Vertex> v5 = {pxminus, pyplus, pzminus};
        test_bd.add_simplex(v5);

        std::vector<Vertex> v3 = {pxplus, pyminus, pzminus};
        test_bd.add_simplex(v3);

        std::vector<Vertex> v7 = {pxminus, pyminus, pzminus};
        test_bd.add_simplex(v7);
        EXPECT_TRUE(actual_bd==test_bd);

// Adjacent Simplex testing
        SComplex actual_adj = SComplex(vertices4);
        actual_adj.add_simplex(s2);
        actual_adj.add_simplex(s1);
        std::vector<Simplex> adj_list = sc.adj_simplex_list(initial);
//      TCERC::Vertex p1=initial.vertex_list()[0];
//      TCERC::Vertex p2=initial.vertex_list()[1];
//      TCERC::Vertex p3=initial.vertex_list()[2];
//      TCERC::Vertex p4=initial.vertex_list()[3];
//     	iRRAM::cout << "origin Point " << p1[0] << " " << p1[1] << " " << p1[2] << endl;
//     	iRRAM::cout << "origin Point " << p2[0] << " " << p2[1] << " " << p2[2] << endl;
//     	iRRAM::cout << "origin Point " << p3[0] << " " << p3[1] << " " << p3[2] << endl;
// 	iRRAM::cout << "origin Point " << p4[0] << " " << p4[1] << " " << p4[2] << endl << endl;
        SComplex test_adj = SComplex(3, 3);
        for (int i=0; i < adj_list.size(); i++) {
        	test_adj.add_simplex(adj_list[i]);
//        	p1=adj_list[i].vertex_list()[0];
//        	p2=adj_list[i].vertex_list()[1];
//        	p3=adj_list[i].vertex_list()[2];
//     		iRRAM::cout << "Locating Point " << p1[0] << " " << p1[1] << " " << p1[2] << endl;
//     		iRRAM::cout << "Locating Point " << p2[0] << " " << p2[1] << " " << p2[2] << endl;
//     		iRRAM::cout << "Locating Point " << p3[0] << " " << p3[1] << " " << p3[2] << endl;
//     		iRRAM::cout << "Locating Point " << p4[0] << " " << p4[1] << " " << p4[2] << endl << endl;
        }
        //iRRAM::cout << "Locating Point " << x << " " << y << " " << z << endl;
//        std::cout << actual_adj.num_simplex() << std::endl;
//        std::cout << test_adj.num_simplex() << std::endl;
        EXPECT_TRUE(actual_adj == test_adj);
}

TEST(TestSComplex, Create_and_Point_Locator)
{
        std::vector<REAL> origin = {0.0, 0.0, 0.0};
        std::vector<REAL> xplus = {1.0, 0.0, 0.0};
        std::vector<REAL> xminus = {-1.0, 0.0, 0.0};
        std::vector<REAL> yplus = {0.0, 1.0, 0.0};
        std::vector<REAL> yminus = {0.0, -1.0, 0.0};
        std::vector<REAL> zplus = {0.0, 0.0, 1.0};
        std::vector<REAL> zminus = {0.0, 0.0, -1.0};
        Vertex porigin = Vertex(origin);
        Vertex pxplus = Vertex(xplus);
        Vertex pxminus = Vertex(xminus);
        Vertex pyplus = Vertex(yplus);
        Vertex pyminus = Vertex(yminus);
        Vertex pzplus = Vertex(zplus);
        Vertex pzminus = Vertex(zminus);

        std::vector<Vertex> vertices_initial = {porigin, pxplus, pyplus, pzplus};
        Simplex initial = Simplex(vertices_initial);

        SComplex sc = SComplex(initial);
        std::vector<Vertex> vertices1 = {porigin, pxplus, pyplus, pzminus};
        sc.add_simplex(vertices1);
        
        std::vector<Vertex> vertices2 = {porigin, pxplus, pyminus, pzplus};
        sc.add_simplex(vertices2);
        
        std::vector<Vertex> vertices3 = {porigin, pxplus, pyminus, pzminus};
        sc.add_simplex(vertices3);
        
        std::vector<Vertex> vertices4 = {porigin, pxminus, pyplus, pzplus};
        sc.add_simplex(vertices4);
        
        std::vector<Vertex> vertices5 = {porigin, pxminus, pyplus, pzminus};
        sc.add_simplex(vertices5);
        
        std::vector<Vertex> vertices6 = {porigin, pxminus, pyminus, pzplus};
        sc.add_simplex(vertices6);
        
        std::vector<Vertex> vertices7 = {porigin, pxminus, pyminus, pzminus};
        sc.add_simplex(vertices7);
        int n = 50;
        for (int i = 0; i < 20; i++) {
                time_t t;
                srand(0x4D12EC28 * i * i * i);
                REAL x = ((float)rand()) / RAND_MAX * 0.9;
                REAL y = ((float)rand()) / RAND_MAX * 0.9 * (0.9 - x);
                REAL z = ((float)rand()) / RAND_MAX * 0.9 * (0.9 - x - y);
		x = max(x, REAL(0.001));
		y = max(y, REAL(0.001));
		z = max(z, REAL(0.001));

                int location = 0;
                if (rand() % 2) {
                        x = -x;
                        location += 4;
                }
                if (rand() % 2) {
                        y = -y;
                        location += 2;
                }
                if (rand() % 2) {
                        z = -z;
                        location += 1;
                }
		//printf("Locating Point [[%lf _ %lf _ %f]] += %lf\n", x.as_double(), y.as_double(), z.as_double(), (x+y+z).as_double());
                //iRRAM::cout << "Locating Point " << x << " " << y << " " << z << " += " << x + y + z << endl;
		//printf("Expecting %d\n", location);
		//iRRAM::cout << "Expecting " << location << endl;
                std::vector<REAL> test_locator = {x, y, z};
                Point test_point = Point(test_locator);

                Simplex location_result = sc.point_locator(test_point, n);
		EXPECT_TRUE(location_result.is_member_vertex((location & 0x4) ? pxminus : pxplus));
		//printf("Expecting %d %s\n", location, (location & 0x2) ? "y minus" : "y plus");
                EXPECT_TRUE(location_result.is_member_vertex((location & 0x2) ? pyminus : pyplus));
                EXPECT_TRUE(location_result.is_member_vertex((location & 0x1) ? pzminus : pzplus));
        }

        REAL small = 0.8;
        for (int i = 0; i < n; i++) {
                small /= 2.0;
        }
        REAL x = (float)1.0 + small;
        std::vector<REAL> test_outside_locator = {x, REAL(0.0), REAL(0.0)};
        Point test_outside_point = Point(test_outside_locator);
        Simplex location_result2 = sc.point_locator(test_outside_point, n);
        EXPECT_TRUE(location_result2.is_member_vertex(pxplus));
        EXPECT_TRUE(location_result2.is_member_vertex(pyplus));
        EXPECT_TRUE(location_result2.is_member_vertex(pzplus));
				
}

TEST(TestSComplex, SpecialSimplex)
{

	int n = 10;

	std::vector<REAL> v1 = {-1.0, -1.0, -1.0};
	std::vector<REAL> v2 = {1.0, -1.0, -1.0};
	std::vector<REAL> v3 = {0.0, 1.0, -1.0};
	std::vector<REAL> v4 = {0.0, 0.0, 1.0};

	Vertex p1 = Vertex(v1);
	Vertex p2 = Vertex(v2);
	Vertex p3 = Vertex(v3);
	Vertex p4 = Vertex(v4);
	
	std::vector<Vertex> points = {p1, p2, p3, p4};

	Simplex special_Simplex = Simplex(points);

	SComplex single_SComplex = SComplex(special_Simplex);
	
	std::vector<REAL> test_v = {0.1, 0.1, -0.5};
	Point test_p = Point(test_v);

	Simplex location_result = single_SComplex.point_locator(test_p, n);

}
