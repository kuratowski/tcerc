#include "test_main.h"

#include "subspace.h"
#include "simplex.h"

using namespace TCERC;
using namespace iRRAM;
using namespace std;

TEST(TestSubspace, TestCreate)
{
	vector<REAL> e1 = {1, 0, 0, 0};
	vector<REAL> e2 = {0, 1, 0, 0};
	vector<REAL> e3 = {0, 0, 1, 0};
	vector<Vector> basis = {Vector(e1), Vector(e2), Vector(e3)};
	Subspace sbsp(basis);
	EXPECT_EQ(sbsp.dimension(), 3);
	EXPECT_EQ(sbsp.ambient_dimension(), 4);

//	Vertex p1({1, 0, 0, 0}), p2({0, 1, 0, 0}), p3({0, 0, 1, 0});
//	Simplex s({p1, p2, p3});

//	Subspace s_sbsp(s);
//       std::cout<<"Fuck"<<std::endl;

//	EXPECT_EQ(s_sbsp.dimension(), 2);
//	EXPECT_EQ(s_sbsp.ambient_dimension(), 4);
}

TEST(TestSubspace, TestProject)
{
	vector<REAL> e1 = {1, 2, 4, 0};
	vector<REAL> e2 = {1, 1, 0, 0};
	vector<REAL> e3 = {0, 0, 3, 0};
	
	vector<Vector> basis = {Vector(e1), Vector(e2), Vector(e3)};
	Subspace sbsp(basis);

	vector<REAL> trans = {3.2, 4.8, 1.0};
	vector<REAL> proj_coord = {2, 3, 7, 0};
	Vector proj = Vector(proj_coord);

        iRRAM::REAL ub1("8.0001");
        iRRAM::REAL lb1("7.9999");
        iRRAM::REAL ub2("11.2001");
        iRRAM::REAL lb2("11.1999");
        iRRAM::REAL ub3("15.8001");
        iRRAM::REAL lb3("15.7999");
        iRRAM::REAL ub4("0.0001");
        iRRAM::REAL lb4("-0.0001");
        iRRAM::REAL ub5("1.0001");
        iRRAM::REAL lb5("0.9999");
        iRRAM::REAL ub6("1.0001");
        iRRAM::REAL lb6("0.9999");
        iRRAM::REAL ub7("1.0001");
        iRRAM::REAL lb7("0.9999");

	Vector trans_result = sbsp.embed(trans);
	vector<REAL> proj_result = sbsp.project(proj);
        
	EXPECT_EQ(choose(trans_result[0] < ub1, trans_result[0] > ub1), 1);
       EXPECT_EQ(choose(trans_result[0] < lb1, trans_result[0] > lb1), 2);
       EXPECT_EQ(choose(trans_result[1] < ub2, trans_result[1] > ub2), 1);
       EXPECT_EQ(choose(trans_result[1] < lb2, trans_result[1] > lb2), 2);
       EXPECT_EQ(choose(trans_result[2] < ub3, trans_result[2] > ub3), 1);
       EXPECT_EQ(choose(trans_result[2] < lb3, trans_result[2] > lb3), 2);
       EXPECT_EQ(choose(trans_result[3] < ub4, trans_result[3] > ub4), 1);
       EXPECT_EQ(choose(trans_result[3] < lb4, trans_result[3] > lb4), 2);
       EXPECT_EQ(choose(proj_result[0] < ub5, proj_result[0] > ub5), 1);
       EXPECT_EQ(choose(proj_result[0] < lb5, proj_result[0] > lb5), 2);
       EXPECT_EQ(choose(proj_result[1] < ub6, proj_result[1] > ub6), 1);
       EXPECT_EQ(choose(proj_result[1] < lb6, proj_result[1] > lb6), 2);
       EXPECT_EQ(choose(proj_result[2] < ub7, proj_result[2] > ub7), 1);
       EXPECT_EQ(choose(proj_result[2] < lb7, proj_result[2] > lb7), 2);
}
