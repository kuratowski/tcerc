#include "test_main.h"
#include "simplex.h"
#include "utils.h"
#include <cstring>

using namespace TCERC;

TEST(TestLinsolver, Create)
{
/*
	std::vector<iRRAM::REAL> val1 = {1.0, 1.0, -1.0};
	std::vector<iRRAM::REAL> val2 = {2.0, -1.0, 3.0};
	std::vector<iRRAM::REAL> val3 = {1.0, 2.0, 1.0};
	std::vector<iRRAM::REAL> ans = {0.0, 9.0, 8.0};
*/
	std::vector<iRRAM::REAL> val1 = {1.0, 3.0};
	std::vector<iRRAM::REAL> val2 = {0.0, 1.0};
	std::vector<iRRAM::REAL> ans = {9.0, 2.0};
	Vector v1(val1);
	Vector v2(val2);
//	Vector v3(val3);
	Vector a(ans);

	std::vector<Vector> matrix = {v1, v2};
	Vector final_ans = linear_solver(matrix, ans);
        iRRAM::REAL ub1("3.0001");
        iRRAM::REAL lb1("2.9999");
        iRRAM::REAL ub2("2.0001");
        iRRAM::REAL lb2("1.9999");
/*        iRRAM::REAL ub3("-0.9999");
        iRRAM::REAL lb3("-1.0001");*/

        EXPECT_EQ(choose(final_ans[0] < ub1, final_ans[0] > ub1), 1);
        EXPECT_EQ(choose(final_ans[0] < lb1, final_ans[0] > lb1), 2);
        
	EXPECT_EQ(choose(final_ans[1] < ub2, final_ans[1] > ub2), 1);
        EXPECT_EQ(choose(final_ans[1] < lb2, final_ans[1] > lb2), 2);
        
//	EXPECT_EQ(choose(final_ans[2] < ub3, final_ans[2] > ub3), 1);
//      EXPECT_EQ(choose(final_ans[2] < lb3, final_ans[2] > lb3), 2);

}

