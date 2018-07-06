#include "test_main.h"
#include "utils.h"

using namespace TCERC;
using namespace std;
using namespace iRRAM;

TEST(TestUtils, Determinant)
{
	vector<vector<REAL>> mat = {{1, 9, 4}, {4, 5, 6}, {-3, -6, 7}};
	REAL d = det(mat);
	EXPECT_EQ(choose(d > -379, d < -378.99), 2);
	EXPECT_EQ(choose(d < -379, d > -379.01), 2);
}
