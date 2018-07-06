#include "test_main.h"

TEST(TestChoose, RealCompare)
{
	iRRAM::REAL a(3), b(4);
	EXPECT_EQ(choose(a < b, b < a), 1);
	EXPECT_EQ(choose(b < a, a < b), 2);
}

TEST(TestChoose, LogisticMap)
{
	iRRAM::REAL x("0.5");
	iRRAM::REAL r(15);
	r /= 4;
	for (int i = 0; i < 100; i++) {
		x = r * x * (iRRAM::REAL(1) - x);
	}
	iRRAM::REAL ub("0.88830");
	iRRAM::REAL lb("0.88829");
	EXPECT_EQ(choose(x < ub, x > ub), 1);
	EXPECT_EQ(choose(x < lb, x > lb), 2);
}

TEST(TestChoose, SmallDivision)
{
	iRRAM::REAL a(1), d(1), ten(10);
	for (int i = 0; i < 30; i++) {
		d *= ten;
	}
	for (int i = 0; i < 20; i++) {
		a /= ten;
	}
	for (int i = 0; i < 20; i++) {
		a *= ten;
	}
	iRRAM::REAL ub("1.0001");
	iRRAM::REAL lb("0.9999");
	EXPECT_EQ(choose(a < ub, a > ub), 1);
	EXPECT_EQ(choose(a < lb, a > lb), 2);
}
