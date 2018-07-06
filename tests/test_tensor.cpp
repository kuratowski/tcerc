#include "test_main.h"

#include "tensor.h"
#include <vector>

using namespace TCERC;
using namespace iRRAM;
using namespace std;

TEST(TestTensor, Contraction)
{
	vector<REAL> data;
	for (int i = 0; i < 16; i++) {
		data.push_back(i);
	}
	vector<REAL> e1 = {1, 0, 0, 0};
	vector<REAL> e2 = {0, 1, 0, 0};
	vector<REAL> e3 = {0, 0, 1, 0};
	vector<REAL> e4 = {0, 0, 0, 1};
	vector<Vector> basis = {Vector(e1), Vector(e2), Vector(e3), Vector(e4)};
	Subspace sbsp(basis);

	EXPECT_EQ(sbsp.dimension(), 4);
	EXPECT_EQ(data.size(), 16);

	Tensor t(sbsp, 1, 1, data);
	t.contraction(0, 1);
	REAL trace = t.coeff({});

	EXPECT_EQ(choose(trace < 30.01, trace > 30), 1);
	EXPECT_EQ(choose(trace > 29.99, trace < 30), 1);


}
TEST(TestTensor, Transpose)
{
	vector<REAL> e1 = {1, 0, 0, 0};
	vector<REAL> e2 = {0, 1, 0, 0};
	vector<REAL> e3 = {0, 0, 1, 0};
	vector<REAL> e4 = {0, 0, 0, 1};
	vector<Vector> basis = {Vector(e1), Vector(e2), Vector(e3), Vector(e4)};
	Subspace sbsp(basis);

	EXPECT_EQ(sbsp.dimension(), 4);
	
	vector<REAL> large_data;
	for (int i = 0; i < 256; i++) {
		large_data.push_back(i);
	}
	Tensor large_t(sbsp, 3, 1, large_data);
	std::vector<int> permutation = {2, 0, 1, 3};
	std::vector<int> test_index = {2, 3, 1, 0};
	EXPECT_EQ(choose(large_t.coeff(test_index) < 180.01, large_t.coeff(test_index) > 180), 1);
	EXPECT_EQ(choose(large_t.coeff(test_index) > 179.99, large_t.coeff(test_index) < 180), 1);

	large_t.transpose(permutation);

	std::vector<int> test_index_permuted = {3, 1, 2, 0};
	EXPECT_EQ(choose(large_t.coeff(test_index_permuted) < 180.01,
			      	large_t.coeff(test_index_permuted) > 180), 1);
	EXPECT_EQ(choose(large_t.coeff(test_index_permuted) > 179.99,
				large_t.coeff(test_index_permuted) < 180), 1);


}

TEST(TestTensor, LiveDemo)
{
	Subspace subspace(2);

	REALMATRIX a(2, 2);
	a(0, 0) = 0.18964764;
	a(0, 1) = 0.78044229;
	a(1, 0) = 0.13703496;
	a(1, 1) = 0.34469586;

	REALMATRIX b(2, 2);
	b(0, 0) = 0.13259441;
	b(0, 1) = 0.60413031;
	b(1, 0) = 0.62412012;
	b(1, 1) = 0.12967686;
	REALMATRIX res = a * b - b * a;

	Tensor t(subspace, res);
	t.contraction(0, 1);
	iRRAM::REAL tr = t.coeff({});

	EXPECT_EQ(choose(tr > 0, tr < 0.0001), 2);
	EXPECT_EQ(choose(tr < 0, tr > -0.0001), 2);
}
