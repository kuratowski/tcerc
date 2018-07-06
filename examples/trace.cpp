#include <iRRAM.h>
#include "tensor.h"

using namespace iRRAM;
using namespace std;
using namespace TCERC;

void compute()
{
	iRRAM::cout << setRwidth(30);

	Subspace sbsp(2);
	REALMATRIX a(2, 2), b(2, 2);

	iRRAM::cout << "Input the entries of A." << endl;
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			iRRAM::cin >> a(i, j);

	iRRAM::cout << "Input the entries of B." << endl;
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			iRRAM::cin >> b(i, j);
	Tensor t_a(sbsp, a);
	Tensor t_b(sbsp, b);

	// Tensor product and contraction is equivalent to matrix product
	Tensor t_ab = t_a * t_b;
	t_ab.contraction(1, 2);
	Tensor t_ba = t_b * t_a;
	t_ba.contraction(1, 2);
	Tensor t = t_ab - t_ba;

	iRRAM::cout << "A" << endl;
	iRRAM::cout << "[[ " << a(0, 0) << ", " << a(0, 1) << "]," << endl;
	iRRAM::cout << " [ " << a(1, 0) << ", " << a(1, 1) << "]]" << endl;
	iRRAM::cout << endl;

	iRRAM::cout << "B" << endl;
	iRRAM::cout << "[[ " << b(0, 0) << ", " << b(0, 1) << "]," << endl;
	iRRAM::cout << " [ " << b(1, 0) << ", " << b(1, 1) << "]]" << endl;
	iRRAM::cout << endl;

	iRRAM::cout << "AB - BA" << endl;
	iRRAM::cout << "[[ " << t.coeff({0, 0}) << ", " << t.coeff({0, 1}) << "]," << endl;
	iRRAM::cout << " [ " << t.coeff({1, 0}) << ", " << t.coeff({1, 1}) << "]]" << endl;
	iRRAM::cout << endl;

	iRRAM::cout << "tr(AB - BA) = " << t.contraction(0, 1).coeff({}) << endl;
}
