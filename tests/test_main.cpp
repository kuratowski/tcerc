#include "test_main.h"

int main(int argc, char **argv)
{
	iRRAM_initialize(argc, argv);
	testing::InitGoogleTest(&argc, argv);
	int ret = RUN_ALL_TESTS();
	iRRAM_finalize();
	return ret;
}

