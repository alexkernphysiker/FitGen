#include <gtest/gtest.h>
#include <genetic_exception.h>
using namespace std;
TEST(GeneticException,Test){
	int c=0;
	SetGeneticErrorOutput([&c](string){c++;});
	{string msg="test1";
		GeneticException ex(msg);
		EXPECT_EQ(0,strcmp(msg.c_str(),ex.what()));
	}
	EXPECT_EQ(1,c);
	SetGeneticErrorOutput([](string){});
	{string msg="test2";
		GeneticException ex(msg);
		EXPECT_EQ(0,strcmp(msg.c_str(),ex.what()));
	}
}