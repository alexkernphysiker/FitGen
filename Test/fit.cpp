#include <gtest/gtest.h>
#include <fit.h>
using namespace Genetic;
using namespace std;
TEST(ParameterFunction,Basic){
	int c=0;
	double res=2;
	ParameterFunction F([&c,&res](ParamSet&,ParamSet&){c++;return res;});
	ParamSet X,P;
	EXPECT_EQ(res,F(X,P));
	EXPECT_EQ(1,c);
}