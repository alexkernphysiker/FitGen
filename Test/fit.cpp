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
typedef FitPoints::Point point;
TEST(point,Base){
	point P;
	EXPECT_EQ(0,P.X.Count());
	EXPECT_EQ(0,P.WX.Count());
	P.y=1;
	P.wy=1;
	P.X<<1;
	P.WX<<0.1;
	{
		point P2(P);
		EXPECT_EQ(P.y,P2.y);
		EXPECT_EQ(P.wy,P2.wy);
		EXPECT_EQ(P.X.Count(),P2.X.Count());
		EXPECT_EQ(P.X[0],P2.X[0]);
		EXPECT_EQ(P.WX.Count(),P2.WX.Count());
		EXPECT_EQ(P.WX[0],P2.WX[0]);
	}
	{
		point P2=P;
		EXPECT_EQ(P.y,P2.y);
		EXPECT_EQ(P.wy,P2.wy);
		EXPECT_EQ(P.X.Count(),P2.X.Count());
		EXPECT_EQ(P.X[0],P2.X[0]);
		EXPECT_EQ(P.WX.Count(),P2.WX.Count());
		EXPECT_EQ(P.WX[0],P2.WX[0]);
	}
}

