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
		point P2;
		P2=P;
		EXPECT_EQ(P.y,P2.y);
		EXPECT_EQ(P.wy,P2.wy);
		EXPECT_EQ(P.X.Count(),P2.X.Count());
		EXPECT_EQ(P.X[0],P2.X[0]);
		EXPECT_EQ(P.WX.Count(),P2.WX.Count());
		EXPECT_EQ(P.WX[0],P2.WX[0]);
	}
}
TEST(FitPoints,Base){
	FitPoints points;
	EXPECT_EQ(0,points.count());
	ASSERT_ANY_THROW(points[-1]);
	ASSERT_ANY_THROW(points[0]);
	point p1;
	p1.X<<0;
	EXPECT_EQ(&points,&(points<<p1));
	EXPECT_EQ(1,points.count());
	ASSERT_ANY_THROW(points[-1]);
	ASSERT_ANY_THROW(points[1]);
	EXPECT_EQ(p1.X[0],points[0].X[0]);
	point p2;
	p2.X<<0;
	EXPECT_EQ(&points,&(points<<p2));
	EXPECT_EQ(2,points.count());
	ASSERT_ANY_THROW(points[-1]);
	ASSERT_ANY_THROW(points[2]);
	EXPECT_EQ(p1.X[0],points[0].X[0]);
	EXPECT_EQ(p2.X[0],points[1].X[0]);
	int c=0;
	for(point&p:points)c++;
	EXPECT_EQ(c,points.count());
}
TEST(FitPoints,Operators){
	auto points=make_shared<FitPoints>();
	point p1;
	EXPECT_EQ(points.get(),(points<<p1).get());
	EXPECT_EQ(points.get(),(points<<make_pair(1.0,1.0)).get());
}
TEST(FitPoints,Select){
	auto points=make_shared<FitPoints>();
	points<<make_pair(0,0)<<make_pair(1,1)<<make_pair(2,2)<<make_pair(3,3)<<make_pair(3,3)<<make_pair(3,1)<<make_pair(1,3);
	auto filter=[](ParamSet&X){return X[0]<2.5;};
	auto y_filter=[](double y){return y<2.5;};
	auto sel1=SelectFitPoints(points,make_shared<Filter>(filter));
	for(point&p:*sel1)EXPECT_EQ(true,filter(p.X));
	auto sel2=SelectFitPoints(points,y_filter);
	for(point&p:*sel2)EXPECT_EQ(true,y_filter(p.y));
	auto sel3=SelectFitPoints(points,make_shared<Filter>(filter),y_filter);
	for(point&p:*sel3){
		EXPECT_EQ(true,filter(p.X));
		EXPECT_EQ(true,y_filter(p.y));
	}
}