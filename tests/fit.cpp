// this file is distributed under 
// MIT license
#include <math.h>
#include <gtest/gtest.h>
#include <math_h/error.h>
#include <Genetic/fit.h>
#include <Genetic/initialconditions.h>
#include <Genetic/paramfunc.h>
#include "engine.h"
using namespace std;
using namespace MathTemplates;
using namespace Genetic;
TEST(ParameterFunction,Basic){
	int c=0;
	double res=2;
	ParameterFunction F([&c,&res](const ParamSet&,const ParamSet&){c++;return res;});
	EXPECT_EQ(res,F(ParamSet(),ParamSet()));
	EXPECT_EQ(1,c);
}

TEST(point,Base){
	{
		Point P({0},0);
		EXPECT_EQ(1,P.X().size());
	}
	Point P({{1,0.1}},{1,1});
	{
		Point P2(P);
		EXPECT_EQ(P.y().val(),P2.y().val());
		EXPECT_EQ(P.y().uncertainty(),P2.y().uncertainty());
		EXPECT_EQ(P.X().size(),P2.X().size());
		EXPECT_EQ(P.X()[0].val(),P2.X()[0].val());
		EXPECT_EQ(P.X()[0].uncertainty(),P2.X()[0].uncertainty());
	}
}
TEST(FitPoints,sizemismatch){
	EXPECT_THROW(FitPoints()<<Point({1},1)<<Point({1,1},1),Exception<FitPoints>);
	EXPECT_NO_THROW(FitPoints()<<Point({1},1)<<Point({1},1));
}
TEST(FitPoints,simple_min_max){
	FitPoints points;
	points<<Point({1},5)<<Point({3},0)<<Point({0},1);
	EXPECT_EQ(1,points.min().size());
	EXPECT_EQ(1,points.max().size());
	EXPECT_EQ(0,points.min()[0]);
	EXPECT_EQ(3,points.max()[0]);
	EXPECT_EQ(0,points.Ymin());
	EXPECT_EQ(5,points.Ymax());
	EXPECT_EQ(1,points.dimensions());
}
TEST(FitPoints,dimensions){
	EXPECT_THROW(FitPoints()<<Point({},1),Exception<FitPoints::Point>);
	EXPECT_EQ(1,(FitPoints()<<Point({1},1)).dimensions());
	EXPECT_EQ(2,(FitPoints()<<Point({1,1},1)).dimensions());
	EXPECT_EQ(3,(FitPoints()<<Point({1,1,1},1)).dimensions());
	EXPECT_EQ(4,(FitPoints()<<Point({1,1,1,1},1)).dimensions());
}
TEST(FitPoints,Base){
	FitPoints points;
	EXPECT_EQ(0,points.size());
	EXPECT_THROW(points[-1],Exception<FitPoints>);
	EXPECT_THROW(points[0],Exception<FitPoints>);
	EXPECT_EQ(&points,&(points<<FitPoints::Point({0},0)));
	EXPECT_EQ(1,points.size());
	EXPECT_THROW(points[-1],Exception<FitPoints>);
	EXPECT_THROW(points[1],Exception<FitPoints>);
	EXPECT_EQ(0,points[0].X()[0].val());
	EXPECT_EQ(&points,&(points<<FitPoints::Point({1},0)));
	EXPECT_EQ(2,points.size());
	EXPECT_THROW(points[-1],Exception<FitPoints>);
	EXPECT_THROW(points[2],Exception<FitPoints>);
	EXPECT_EQ(0,points[0].X()[0].val());
	EXPECT_EQ(1,points[1].X()[0].val());
	int c=0;
	for(const Point&p:points){
		EXPECT_EQ(1,p.X().size());
		c++;
	}
	EXPECT_EQ(c,points.size());
}
TEST(FitPoints,Operators){
	auto points=make_shared<FitPoints>();
	EXPECT_EQ(points.get(),(points<<FitPoints::Point({0},0)).get());
	EXPECT_EQ(points.get(),(points<<point<double>(1.0,1.0)).get());
}
TEST(FitPoints,Select){
	auto points=make_shared<FitPoints>();
	points<<point<double>(0,0)<<point<double>(1,1)<<point<double>(2,2)<<point<double>(3,3)<<point<double>(3,3)<<point<double>(3,1)<<point<double>(1,3);
	auto filter=[](const ParamSet&X){return X[0]<2.5;};
	auto y_filter=[](double y){return y<2.5;};
	auto sel1=SelectFitPoints(points,make_shared<Filter>(filter));
	for(const Point&p:*sel1){
		auto x=p.x();
		EXPECT_EQ(true,filter(x));
	}
	auto sel2=SelectFitPoints(points,y_filter);
	for(const Point&p:*sel2)EXPECT_EQ(true,y_filter(p.y().val()));
	auto sel3=SelectFitPoints(points,make_shared<Filter>(filter),y_filter);
	for(const Point&p:*sel3){
		auto x=p.x();
		EXPECT_EQ(true,filter(x));
		EXPECT_EQ(true,y_filter(p.y().val()));
	}
}
TEST(OptimalityForPoints,Base){
	auto points=make_shared<FitPoints>();
	int func_calls=0;
	auto f=[&func_calls](const ParamSet&,const ParamSet&){func_calls++;return 0.0;};
	int summand_calls=0;
	auto s=[&summand_calls](const FitPoints::Point&,const ParamSet&,const IParamFunc&){summand_calls++;return 1.0;};
	int coef_calls=0;
	auto c=[&coef_calls](const ParamSet&,const IParamFunc&){coef_calls++;return 1.0;};
	OptimalityForPoints S(points,make_shared<ParameterFunction>(f),c,s);
	EXPECT_EQ(0,S(ParamSet()));
	EXPECT_EQ(0,func_calls);
	EXPECT_EQ(points->size(),summand_calls);
	EXPECT_EQ(1,coef_calls);
	for(int count=1;count<5;count++){
		func_calls=summand_calls=coef_calls=0;
		points<<point<double>(0,0);
		EXPECT_EQ(count,S(ParamSet()));
		EXPECT_EQ(0,func_calls);
		EXPECT_EQ(points->size(),summand_calls);
		EXPECT_EQ(1,coef_calls);
	}
}
template<shared_ptr<OptimalityForPoints> OptimalityAlgorithm(shared_ptr<FitPoints>,shared_ptr<IParamFunc>)>
void test_optimality1(double v=INFINITY){
	auto points=make_shared<FitPoints>()
		<<FitPoints::Point({{0,1}},{0,1})
		<<FitPoints::Point({{1,1}},{0,1})
		<<FitPoints::Point({{2,1}},{0,1});
	auto F=make_shared<ParameterFunction>([](const ParamSet&,const ParamSet&){return 0;});
	auto S=OptimalityAlgorithm(points,F);
	EXPECT_NE(nullptr,S.get());
	EXPECT_EQ(0,S->operator()(ParamSet()));
	auto F1=make_shared<ParameterFunction>([](const ParamSet&,const ParamSet&){return 1;});
	auto S1=OptimalityAlgorithm(points,F1);
	EXPECT_NE(nullptr,S1.get());
	EXPECT_EQ(true,S1->operator()(ParamSet())>0);
	if(isfinite(v))
		EXPECT_EQ(v,S1->operator()(ParamSet()));
}
TEST(OptimalityForPoints,SumSquareDiff){
	test_optimality1<SumSquareDiff>(3);
}
TEST(OptimalityForPoints,ChiSquare){
	test_optimality1<ChiSquare>(3);
}
TEST(OptimalityForPoints,ChiSquareWithXError){
	test_optimality1<ChiSquareWithXError>(3);
}
typedef Add<Mul<Arg<0>,Par<0>>,Par<1>> Fit_Func;
typedef Const<1> Fit_Func_err;
auto Points=make_shared<FitPoints>()<<Point({0},{1,1})<<Point({1},{2,1})<<Point({2},{3,1});
auto Init=make_shared<GenerateUniform>()<<make_pair(0,2)<<make_pair(0,2);
TEST(Fit,Basetest){
	Fit<DifferentialMutations<>,SumSquareDiff> fit(Points,make_shared<Fit_Func>());
	fit.Init(20,Init,engine);
	while(!fit.ConcentratedInOnePoint())
		fit.Iterate(engine);
	EXPECT_TRUE(fit.ParamCount()==2);
	EXPECT_TRUE(fit.PopulationSize()==20);
	EXPECT_TRUE(fit.Optimality()==0);
	EXPECT_TRUE(fit.Optimality(fit.PopulationSize()-1)==0);
	EXPECT_EQ(1,fit.Parameters()[0]);
	EXPECT_EQ(1,fit.Parameters()[1]);
}
TEST(FitFunction,Basetest){
	FitFunction<DifferentialMutations<>,Fit_Func,SumSquareDiff> fit(Points);
	fit.Init(20,Init,engine);
	while(!fit.ConcentratedInOnePoint())
		fit.Iterate(engine);
	EXPECT_TRUE(fit.ParamCount()==2);
	EXPECT_TRUE(fit.PopulationSize()==20);
	EXPECT_TRUE(fit.Optimality()==0);
	EXPECT_TRUE(fit.Optimality(fit.PopulationSize()-1)==0);
	EXPECT_EQ(1,fit.Parameters()[0]);
	EXPECT_EQ(1,fit.Parameters()[1]);
}
