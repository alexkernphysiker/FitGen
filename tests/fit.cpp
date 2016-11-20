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
		EXPECT_EQ(P.var_y().val(),P2.y().val());
		EXPECT_EQ(P.y().val(),P2.var_y().val());
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
	EXPECT_THROW(points.dimensions(),Exception<FitPoints>);
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
TEST(FitPoints,GetCurves){
    FitPoints p;
    p<<Point({1,2},3)<<Point({2,3},4)<<Point({3,4},5);
    const auto H1=p.Hist1(1);
    const auto H2=p.Hist1(1,0);
    EXPECT_EQ(3,H1.size());
    EXPECT_EQ(2,H1[0].X().val());
    EXPECT_EQ(3,H1[1].X().val());
    EXPECT_EQ(4,H1[2].X().val());
    EXPECT_EQ(3,H1[0].Y().val());
    EXPECT_EQ(4,H1[1].Y().val());
    EXPECT_EQ(5,H1[2].Y().val());
    EXPECT_EQ(3,H2.size());
    EXPECT_EQ(2,H2[0].X().val());
    EXPECT_EQ(3,H2[1].X().val());
    EXPECT_EQ(4,H2[2].X().val());
    EXPECT_EQ(1,H2[0].Y().val());
    EXPECT_EQ(2,H2[1].Y().val());
    EXPECT_EQ(3,H2[2].Y().val());
    const auto L1=p.Line(1);
    const auto L2=p.Line(1,0);
    EXPECT_EQ(3,L1.size());
    EXPECT_EQ(2,L1[0].X());
    EXPECT_EQ(3,L1[1].X());
    EXPECT_EQ(4,L1[2].X());
    EXPECT_EQ(3,L1[0].Y());
    EXPECT_EQ(4,L1[1].Y());
    EXPECT_EQ(5,L1[2].Y());
    EXPECT_EQ(3,L2.size());
    EXPECT_EQ(2,L2[0].X());
    EXPECT_EQ(3,L2[1].X());
    EXPECT_EQ(4,L2[2].X());
    EXPECT_EQ(1,L2[0].Y());
    EXPECT_EQ(2,L2[1].Y());
    EXPECT_EQ(3,L2[2].Y());
    FitPoints P1(H1),P2(H2),P3(L1),P4(L2);
    EXPECT_EQ(3,P1.size());
    EXPECT_EQ(2,P1[0].X()[0].val());
    EXPECT_EQ(3,P1[1].X()[0].val());
    EXPECT_EQ(4,P1[2].X()[0].val());
    EXPECT_EQ(3,P2.size());
    EXPECT_EQ(2,P2[0].X()[0].val());
    EXPECT_EQ(3,P2[1].X()[0].val());
    EXPECT_EQ(4,P2[2].X()[0].val());
    EXPECT_EQ(3,P3.size());
    EXPECT_EQ(2,P3[0].X()[0].val());
    EXPECT_EQ(3,P3[1].X()[0].val());
    EXPECT_EQ(4,P3[2].X()[0].val());
    EXPECT_EQ(3,P4.size());
    EXPECT_EQ(2,P4[0].X()[0].val());
    EXPECT_EQ(3,P4[1].X()[0].val());
    EXPECT_EQ(4,P4[2].X()[0].val());
}
TEST(FitPoints,Operators){
	auto points=make_shared<FitPoints>();
	EXPECT_EQ(points.get(),(points<<FitPoints::Point({0},0)).get());
	EXPECT_EQ(points.get(),(points<<point<double>(1.0,1.0)).get());
	FitPoints p; p<<FitPoints::Point({2},2);
	EXPECT_EQ(points.get(),(points<<p).get());
	EXPECT_EQ(3,points->size());
	EXPECT_EQ(1,points->dimensions());
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
	EXPECT_EQ(S.Points(),points);
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
