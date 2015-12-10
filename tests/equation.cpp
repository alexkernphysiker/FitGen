// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <equation.h>
#include <initialconditions.h>
#include <math_h/exception_math_h.h>
#include <filter.h>
#include "engine.h"
using namespace Genetic;
using namespace std;
double testfunc(const ParamSet&X){return X[0]*X[0]-X[0];}
double zero(const ParamSet&){return 0;}
auto init=make_shared<GenerateUniform>()<<make_pair(0.5,1);
auto filter=make_shared<And>()<<(make_shared<Above>()<<0.4)<<(make_shared<Below>()<<1.5);
TEST(Equation,OneArg){
	Equation<DifferentialMutations<>> solve(testfunc);
	solve.SetFilter(filter);
	solve.Init(10,init,engine);
	while(!solve.AbsoluteOptimalityExitCondition(0.00001))
		solve.Iterate(engine);
	EXPECT_TRUE(pow(solve[0]-1,2)<0.001);
}
TEST(Equation,TwoArgs){
	Equation<DifferentialMutations<>> solve(testfunc,zero);
	solve.SetFilter(filter);
	solve.Init(10,init,engine);
	while(!solve.AbsoluteOptimalityExitCondition(0.00001))
		solve.Iterate(engine);
	EXPECT_TRUE(pow(solve[0]-1,2)<0.001);
}
TEST(SearchMin,OneArg){
	SearchMin<DifferentialMutations<>> solve(testfunc);
	solve.SetFilter(filter);
	solve.Init(10,init,engine);
	while(!solve.AbsoluteOptimalityExitCondition(0.00001))
		solve.Iterate(engine);
	EXPECT_TRUE(pow(solve[0]-0.5,2)<0.001);
}
