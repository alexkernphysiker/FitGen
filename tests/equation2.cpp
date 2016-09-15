// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <math_h/error.h>
#include <Genetic/equation2.h>
#include <Genetic/initialconditions.h>
#include "engine.h"
using namespace std;
using namespace MathTemplates;
using namespace Genetic;
TEST(InexactEquationSystem,empty){
	InexactEquationSystem A{};
	EXPECT_EQ(A({}),0);
}
TEST(InexactEquationSystem,simple){
	InexactEquationSystem A{
		in_eq([](const ParamSet&P){return P[0];},{0,1})
	};
	EXPECT_EQ(A({0}),0);
	EXPECT_EQ(A({1}),1);
	EXPECT_EQ(A({2}),4);
}
