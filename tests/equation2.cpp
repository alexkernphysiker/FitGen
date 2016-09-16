// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <math_h/error.h>
#include <Genetic/equation2.h>
#include <Genetic/genetic.h>
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
		in_eq([](const ParamSet&P)->double{return P[0];},{0,1})
	};
	EXPECT_EQ(A({0.0}),0);
	EXPECT_EQ(A({0.5}),0.25);
	EXPECT_EQ(A({1.0}),1);
	EXPECT_EQ(A({2.0}),4);
}
TEST(InexactEquationSystem,twoparams){
	InexactEquationSystem A{
		in_eq([](const ParamSet&P)->double{return P[0]+P[1];},{0,1})
	};
	EXPECT_EQ(A({0.0,0.0}),0);
	EXPECT_EQ(A({0.5,0.0}),0.25);
	EXPECT_EQ(A({1.0,0.0}),1);
	EXPECT_EQ(A({2.0,0.0}),4);
	EXPECT_EQ(A({0.0,0.0}),0);
	EXPECT_EQ(A({0.0,0.5}),0.25);
	EXPECT_EQ(A({0.0,1.0}),1);
	EXPECT_EQ(A({0.0,2.0}),4);
}
TEST(InexactEquationSystem,two_eq){
	InexactEquationSystem A{
		in_eq([](const ParamSet&P)->double{return P[0]+P[1];},{0,1}),
		in_eq([](const ParamSet&P)->double{return P[0]-P[1];},{0,1})
	};
	EXPECT_EQ(A({0.0,0.0}),0);
	EXPECT_EQ(A({0.5,0.0}),0.5);
	EXPECT_EQ(A({1.0,0.0}),2);
	EXPECT_EQ(A({2.0,0.0}),8);
	EXPECT_EQ(A({0.0,0.0}),0);
	EXPECT_EQ(A({0.0,0.5}),0.5);
	EXPECT_EQ(A({0.0,1.0}),2);
	EXPECT_EQ(A({0.0,2.0}),8);
}
TEST(InexactEquationSolver,Integrationtest){
	InexactEquationSolver<DifferentialMutations<>> eq{
		in_eq([](const ParamSet&P)->double{return P[0]+P[1];},{0,1}),
		in_eq([](const ParamSet&P)->double{return P[0]-P[1];},{0,1})
	};
	eq.Init(100,make_shared<GenerateUniform>()<<make_pair(-20,20)<<make_pair(-20,20),engine);
	Find(eq,engine);
	EXPECT_TRUE(pow(eq[0],2)<0.0000001);
	EXPECT_TRUE(pow(eq[1],2)<0.0000001);
}
