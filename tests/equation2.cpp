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
TEST(InexactEquationSystem,empty2){
	InexactEquationSystem A(list<InexactEquation>{});
	EXPECT_EQ(A({}),0);
}
TEST(InexactEquationSystem,simple){
	InexactEquationSystem A{
	    {.left=[](const ParamSet&P)->double{return P[0];},.right={0,1}}
	};
	EXPECT_EQ(A({0.0}),0);
	EXPECT_EQ(A({0.5}),0.25);
	EXPECT_EQ(A({1.0}),1);
	EXPECT_EQ(A({2.0}),4);
}
TEST(InexactEquationSystem,simple2){
	InexactEquationSystem A(list<InexactEquation>{
	    {.left=[](const ParamSet&P)->double{return P[0];},.right={0,1}}
	});
	EXPECT_EQ(A({0.0}),0);
	EXPECT_EQ(A({0.5}),0.25);
	EXPECT_EQ(A({1.0}),1);
	EXPECT_EQ(A({2.0}),4);
}
TEST(InexactEquationSystem,twoparams){
	InexactEquationSystem A{
	    {.left=[](const ParamSet&P)->double{return P[0]+P[1];},.right={0,1}}
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
TEST(InexactEquationSystem,twoparams2){
	InexactEquationSystem A(list<InexactEquation>{
	    {.left=[](const ParamSet&P)->double{return P[0]+P[1];},.right={0,1}}
	});
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
	    {.left=[](const ParamSet&P)->double{return P[0]+P[1];},.right={0,1}},
	    {.left=[](const ParamSet&P)->double{return P[0]-P[1];},.right={0,1}}
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
TEST(InexactEquationSystem,two_eq2){
	InexactEquationSystem A(list<InexactEquation>{
	    {.left=[](const ParamSet&P)->double{return P[0]+P[1];},.right={0,1}},
	    {.left=[](const ParamSet&P)->double{return P[0]-P[1];},.right={0,1}}
	});
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
	InexactEquationSolver<DifferentialMutations<>> test{
	    {.left=[](const ParamSet&P)->double{return P[0]+P[1];},.right={0,1}},
	    {.left=[](const ParamSet&P)->double{return P[0]-P[1];},.right={0,1}}
	};
	test.Init(100,make_shared<GenerateUniform>()<<make_pair(-20,20)<<make_pair(-20,20),engine);
	Find(test,engine);
	EXPECT_TRUE(pow(test.Parameters()[0],2)<0.0000001);
	EXPECT_TRUE(pow(test.Parameters()[1],2)<0.0000001);
	for(const auto&eq:test.equations()){
	    EXPECT_EQ(0,eq.right.val());
	    EXPECT_EQ(1,eq.right.uncertainty());
	    EXPECT_TRUE(pow(eq.left(test.Parameters()).val(),2)<0.0000001);
	}
}
