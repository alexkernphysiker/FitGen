// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <math_h/error.h>
#include <Genetic/equation.h>
#include <Genetic/genetic.h>
#include <Genetic/initialconditions.h>
#include "engine.h"
using namespace std;
using namespace MathTemplates;
using namespace Genetic;
TEST(EquationSystem,empty){
    EquationSystem A{};
    EXPECT_EQ(A({}),0);
}
TEST(EquationSystem,empty2){
    EquationSystem A(list<Equation>{});
    EXPECT_EQ(A({}),0);
}
TEST(EquationSystem,simple){
    EquationSystem A{
	equation([](const ParamSet&P)->double{return P[0];},0)
    };
    EXPECT_EQ(A({0.0}),0);
    EXPECT_EQ(A({0.5}),0.25);
    EXPECT_EQ(A({1.0}),1);
    EXPECT_EQ(A({2.0}),4);
}
TEST(EquationSystem,simple2){
    EquationSystem A(list<Equation>{
	equation([](const ParamSet&P)->double{return P[0];},0)
    });
    EXPECT_EQ(A({0.0}),0);
    EXPECT_EQ(A({0.5}),0.25);
    EXPECT_EQ(A({1.0}),1);
    EXPECT_EQ(A({2.0}),4);
}
TEST(EquationSystem,twoparams){
    EquationSystem A{
	equation([](const ParamSet&P)->double{return P[0]+P[1];},0)
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
TEST(EquationSystem,twoparams2){
    EquationSystem A(list<Equation>{
	equation([](const ParamSet&P)->double{return P[0]+P[1];},0)
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
TEST(EquationSystem,two_eq){
    EquationSystem A{
	equation([](const ParamSet&P)->double{return P[0]+P[1];},0),
	equation([](const ParamSet&P)->double{return P[0]-P[1];},0)
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
TEST(EquationSystem,two_eq2){
    EquationSystem A(list<Equation>{
	equation([](const ParamSet&P)->double{return P[0]+P[1];},0),
	equation([](const ParamSet&P)->double{return P[0]-P[1];},0)
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
TEST(EquationSolver,Integrationtest){
    EquationSolver<DifferentialMutations<>> test{
	equation([](const ParamSet&P)->double{return P[0]+P[1];},0),
	equation([](const ParamSet&P)->double{return P[0]-P[1];},0)
    };
    test.Init(100,make_shared<GenerateUniform>()<<make_pair(-20,20)<<make_pair(-20,20),engine);
    Find(test,engine);
    EXPECT_TRUE(pow(test.Parameters()[0],2)<0.0000001);
    EXPECT_TRUE(pow(test.Parameters()[1],2)<0.0000001);
    for(const auto&eq:test.equations()){
	EXPECT_EQ(0,eq.right);
	EXPECT_TRUE(pow(eq.left(test.Parameters()),2)<0.0000001);
    }
}
