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
	eq([](const ParamSet&P)->double{return P[0];},0)
    };
    EXPECT_EQ(A({0.0}),0);
    EXPECT_EQ(A({0.5}),0.25);
    EXPECT_EQ(A({1.0}),1);
    EXPECT_EQ(A({2.0}),4);
}
TEST(EquationSystem,simple2){
    EquationSystem A(list<Equation>{
	eq([](const ParamSet&P)->double{return P[0];},0)
    });
    EXPECT_EQ(A({0.0}),0);
    EXPECT_EQ(A({0.5}),0.25);
    EXPECT_EQ(A({1.0}),1);
    EXPECT_EQ(A({2.0}),4);
}
TEST(EquationSystem,twoparams){
    EquationSystem A{
	eq([](const ParamSet&P)->double{return P[0]+P[1];},0)
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
	eq([](const ParamSet&P)->double{return P[0]+P[1];},0)
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
	eq([](const ParamSet&P)->double{return P[0]+P[1];},0),
	eq([](const ParamSet&P)->double{return P[0]-P[1];},0)
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
	eq([](const ParamSet&P)->double{return P[0]+P[1];},0),
	eq([](const ParamSet&P)->double{return P[0]-P[1];},0)
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
    EquationSolver<DifferentialMutations<>> equation{
	eq([](const ParamSet&P)->double{return P[0]+P[1];},0),
	eq([](const ParamSet&P)->double{return P[0]-P[1];},0)
    };
    equation.Init(100,make_shared<GenerateUniform>()<<make_pair(-20,20)<<make_pair(-20,20),engine);
    Find(equation,engine);
    EXPECT_TRUE(pow(equation[0],2)<0.0000001);
    EXPECT_TRUE(pow(equation[1],2)<0.0000001);
}
