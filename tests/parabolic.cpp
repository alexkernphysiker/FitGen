// this file is distributed under 
// MIT license
#include <math.h>
#include <gtest/gtest.h>
#include <math_h/error.h>
#include <Genetic/parabolic.h>
#include <Genetic/initialconditions.h>
#include <Genetic/paramfunc.h>
#include "engine.h"
using namespace std;
using namespace MathTemplates;
using namespace Genetic;
class ParabolicTest:public virtual ParabolicErrorEstimationFromChisq{
public:
	ParabolicTest(const double test_param=1):AbstractGenetic(make_shared<OptimalityFunction>([test_param](const ParamSet&P){
		double res=0;
		for(int i=0,n=P.size();i<n;i++)res+=pow(P[i],2);
		return res*test_param;
	})),ParabolicErrorEstimationFromChisq(){}
	virtual ~ParabolicTest(){}
};
TEST(ParabolicErrorEstimationFromChisq,Base){
	ParabolicTest gen;
	gen.SetUncertaintyCalcDeltas({0.01})
	.Init(1,make_shared<InitialDistributions>()<<make_shared<RandomValueGenerator<double>>(-0.0001,0.0001),engine);
	EXPECT_TRUE(pow(gen.ParametersWithUncertainties()[0].uncertainty()-1.0,2)<0.0001);
}
TEST(ParabolicErrorEstimationFromChisq,NegDelta){
	ParabolicTest gen;
	gen.SetUncertaintyCalcDeltas({-0.01})
	.Init(1,make_shared<InitialDistributions>()<<make_shared<RandomValueGenerator<double>>(-0.0001,0.0001),engine);
	EXPECT_THROW(gen.ParametersWithUncertainties()[0],Exception<ParabolicErrorEstimationFromChisq>);
}
TEST(ParabolicErrorEstimationFromChisq,NegZero){
	ParabolicTest gen;
	gen.SetUncertaintyCalcDeltas({0})
	.Init(1,make_shared<InitialDistributions>()<<make_shared<RandomValueGenerator<double>>(-0.0001,0.0001),engine);
	EXPECT_THROW(gen.ParametersWithUncertainties()[0],Exception<ParabolicErrorEstimationFromChisq>);
}
TEST(ParabolicErrorEstimationFromChisq,BaseTest){
	for(int count=1;count<10;count++){
		ParabolicTest gen;
		auto init=make_shared<InitialDistributions>();
		for(int i=0;i<count;i++)
			init<<make_shared<RandomValueGenerator<double>>(-0.001,0.001);
		gen.SetUncertaintyCalcDeltas(parEq(count,0.01)).Init(1,init,engine);
		ASSERT_EQ(count,gen.ParametersWithUncertainties().size());
		for(size_t i=0;i<gen.ParamCount();i++){
			EXPECT_EQ(gen.Parameters()[i],gen.ParametersWithUncertainties()[i].val());
			EXPECT_TRUE(pow(gen.ParametersWithUncertainties()[i].uncertainty()-1.0,2)<0.0001);
		}
	}
}
TEST(ParabolicErrorEstimationFromChisq,MoreInteresting){
	for(int count=1;count<10;count++){
		ParabolicTest gen(2);
		auto init=make_shared<InitialDistributions>();
		for(int i=0;i<count;i++)
			init<<make_shared<RandomValueGenerator<double>>(-0.001,0.001);
		gen.SetUncertaintyCalcDeltas(parEq(count,0.01)).Init(1,init,engine);
		ASSERT_EQ(count,gen.ParametersWithUncertainties().size());
		for(size_t i=0;i<gen.ParamCount();i++){
			EXPECT_EQ(gen.Parameters()[i],gen.ParametersWithUncertainties()[i].val());
			EXPECT_TRUE(pow(gen.ParametersWithUncertainties()[i].uncertainty()-sqrt(0.5),2)<0.0001);
		}
	}
	for(int count=1;count<10;count++){
		ParabolicTest gen(0.5);
		auto init=make_shared<InitialDistributions>();
		for(int i=0;i<count;i++)
			init<<make_shared<RandomValueGenerator<double>>(-0.001,0.001);
		gen.SetUncertaintyCalcDeltas(parEq(count,0.01)).Init(1,init,engine);
		ASSERT_EQ(count,gen.ParametersWithUncertainties().size());
		for(size_t i=0;i<gen.ParamCount();i++){
			EXPECT_EQ(gen.Parameters()[i],gen.ParametersWithUncertainties()[i].val());
			EXPECT_TRUE(pow(gen.ParametersWithUncertainties()[i].uncertainty()-sqrt(2.0),2)<0.0001);
		}
	}
}
TEST(ParabolicErrorEstimationFromChisq,Infinite){
	for(int count=1;count<10;count++){
		ParabolicTest gen(0.);
		auto init=make_shared<InitialDistributions>();
		for(int i=0;i<count;i++)
			init<<make_shared<RandomValueGenerator<double>>(-0.001,0.001);
		gen.SetUncertaintyCalcDeltas(parEq(count,0.01)).Init(1,init,engine);
		ASSERT_EQ(count,gen.ParametersWithUncertainties().size());
		for(size_t i=0;i<gen.ParamCount();i++){
			EXPECT_EQ(gen.Parameters()[i],gen.ParametersWithUncertainties()[i].val());
			EXPECT_FALSE(isfinite(gen.ParametersWithUncertainties()[i].uncertainty()));
		}
	}
	for(int count=1;count<10;count++){
		ParabolicTest gen(-1);
		auto init=make_shared<InitialDistributions>();
		for(int i=0;i<count;i++)
			init<<make_shared<RandomValueGenerator<double>>(-0.001,0.001);
		gen.SetUncertaintyCalcDeltas(parEq(count,0.01)).Init(1,init,engine);
		ASSERT_EQ(count,gen.ParametersWithUncertainties().size());
		for(size_t i=0;i<gen.ParamCount();i++){
			EXPECT_EQ(gen.Parameters()[i],gen.ParametersWithUncertainties()[i].val());
			EXPECT_FALSE(isfinite(gen.ParametersWithUncertainties()[i].uncertainty()));
		}
	}
}
