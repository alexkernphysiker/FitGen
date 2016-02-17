// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <math_h/error.h>
#include <math_h/hist.h>
#include <Genetic/genetic.h>
#include <Genetic/initialconditions.h>
#include "engine.h"
using namespace std;
using namespace MathTemplates;
using namespace Genetic;
template<class GENETIC>
class TestClass:public virtual GENETIC{
public:
	TestClass():AbstractGenetic(make_shared<OptimalityFunction>([](const ParamSet&){return 0;})),GENETIC(){}
	virtual ~TestClass(){}
	void MAKE_TEST(ParamSet&P,RANDOM&R){
		GENETIC::mutations(P,R);
	}
};
TEST(DifferentialMutations, Throws){
	TestClass<DifferentialMutations<>> gen;
	double c=gen.MutationCoefficient();
	EXPECT_THROW(gen.SetMutationCoefficient(-1),Exception<DifferentialMutations<>>);
	EXPECT_EQ(c,gen.MutationCoefficient());
	EXPECT_NO_THROW(gen.SetMutationCoefficient(0));
	EXPECT_EQ(0,gen.MutationCoefficient());
	EXPECT_NO_THROW(gen.SetMutationCoefficient(0.5));
	EXPECT_EQ(0.5,gen.MutationCoefficient());
	EXPECT_NO_THROW(gen.SetMutationCoefficient(1));
	EXPECT_EQ(1,gen.MutationCoefficient());
	EXPECT_NO_THROW(gen.SetMutationCoefficient(2));
	EXPECT_EQ(2,gen.MutationCoefficient());
}
TEST(DifferentialMutations,Zeros){
	for(size_t count=0;count<10;count++){
		TestClass<DifferentialMutations<>> gen;
		auto init=make_shared<InitialDistributions>();
		for(size_t i=0; i<count;i++)
			init<<make_shared<RandomValueGenerator<double>>(0,0.001);
		gen.Init(10,init,engine);
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P,engine);
		EXPECT_TRUE(P.size()==count);
		for(double p:P)
			EXPECT_TRUE(pow(p,2)<0.0001);
	}
}
TEST(DifferentialMutations,Upper){
	for(size_t count=0;count<5;count++){
		TestClass<DifferentialMutations<>> gen;
		auto init=make_shared<InitialDistributions>();
		for(size_t i=0; i<count;i++)
			init<<make_shared<RandomValueGenerator<double>>(-0.5,0.5);
		gen.Init(5,init,engine);
		for(gen.SetMutationCoefficient(0);
			gen.MutationCoefficient()<=1;
			gen.SetMutationCoefficient(gen.MutationCoefficient()+0.1)
		)for(size_t i=0;i<50;i++){
			ParamSet P=parZeros(count);
			gen.MAKE_TEST(P,engine);
			EXPECT_TRUE(P.size()==count);
			for(double p:P){
				EXPECT_TRUE(p<=gen.MutationCoefficient());
				EXPECT_TRUE((-p)<=gen.MutationCoefficient());
			}
		}
	}
}
TEST(Crossing,Throws){
	TestClass<Crossing<>> gen;
	double c=gen.CrossingProbability();
	EXPECT_THROW(gen.SetCrossingProbability(-1),Exception<Crossing<>>);
	EXPECT_EQ(c,gen.CrossingProbability());
	EXPECT_NO_THROW(gen.SetCrossingProbability(0));
	EXPECT_EQ(0,gen.CrossingProbability());
	EXPECT_NO_THROW(gen.SetCrossingProbability(0.5));
	EXPECT_EQ(0.5,gen.CrossingProbability());
	EXPECT_NO_THROW(gen.SetCrossingProbability(1));
	EXPECT_EQ(1,gen.CrossingProbability());
	EXPECT_THROW(gen.SetCrossingProbability(2),Exception<Crossing<>>);
	EXPECT_EQ(1,gen.CrossingProbability());
}
TEST(Crossing,No){
	for(size_t count=0;count<10;count++){
		TestClass<Crossing<>> gen;
		auto init=make_shared<InitialDistributions>();
		for(size_t i=0; i<count;i++)
			init<<make_shared<RandomValueGenerator<double>>(-0.5,0.5);
		gen.Init(10,init,engine);
		gen.SetCrossingProbability(0);
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P,engine);
		EXPECT_TRUE(P.size()==count);
		for(double p:P)
			EXPECT_TRUE(p==0);
	}
}
TEST(Crossing,Yes){
	for(size_t count=0;count<10;count++){
		TestClass<Crossing<>> gen;
		auto init=make_shared<InitialDistributions>();
		for(size_t i=0; i<count;i++)
			init<<make_shared<Distrib>(0.9,1.0);
		gen.Init(10,init,engine);
		gen.SetCrossingProbability(1);
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P,engine);
		EXPECT_TRUE(P.size()==count);
		for(double p:P)
			EXPECT_TRUE((p==0)||((p>=0.89)&&(p<=1.01)));
	}
}
TEST(AbsoluteMutations,Throws){
	TestClass<AbsoluteMutations<>> gen;
	double c=gen.AbsoluteMutationsProbability();
	EXPECT_THROW(gen.SetAbsoluteMutationsProbability(-1),Exception<AbsoluteMutations<>>);
	EXPECT_EQ(c,gen.AbsoluteMutationsProbability());
	EXPECT_NO_THROW(gen.SetAbsoluteMutationsProbability(0));
	EXPECT_EQ(0,gen.AbsoluteMutationsProbability());
	EXPECT_NO_THROW(gen.SetAbsoluteMutationsProbability(0.5));
	EXPECT_EQ(0.5,gen.AbsoluteMutationsProbability());
	EXPECT_NO_THROW(gen.SetAbsoluteMutationsProbability(1));
	EXPECT_EQ(1,gen.AbsoluteMutationsProbability());
	EXPECT_THROW(gen.SetAbsoluteMutationsProbability(2),Exception<AbsoluteMutations<>>);
	EXPECT_EQ(1,gen.AbsoluteMutationsProbability());
	EXPECT_EQ(0,gen.AbsoluteMutationCoefficients().size());
	EXPECT_THROW(gen.SetAbsoluteMutationCoefficients({-1}),Exception<AbsoluteMutations<>>);
	EXPECT_NO_THROW(gen.SetAbsoluteMutationCoefficients({0}));
	EXPECT_EQ(0,gen.AbsoluteMutationCoefficients()[0]);
	EXPECT_NO_THROW(gen.SetAbsoluteMutationCoefficients({0.5}));
	EXPECT_EQ(0.5,gen.AbsoluteMutationCoefficients()[0]);
	EXPECT_NO_THROW(gen.SetAbsoluteMutationCoefficients({1}));
	EXPECT_EQ(1,gen.AbsoluteMutationCoefficients()[0]);
	EXPECT_NO_THROW(gen.SetAbsoluteMutationCoefficients({2}));
	EXPECT_EQ(2,gen.AbsoluteMutationCoefficients()[0]);
}
TEST(AbsoluteMutations,Size){
	for(size_t count=0;count<10;count++){
		TestClass<Crossing<>> gen;
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P,engine);
		EXPECT_TRUE(P.size()==count);
	}
}
TEST(RelativeMutations,Throws){
	TestClass<RelativeMutations<>> gen;
	double c=gen.RelativeMutationsProbability();
	EXPECT_THROW(gen.SetRelativeMutationsProbability(-1),Exception<RelativeMutations<>>);
	EXPECT_EQ(c,gen.RelativeMutationsProbability());
	EXPECT_NO_THROW(gen.SetRelativeMutationsProbability(0));
	EXPECT_EQ(0,gen.RelativeMutationsProbability());
	EXPECT_NO_THROW(gen.SetRelativeMutationsProbability(0.5));
	EXPECT_EQ(0.5,gen.RelativeMutationsProbability());
	EXPECT_NO_THROW(gen.SetRelativeMutationsProbability(1));
	EXPECT_EQ(1,gen.RelativeMutationsProbability());
	EXPECT_THROW(gen.SetRelativeMutationsProbability(2),Exception<RelativeMutations<>>);
	EXPECT_EQ(1,gen.RelativeMutationsProbability());
	EXPECT_EQ(0,gen.RelativeMutationCoefficients().size());
	EXPECT_THROW(gen.SetRelativeMutationCoefficients({-1}),Exception<RelativeMutations<>>);
	EXPECT_NO_THROW(gen.SetRelativeMutationCoefficients({0}));
	EXPECT_EQ(0,gen.RelativeMutationCoefficients()[0]);
	EXPECT_NO_THROW(gen.SetRelativeMutationCoefficients({0.5}));
	EXPECT_EQ(0.5,gen.RelativeMutationCoefficients()[0]);
	EXPECT_NO_THROW(gen.SetRelativeMutationCoefficients({1}));
	EXPECT_EQ(1,gen.RelativeMutationCoefficients()[0]);
	EXPECT_NO_THROW(gen.SetRelativeMutationCoefficients({2}));
	EXPECT_EQ(2,gen.RelativeMutationCoefficients()[0]);
}
TEST(RelativeMutations,Size){
	for(size_t count=0;count<10;count++){
		TestClass<Crossing<>> gen;
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P,engine);
		EXPECT_TRUE(P.size()==count);
	}
}
class TestMutations: public virtual AbstractGenetic{
public:
	TestMutations():AbstractGenetic(){}
	virtual ~TestMutations(){}
protected:
	virtual void mutations(ParamSet &C,RANDOM&)override{
		C=parOnes(C.size());
	}
};
TEST(ExactCopying,Throws){
	TestClass<ExactCopying<TestMutations>> gen;
	double c=gen.ExactCopyingProbability();
	EXPECT_THROW(gen.SetExactCopyingProbability(-1),Exception<ExactCopying<TestMutations>>);
	EXPECT_EQ(c,gen.ExactCopyingProbability());
	EXPECT_NO_THROW(gen.SetExactCopyingProbability(0));
	EXPECT_EQ(0,gen.ExactCopyingProbability());
	EXPECT_NO_THROW(gen.SetExactCopyingProbability(0.5));
	EXPECT_EQ(0.5,gen.ExactCopyingProbability());
	EXPECT_NO_THROW(gen.SetExactCopyingProbability(1));
	EXPECT_EQ(1,gen.ExactCopyingProbability());
	EXPECT_THROW(gen.SetExactCopyingProbability(2),Exception<ExactCopying<TestMutations>>);
	EXPECT_EQ(1,gen.ExactCopyingProbability());
}
TEST(ExactCopying,Size){
	for(size_t count=0;count<10;count++){
		TestClass<ExactCopying<TestMutations>> gen;
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P,engine);
		EXPECT_TRUE(P.size()==count);
	}
}

TEST(ExactCopying,Check){
		double S=0;
		TestClass<ExactCopying<TestMutations>> gen;
		for(double P=0;P<=1;P+=0.1){
			gen.SetExactCopyingProbability(P);
			Distribution1D<double> D(BinsByCount(2,-0.5,1.5));
			for(int i=0;i<1000;i++){
				ParamSet P{0};
				gen.MAKE_TEST(P,engine);
				D<<P[0];
			}
			double P_exp=D[0].Y().val()/(D[0].Y().val()+D[1].Y().val());
			double dP_exp=sqrt(D[0].Y().val())/(D[0].Y().val()+D[1].Y().val())+sqrt(D[1].Y().val())/pow(D[0].Y().val()+D[1].Y().val(),2);
			S+=pow((P-P_exp)/dP_exp,2);
		}
		S/=11.0;
		EXPECT_TRUE(S<=2.0);
}
