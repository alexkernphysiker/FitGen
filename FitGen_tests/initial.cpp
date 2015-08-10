// this file is distributed under 
// GPL v 3.0 license
#include <gtest/gtest.h>
#include <initialconditions.h>
#include <genetic_exception.h>
#include <math_h/functions.h>
using namespace Genetic;
using namespace std;
const int n=5;
RANDOM engine;
shared_ptr<Distrib> Distrs[]={
	make_shared<Distrib>([](double x){return Gaussian<double>(x,0,1);},0,10,10),
	make_shared<Distrib>([](double x){return Gaussian<double>(x,1,1);},0,10,10),
	make_shared<Distrib>([](double x){return Gaussian<double>(x,2,1);},0,10,10),
	make_shared<Distrib>([](double x){return Gaussian<double>(x,3,1);},0,10,10),
	make_shared<Distrib>([](double x){return Gaussian<double>(x,4,1);},0,10,10)
};
TEST(InitialDistributions,Create){
	InitialDistributions I;
	EXPECT_EQ(0,I.Count());
}
TEST(InitialDistributions,Add){
	for(int count=1;count<n;count++){
		InitialDistributions I;
		for(int i=0;i<count;i++)
			EXPECT_EQ(&I,&(I<<Distrs[i]));
		EXPECT_EQ(count,I.Count());
		for(int i=0;i<count;i++)
			EXPECT_EQ(Distrs[i].get(),&I[i]);
		EXPECT_THROW(I[count](engine),GeneticException);
		EXPECT_THROW(I[-1](engine),GeneticException);
	}
	for(int count=1;count<n;count++){
		auto I=make_shared<InitialDistributions>();
		for(int i=0;i<count;i++)
			EXPECT_EQ(I.get(),(I<<Distrs[i]).get());
		EXPECT_EQ(count,I->Count());
		for(int i=0;i<count;i++)
			EXPECT_EQ(Distrs[i].get(),&(I->operator[](i)));
	}
}
TEST(InitialDistributions,Generate){
	for(int count=1;count<5;count++){
		InitialDistributions I;
		for(int i=0;i<count;i++)I<<Distrs[i];
		ParamSet P=I.Generate(engine);
		EXPECT_EQ(count,P.Count());
		for(int i=0;i<count;i++)
			EXPECT_EQ(true,(P[i]>=0)&&(P[i]<=10));
	}
}


TEST(GenerateUniform,Create){
	GenerateUniform I;
	EXPECT_EQ(0,I.Count());
}
TEST(GenerateUniform,Add){
	for(int count=1;count<n;count++){
		GenerateUniform I;
		for(int i=0;i<count;i++)
			EXPECT_EQ(&I,&(I.Add(i,2*n-i)));
		EXPECT_EQ(count,I.Count());
		for(int i=0;i<count;i++){
			EXPECT_EQ(i,I.Min(i));
			EXPECT_EQ(2*n-i,I.Max(i));
		}
		EXPECT_THROW(I.Min(count),GeneticException);
		EXPECT_THROW(I.Max(count),GeneticException);
		EXPECT_THROW(I.Min(-1),GeneticException);
		EXPECT_THROW(I.Max(-1),GeneticException);
	}
	{
		GenerateUniform I;
		EXPECT_THROW(I.Add(2,1),GeneticException);
	}
}
TEST(GenerateUniform,AddSharedPtr){
	for(int count=1;count<n;count++){
		auto I=make_shared<GenerateUniform>();
		for(int i=0;i<count;i++)
			EXPECT_EQ(I.get(),(I<<make_pair<double,double>(i,2*n-i)).get());
		EXPECT_EQ(count,I->Count());
		for(int i=0;i<count;i++){
			EXPECT_EQ(i,I->Min(i));
			EXPECT_EQ(2*n-i,I->Max(i));
		}
	}
}
TEST(GenerateUniform,Generate){
	for(int count=1;count<n;count++){
		GenerateUniform I;
		for(int i=0;i<count;i++)
			I.Add(i,2*n-i);
		for(int i=0;i<100;i++){
			ParamSet P=I.Generate(engine);
			EXPECT_EQ(count,P.Count());
			for(int c=0;c<count;c++)
				EXPECT_TRUE((P[c]>=double(c))&&(P[c]<=double(2*n-c)));
		}
	}
}

TEST(GenerateByGauss,Create){
	GenerateByGauss I;
	EXPECT_EQ(0,I.Count());
}
TEST(GenerateByGauss,Add){
	for(int count=1;count<n;count++){
		GenerateByGauss I;
		for(int i=0;i<count;i++)
			EXPECT_EQ(&I,&(I.Add(i,2*n-i)));
		EXPECT_EQ(count,I.Count());
		for(int i=0;i<count;i++){
			EXPECT_EQ(i,I.Mean(i));
			EXPECT_EQ(2*n-i,I.Sigma(i));
		}
		EXPECT_THROW(I.Mean(count),GeneticException);
		EXPECT_THROW(I.Sigma(count),GeneticException);
		EXPECT_THROW(I.Mean(-1),GeneticException);
		EXPECT_THROW(I.Sigma(-1),GeneticException);
	}
}
TEST(GenerateByGauss,AddSharedPtr){
	for(int count=1;count<n;count++){
		auto I=make_shared<GenerateByGauss>();
		for(int i=0;i<count;i++)
			EXPECT_EQ(I.get(),(I<<make_pair<double,double>(i,2*n-i)).get());
		EXPECT_EQ(count,I->Count());
		for(int i=0;i<count;i++){
			EXPECT_EQ(i,I->Mean(i));
			EXPECT_EQ(2*n-i,I->Sigma(i));
		}
	}
}
TEST(GenerateByGauss,Generate){
	for(int count=1;count<n;count++){
		GenerateByGauss I;
		for(int i=0;i<count;i++)
			I.Add(i,2*n-i);
		for(int i=0;i<100;i++){
			ParamSet P=I.Generate(engine);
			EXPECT_EQ(count,P.Count());
			for(int c=0;c<count;c++)
				EXPECT_TRUE(isfinite(P[c]));
		}
	}
}
