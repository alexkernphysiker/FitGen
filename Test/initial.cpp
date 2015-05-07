#include <gtest/gtest.h>
#include <initialconditions.h>
#include <genetic_exception.h>
#include <math_h/functions.h>
using namespace Genetic;
using namespace std;
const int n=5;
generator Funcs[]={[](){return 0;},[](){return 1;},[](){return 2;},[](){return 3;},[](){return 4;}};
shared_ptr<Distrib> Distrs[]={
	make_shared<Distrib>([](double x){return Gaussian<double>(x,0,1);},0,10,0.1),
	make_shared<Distrib>([](double x){return Gaussian<double>(x,1,1);},0,10,0.1),
	make_shared<Distrib>([](double x){return Gaussian<double>(x,2,1);},0,10,0.1),
	make_shared<Distrib>([](double x){return Gaussian<double>(x,3,1);},0,10,0.1),
	make_shared<Distrib>([](double x){return Gaussian<double>(x,4,1);},0,10,0.1)
};
TEST(Initialiser,Create){
	Initialiser I;
	EXPECT_EQ(0,I.Count());
}
TEST(Initialiser,Add){
	for(int count=1;count<n;count++){
		Initialiser I;
		for(int i=0;i<count;i++)
			EXPECT_EQ(&I,&(I<<Funcs[i]));
		EXPECT_EQ(count,I.Count());
		for(int i=0;i<count;i++)
			EXPECT_EQ(Funcs[i](),I[i]());
		int c=0;
		for(auto i:I)c++;
		EXPECT_EQ(c,count);
		ASSERT_ANY_THROW(I[count]());
		ASSERT_ANY_THROW(I[-1]());
	}
	for(int count=1;count<n;count++){
		auto I=make_shared<Initialiser>();
		for(int i=0;i<count;i++)
			EXPECT_EQ(I.get(),(I<<Funcs[i]).get());
		EXPECT_EQ(count,I->Count());
		for(int i=0;i<count;i++)
			EXPECT_EQ(Funcs[i](),I->operator[](i)());
	}
}
TEST(Initialiser,Generate){
	for(int count=1;count<5;count++){
		Initialiser I;
		for(int i=0;i<count;i++)I<<Funcs[i];
		ParamSet P=I.Generate();
		EXPECT_EQ(count,P.Count());
		for(int i=0;i<count;i++)
			EXPECT_EQ(Funcs[i](),P[i]);
	}
}
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
		ASSERT_ANY_THROW(I[count]());
		ASSERT_ANY_THROW(I[-1]());
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
		ASSERT_ANY_THROW(I.Min(count));
		ASSERT_ANY_THROW(I.Max(count));
		ASSERT_ANY_THROW(I.Min(-1));
		ASSERT_ANY_THROW(I.Max(-1));
	}
	{
		GenerateUniform I;
		ASSERT_ANY_THROW(I.Add(2,1));
	}
}
TEST(GenerateUniform,Generate){
	for(int count=1;count<n;count++){
		GenerateUniform I;
		for(int i=0;i<count;i++)
			I.Add(i,2*n-i);
		for(int i=0;i<100;i++){
			ParamSet P=I.Generate();
			EXPECT_EQ(count,P.Count());
			for(int c=0;c<count;c++)
				ASSERT_TRUE((P[c]>=double(c))&&(P[c]<=double(2*n-c)));
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
		ASSERT_ANY_THROW(I.Mean(count));
		ASSERT_ANY_THROW(I.Sigma(count));
		ASSERT_ANY_THROW(I.Mean(-1));
		ASSERT_ANY_THROW(I.Sigma(-1));
	}
}