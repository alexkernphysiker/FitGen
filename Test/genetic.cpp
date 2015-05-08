#include <gtest/gtest.h>
#include <genetic.h>
#include <initialconditions.h>
#include <genetic_exception.h>
using namespace Genetic;
using namespace std;
template<class GENETIC>
class TestClass:public virtual GENETIC{
public:
	TestClass():AbstractGenetic(make_shared<OptimalityFunction>([](ParamSet&){return 0;})),GENETIC(){}
	virtual ~TestClass(){}
	void MAKE_TEST(ParamSet&P){
		GENETIC::mutations(P);
	}
};
TEST(DifferentialMutations, Throws){
	TestClass<DifferentialMutations<>> gen;
	double c=gen.MutationCoefficient();
	ASSERT_THROW(gen.SetMutationCoefficient(-1),GeneticException);
	EXPECT_EQ(c,gen.MutationCoefficient());
	ASSERT_NO_THROW(gen.SetMutationCoefficient(0));
	EXPECT_EQ(0,gen.MutationCoefficient());
	ASSERT_NO_THROW(gen.SetMutationCoefficient(0.5));
	EXPECT_EQ(0.5,gen.MutationCoefficient());
	ASSERT_NO_THROW(gen.SetMutationCoefficient(1));
	EXPECT_EQ(1,gen.MutationCoefficient());
	ASSERT_NO_THROW(gen.SetMutationCoefficient(2));
	EXPECT_EQ(2,gen.MutationCoefficient());
}
TEST(DifferentialMutations,Zeros){
	for(int count=0;count<10;count++){
		TestClass<DifferentialMutations<>> gen;
		auto init=make_shared<Initialiser>();
		for(int i=0; i<count;i++)
			init<<[](){return 0;};
		gen.Init(10,init);
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P);
		ASSERT_TRUE(P.Count()==count);
		for(double p:P)
			ASSERT_TRUE(p==0);
	}
}
TEST(DifferentialMutations,Upper){
	for(int count=0;count<10;count++){
		TestClass<DifferentialMutations<>> gen;
		auto init=make_shared<Initialiser>();
		for(int i=0; i<count;i++)
			init<<[](){return RandomUniformlyR(-1,1);};
		gen.Init(10,init);
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P);
		ASSERT_TRUE(P.Count()==count);
		for(double p:P)
			ASSERT_TRUE(p*p<=gen.MutationCoefficient());
	}
}
TEST(Crossing,Throws){
	TestClass<Crossing<>> gen;
	double c=gen.CrossingProbability();
	ASSERT_THROW(gen.SetCrossingProbability(-1),GeneticException);
	EXPECT_EQ(c,gen.CrossingProbability());
	ASSERT_NO_THROW(gen.SetCrossingProbability(0));
	EXPECT_EQ(0,gen.CrossingProbability());
	ASSERT_NO_THROW(gen.SetCrossingProbability(0.5));
	EXPECT_EQ(0.5,gen.CrossingProbability());
	ASSERT_NO_THROW(gen.SetCrossingProbability(1));
	EXPECT_EQ(1,gen.CrossingProbability());
	ASSERT_THROW(gen.SetCrossingProbability(2),GeneticException);
	EXPECT_EQ(1,gen.CrossingProbability());
}
TEST(Crossing,No){
	for(int count=0;count<10;count++){
		TestClass<Crossing<>> gen;
		auto init=make_shared<Initialiser>();
		for(int i=0; i<count;i++)
			init<<[](){return 1;};
		gen.Init(10,init);
		gen.SetCrossingProbability(0);
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P);
		ASSERT_TRUE(P.Count()==count);
		for(double p:P)
			ASSERT_TRUE(p==0);
	}
}
TEST(Crossing,Yes){
	for(int count=0;count<10;count++){
		TestClass<Crossing<>> gen;
		auto init=make_shared<Initialiser>();
		for(int i=0; i<count;i++)
			init<<[](){return 1;};
		gen.Init(10,init);
		gen.SetCrossingProbability(1);
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P);
		ASSERT_TRUE(P.Count()==count);
		for(double p:P)
			ASSERT_TRUE((p==0)||(p==1));
	}
}
TEST(AbsoluteMutations,Throws){
	TestClass<AbsoluteMutations<>> gen;
	double c=gen.AbsoluteMutationsProbability();
	ASSERT_THROW(gen.SetAbsoluteMutationsProbability(-1),GeneticException);
	EXPECT_EQ(c,gen.AbsoluteMutationsProbability());
	ASSERT_NO_THROW(gen.SetAbsoluteMutationsProbability(0));
	EXPECT_EQ(0,gen.AbsoluteMutationsProbability());
	ASSERT_NO_THROW(gen.SetAbsoluteMutationsProbability(0.5));
	EXPECT_EQ(0.5,gen.AbsoluteMutationsProbability());
	ASSERT_NO_THROW(gen.SetAbsoluteMutationsProbability(1));
	EXPECT_EQ(1,gen.AbsoluteMutationsProbability());
	ASSERT_THROW(gen.SetAbsoluteMutationsProbability(2),GeneticException);
	EXPECT_EQ(1,gen.AbsoluteMutationsProbability());
	EXPECT_EQ(0,gen.AbsoluteMutationCoefficients().Count());
	ASSERT_THROW(gen.SetAbsoluteMutationCoefficients(ParamSet(-1)),GeneticException);
	ASSERT_NO_THROW(gen.SetAbsoluteMutationCoefficients(0));
	EXPECT_EQ(0,gen.AbsoluteMutationCoefficients()[0]);
	ASSERT_NO_THROW(gen.SetAbsoluteMutationCoefficients(0.5));
	EXPECT_EQ(0.5,gen.AbsoluteMutationCoefficients()[0]);
	ASSERT_NO_THROW(gen.SetAbsoluteMutationCoefficients(1));
	EXPECT_EQ(1,gen.AbsoluteMutationCoefficients()[0]);
	ASSERT_NO_THROW(gen.SetAbsoluteMutationCoefficients(2));
	EXPECT_EQ(2,gen.AbsoluteMutationCoefficients()[0]);
}
TEST(AbsoluteMutations,Size){
	for(int count=0;count<10;count++){
		TestClass<Crossing<>> gen;
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P);
		ASSERT_TRUE(P.Count()==count);
	}
}
TEST(RelativeMutations,Throws){
	TestClass<RelativeMutations<>> gen;
	double c=gen.RelativeMutationsProbability();
	ASSERT_THROW(gen.SetRelativeMutationsProbability(-1),GeneticException);
	EXPECT_EQ(c,gen.RelativeMutationsProbability());
	ASSERT_NO_THROW(gen.SetRelativeMutationsProbability(0));
	EXPECT_EQ(0,gen.RelativeMutationsProbability());
	ASSERT_NO_THROW(gen.SetRelativeMutationsProbability(0.5));
	EXPECT_EQ(0.5,gen.RelativeMutationsProbability());
	ASSERT_NO_THROW(gen.SetRelativeMutationsProbability(1));
	EXPECT_EQ(1,gen.RelativeMutationsProbability());
	ASSERT_THROW(gen.SetRelativeMutationsProbability(2),GeneticException);
	EXPECT_EQ(1,gen.RelativeMutationsProbability());
	EXPECT_EQ(0,gen.RelativeMutationCoefficients().Count());
	ASSERT_THROW(gen.SetRelativeMutationCoefficients(ParamSet(-1)),GeneticException);
	ASSERT_NO_THROW(gen.SetRelativeMutationCoefficients(0));
	EXPECT_EQ(0,gen.RelativeMutationCoefficients()[0]);
	ASSERT_NO_THROW(gen.SetRelativeMutationCoefficients(0.5));
	EXPECT_EQ(0.5,gen.RelativeMutationCoefficients()[0]);
	ASSERT_NO_THROW(gen.SetRelativeMutationCoefficients(1));
	EXPECT_EQ(1,gen.RelativeMutationCoefficients()[0]);
	ASSERT_NO_THROW(gen.SetRelativeMutationCoefficients(2));
	EXPECT_EQ(2,gen.RelativeMutationCoefficients()[0]);
}
TEST(RelativeMutations,Size){
	for(int count=0;count<10;count++){
		TestClass<Crossing<>> gen;
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P);
		ASSERT_TRUE(P.Count()==count);
	}
}
TEST(ExactCopying,Throws){
	TestClass<ExactCopying<AbstractGenetic>> gen;
	double c=gen.ExactCopyingProbability();
	ASSERT_THROW(gen.SetExactCopyingProbability(-1),GeneticException);
	EXPECT_EQ(c,gen.ExactCopyingProbability());
	ASSERT_NO_THROW(gen.SetExactCopyingProbability(0));
	EXPECT_EQ(0,gen.ExactCopyingProbability());
	ASSERT_NO_THROW(gen.SetExactCopyingProbability(0.5));
	EXPECT_EQ(0.5,gen.ExactCopyingProbability());
	ASSERT_NO_THROW(gen.SetExactCopyingProbability(1));
	EXPECT_EQ(1,gen.ExactCopyingProbability());
	ASSERT_THROW(gen.SetExactCopyingProbability(2),GeneticException);
	EXPECT_EQ(1,gen.ExactCopyingProbability());
}

TEST(ExactCopying,Size){
	for(int count=0;count<10;count++){
		TestClass<Crossing<>> gen;
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P);
		ASSERT_TRUE(P.Count()==count);
	}
}
