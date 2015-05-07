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
TEST(AbsoluteMutations,Size){
	for(int count=0;count<10;count++){
		TestClass<Crossing<>> gen;
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P);
		ASSERT_TRUE(P.Count()==count);
	}
}
TEST(RelativeMutations,Size){
	for(int count=0;count<10;count++){
		TestClass<Crossing<>> gen;
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P);
		ASSERT_TRUE(P.Count()==count);
	}
}
TEST(ExactCopying,Size){
	for(int count=0;count<10;count++){
		TestClass<Crossing<>> gen;
		ParamSet P=parZeros(count);
		gen.MAKE_TEST(P);
		ASSERT_TRUE(P.Count()==count);
	}
}
