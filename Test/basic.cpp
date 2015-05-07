#include <gtest/gtest.h>
#include <abstract.h>
#include <genetic_exception.h>
#include <math_h/randomfunc.h>
#include "testgen.h"
using namespace Genetic;
using namespace std;
TEST(OptimalityFunction,BaseTest){
	auto f=[](ParamSet&p){return pow(p[0],2);};
	OptimalityFunction F(f);
	for(double x=-1;x<=1;x+=0.1){
		ParamSet X(x);
		EXPECT_EQ(f(X),F(X));
	}
}
class Init:public IInitialConditions{
private:
	function<ParamSet()> m_func;
public:
	Init(function<ParamSet()>f){
		m_func=f;
	}
	virtual ~Init(){}
	virtual ParamSet Generate(){
		return m_func();
	}
};
auto optimality=make_shared<OptimalityFunction>([](ParamSet&p){return pow(p[0],2);});
auto initial=make_shared<Init>([](){return ParamSet(0);});
void test_init(unsigned int threads,int population){
	GeneticTest gen(optimality);
	EXPECT_EQ(optimality.get(),gen.OptimalityCalculator().get());
	ASSERT_THROW(gen.SetThreadCount(0),GeneticException);
	gen.SetThreadCount(threads);
	EXPECT_EQ(threads,gen.ThreadCount());
	EXPECT_EQ(0,gen.PopulationSize());
	ASSERT_THROW(gen.Init(0,initial),GeneticException);
	ASSERT_THROW(gen.ParamCount(),GeneticException);
	EXPECT_EQ(0,gen.PopulationSize());
	ASSERT_NO_THROW(gen.Init(population,initial));
	EXPECT_EQ(population,gen.PopulationSize());
	ASSERT_THROW(gen.Init(population,initial),GeneticException);
	EXPECT_EQ(population,gen.PopulationSize());
	EXPECT_EQ(1,gen.ParamCount());
}
TEST(AbstractGenetic,InitSync){test_init(1,10);}
TEST(AbstractGenetic,InitAsync2){test_init(2,10);}
TEST(AbstractGenetic,InitAsync3){test_init(3,10);}
TEST(AbstractGenetic,InitAsync4){test_init(4,10);}
TEST(AbstractGenetic,InitSync_){test_init(1,5);}
TEST(AbstractGenetic,InitAsync2_){test_init(2,5);}
TEST(AbstractGenetic,InitAsync3_){test_init(3,5);}
TEST(AbstractGenetic,InitAsync4_){test_init(4,5);}
TEST(AbstractGenetic,InitSync__){test_init(1,2);}
TEST(AbstractGenetic,InitAsync2__){test_init(2,2);}
TEST(AbstractGenetic,InitAsync3__){test_init(3,2);}
TEST(AbstractGenetic,InitAsync4__){test_init(4,2);}
TEST(AbstractGenetic,Throwing){
	GeneticTest gen(optimality);
	ASSERT_THROW(gen.Parameters(0),GeneticException);
	ASSERT_THROW(gen.ConcentratedInOnePoint(),GeneticException);
	ASSERT_THROW(gen.AbsoluteOptimalityExitCondition(1),GeneticException);
	ASSERT_THROW(gen.RelativeOptimalityExitCondition(1),GeneticException);
	ASSERT_THROW(gen.Iterate(),GeneticException);
	ASSERT_THROW(gen.Optimality(0),GeneticException);
	ASSERT_NO_THROW(gen.Init(2,initial));
	ASSERT_THROW(gen.AbsoluteOptimalityExitCondition(-1),GeneticException);
	ASSERT_THROW(gen.RelativeOptimalityExitCondition(-1),GeneticException);
	ASSERT_THROW(gen.Optimality(-1),GeneticException);
	ASSERT_NO_THROW(gen.Optimality(0));
	ASSERT_NO_THROW(gen.Optimality(1));
	ASSERT_THROW(gen.Optimality(2),GeneticException);
	ASSERT_THROW(gen[-1],GeneticException);
	ASSERT_NO_THROW(gen[0]);
	ASSERT_THROW(gen[1],GeneticException);
}
void test_iterate(unsigned int threads,int population,unsigned int iterations){
	GeneticTest gen(optimality);
	gen.SetThreadCount(threads);
	gen.Init(population,initial);
	for(unsigned int i=0;i<iterations;i++){
		gen.Iterate();
		EXPECT_EQ(population,gen.PopulationSize());
		EXPECT_EQ(1,gen.ParamCount());
		EXPECT_EQ(i+1,gen.iteration_count());
		EXPECT_EQ(true,gen.ConcentratedInOnePoint());
		EXPECT_EQ(0,gen[0]);
		for(int i=0;i<population;i++)
			EXPECT_EQ(0,gen.Parameters(i)[0]);
		EXPECT_EQ(true,gen.ConcentratedInOnePoint());
		EXPECT_EQ(true,gen.AbsoluteOptimalityExitCondition(0));
		EXPECT_EQ(true,gen.RelativeOptimalityExitCondition(0));
		int c=0;for(auto p:gen)c++;
		EXPECT_EQ(c,gen.ParamCount());
		EXPECT_EQ(0,gen.ParamAverage()[0]);
		EXPECT_EQ(0,gen.ParamDispersion()[0]);
		EXPECT_EQ(0,gen.ParamMaxDeviation()[0]);
	}
}
TEST(AbstractGenetic,ItSync){test_iterate(1,10,100);}
TEST(AbstractGenetic,ItAsync2){test_iterate(2,10,100);}
TEST(AbstractGenetic,ItAsync3){test_iterate(3,10,100);}
TEST(AbstractGenetic,ItAsync4){test_iterate(4,10,100);}
TEST(AbstractGenetic,ItSync_){test_iterate(1,5,5);}
TEST(AbstractGenetic,ItAsync2_){test_iterate(2,5,5);}
TEST(AbstractGenetic,ItAsync3_){test_iterate(3,5,5);}
TEST(AbstractGenetic,ItAsync4_){test_iterate(4,5,5);}
TEST(AbstractGenetic,ItSync__){test_iterate(1,2,5);}
TEST(AbstractGenetic,ItAsync2__){test_iterate(2,2,5);}
TEST(AbstractGenetic,ItAsync3__){test_iterate(3,2,5);}
TEST(AbstractGenetic,ItAsync4__){test_iterate(4,2,5);}
class GeneticTestWithMutations:public AbstractGenetic{
public:
	GeneticTestWithMutations(shared_ptr<IOptimalityFunction> optimality):AbstractGenetic(optimality){}
	virtual ~GeneticTestWithMutations(){}
protected:
	virtual void mutations(ParamSet&P)override{P.Set(0,RandomUniformlyR(0.0,1.0));}
};
TEST(AbstractGenetic,FilterSettingFunc){
	auto initial_uniform=make_shared<Init>([](){return ParamSet(RandomUniformlyR(0.0,1.0));});
	auto filter=[](ParamSet&P){return P[0]>0.5;};
	GeneticTestWithMutations gen(optimality);
	gen.SetFilter(filter);
	gen.Init(100,initial_uniform);
	for(int i=0,n=gen.PopulationSize();i<n;i++)
		EXPECT_EQ(true,filter(gen.Parameters(i)));
	for(int j=0;j<10;j++){
		gen.Iterate();
		for(int i=0,n=gen.PopulationSize();i<n;i++)
			EXPECT_EQ(true,filter(gen.Parameters(i)));
	}
	gen.RemoveFilter();
	for(int j=0;j<10;j++)gen.Iterate();
	for(int i=0,n=gen.PopulationSize();i<n;i++)
		EXPECT_EQ(false,filter(gen.Parameters(i)));
}
TEST(AbstractGenetic,FilterSetting){
	auto initial_uniform=make_shared<Init>([](){return ParamSet(RandomUniformlyR(0.0,1.0));});
	auto filter=[](ParamSet&P){return P[0]>0.5;};
	GeneticTestWithMutations gen(optimality);
	gen.SetFilter(make_shared<Filter>(filter));
	gen.Init(100,initial_uniform);
	for(int i=0,n=gen.PopulationSize();i<n;i++)
		EXPECT_EQ(true,filter(gen.Parameters(i)));
	for(int j=0;j<10;j++){
		gen.Iterate();
		for(int i=0,n=gen.PopulationSize();i<n;i++)
			EXPECT_EQ(true,filter(gen.Parameters(i)));
	}
	gen.RemoveFilter();
	for(int j=0;j<10;j++)gen.Iterate();
	for(int i=0,n=gen.PopulationSize();i<n;i++)
		EXPECT_EQ(false,filter(gen.Parameters(i)));
}
