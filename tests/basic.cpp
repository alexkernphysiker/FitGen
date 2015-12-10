// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <abstract.h>
#include <math_h/exception_math_h.h>
#include <math_h/randomfunc.h>
#include <initialconditions.h>
#include "engine.h"
using namespace Genetic;
using namespace std;
TEST(OptimalityFunction,BaseTest){
	auto f=[](const ParamSet&p){return pow(p[0],2);};
	OptimalityFunction F(f);
	for(double x=-1;x<=1;x+=0.1){
		EXPECT_EQ(f(ParamSet(x)),F(ParamSet(x)));
	}
}
class GeneticTest:public Genetic::AbstractGenetic{
public:
	GeneticTest(std::shared_ptr<Genetic::IOptimalityFunction> optimality):AbstractGenetic(optimality){}
	virtual ~GeneticTest(){}
};
auto optimality=make_shared<OptimalityFunction>([](const ParamSet&p){return pow(p[0],2);});
auto initial=make_shared<InitialDistributions>()<<make_shared<RandomValueGenerator<double>>(0,0.001);
void test_init(size_t threads,size_t population){
	GeneticTest gen(optimality);
	EXPECT_EQ(optimality.get(),gen.OptimalityCalculator().get());
	EXPECT_THROW(gen.SetThreadCount(0),math_h_error<AbstractGenetic>);
	gen.SetThreadCount(threads);
	EXPECT_EQ(threads,gen.ThreadCount());
	EXPECT_EQ(0,gen.PopulationSize());
	EXPECT_THROW(gen.Init(0,initial,engine),math_h_error<AbstractGenetic>);
	EXPECT_THROW(gen.ParamCount(),math_h_error<AbstractGenetic>);
	EXPECT_EQ(0,gen.PopulationSize());
	EXPECT_NO_THROW(gen.Init(population,initial,engine));
	EXPECT_EQ(population,gen.PopulationSize());
	EXPECT_THROW(gen.Init(population,initial,engine),math_h_error<AbstractGenetic>);
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
	EXPECT_THROW(gen.Parameters(0),math_h_error<AbstractGenetic>);
	EXPECT_THROW(gen.ConcentratedInOnePoint(),math_h_error<AbstractGenetic>);
	EXPECT_THROW(gen.AbsoluteOptimalityExitCondition(1),math_h_error<AbstractGenetic>);
	EXPECT_THROW(gen.RelativeOptimalityExitCondition(1),math_h_error<AbstractGenetic>);
	EXPECT_THROW(gen.Iterate(engine),math_h_error<AbstractGenetic>);
	EXPECT_THROW(gen.ParametersDispersionExitCondition(ParamSet(0)),math_h_error<AbstractGenetic>);
	EXPECT_THROW(gen.RelativeParametersDispersionExitCondition(ParamSet(0)),math_h_error<AbstractGenetic>);
	EXPECT_THROW(gen.Optimality(0),math_h_error<AbstractGenetic>);
	EXPECT_NO_THROW(gen.Init(2,initial,engine));
	EXPECT_THROW(gen.AbsoluteOptimalityExitCondition(-1),math_h_error<AbstractGenetic>);
	EXPECT_THROW(gen.RelativeOptimalityExitCondition(-1),math_h_error<AbstractGenetic>);
	EXPECT_FALSE(gen.ParametersDispersionExitCondition(parZeros(1)));
	EXPECT_FALSE(gen.ParametersDispersionExitCondition(parOnes(1)));
	EXPECT_FALSE(gen.RelativeParametersDispersionExitCondition(parZeros(1)));
	EXPECT_FALSE(gen.RelativeParametersDispersionExitCondition(parOnes(1)));
	EXPECT_THROW(gen.Optimality(-1),math_h_error<AbstractGenetic>);
	EXPECT_NO_THROW(gen.Optimality(0));
	EXPECT_NO_THROW(gen.Optimality(1));
	EXPECT_THROW(gen.Optimality(2),math_h_error<AbstractGenetic>);
	EXPECT_THROW(gen[-1],math_h_error<AbstractGenetic>);
	EXPECT_NO_THROW(gen[0]);
	EXPECT_THROW(gen[1],math_h_error<AbstractGenetic>);
	gen.Iterate(engine);
	EXPECT_THROW(gen.ParametersDispersionExitCondition(parEq(1,-1)),math_h_error<AbstractGenetic>);
	EXPECT_THROW(gen.RelativeParametersDispersionExitCondition(parEq(1,-1)),math_h_error<AbstractGenetic>);
}
#define EXPECT_CLOSE(A,B) EXPECT_TRUE(pow((A)-(B),2)<0.0001);
void test_iterate(size_t threads,size_t population,size_t iterations){
	GeneticTest gen(optimality);
	gen.SetThreadCount(threads);
	gen.Init(population,initial,engine);
	EXPECT_FALSE(gen.ParametersDispersionExitCondition(parOnes(1)));
	EXPECT_FALSE(gen.RelativeParametersDispersionExitCondition(parOnes(1)));
	for(size_t i=0;i<iterations;i++){
		gen.Iterate(engine);
		EXPECT_EQ(population,gen.PopulationSize());
		EXPECT_EQ(1,gen.ParamCount());
		EXPECT_EQ(i+1,gen.iteration_count());
		EXPECT_CLOSE(gen[0],0);
		for(size_t i=0;i<population;i++){
			EXPECT_CLOSE(0,gen.Parameters(i)[0]);
			if(i>0)
				EXPECT_TRUE(gen.Optimality(i-1)<=gen.Optimality(i));
		}
		EXPECT_EQ(true,gen.AbsoluteOptimalityExitCondition(0.001));
		size_t c=0;for(auto p:gen)c++;
		EXPECT_EQ(c,gen.ParamCount());
		EXPECT_CLOSE(0,gen.ParamAverage()[0]);
		EXPECT_CLOSE(0,gen.ParamDispersion()[0]);
		EXPECT_CLOSE(0,gen.ParamMaxDeviation()[0]);
		EXPECT_TRUE(gen.ParametersDispersionExitCondition(parEq(c,INFINITY)));
		EXPECT_TRUE(gen.ParametersDispersionExitCondition(parOnes(c)));
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
private:
	RandomValueGenerator<double> G;
public:
	GeneticTestWithMutations(shared_ptr<IOptimalityFunction> optimality):AbstractGenetic(optimality),G(-1,1){}
	virtual ~GeneticTestWithMutations(){}
protected:
	virtual void mutations(ParamSet&P,RANDOM&R)override{P.Set(0,G(R));}
};
TEST(AbstractGenetic,FilterSettingFunc){
	auto initial_uniform=make_shared<InitialDistributions>()<<make_shared<RandomValueGenerator<double>>(0,1);
	auto filter=[](const ParamSet&P){return P[0]>0.5;};
	GeneticTestWithMutations gen(optimality);
	gen.SetFilter(filter);
	gen.Init(100,initial_uniform,engine);
	for(size_t i=0,n=gen.PopulationSize();i<n;i++)
		EXPECT_EQ(true,filter(gen.Parameters(i)));
	for(size_t j=0;j<10;j++){
		gen.Iterate(engine);
		for(size_t i=0,n=gen.PopulationSize();i<n;i++)
			EXPECT_EQ(true,filter(gen.Parameters(i)));
	}
	gen.RemoveFilter();
	for(size_t j=0;j<10;j++)gen.Iterate(engine);
	for(size_t i=0,n=gen.PopulationSize();i<n;i++)
		EXPECT_EQ(false,filter(gen.Parameters(i)));
}
TEST(AbstractGenetic,FilterSetting){
	auto initial_uniform=make_shared<InitialDistributions>()<<make_shared<RandomValueGenerator<double>>(-1,1);
	auto filter=[](const ParamSet&P){return P[0]>0.5;};
	GeneticTestWithMutations gen(optimality);
	gen.SetFilter(make_shared<Filter>(filter));
	gen.Init(100,initial_uniform,engine);
	for(size_t i=0,n=gen.PopulationSize();i<n;i++)
		EXPECT_EQ(true,filter(gen.Parameters(i)));
	for(size_t j=0;j<10;j++){
		gen.Iterate(engine);
		for(int i=0,n=gen.PopulationSize();i<n;i++)
			EXPECT_EQ(true,filter(gen.Parameters(i)));
	}
	gen.RemoveFilter();
	for(size_t j=0;j<10;j++)gen.Iterate(engine);
	for(size_t i=0,n=gen.PopulationSize();i<n;i++)
		EXPECT_EQ(false,filter(gen.Parameters(i)));
}
TEST(AbstractGenetic,Infinite){
	auto initial_uniform=make_shared<InitialDistributions>()<<make_shared<RandomValueGenerator<double>>(0,1);
	auto opt=make_shared<OptimalityFunction>([](const ParamSet&P){
		double res=P[0];
		res*=res;
		if(res<0.00001)return double(INFINITY);
		return res;
	});
	GeneticTest gen(opt);
	gen.Init(100,initial_uniform,engine);
	EXPECT_FALSE(gen.ParametersDispersionExitCondition(ParamSet(0)));
	for(size_t i=0,n=gen.PopulationSize();i<n;i++)
		EXPECT_TRUE(isfinite(gen.Parameters(i)[0]));
	for(size_t j=0;j<30;j++){
		gen.Iterate(engine);
		for(size_t i=0,n=gen.PopulationSize();i<n;i++){
			EXPECT_TRUE(isfinite(gen.Optimality(i)));
			if(i>0)
				EXPECT_TRUE(gen.Optimality(i-1)<=gen.Optimality(i));
		}
	}
}
