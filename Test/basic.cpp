#include <gtest/gtest.h>
#include <abstract.h>
#include <math_h/randomfunc.h>
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

class GeneticTest:public AbstractGenetic{
public:
	GeneticTest(shared_ptr<IOptimalityFunction> optimality):AbstractGenetic(optimality){}
	virtual ~GeneticTest(){}
};
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
auto initial_uniform=make_shared<Init>([](){return ParamSet(RandomUniformlyR(0.0,1.0));});
void test_init(unsigned int threads,int population){
	GeneticTest gen(optimality);
	EXPECT_EQ(optimality.get(),gen.OptimalityCalculator().get());
	ASSERT_ANY_THROW(gen.SetThreadCount(0));
	gen.SetThreadCount(threads);
	EXPECT_EQ(threads,gen.ThreadCount());
	EXPECT_EQ(0,gen.PopulationSize());
	ASSERT_ANY_THROW(gen.Init(0,initial));
	ASSERT_ANY_THROW(gen.ParamCount());
	EXPECT_EQ(0,gen.PopulationSize());
	ASSERT_NO_THROW(gen.Init(population,initial));
	EXPECT_EQ(population,gen.PopulationSize());
	ASSERT_ANY_THROW(gen.Init(population,initial));
	EXPECT_EQ(population,gen.PopulationSize());
	EXPECT_EQ(1,gen.ParamCount());
}
TEST(AbstractGenetic,InitSync){test_init(1,10);}
TEST(AbstractGenetic,InitAsync){test_init(2,10);}
