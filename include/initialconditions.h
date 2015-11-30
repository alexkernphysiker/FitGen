// this file is distributed under 
// GPL v 3.0 license
#ifndef ___DnoQpUEN
#define ___DnoQpUEN
#include <random>
#include <functional>
#include "../math_h/include/randomfunc.h"
#include "abstract.h"
namespace Genetic{
	using namespace std;
	typedef RandomValueGenerator<double> Distrib;
	class InitialDistributions:public IInitialConditions{
	public:
		InitialDistributions();
		virtual ~InitialDistributions();
		virtual ParamSet Generate(RANDOM&R)override;
		InitialDistributions &operator<<(shared_ptr<Distrib> distr);
		size_t Count();
		Distrib &operator[](size_t i);
	private:
		vector<shared_ptr<Distrib>> ParamDistr;
	};
	shared_ptr<InitialDistributions> operator<<(shared_ptr<InitialDistributions> init,shared_ptr<Distrib> func);
	
	class GenerateUniform:public IInitialConditions{
	public:
		GenerateUniform();
		virtual ~GenerateUniform();
		size_t Count();
		double Min(size_t i);
		double Max(size_t i);
		GenerateUniform &Add(double min,double max);
		virtual ParamSet Generate(RANDOM&R)override;
	private:
		vector<double> m_min;
		vector<double> m_max;
	};
	class GenerateByGauss:public IInitialConditions{
	public:
		GenerateByGauss();
		virtual ~GenerateByGauss();
		size_t Count();
		double Mean(size_t i);
		double Sigma(size_t i);
		GenerateByGauss &Add(double mean,double sig);
		virtual ParamSet Generate(RANDOM&R)override;
	private:
		vector<double> m_mean;
		vector<double> m_sig;
	};
	inline shared_ptr<GenerateByGauss> operator<<(shared_ptr<GenerateByGauss> G, pair<double,double> value){
		G->Add(value.first,value.second);
		return G;
	}
	inline shared_ptr<GenerateUniform> operator<<(shared_ptr<GenerateUniform> G, pair<double,double> value){
		G->Add(value.first,value.second);
		return G;
	}
}
#endif