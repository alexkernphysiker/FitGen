// this file is distributed under 
// MIT license
#ifndef ___DnoQpUEN
#define ___DnoQpUEN
#include <random>
#include <functional>
#include "../math_h/randomfunc.h"
#include "abstract.h"
namespace Genetic{
	using namespace std;
	using namespace MathTemplates;
	typedef RandomValueGenerator<double> Distrib;
	class InitialDistributions:public IInitialConditions{
	public:
		InitialDistributions();
		virtual ~InitialDistributions();
		virtual ParamSet Generate(RANDOM&R)const override;
		InitialDistributions &operator<<(const shared_ptr<Distrib> distr);
		const size_t Count()const ;
		const Distrib &operator[](const size_t i)const ;
	private:
		vector<shared_ptr<Distrib>> ParamDistr;
	};
	shared_ptr<InitialDistributions> operator<<(shared_ptr<InitialDistributions> init,const shared_ptr<Distrib> func);
	
	class GenerateUniform:public IInitialConditions{
	public:
		GenerateUniform();
		virtual ~GenerateUniform();
		const size_t Count()const ;
		const double Min(const size_t i)const ;
		const double Max(const size_t i)const ;
		GenerateUniform &Add(const double min,const double max);
		virtual ParamSet Generate(RANDOM&R)const override;
	private:
		vector<double> m_min;
		vector<double> m_max;
	};
	class GenerateByGauss:public IInitialConditions{
	public:
		GenerateByGauss();
		virtual ~GenerateByGauss();
		const size_t Count()const ;
		const double Mean(const size_t i)const ;
		const double Sigma(const size_t i)const ;
		GenerateByGauss &Add(const double mean,const double sig);
		virtual ParamSet Generate(RANDOM&R)const override;
	private:
		vector<double> m_mean;
		vector<double> m_sig;
	};
	inline shared_ptr<GenerateByGauss> operator<<(shared_ptr<GenerateByGauss> G,const pair<double,double>&&value){
		G->Add(value.first,value.second);
		return G;
	}
	inline shared_ptr<GenerateUniform> operator<<(shared_ptr<GenerateUniform> G,const pair<double,double>&&value){
		G->Add(value.first,value.second);
		return G;
	}
}
#endif