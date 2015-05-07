#ifndef ___DnoQpUEN
#define ___DnoQpUEN
#include <functional>
#include "abstract.h"
#include "math_h/randomfunc.h"
namespace Genetic{
	using namespace std;
	typedef function<double()> generator;
	typedef RandomValueGenerator<double> Distrib;
	class Initialiser:public IInitialConditions{
	public:
		Initialiser();
		virtual ~Initialiser();
		virtual ParamSet Generate()override;
		Initialiser &operator<<(generator gen);
		int Count();
		generator operator[](int i);
		typedef vector<generator>::iterator iterator;
		typedef vector<generator>::const_iterator const_iterator;
		iterator begin();
		const_iterator cbegin()const;
		iterator end();
		const_iterator cend() const;
	private:
		vector<generator> generators;
	};
	shared_ptr<Initialiser> operator<<(shared_ptr<Initialiser> init, generator gen);
	
	class InitialDistributions:public IInitialConditions{
	public:
		InitialDistributions();
		virtual ~InitialDistributions();
		virtual ParamSet Generate()override;
		InitialDistributions &operator<<(shared_ptr<Distrib> distr);
		int Count();
		Distrib &operator[](int i);
	private:
		vector<shared_ptr<Distrib>> ParamDistr;
	};
	shared_ptr<InitialDistributions> operator<<(shared_ptr<InitialDistributions> init,shared_ptr<Distrib> func);
	
	class GenerateUniform:public IInitialConditions{
	public:
		GenerateUniform();
		virtual ~GenerateUniform();
		int Count();
		double Min(int i);
		double Max(int i);
		GenerateUniform &Add(double min,double max);
		virtual ParamSet Generate()override;
	private:
		vector<double> m_min;
		vector<double> m_max;
	};
	class GenerateByGauss:public IInitialConditions{
	public:
		GenerateByGauss();
		virtual ~GenerateByGauss();
		int Count();
		double Mean(int i);
		double Sigma(int i);
		GenerateByGauss &Add(double mean,double sig);
		virtual ParamSet Generate()override;
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