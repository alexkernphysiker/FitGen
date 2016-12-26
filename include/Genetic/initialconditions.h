// this file is distributed under 
// LGPL license
#ifndef ___DnoQpUEN
#define ___DnoQpUEN
#include <random>
#include <functional>
#include <math_h/randomfunc.h>
#include "abstract.h"
namespace Genetic{
	typedef MathTemplates::IFunction<double,RANDOM&> Distrib;
	typedef MathTemplates::RandomValueTableDistr<double,RANDOM&> DistribTable;
	typedef MathTemplates::RandomGauss<double,RANDOM&> DistribGauss;
	typedef MathTemplates::RandomUniform<double,RANDOM&> DistribUniform;
	class FixParam:public Distrib{
	public:
	    FixParam(const double&x);
	    FixParam(const FixParam&source);
	    virtual ~FixParam();
	    virtual double operator ()(RANDOM&)const override;
	private:
	    double value;
	};
	class InitialDistributions:public IInitialConditions{
	public:
		InitialDistributions();
		virtual ~InitialDistributions();
		virtual ParamSet Generate(RANDOM&R)const override;
		InitialDistributions &operator<<(const std::shared_ptr<Distrib> distr);
		const size_t Count()const ;
		const Distrib&operator[](const size_t i)const ;
	private:
		std::vector<std::shared_ptr<Distrib>> ParamDistr;
	};
	std::shared_ptr<InitialDistributions> operator<<(std::shared_ptr<InitialDistributions> init,const std::shared_ptr<Distrib> func);
}
#endif
