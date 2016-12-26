// this file is distributed under 
// LGPL license
#ifndef ______searchingminimum____H________
#	define ______searchingminimum____H________
#include <utility>
#include <list>
#include "abstract.h"
namespace Genetic{
    class OptimalityFunction:public IOptimalityFunction{
    public:
	OptimalityFunction(const std::function<double(const ParamSet&)> f);
	virtual ~OptimalityFunction();
	virtual double operator()(const ParamSet&P)const override;
    private:
	std::function<double(const ParamSet&)> func;
    };
    template<class MUTATION_TYPE>
    class SearchMin:public virtual MUTATION_TYPE{
    public:
	SearchMin(const std::function<double(const ParamSet&)> f)
	:AbstractGenetic(std::make_shared<OptimalityFunction>(f)){}
	virtual ~SearchMin(){}
    };
}
#endif

