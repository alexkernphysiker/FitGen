// this file is distributed under 
// MIT license
#ifndef GIQYSBLPGFAJYVCC
#define GIQYSBLPGFAJYVCC
#include "abstract.h"
#include "genetic.h"
namespace Genetic{
	template<class GENETIC>
	class Equation:public virtual GENETIC{
	public:
		Equation(const std::function<double(const ParamSet&)> f)
		:AbstractGenetic(std::make_shared<OptimalityFunction>([f](const ParamSet&P){
			return pow(f(P),2);
		})),GENETIC(){}
		Equation(const std::function<double(const ParamSet&)> left,const std::function<double(const ParamSet&)> right)
		:AbstractGenetic(std::make_shared<OptimalityFunction>([left,right](const ParamSet&P){
			return pow(left(P)-right(P),2);
		})),GENETIC(){}
	};
	template<class GENETIC>
	class SearchMin:public GENETIC{
	public:
		SearchMin(const std::function<double(const ParamSet&)> f)
		:AbstractGenetic(std::make_shared<OptimalityFunction>(f)),GENETIC(){}
	};
}
#endif
