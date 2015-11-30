// this file is distributed under 
// GPL v 3.0 license
#ifndef GIQYSBLPGFAJYVCC
#define GIQYSBLPGFAJYVCC
#include "abstract.h"
#include "genetic.h"
namespace Genetic{
	using namespace std;
	template<class GENETIC>
	class Equation:public virtual GENETIC{
	public:
		Equation(function<double(const ParamSet&)> f)
		:AbstractGenetic(make_shared<OptimalityFunction>([f](const ParamSet&P){
			return pow(f(P),2);
		})),GENETIC(){}
		Equation(function<double(const ParamSet&)> left,function<double(const ParamSet&)> right)
		:AbstractGenetic(make_shared<OptimalityFunction>([left,right](const ParamSet&P){
			return pow(left(P)-right(P),2);
		})),GENETIC(){}
	};
	template<class GENETIC>
	class SearchMin:public GENETIC{
	public:
		SearchMin(function<double(const ParamSet&)> f)
		:AbstractGenetic(make_shared<OptimalityFunction>(f)),GENETIC(){}
	};
}
#endif
