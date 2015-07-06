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
		Equation(function<double(ParamSet&&)> f)
		:AbstractGenetic(make_shared<OptimalityFunction>([f](ParamSet&&P){
			return pow(f(static_cast<ParamSet&&>(P)),2);
		})),GENETIC(){}
		Equation(function<double(ParamSet&&)> left,function<double(ParamSet&&)> right)
		:AbstractGenetic(make_shared<OptimalityFunction>([left,right](ParamSet&&P){
			return pow(left(static_cast<ParamSet&&>(P))-right(static_cast<ParamSet&&>(P)),2);
		})),GENETIC(){}
	};
	template<class GENETIC>
	class SearchMin:public GENETIC{
	public:
		SearchMin(function<double(ParamSet&&)> f):AbstractGenetic(make_shared<OptimalityFunction>(f)),GENETIC(){}
	};
}
#endif
