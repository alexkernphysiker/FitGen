#ifndef GIQYSBLPGFAJYVCC
#define GIQYSBLPGFAJYVCC
#include "abstract.h"
#include "genetic.h"
namespace Genetic{
	using namespace std;
	template<class GENETIC>
	class Equation:public GENETIC{
	public:
		Equation(function<double(ParamSet&)> f)
			:GENETIC(make_shared<OptimalityFunction>([f](ParamSet&P){return pow(f(P),2);})){}
		Equation(function<double(ParamSet&)> left,function<double(ParamSet&)> right)
			:GENETIC(make_shared<OptimalityFunction>([left,right](ParamSet&P){return pow(left(P)-right(P),2);})){}
	};
	template<class GENETIC>
	class SearchMin:public GENETIC{
	public:
		SearchMin(function<double(ParamSet&)> f)
		:GENETIC(make_shared<OptimalityFunction>([f](ParamSet&P){return f(P);})){}
	};
}
#endif
