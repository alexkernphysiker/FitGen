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
		:GENETIC(make_shared<OptimalityFunction>([f](ParamSet&P){
			double res=f(P);
			if(res>=0)
				return res;
			else
				return -res;
		})){}
	};
	template<class GENETIC>
	class SearchMin:public GENETIC{
	public:
		SearchMin(function<double(ParamSet&)> f)
			:GENETIC(make_shared<OptimalityFunction>(f)){}
	};
}
#endif
