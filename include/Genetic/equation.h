// this file is distributed under 
// MIT license
#ifndef GIQYSBLPGFAJYVCC
#define GIQYSBLPGFAJYVCC
#include "abstract.h"
#include "genetic.h"
namespace Genetic{
	template<class GENETIC>
	class AbstractEq:public virtual GENETIC{
	public:
		const double operator [](const size_t i)const{return AbstractGenetic::Parameters()[i];}
		typedef ParamSet::const_iterator const_iterator;
		const_iterator begin()const{return AbstractGenetic::Parameters().begin();}
		const_iterator cbegin()const{return AbstractGenetic::Parameters().cbegin();}
		const_iterator end() const{return AbstractGenetic::Parameters().end();}
		const_iterator cend() const{return AbstractGenetic::Parameters().cend();}
	};
	template<class GENETIC>
	class Equation:public virtual AbstractEq<GENETIC>{
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
	class SearchMin:public AbstractEq<GENETIC>{
	public:
		SearchMin(const std::function<double(const ParamSet&)> f)
		:AbstractGenetic(std::make_shared<OptimalityFunction>(f)),GENETIC(){}
	};
	template<class GENETIC>
	class SearchMax:public AbstractEq<GENETIC>{
	public:
		SearchMax(const std::function<double(const ParamSet&)> f)
		:AbstractGenetic(std::make_shared<OptimalityFunction>([f](const ParamSet&P){return -f(P);})),GENETIC(){}
	};
}
#endif
