// this file is distributed under 
// MIT license
#ifndef GIQYSBLPGFAJYVCC
#define GIQYSBLPGFAJYVCC
#include "abstract.h"
#include "genetic.h"
namespace Genetic{
    template<class MUTATION_TYPE>
    	class AbstractEq:public virtual MUTATION_TYPE{
	protected:
		AbstractEq(){}
	public:
		virtual ~AbstractEq(){}
		const double&operator [](const size_t i)const{return AbstractGenetic::Parameters()[i];}
		typedef ParamSet::const_iterator const_iterator;
		const_iterator begin()const{return AbstractGenetic::Parameters().begin();}
		const_iterator cbegin()const{return AbstractGenetic::Parameters().cbegin();}
		const_iterator end() const{return AbstractGenetic::Parameters().end();}
		const_iterator cend() const{return AbstractGenetic::Parameters().cend();}
	};
	template<class MUTATION_TYPE>
	class Equation:public virtual AbstractEq<MUTATION_TYPE>{
	public:
		Equation(const std::function<double(const ParamSet&)> f)
		:AbstractGenetic(std::make_shared<OptimalityFunction>([f](const ParamSet&P){
			return pow(f(P),2);
		})){}
		Equation(const std::function<double(const ParamSet&)> left,const std::function<double(const ParamSet&)> right)
		:AbstractGenetic(std::make_shared<OptimalityFunction>([left,right](const ParamSet&P){
			return pow(left(P)-right(P),2);
		})){}
	};
	template<class MUTATION_TYPE>
	class SearchMin:public AbstractEq<MUTATION_TYPE>{
	public:
		SearchMin(const std::function<double(const ParamSet&)> f)
		:AbstractGenetic(std::make_shared<OptimalityFunction>(f)){}
	};
	template<class GENETIC>
	class SearchMax:public AbstractEq<GENETIC>{
	public:
		SearchMax(const std::function<double(const ParamSet&)> f)
		:AbstractGenetic(std::make_shared<OptimalityFunction>([f](const ParamSet&P){return -f(P);})){}
	};
}
#endif
