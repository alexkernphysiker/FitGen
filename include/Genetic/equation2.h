// this file is distributed under 
// MIT license
#ifndef ______EQUATIONS_WITH_UNCERTAINTIES____H________
#	define ______EQUATIONS_WITH_UNCERTAINTIES____H________
#include <utility>
#include <list>
#include <math_h/sigma.h>
#include "abstract.h"
#include "parabolic.h"
namespace Genetic{
	typedef std::pair<std::function<double(const ParamSet&)>,MathTemplates::value<double>> InexactEquation;
	inline const InexactEquation in_eq(const std::function<double(const ParamSet&)>l,const MathTemplates::value<double>&r){return std::make_pair(l,r);}
	inline const InexactEquation in_eq(const std::function<double(const ParamSet&)>l,const MathTemplates::value<double>&&r){return std::make_pair(l,r);}
	class InexactEquationSystem:public virtual IOptimalityFunction{
	public:
		InexactEquationSystem(const std::initializer_list<InexactEquation>&source);
		InexactEquationSystem(const std::list<InexactEquation>&source);
		virtual ~InexactEquationSystem();
		virtual double operator()(const ParamSet&P)const override;
	private:
		std::list<InexactEquation> f_data;
	};
	template<class GENETIC>
	class InexactEquationSolver:public virtual GENETIC,public virtual ParabolicErrorEstimationFromChisq{
	public:
		InexactEquationSolver(const std::initializer_list<InexactEquation>&source):AbstractGenetic(std::make_shared<InexactEquationSystem>(source)){}
		InexactEquationSolver(const std::list<InexactEquation>&source):AbstractGenetic(std::make_shared<InexactEquationSystem>(source)){}
		virtual ~InexactEquationSolver(){}
		const double&operator [](const size_t i)const{return AbstractGenetic::Parameters()[i];}
		typedef ParamSet::const_iterator const_iterator;
		const_iterator begin()const{return AbstractGenetic::Parameters().begin();}
		const_iterator cbegin()const{return AbstractGenetic::Parameters().cbegin();}
		const_iterator end() const{return AbstractGenetic::Parameters().end();}
		const_iterator cend() const{return AbstractGenetic::Parameters().cend();}
		
	};
}
#endif
