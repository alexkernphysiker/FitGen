// this file is distributed under 
// LGPL license
#ifndef ______EQUATIONS_WITH_UNCERTAINTIES____H________
#	define ______EQUATIONS_WITH_UNCERTAINTIES____H________
#include <utility>
#include <list>
#include <math_h/sigma.h>
#include "abstract.h"
namespace Genetic{
	struct InexactEquation{
	    std::function<MathTemplates::value<double>(const ParamSet&)> left;
	    MathTemplates::value<double> right;
	};
	class InexactEquationSystem:public virtual IOptimalityFunction{
	public:
		InexactEquationSystem(const std::initializer_list<InexactEquation>&source);
		InexactEquationSystem(const std::list<InexactEquation>&source);
		virtual ~InexactEquationSystem();
		virtual double operator()(const ParamSet&P)const override;
		const std::list<InexactEquation>&equations()const;
	private:
		std::list<InexactEquation> f_data;
	};
	template<class MUTATION_TYPE>
	class InexactEquationSolver:public virtual MUTATION_TYPE{
	public:
	    InexactEquationSolver(const std::initializer_list<InexactEquation>&source)
	    :AbstractGenetic(std::make_shared<InexactEquationSystem>(source)){}
	    InexactEquationSolver(const std::list<InexactEquation>&source)
	    :AbstractGenetic(std::make_shared<InexactEquationSystem>(source)){}
	    virtual ~InexactEquationSolver(){}
	    const std::list<InexactEquation>&equations()const{
		return std::dynamic_pointer_cast<InexactEquationSystem>(
		    AbstractGenetic::OptimalityCalculator()
		)->equations();
	    }
	};
}
#endif
