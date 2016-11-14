// this file is distributed under 
// MIT license
#ifndef ______EQUATIONS_WITHOUT_UNCERTAINTIES____H________
#	define ______EQUATIONS_WITHOUT_UNCERTAINTIES____H________
#include <utility>
#include <list>
#include "abstract.h"
namespace Genetic{
    typedef std::pair<std::function<double(const ParamSet&)>,double> Equation;
    inline const Equation eq(const std::function<double(const ParamSet&)>l,const double&r){return std::make_pair(l,r);}
    inline const Equation eq(const std::function<double(const ParamSet&)>l,const double&&r){return std::make_pair(l,r);}
    class EquationSystem:public virtual IOptimalityFunction{
    public:
	EquationSystem(const std::initializer_list<Equation>&source);
	EquationSystem(const std::list<Equation>&source);
	virtual ~EquationSystem();
	virtual double operator()(const ParamSet&P)const override;
    private:
	std::list<Equation> f_data;
    };
    template<class MUTATION_TYPE>
    class EquationSolver:public virtual MUTATION_TYPE{
    public:
	EquationSolver(const std::initializer_list<Equation>&source):AbstractGenetic(std::make_shared<EquationSystem>(source)){}
	EquationSolver(const std::list<Equation>&source):AbstractGenetic(std::make_shared<EquationSystem>(source)){}
	virtual ~EquationSolver(){}
	const double&operator [](const size_t i)const{return AbstractGenetic::Parameters()[i];}
	typedef ParamSet::const_iterator const_iterator;
	const_iterator begin()const{return AbstractGenetic::Parameters().begin();}
	const_iterator cbegin()const{return AbstractGenetic::Parameters().cbegin();}
	const_iterator end() const{return AbstractGenetic::Parameters().end();}
	const_iterator cend() const{return AbstractGenetic::Parameters().cend();}
    };
}
#endif
