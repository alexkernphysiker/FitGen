// this file is distributed under
// LGPL license
#ifndef ______EQUATIONS_WITHOUT_UNCERTAINTIES____H________
#	define ______EQUATIONS_WITHOUT_UNCERTAINTIES____H________
#include <utility>
#include <list>
#include "abstract.h"
namespace Genetic
{
struct Equation {
    std::function<double(const ParamSet &)> left, right;
};
class EquationSystem: public virtual IOptimalityFunction
{
public:
    EquationSystem(const std::initializer_list<Equation> &source);
    EquationSystem(const std::list<Equation> &source);
    virtual ~EquationSystem();
    virtual double operator()(const ParamSet &P)const override;
    const std::list<Equation> &equations()const;
private:
    std::list<Equation> f_data;
};
template<class MUTATION_TYPE>
class EquationSolver: public virtual MUTATION_TYPE
{
public:
    EquationSolver(const std::initializer_list<Equation> &source)
        : AbstractGenetic(std::make_shared<EquationSystem>(source)) {}
    EquationSolver(const std::list<Equation> &source)
        : AbstractGenetic(std::make_shared<EquationSystem>(source)) {}
    virtual ~EquationSolver() {}
    const std::list<Equation> &equations()const
    {
        return std::dynamic_pointer_cast<EquationSystem>(AbstractGenetic::OptimalityCalculator())->equations();
    }
};
}
#endif
