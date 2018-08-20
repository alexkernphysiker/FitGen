// this file is distributed under
// LGPLv3 license
#ifndef ______EQUATIONS_WITHOUT_UNCERTAINTIES____H________
#	define ______EQUATIONS_WITHOUT_UNCERTAINTIES____H________
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <utility>
#include <vector>
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
    EquationSystem(const std::vector<Equation> &source);
    virtual ~EquationSystem();
    virtual double operator()(const ParamSet &P)const override;
    const std::vector<Equation> &equations()const;
private:
    std::vector<Equation> f_data;
};
template<class MUTATION_TYPE>
class EquationSolver: public virtual MUTATION_TYPE
{
public:
    EquationSolver(const std::initializer_list<Equation> &source)
        : AbstractGenetic(std::make_shared<EquationSystem>(source)) {}
    EquationSolver(const std::vector<Equation> &source)
        : AbstractGenetic(std::make_shared<EquationSystem>(source)) {}
    EquationSolver(const EquationSolver&source):AbstractGenetic(source),MUTATION_TYPE(source){}
    virtual ~EquationSolver() {}
    const std::vector<Equation> &equations()const
    {
        return std::dynamic_pointer_cast<EquationSystem>(AbstractGenetic::OptimalityCalculator())->equations();
    }
};
}
#endif
