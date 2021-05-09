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
#include "genetic.h"
namespace Genetic
{
struct Equation {
    const std::function<double(const ParamSet &)> left;
    const std::function<double(const ParamSet &)> right;
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
template<class MUTATION_TYPE,class...Parrents>
class EquationSolver: public virtual MUTATION_TYPE, public virtual Parrents...
{
    static_assert(std::is_base_of<EmptyMutation,MUTATION_TYPE>::value,"Mutation algorithm must be a class derived from EmptyMutation");
public:
    EquationSolver(const std::initializer_list<Equation> &source)
        : AbstractGenetic(std::make_shared<EquationSystem>(source)) {}
    EquationSolver(const std::vector<Equation> &source)
        : AbstractGenetic(std::make_shared<EquationSystem>(source)) {}
    EquationSolver(EquationSolver&&source):AbstractGenetic(std::move(source)),MUTATION_TYPE(std::move(source)),Parrents(std::move(source))...{}
    virtual ~EquationSolver() {}
    const std::vector<Equation> &equations()const
    {
        return std::dynamic_pointer_cast<EquationSystem>(AbstractGenetic::OptimalityCalculator())->equations();
    }
};
}
#endif
