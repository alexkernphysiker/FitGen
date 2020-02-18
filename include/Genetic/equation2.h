// this file is distributed under
// LGPLv3 license
#ifndef ______EQUATIONS_WITH_UNCERTAINTIES____H________
#	define ______EQUATIONS_WITH_UNCERTAINTIES____H________
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <utility>
#include <list>
#include <vector>
#include <math_h/sigma.h>
#include "abstract.h"
#include "genetic.h"
namespace Genetic
{
struct InexactEquation {
    std::function<MathTemplates::value<double>(const ParamSet &)> left;
    MathTemplates::value<double> right;
};
class InexactEquationSystem: public virtual IOptimalityFunction
{
public:
    InexactEquationSystem(const std::initializer_list<InexactEquation> &source);
    InexactEquationSystem(const std::vector<InexactEquation> &source);
    virtual ~InexactEquationSystem();
    virtual double operator()(const ParamSet &P)const override;
    const std::vector<InexactEquation> &equations()const;
private:
    std::vector<InexactEquation> f_data;
};
template<class MUTATION_TYPE,class...Parrents>
class InexactEquationSolver: public virtual MUTATION_TYPE,public virtual Parrents...
{
    static_assert(std::is_base_of<EmptyMutation,MUTATION_TYPE>::value,"Mutation algorithm must be a class derived from EmptyMutation");
public:
    InexactEquationSolver(const std::initializer_list<InexactEquation> &source)
        : AbstractGenetic(std::make_shared<InexactEquationSystem>(source)){}
    InexactEquationSolver(const std::vector<InexactEquation> &source)
        : AbstractGenetic(std::make_shared<InexactEquationSystem>(source)){}
    InexactEquationSolver(InexactEquationSolver &&source)
        : AbstractGenetic(std::move(source)),MUTATION_TYPE(std::move(source)),Parrents(std::move(source))...{}
    virtual ~InexactEquationSolver() {}
    const std::vector<InexactEquation> &equations()const
    {
        return std::dynamic_pointer_cast<InexactEquationSystem>(
                   AbstractGenetic::OptimalityCalculator()
               )->equations();
    }
};
}
#endif
