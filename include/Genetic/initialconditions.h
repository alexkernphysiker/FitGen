// this file is distributed under
// LGPLv3 license
#ifndef ___DnoQpUEN
#define ___DnoQpUEN
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <math_h/randomfunc.h>
#include "abstract.h"
namespace Genetic
{
typedef MathTemplates::RandomValueGenerator<> Distrib;
typedef MathTemplates::RandomValueTableDistr<> DistribTable;
typedef MathTemplates::RandomGauss<> DistribGauss;
typedef MathTemplates::RandomUniform<> DistribUniform;
class FixParam: public Distrib
{
public:
    FixParam(const double &x);
    FixParam(const FixParam &source);
    virtual ~FixParam();
    virtual double operator()()const override;
private:
    double value;
};
class InitialDistributions: public IInitialConditions
{
public:
    InitialDistributions();
    virtual ~InitialDistributions();
    virtual ParamSet Generate()const override;
    InitialDistributions &operator<<(const std::shared_ptr<Distrib> distr);
    const size_t Count()const ;
    const Distrib &operator[](const size_t i)const ;
private:
    std::vector<std::shared_ptr<Distrib>> ParamDistr;
};
std::shared_ptr<InitialDistributions> operator<<(std::shared_ptr<InitialDistributions> init, const std::shared_ptr<Distrib> func);
}
#endif
