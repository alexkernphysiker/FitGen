// this file is distributed under
// LGPLv3 license
#ifndef ____WrKDhKHP___
#define ____WrKDhKHP___
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#ifdef using_multithread
#include <mutex>
#endif
#include <vector>
#include <memory>
#include <utility>
#include <math.h>
#include <math_h/sigma.h>
#include <math_h/randomfunc.h>
#include <math_h/tabledata.h>
#include "paramset.h"
namespace Genetic
{
class IInitialConditions
{
public:
    virtual ~IInitialConditions() {}
    virtual ParamSet Generate()const = 0;
};
class IParamCheck
{
public:
    virtual ~IParamCheck() {}
    virtual bool operator()(const ParamSet &P)const = 0;
};
class IOptimalityFunction
{
public:
    virtual ~IOptimalityFunction() {}
    virtual double operator()(const ParamSet &P)const = 0;
};
class AbstractGenetic
{
protected:
    AbstractGenetic();
    AbstractGenetic(const AbstractGenetic&source);
    AbstractGenetic(const std::shared_ptr<IOptimalityFunction> optimality);
public:
    virtual ~AbstractGenetic();
    std::shared_ptr<IOptimalityFunction> OptimalityCalculator()const;

    AbstractGenetic &SetFilter(const std::shared_ptr<IParamCheck> filter);
    AbstractGenetic &SetFilter(const std::function<bool(const ParamSet &)>);
    AbstractGenetic &RemoveFilter();
#ifdef using_multithread
    AbstractGenetic &SetThreadCount(const size_t threads_count);
    const size_t ThreadCount()const;
#endif
    AbstractGenetic &Init(const size_t population_size, const std::shared_ptr<IInitialConditions> initial_conditions);
    void Iterate();

    const unsigned long long int &iteration_count()const;
    const size_t PopulationSize()const;
    const size_t ParamCount()const;
    const double &Optimality(const size_t point_index = 0)const;
    const ParamSet &Parameters(const size_t point_index = 0)const;
    const std::vector<MathTemplates::value<double>> &ParametersStatistics()const;

    const bool ConcentratedInOnePoint()const;
    const bool AbsoluteOptimalityExitCondition(const double &accuracy)const;
    const bool RelativeOptimalityExitCondition(const double &accuracy)const;
    const bool ParametersDispersionExitCondition(const ParamSet &max_disp)const;
    const bool RelativeParametersDispersionExitCondition(const ParamSet &max_disp)const;
protected:
    virtual void mutations(ParamSet &)const;//contains empty implementation for templates from genetic.h could work correctly
    virtual void HandleIteration();
private:
    std::shared_ptr<IOptimalityFunction> m_optimality;
    std::shared_ptr<IParamCheck> m_filter;
    MathTemplates::SortedPoints<double,ParamSet> m_population;
    std::vector<MathTemplates::value<double>> m_stat;
    unsigned long long int m_itercount;
#ifdef using_multithread
    std::mutex m_mutex;
    size_t threads;
#endif
};

class Filter: public IParamCheck
{
public:
    Filter(const std::function<bool(const ParamSet &)> c);
    virtual ~Filter();
    virtual bool operator()(const ParamSet &P)const override;
private:
    std::function<bool(const ParamSet &)> condition;
};
}
#endif
