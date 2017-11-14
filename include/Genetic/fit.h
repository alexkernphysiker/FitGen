// this file is distributed under
// MIT license
#ifndef ____lrEPWamH___
#define ____lrEPWamH___
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <functional>
#include <math_h/tabledata.h>
#include "abstract.h"
#include "genetic.h"
namespace Genetic
{
class IParamFunc
{
public:
    virtual ~IParamFunc() {}
    virtual double operator()(const ParamSet &X, const ParamSet &P)const = 0;
};
typedef std::function<double(const ParamSet &, const ParamSet &)> paramFunc;
class ParameterFunction: public IParamFunc
{
public:
    ParameterFunction(const paramFunc f);
    virtual ~ParameterFunction();
    virtual double operator()(const ParamSet &X, const ParamSet &P)const override;
private:
    paramFunc func;
};
typedef MathTemplates::point<ParamSet,MathTemplates::value<>> Point;
typedef MathTemplates::Points<ParamSet,MathTemplates::value<>> FitPoints;
std::shared_ptr<FitPoints> operator<<(std::shared_ptr<FitPoints> src, const Point &p);
std::shared_ptr<FitPoints> operator<<(std::shared_ptr<FitPoints> src, const MathTemplates::point<MathTemplates::value<>> &p);
std::shared_ptr<FitPoints> operator<<(std::shared_ptr<FitPoints> src, const FitPoints&data);
std::shared_ptr<FitPoints> operator<<(std::shared_ptr<FitPoints> src, const MathTemplates::Points<> &h);
std::shared_ptr<FitPoints> operator<<(std::shared_ptr<FitPoints> src, const MathTemplates::Points<MathTemplates::value<>> &h);
inline std::shared_ptr<FitPoints> operator<<(std::shared_ptr<FitPoints> src, const MathTemplates::SortedPoints<> &h){return src<<h();}
inline std::shared_ptr<FitPoints> operator<<(std::shared_ptr<FitPoints> src, const MathTemplates::SortedPoints<MathTemplates::value<>> &h){return src<<h();}

class OptimalityForPoints: public IOptimalityFunction
{
public:
    typedef std::function<double(const ParamSet &, const IParamFunc &)> Coefficient;
    typedef std::function<double(const Point &, const ParamSet &, const IParamFunc &)> Summand;
    OptimalityForPoints(const std::shared_ptr<FitPoints> p, const  std::shared_ptr<IParamFunc> f, const Coefficient c, const Summand s);
    virtual ~OptimalityForPoints();
    virtual double operator()(const ParamSet &P)const override;
    std::shared_ptr<FitPoints> Points()const;
protected:
    std::shared_ptr<FitPoints> points;
    std::shared_ptr<IParamFunc> func;
    Coefficient C;
    Summand S;
};
std::shared_ptr<OptimalityForPoints> SumSquareDiff(const std::shared_ptr<FitPoints> points, const std::shared_ptr<IParamFunc> f);
std::shared_ptr<OptimalityForPoints> ChiSquare(const std::shared_ptr<FitPoints> points, const std::shared_ptr<IParamFunc> f);
template <
    class MUTATION_TYPE,
    std::shared_ptr<OptimalityForPoints> OptimalityAlgorithm(const std::shared_ptr<FitPoints>, const std::shared_ptr<IParamFunc>) = ChiSquare
    >
class Fit: public virtual MUTATION_TYPE
{
private:
    std::shared_ptr<IParamFunc> m_func;
protected:
    Fit(std::shared_ptr<IParamFunc> f)
    {
        m_func = f;
    }
public:
    Fit(
        const std::shared_ptr<FitPoints> points,
        const std::shared_ptr<IParamFunc> f
    ): AbstractGenetic(OptimalityAlgorithm(points, f))
    {
        m_func = f;
    }
    Fit(const std::shared_ptr<FitPoints> points, const paramFunc f): Fit(points, std::make_shared<ParameterFunction>(f)) {}
    virtual ~Fit() {}
    double operator()(const ParamSet &X)const
    {
        return m_func->operator()(X, AbstractGenetic::Parameters());
    }
    std::shared_ptr<IParamFunc> Func()const
    {
        return m_func;
    }
    const FitPoints& Points()const
    {
        return *std::dynamic_pointer_cast<OptimalityForPoints>(AbstractGenetic::OptimalityCalculator())->Points();
    }
    const MathTemplates::Points<double,MathTemplates::value<>> PointsProjection(const size_t index)
    {
	MathTemplates::Points<double,MathTemplates::value<>> res;
	for(const auto&p:Points())
	    res.push_back(make_point(p.X()[index],p.Y()));
	return res;
    }
};
template<class GENETIC, class FUNC, std::shared_ptr<OptimalityForPoints> OptimalityAlgorithm(const std::shared_ptr<FitPoints>, const std::shared_ptr<IParamFunc>) = ChiSquare>
class FitFunction: public virtual Fit<GENETIC, OptimalityAlgorithm>
{
public:
    typedef FUNC functype;
    FitFunction(const std::shared_ptr<FitPoints> points):
        AbstractGenetic(OptimalityAlgorithm(points, std::make_shared<FUNC>())),
        Fit<GENETIC, OptimalityAlgorithm>(std::make_shared<FUNC>()) {}
    virtual ~FitFunction() {}
};
}
#endif
