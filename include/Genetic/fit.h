// this file is distributed under
// LGPLv3 license
#ifndef ____lrEPWamH___
#define ____lrEPWamH___
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <functional>
#include <math_h/tabledata.h>
#include "abstract.h"
#include "genetic.h"
#include "uncertainties.h"
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
typedef MathTemplates::point<double,MathTemplates::value<>> Point1D;
typedef MathTemplates::Points<double,MathTemplates::value<>> FitPoints1D;
typedef MathTemplates::SortedPoints<double,MathTemplates::value<>> FitPoints1DSorted;
class OptimalityForPoints: public IOptimalityFunction
{
public:
    typedef std::function<double(const ParamSet &, const IParamFunc &)> Coefficient;
    typedef std::function<double(const Point &, const ParamSet &, const IParamFunc &)> Summand;
    OptimalityForPoints(const FitPoints&p, const  std::shared_ptr<IParamFunc> f, const Coefficient c, const Summand s);
    virtual ~OptimalityForPoints();
    virtual double operator()(const ParamSet &P)const override;
    const FitPoints&Points()const;
protected:
    FitPoints points;
    Coefficient C;
    Summand S;
    std::shared_ptr<IParamFunc> func;
};
std::shared_ptr<OptimalityForPoints> SumSquareDiff(const FitPoints&points, const std::shared_ptr<IParamFunc> f);
std::shared_ptr<OptimalityForPoints> ChiSquare(const FitPoints&points, const std::shared_ptr<IParamFunc> f);
FitPoints ConvertPoints(const FitPoints1D&source);
FitPoints ConvertPoints(const FitPoints1DSorted&source);
template <
    class MUTATION_TYPE,
    std::shared_ptr<OptimalityForPoints> OptimalityAlgorithm(const FitPoints&, const std::shared_ptr<IParamFunc>) = ChiSquare
    >
class Fit: public virtual MUTATION_TYPE
{
private:
    std::shared_ptr<IParamFunc> m_func;
protected:
    Fit(std::shared_ptr<IParamFunc> f){m_func = f;}
public:
    Fit(const FitPoints&points,const std::shared_ptr<IParamFunc> f):AbstractGenetic(OptimalityAlgorithm(points, f)){m_func = f;}
    inline Fit(const FitPoints&points, const paramFunc f): Fit(points, std::make_shared<ParameterFunction>(f)) {}
    template<class Source>
    inline Fit(const Source&points,const std::shared_ptr<IParamFunc> f):Fit(ConvertPoints(points),f){}
    template<class Source>
    inline Fit(const Source&points,const paramFunc f):Fit(ConvertPoints(points),f){}
    virtual ~Fit() {}
    inline std::shared_ptr<IParamFunc> Func()const{return m_func;}
    inline double func(const ParamSet&X,const ParamSet&P)const{return m_func->operator()(X,P);}
    inline double operator()(const ParamSet &X)const
    {
        return func(X, AbstractGenetic::Parameters());
    }
    inline const FitPoints& Points()const
    {
        return std::dynamic_pointer_cast<OptimalityForPoints>(AbstractGenetic::OptimalityCalculator())->Points();
    }
    const MathTemplates::Points<double,MathTemplates::value<>> PointsProjection(const size_t index)
    {
	MathTemplates::Points<double,MathTemplates::value<>> res;
	for(const auto&p:Points())
	    res.push_back(make_point(p.X()[index],p.Y()));
	return res;
    }
};
template<class GENETIC, class FUNC, std::shared_ptr<OptimalityForPoints> OptimalityAlgorithm(const FitPoints&, const std::shared_ptr<IParamFunc>) = ChiSquare>
class FitFunction: public virtual Fit<GENETIC, OptimalityAlgorithm>
{
public:
    typedef FUNC functype;
    FitFunction(const FitPoints& points):
        AbstractGenetic(OptimalityAlgorithm(points, std::make_shared<FUNC>())),
        Fit<GENETIC, OptimalityAlgorithm>(std::make_shared<FUNC>()){}
    template<class Source>
    inline FitFunction(const Source&points):FitFunction(ConvertPoints(points)){}
    virtual ~FitFunction() {}
};

template<class MUTATION_TYPE>
class Fit2: public virtual Fit<MUTATION_TYPE,ChiSquare>,public virtual ParabolicErrorEstimationFromChisq
{
public:
    Fit2(const FitPoints&points,const std::shared_ptr<IParamFunc> f):AbstractGenetic(ChiSquare(points, f))
	,Fit<MUTATION_TYPE,ChiSquare>(f),ParabolicErrorEstimationFromChisq(){}
    inline Fit2(const FitPoints&points, const paramFunc f): Fit2(points, std::make_shared<ParameterFunction>(f)) {}
    template<class Source>
    inline Fit2(const Source&points,const std::shared_ptr<IParamFunc> f):Fit2(ConvertPoints(points),f){}
    template<class Source>
    inline Fit2(const Source&points,const paramFunc f):Fit2(ConvertPoints(points),f){}
    virtual ~Fit2() {}
    MathTemplates::value<> FuncWithUncertainties(const ParamSet &X)const
    {
        const std::function<double(const std::vector<double>&)> F=[&X,this](const std::vector<double>&P){
	    return Fit<MUTATION_TYPE,ChiSquare>::Func()->operator()(X,ParamSet{P});
	};
	return MathTemplates::FUNC(F,ParabolicErrorEstimationFromChisq::ParametersWithUncertainties());
    }
};
template<class GENETIC, class FUNC>
class FitFunction2: public virtual Fit2<GENETIC>
{
public:
    typedef FUNC functype;
    FitFunction2(const FitPoints& points):
	AbstractGenetic(ChiSquare(points, std::make_shared<FUNC>())),
	Fit<GENETIC,ChiSquare>(points,std::make_shared<FUNC>()),
	Fit2<GENETIC>(points,std::make_shared<FUNC>()){}
    template<class Source>
    inline FitFunction2(const Source&points):FitFunction2(ConvertPoints(points)){}
    virtual ~FitFunction2() {}
};
}
#endif
