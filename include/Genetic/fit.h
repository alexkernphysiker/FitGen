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
class FunctionContainer{
protected:
    FunctionContainer(std::shared_ptr<IParamFunc> f);
public:
    FunctionContainer(const FunctionContainer&source);
    virtual ~FunctionContainer();
    std::shared_ptr<IParamFunc> Func()const;
    double func(const ParamSet&X,const ParamSet&P)const;
    double operator()(const ParamSet &X)const;
private:
    std::shared_ptr<IParamFunc> m_func;
};
template <
    class MUTATION_TYPE,
    std::shared_ptr<OptimalityForPoints> OptimalityAlgorithm(const FitPoints&, const std::shared_ptr<IParamFunc>) = ChiSquare,
    class... Parrents
>
class Fit: public virtual FunctionContainer,public virtual MUTATION_TYPE,public virtual Parrents...
{
    static_assert(std::is_base_of<EmptyMutation,MUTATION_TYPE>::value,"Mutation algorithm must be a class derived from EmptyMutation");
public:
    Fit(const FitPoints&points,const std::shared_ptr<IParamFunc> f):FunctionContainer(f),AbstractGenetic(OptimalityAlgorithm(points, f)){}
    Fit(const Fit&source):FunctionContainer(source),AbstractGenetic(source),MUTATION_TYPE(source),Parrents(source)...{}
    inline Fit(const FitPoints&points, const paramFunc f):Fit(points, std::make_shared<ParameterFunction>(f)) {}
    template<class Source>
    inline Fit(const Source&points,const std::shared_ptr<IParamFunc> f):Fit(ConvertPoints(points),f){}
    template<class Source>
    inline Fit(const Source&points,const paramFunc f):Fit(ConvertPoints(points),f){}
    virtual ~Fit() {}
    double operator()(const ParamSet &X)const{
	return FunctionContainer::func(X, AbstractGenetic::Parameters());
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
template<class GENETIC, class FUNC, std::shared_ptr<OptimalityForPoints> OptimalityAlgorithm(const FitPoints&, const std::shared_ptr<IParamFunc>) = ChiSquare,class...Parrents>
class FitFunction: public Fit<GENETIC, OptimalityAlgorithm,Parrents...>
{
    static_assert(std::is_base_of<EmptyMutation,GENETIC>::value,"Mutation algorithm must be a class derived from EmptyMutation");
public:
    typedef FUNC functype;
    FitFunction(const FitPoints& points):
	FunctionContainer(std::make_shared<FUNC>()),
        AbstractGenetic(OptimalityAlgorithm(points,FunctionContainer::Func())),
        Fit<GENETIC, OptimalityAlgorithm,Parrents...>(points,FunctionContainer::Func()){}
    FitFunction(const FitFunction&source):FunctionContainer(source),AbstractGenetic(source),Fit<GENETIC, OptimalityAlgorithm,Parrents...>(source),Parrents(source)...{}
    template<class Source>
    inline FitFunction(const Source&points):FitFunction(ConvertPoints(points)){}
    virtual ~FitFunction(){}
};

};
#endif
