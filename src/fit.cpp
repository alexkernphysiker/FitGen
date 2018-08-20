// this file is distributed under
// LGPLv3 license
#include <math.h>
#include <climits>
#include <math_h/error.h>
#include <Genetic/fit.h>
namespace Genetic
{
using namespace std;
using namespace MathTemplates;
ParameterFunction::ParameterFunction(const function<double(const ParamSet &, const ParamSet &)> f)
{
    func = f;
}
ParameterFunction::~ParameterFunction() {}
double ParameterFunction::operator()(const ParamSet &X, const ParamSet &P)const
{
    return func(X, P);
}
FitPoints ConvertPoints(const FitPoints1D&source){
    FitPoints dest;
    for(const auto&p:source)dest.push_back(make_point(ParamSet{p.X()},p.Y()));
    return dest;
}
FitPoints ConvertPoints(const FitPoints1DSorted&source){
    FitPoints dest;
    for(const auto&p:source)dest.push_back(make_point(ParamSet{p.X()},p.Y()));
    return dest;
}
OptimalityForPoints::OptimalityForPoints(const FitPoints &p, const std::shared_ptr< Genetic::IParamFunc > f, const OptimalityForPoints::Coefficient c, const OptimalityForPoints::Summand s):points(p)
{
    func = f;
    C = c;
    S = s;
}
OptimalityForPoints::~OptimalityForPoints() {}
const FitPoints&OptimalityForPoints::Points()const
{
    return points;
}
double OptimalityForPoints::operator()(const ParamSet &P)const
{
    double res = 0;
    for (const auto&p : points)
        res += S(p, P, *func);
    return res * C(P, *func);
}
shared_ptr<OptimalityForPoints> SumSquareDiff(const FitPoints&points, const shared_ptr<IParamFunc> f)
{
    OptimalityForPoints::Coefficient c = [](const ParamSet &, const IParamFunc &) {
        return 1.0;
    };
    OptimalityForPoints::Summand s = [](const Point & p, const ParamSet & P, const IParamFunc & F) {
        return pow(p.Y().val() - F(p.X(), P), 2);
    };
    return make_shared<OptimalityForPoints>(points, f, c, s);
}

shared_ptr< OptimalityForPoints > ChiSquare(const FitPoints &points, const shared_ptr<IParamFunc > f)
{
    OptimalityForPoints::Coefficient c = [points](const ParamSet &, const IParamFunc &) {
        return 1.0;
    };
    OptimalityForPoints::Summand s = [](const Point & p, const ParamSet & P, const IParamFunc & F) {
        return p.Y().NumCompare(F(p.X(), P));
    };
    return make_shared<OptimalityForPoints>(points, f, c, s);
}


FunctionContainer::FunctionContainer(shared_ptr<IParamFunc> f):m_func(f){}
FunctionContainer::FunctionContainer(const FunctionContainer&source):m_func(source.m_func){}
FunctionContainer::~FunctionContainer(){}
shared_ptr<IParamFunc> FunctionContainer::Func()const{return m_func;}
double FunctionContainer::func(const ParamSet&X,const ParamSet&P)const{return m_func->operator()(X,P);}

}
