// this file is distributed under
// LGPL license
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

shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints>src, const Point &p)
{
    src->push_back(p);
    return src;
}
shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints> src, const FitPoints &data)
{
    for (const Point &p : data)src << p;
    return src;
}
OptimalityForPoints::OptimalityForPoints(
    const std::shared_ptr< FitPoints > p,
    const shared_ptr<IParamFunc> f,
    const OptimalityForPoints::Coefficient c,
    const OptimalityForPoints::Summand s
)
{
    points = p;
    func = f;
    C = c;
    S = s;
}
OptimalityForPoints::~OptimalityForPoints() {}
shared_ptr<FitPoints> OptimalityForPoints::Points()const
{
    return points;
}
double OptimalityForPoints::operator()(const ParamSet &P)const
{
    double res = 0;
    for (const auto&p : *points)
        res += S(p, P, *func);
    return res * C(P, *func);
}
shared_ptr<OptimalityForPoints> SumSquareDiff(const shared_ptr<FitPoints> points, const shared_ptr<IParamFunc> f)
{
    OptimalityForPoints::Coefficient c = [](const ParamSet &, const IParamFunc &) {
        return 1.0;
    };
    OptimalityForPoints::Summand s = [](const Point & p, const ParamSet & P, const IParamFunc & F) {
        return pow(p.Y().val() - F(p.X(), P), 2);
    };
    return make_shared<OptimalityForPoints>(points, f, c, s);
}

shared_ptr<OptimalityForPoints> ChiSquare(const shared_ptr<FitPoints> points, const shared_ptr<IParamFunc> f)
{
    OptimalityForPoints::Coefficient c = [points](const ParamSet &, const IParamFunc &) {
        return 1.0;
    };
    OptimalityForPoints::Summand s = [](const Point & p, const ParamSet & P, const IParamFunc & F) {
        return p.Y().NumCompare(F(p.X(), P));
    };
    return make_shared<OptimalityForPoints>(points, f, c, s);
}
}
