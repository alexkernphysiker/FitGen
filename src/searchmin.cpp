// this file is distributed under
// LGPL license
#ifdef using_multithread
#include <thread>
#endif
#include <math.h>
#include <math_h/error.h>
#include <math_h/interpolate.h>
#include <Genetic/searchmin.h>
namespace Genetic
{
using namespace std;
using namespace MathTemplates;
OptimalityFunction::OptimalityFunction(const function<double(const ParamSet &)> f)
{
    func = f;
}
OptimalityFunction::~OptimalityFunction() {}
double OptimalityFunction::operator()(const ParamSet &P)const
{
    return func(P);
}
};
