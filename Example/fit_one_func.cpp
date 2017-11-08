// this file is distributed under
// LGPL license
#include <iostream>
#include <fstream>
#include <random>
#include <gnuplot_wrap.h>
#include <math_h/functions.h>
#include <Genetic/fit.h>
#include <Genetic/initialconditions.h>
using namespace std;
using namespace Genetic;
using namespace GnuplotWrap;
using namespace MathTemplates;
int main()
{
    RANDOM random_engine;

    //Generating the distribution for fitting
    double left = 0, right = 10, sigma = 1.0;
    size_t count = 500;
    Distribution1D<> distribution(BinsByStep(left, 1.0, right));
    DistribGauss Distr((right + left) / 2.0, sigma);
    for (size_t i = 0; i < count; i++)
        distribution.Fill(Distr(random_engine));

    //Fitting generated distribution
    Fit<DifferentialMutations<>> fit(make_shared<FitPoints>() << distribution,
    [](const ParamSet & X, const ParamSet & P) {
        return Gaussian(X[0], P[0], P[1]) * P[2];
    });
    fit.SetFilter(
    [](const ParamSet & P) {
        return (P[1] > 0) && (P[2] > 0);
    }
    );
    fit.Init(30,
             make_shared<InitialDistributions>()
             << make_shared<DistribUniform>(left, right)
             << make_shared<DistribUniform>(0, right - left)
             << make_shared<DistribUniform>(0, 2.0 * count),
             random_engine
            );
    while (!fit.AbsoluteOptimalityExitCondition(0.000001))
        fit.Iterate(random_engine);

    //Output results
    cout << "Chi^2 = " << fit.Optimality() << endl;
    cout << "Chi^2 divided by degrees of freedom = "
         << fit.Optimality() / (fit.Points().size() - fit.ParamCount()) << endl;
    cout << endl;

    cout << "Fit parameters" << endl;
    for (const auto &P : fit.Parameters())cout << P << "\t";
    cout << endl;

    //plotting results
    Plot("FitGen-example1").Hist(distribution)
    .Line(SortedPoints<>([&fit](double x) {return fit({x});}, ChainWithStep(0.0, 0.01, 10.0)));
    return 0;
}
