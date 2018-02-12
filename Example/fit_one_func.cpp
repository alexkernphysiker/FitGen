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
    //Generating the distribution for fitting
    double left = 0, right = 10, sigma = 1.0;
    size_t count = 500;
    Distribution1D<> distribution(BinsByStep(left, 1.0, right));
    DistribGauss Distr((right + left) / 2.0, sigma);
    for (size_t i = 0; i < count; i++)
        distribution.Fill(Distr());

    //Fitting generated distribution
#ifdef ____middle_version_of_math_h_____
    //in c++14+ version you can convert histogram to proper type
    Fit<DifferentialMutations<>> fit(
	distribution.removeXerorbars(),
	[](const ParamSet & X, const ParamSet & P) {
	    return Gaussian(X[0], P[0], P[1]) * P[2];
	}
    );
#else
    //in c++11 version you have to convert histogram points one-by-one
    FitPoints points;
    for(const auto&P:distribution)points.push_back(make_point(ParamSet{P.X().val()},P.Y()));
    Fit<DifferentialMutations<>> fit(
	points,
	[](const ParamSet & X, const ParamSet & P) {
	    return Gaussian(X[0], P[0], P[1]) * P[2];
	}
    );
#endif
    fit.SetFilter(
    [](const ParamSet & P) {
        return (P[1] > 0) && (P[2] > 0);
    }
    );
    fit.Init(30,
             make_shared<InitialDistributions>()
             << make_shared<DistribUniform>(left, right)
             << make_shared<DistribUniform>(0, right - left)
             << make_shared<DistribUniform>(0, 2.0 * count)
            );
    while (!fit.AbsoluteOptimalityExitCondition(0.000001))
        fit.Iterate();

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
