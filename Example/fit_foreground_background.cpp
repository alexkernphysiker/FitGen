#include <iostream>
#include <fstream>
#include <gnuplot_wrap.h>
#include <Genetic/fit.h>
#include <Genetic/parabolic.h>
#include <Genetic/filter.h>
#include <Genetic/initialconditions.h>
using namespace std;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
typedef Mul<Par<0>, Func3<BreitWigner, Arg<0>, Par<2>, Par<1>>> Foreground;
typedef PolynomFunc<Arg<0>, Foreground::ParamCount, 4> Background;
typedef Add<Foreground, Background> TotalFunc;
int main()
{
    RANDOM random_engine;
    //Fitting
    FitFunction<DifferentialMutations<Uncertainty>, TotalFunc> fit(
        make_shared<FitPoints>()
    << Point({{ -67.5, 2.5}}, {179.4, 12.5})
    << Point({{ -62.5, 2.5}}, {213.1, 13.0})
    << Point({{ -57.5, 2.5}}, {221.6, 12.0})
    << Point({{ -52.5, 2.5}}, {220.1, 11.5})
    << Point({{ -47.5, 2.5}}, {214.4, 10.5})
    << Point({{ -42.5, 2.5}}, {241.0, 10.5})
    << Point({{ -37.5, 2.5}}, {237.3, 10.5})
    << Point({{ -32.5, 2.5}}, {240.1, 10.0})
    << Point({{ -27.5, 2.5}}, {268.7, 10.0})
    << Point({{ -22.5, 2.5}}, {310.3, 10.0})
    << Point({{ -17.5, 2.5}}, {326.5, 11.0})
    << Point({{ -12.5, 2.5}}, {317.2, 10.5})
    << Point({{ -07.5, 2.5}}, {360.4, 12.0})
    << Point({{ -02.5, 2.5}}, {371.7, 11.5})
    << Point({{ 02.5, 2.5}}, {406.6, 12.0})
    << Point({{ 07.5, 2.5}}, {413.2, 12.0})
    << Point({{ 12.5, 2.5}}, {455.1, 13.0})
    << Point({{ 17.5, 2.5}}, {495.3, 14.0})
    << Point({{ 22.5, 2.5}}, {497.3, 14.5})
    << Point({{ 27.5, 2.5}}, {511.4, 15.0})
    );
    fit.SetFilter([](const ParamSet & P) {
        return (P[0] > 0) && (P[1] > 0);
    });
    auto initial = make_shared<InitialDistributions>()
                   << make_shared<DistribGauss>(100., 100.)
                   << make_shared<DistribUniform>(0., 50.)
                   << make_shared<DistribGauss>(-20., 1.)
                   << make_shared<DistribGauss>(400., 100.)
                   << make_shared<DistribGauss>(5., 2.)
                   << make_shared<DistribGauss>(0., 0.5);
    while (initial->Count() < TotalFunc::ParamCount)
        initial << make_shared<DistribGauss>(0., 0.01);
    fit.Init(TotalFunc::ParamCount * 10, initial, random_engine);

    while (!fit.AbsoluteOptimalityExitCondition(0.00001)) {
        fit.Iterate(random_engine);
        cout << fit.iteration_count() << " iterations; "
             << fit.Optimality() << " < Chi^2 < "
             << fit.Optimality(fit.PopulationSize() - 1)
             << "           \r";
    }
    cout << endl;

    //Output results
    cout << "Chi^2 divided by degrees of freedom = "
         << fit.Optimality() / (fit.Points().size() - fit.ParamCount())
         << endl;
    cout << "Fit parameters with uncertainties" << endl;
    fit.SetUncertaintyCalcDeltas(parEq(fit.ParamCount(), 0.01));
    for (const auto &P : fit.ParametersWithUncertainties())
        cout << P << endl;

    //Plotting total fit and background
    const auto &P = fit.Parameters();
    const auto chain = ChainWithStep(-70.0, 0.1, 30.0);
    const SortedPoints<>
    totalfit([&fit](double x) {
        return fit({x});
    }, chain),
    background([&P](double x) {
        return Background()({x}, P);
    }, chain);
    Plot("FitGen-example2").YUncertainties(fit.PointsProjection(0),"points").Line(totalfit, "fit")
    .Line(background, "background") << "set key on";
    return 0;
}
