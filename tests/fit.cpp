// this file is distributed under
// LGPLv3 license
#include <math.h>
#include <gtest/gtest.h>
#include <math_h/error.h>
#include <Genetic/fit.h>
#include <Genetic/uncertainties.h>
#include <Genetic/initialconditions.h>
#include <Genetic/paramfunc.h>
using namespace std;
using namespace MathTemplates;
using namespace Genetic;
const double epsilon=0.000001;
#define ALMOST_EQ(a,b) EXPECT_TRUE(pow(a-b,2)<epsilon)
TEST(ParameterFunction, Basic)
{
    int c = 0;
    double res = 2;
    ParameterFunction F([&c, &res](const ParamSet &, const ParamSet &) {
        c++;
        return res;
    });
    EXPECT_EQ(res, F(ParamSet(), ParamSet()));
    EXPECT_EQ(1, c);
}

TEST(point, Base)
{
    {
        Point P({0}, 0);
        EXPECT_EQ(1, P.X().size());
    }
    Point P({{1, 0.1}}, {1, 1});
    {
        Point P2(P);
        EXPECT_EQ(P.Y().val(), P2.Y().val());
        EXPECT_EQ(P.Y().uncertainty(), P2.Y().uncertainty());
        EXPECT_EQ(P.X().size(), P2.X().size());
        EXPECT_EQ(P.X()[0], P2.X()[0]);
    }
}
TEST(OptimalityForPoints, Base)
{
    FitPoints points;
    int func_calls = 0;
    auto f = [&func_calls](const ParamSet &, const ParamSet &) {
        func_calls++;
        return 0.0;
    };
    int summand_calls = 0;
    auto s = [&summand_calls](const Point &, const ParamSet &, const IParamFunc &) {
        summand_calls++;
        return 1.0;
    };
    int coef_calls = 0;
    auto c = [&coef_calls](const ParamSet &, const IParamFunc &) {
        coef_calls++;
        return 1.0;
    };
    for (int count = 1; count < 5; count++) {
        func_calls = summand_calls = coef_calls = 0;
        points.push_back(make_point(ParamSet{0}, 0));
	OptimalityForPoints S(points, make_shared<ParameterFunction>(f), c, s);
        EXPECT_EQ(count, S(ParamSet()));
        EXPECT_EQ(0, func_calls);
        EXPECT_EQ(points.size(), summand_calls);
        EXPECT_EQ(1, coef_calls);
    }
}
template<shared_ptr<OptimalityForPoints> OptimalityAlgorithm(const FitPoints&, shared_ptr<IParamFunc>)>
void test_optimality1(double v = INFINITY)
{
    FitPoints points{
	Point({0}, {0, 1}),
	Point({1}, {0, 1}),
	Point({2}, {0, 1})
    };
    auto F = make_shared<ParameterFunction>([](const ParamSet &, const ParamSet &) {
        return 0;
    });
    auto S = OptimalityAlgorithm(points, F);
    EXPECT_NE(nullptr, S.get());
    EXPECT_EQ(0, S->operator()(ParamSet()));
    auto F1 = make_shared<ParameterFunction>([](const ParamSet &, const ParamSet &) {
        return 1;
    });
    auto S1 = OptimalityAlgorithm(points, F1);
    EXPECT_NE(nullptr, S1.get());
    EXPECT_EQ(true, S1->operator()(ParamSet()) > 0);
    if (isfinite(v)) {
        EXPECT_EQ(v, S1->operator()(ParamSet()));
    }
}
TEST(OptimalityForPoints, SumSquareDiff)
{
    test_optimality1<SumSquareDiff>(3.);
}
TEST(OptimalityForPoints, ChiSquare)
{
    test_optimality1<ChiSquare>(3.);
}
typedef Add<Mul<Arg<0>, Par<0>>, Par<1>> Fit_Func;
typedef Const<1> Fit_Func_err;
FitPoints TestPoints{Point({0}, {1, 1}),Point({1}, {2, 1}),Point({2}, {3, 1})};
auto Init = make_shared<InitialDistributions>() << make_shared<DistribUniform>(0, 2) << make_shared<DistribUniform>(0, 2);
TEST(Fit, Basetest)
{
    Fit<DifferentialMutations<>, SumSquareDiff> fit(TestPoints, make_shared<Fit_Func>());
    fit.Init(30, Init);
    while (!fit.AbsoluteOptimalityExitCondition(epsilon))
        fit.Iterate();
    EXPECT_TRUE(fit.ParamCount() == 2);
    EXPECT_TRUE(fit.PopulationSize() == 30);
    ALMOST_EQ(fit.Optimality(),0.0);
    ALMOST_EQ(fit.Optimality(fit.PopulationSize() - 1),0.0);
    ALMOST_EQ(1.0, fit.Parameters()[0]);
    ALMOST_EQ(1.0, fit.Parameters()[1]);
    const auto fit_copy=fit;
    EXPECT_TRUE(fit_copy.ParamCount() == 2);
    EXPECT_TRUE(fit_copy.PopulationSize() == 30);
    ALMOST_EQ(fit_copy.Optimality(),0.0);
    ALMOST_EQ(fit_copy.Optimality(fit.PopulationSize() - 1),0.0);
    ALMOST_EQ(1.0, fit_copy.Parameters()[0]);
    ALMOST_EQ(1.0, fit_copy.Parameters()[1]);
}
TEST(Fit, Basetest2)
{
    Fit<DifferentialMutations<>,ChiSquare,UncertaintiesEstimation> fit(TestPoints, make_shared<Fit_Func>());
    fit.Init(30, Init);
    while (!fit.AbsoluteOptimalityExitCondition(epsilon))
        fit.Iterate();
    EXPECT_TRUE(fit.ParamCount() == 2);
    EXPECT_TRUE(fit.PopulationSize() == 30);
    ALMOST_EQ(fit.Optimality(),0.0);
    ALMOST_EQ(fit.Optimality(fit.PopulationSize() - 1),0.0);
    ALMOST_EQ(1.0, fit.Parameters()[0]);
    ALMOST_EQ(1.0, fit.Parameters()[1]);
    fit.SetUncertaintyCalcDeltas({0.01,0.01});
    EXPECT_TRUE(fit.ParametersWithUncertainties()[0].Contains(fit.Parameters()[0]));
    EXPECT_TRUE(fit.ParametersWithUncertainties()[1].Contains(fit.Parameters()[1]));
    const auto fit_copy=fit;
    EXPECT_TRUE(fit_copy.ParamCount() == 2);
    EXPECT_TRUE(fit_copy.PopulationSize() == 30);
    ALMOST_EQ(fit_copy.Optimality(),0.0);
    ALMOST_EQ(fit_copy.Optimality(fit.PopulationSize() - 1),0.0);
    ALMOST_EQ(1.0, fit_copy.Parameters()[0]);
    ALMOST_EQ(1.0, fit_copy.Parameters()[1]);
    EXPECT_TRUE(fit_copy.ParametersWithUncertainties()[0].Contains(fit.Parameters()[0]));
    EXPECT_TRUE(fit_copy.ParametersWithUncertainties()[1].Contains(fit.Parameters()[1]));
}
TEST(FitFunction, Basetest)
{
    FitFunction<DifferentialMutations<>, Fit_Func, SumSquareDiff> fit(TestPoints);
    fit.Init(30, Init);
    while (!fit.AbsoluteOptimalityExitCondition(epsilon))
        fit.Iterate();
    EXPECT_TRUE(fit.ParamCount() == 2);
    EXPECT_TRUE(fit.PopulationSize() == 30);
    ALMOST_EQ(fit.Optimality(),0.0);
    ALMOST_EQ(fit.Optimality(fit.PopulationSize() - 1),0.0);
    ALMOST_EQ(1.0, fit.Parameters()[0]);
    ALMOST_EQ(1.0, fit.Parameters()[1]);
    const auto fit_copy=fit;
    EXPECT_TRUE(fit_copy.ParamCount() == 2);
    EXPECT_TRUE(fit_copy.PopulationSize() == 30);
    ALMOST_EQ(fit_copy.Optimality(),0.0);
    ALMOST_EQ(fit_copy.Optimality(fit.PopulationSize() - 1),0.0);
    ALMOST_EQ(1.0, fit_copy.Parameters()[0]);
    ALMOST_EQ(1.0, fit_copy.Parameters()[1]);
}
TEST(FitFunction, Basetest2)
{
    FitFunction<DifferentialMutations<>, Fit_Func,ChiSquare,UncertaintiesEstimation> fit(TestPoints);
    fit.Init(30, Init);
    while (!fit.AbsoluteOptimalityExitCondition(epsilon))
        fit.Iterate();
    EXPECT_TRUE(fit.ParamCount() == 2);
    EXPECT_TRUE(fit.PopulationSize() == 30);
    ALMOST_EQ(fit.Optimality(),0.0);
    ALMOST_EQ(fit.Optimality(fit.PopulationSize() - 1),0.0);
    ALMOST_EQ(1.0, fit.Parameters()[0]);
    ALMOST_EQ(1.0, fit.Parameters()[1]);
    fit.SetUncertaintyCalcDeltas({0.01,0.01});
    EXPECT_TRUE(fit.ParametersWithUncertainties()[0].Contains(fit.Parameters()[0]));
    EXPECT_TRUE(fit.ParametersWithUncertainties()[1].Contains(fit.Parameters()[1]));
    const auto fit_copy=fit;
    EXPECT_TRUE(fit_copy.ParamCount() == 2);
    EXPECT_TRUE(fit_copy.PopulationSize() == 30);
    ALMOST_EQ(fit_copy.Optimality(),0.0);
    ALMOST_EQ(fit_copy.Optimality(fit.PopulationSize() - 1),0.0);
    ALMOST_EQ(1.0, fit_copy.Parameters()[0]);
    ALMOST_EQ(1.0, fit_copy.Parameters()[1]);
    EXPECT_TRUE(fit_copy.ParametersWithUncertainties()[0].Contains(fit.Parameters()[0]));
    EXPECT_TRUE(fit_copy.ParametersWithUncertainties()[1].Contains(fit.Parameters()[1]));
}
TEST(FitFunction, Basetest3)
{
    FitFunction<DifferentialMutations<>, Fit_Func,ChiSquare,FunctionUncertaintiesEstimation> fit(TestPoints);
    fit.Init(30, Init);
    while (!fit.AbsoluteOptimalityExitCondition(epsilon))
        fit.Iterate();
    EXPECT_TRUE(fit.ParamCount() == 2);
    EXPECT_TRUE(fit.PopulationSize() == 30);
    ALMOST_EQ(fit.Optimality(),0.0);
    ALMOST_EQ(fit.Optimality(fit.PopulationSize() - 1),0.0);
    ALMOST_EQ(1.0, fit.Parameters()[0]);
    ALMOST_EQ(1.0, fit.Parameters()[1]);
    fit.SetUncertaintyCalcDeltas({0.01,0.01});
    EXPECT_TRUE(fit.ParametersWithUncertainties()[0].Contains(fit.Parameters()[0]));
    EXPECT_TRUE(fit.ParametersWithUncertainties()[1].Contains(fit.Parameters()[1]));
    for(double x=0;x<=2;x+=0.1)EXPECT_TRUE(fit.FuncWithUncertainties({x}).Contains(fit({x})));
    const auto fit_copy=fit;
    EXPECT_TRUE(fit_copy.ParamCount() == 2);
    EXPECT_TRUE(fit_copy.PopulationSize() == 30);
    ALMOST_EQ(fit_copy.Optimality(),0.0);
    ALMOST_EQ(fit_copy.Optimality(fit.PopulationSize() - 1),0.0);
    ALMOST_EQ(1.0, fit_copy.Parameters()[0]);
    ALMOST_EQ(1.0, fit_copy.Parameters()[1]);
    EXPECT_TRUE(fit_copy.ParametersWithUncertainties()[0].Contains(fit.Parameters()[0]));
    EXPECT_TRUE(fit_copy.ParametersWithUncertainties()[1].Contains(fit.Parameters()[1]));
    for(double x=0;x<=2;x+=0.1)EXPECT_TRUE(fit_copy.FuncWithUncertainties({x}).Contains(fit_copy({x})));
}
