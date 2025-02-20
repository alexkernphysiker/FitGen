// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/error.h>
#include <Genetic/equation2.h>
#include <Genetic/genetic.h>
#include <Genetic/uncertainties.h>
#include <Genetic/initialconditions.h>
using namespace std;
using namespace MathTemplates;
using namespace Genetic;
TEST(InexactEquationSystem, empty)
{
    InexactEquationSystem A{};
    EXPECT_EQ(A({}), 0);
}
TEST(InexactEquationSystem, empty2)
{
    InexactEquationSystem A(vector<InexactEquation> {});
    EXPECT_EQ(A({}), 0);
}
TEST(InexactEquationSystem, simple)
{
    InexactEquationSystem A{
        {
        .left = [](const ParamSet& P)->double {return P[0];},
        .right = {0, 1}
    }
    };
    EXPECT_EQ(A({ 0.0 }), 0);
    EXPECT_EQ(A({ 0.5 }), 0.25);
    EXPECT_EQ(A({ 1.0 }), 1);
    EXPECT_EQ(A({ 2.0 }), 4);
}
TEST(InexactEquationSystem, simple2)
{
    InexactEquationSystem A(vector<InexactEquation> {
        {
            .left = [](const ParamSet& P)->double {return P[0];},
                .right = { 0, 1 }
        }
    });
    EXPECT_EQ(A({ 0.0 }), 0);
    EXPECT_EQ(A({ 0.5 }), 0.25);
    EXPECT_EQ(A({ 1.0 }), 1);
    EXPECT_EQ(A({ 2.0 }), 4);
}
TEST(InexactEquationSystem, twoparams)
{
    InexactEquationSystem A{
        {
        .left = [](const ParamSet& P)->double {return P[0] + P[1];},
        .right = {0, 1}
    }
    };
    EXPECT_EQ(A({ 0.0, 0.0 }), 0);
    EXPECT_EQ(A({ 0.5, 0.0 }), 0.25);
    EXPECT_EQ(A({ 1.0, 0.0 }), 1);
    EXPECT_EQ(A({ 2.0, 0.0 }), 4);
    EXPECT_EQ(A({ 0.0, 0.0 }), 0);
    EXPECT_EQ(A({ 0.0, 0.5 }), 0.25);
    EXPECT_EQ(A({ 0.0, 1.0 }), 1);
    EXPECT_EQ(A({ 0.0, 2.0 }), 4);
}
TEST(InexactEquationSystem, twoparams2)
{
    InexactEquationSystem A(vector<InexactEquation> {
        {
            .left = [](const ParamSet& P)->double {return P[0] + P[1];},
                .right = { 0, 1 }
        }
    });
    EXPECT_EQ(A({ 0.0, 0.0 }), 0);
    EXPECT_EQ(A({ 0.5, 0.0 }), 0.25);
    EXPECT_EQ(A({ 1.0, 0.0 }), 1);
    EXPECT_EQ(A({ 2.0, 0.0 }), 4);
    EXPECT_EQ(A({ 0.0, 0.0 }), 0);
    EXPECT_EQ(A({ 0.0, 0.5 }), 0.25);
    EXPECT_EQ(A({ 0.0, 1.0 }), 1);
    EXPECT_EQ(A({ 0.0, 2.0 }), 4);
}
TEST(InexactEquationSystem, two_eq)
{
    InexactEquationSystem A{
        {
        .left = [](const ParamSet& P)->double {return P[0] + P[1];},
        .right = {0, 1}
    },{
        .left = [](const ParamSet& P)->double {return P[0] - P[1];},
        .right = {0, 1}
    }
    };
    EXPECT_EQ(A({ 0.0, 0.0 }), 0);
    EXPECT_EQ(A({ 0.5, 0.0 }), 0.5);
    EXPECT_EQ(A({ 1.0, 0.0 }), 2);
    EXPECT_EQ(A({ 2.0, 0.0 }), 8);
    EXPECT_EQ(A({ 0.0, 0.0 }), 0);
    EXPECT_EQ(A({ 0.0, 0.5 }), 0.5);
    EXPECT_EQ(A({ 0.0, 1.0 }), 2);
    EXPECT_EQ(A({ 0.0, 2.0 }), 8);
}
TEST(InexactEquationSystem, two_eq2)
{
    InexactEquationSystem A(vector<InexactEquation> {
        {
            .left = [](const ParamSet& P)->double {return P[0] + P[1];},
                .right = { 0, 1 }
        }, {
            .left = [](const ParamSet& P)->double {return P[0] - P[1];},
            .right = {0, 1}
        }
    });
    EXPECT_EQ(A({ 0.0, 0.0 }), 0);
    EXPECT_EQ(A({ 0.5, 0.0 }), 0.5);
    EXPECT_EQ(A({ 1.0, 0.0 }), 2);
    EXPECT_EQ(A({ 2.0, 0.0 }), 8);
    EXPECT_EQ(A({ 0.0, 0.0 }), 0);
    EXPECT_EQ(A({ 0.0, 0.5 }), 0.5);
    EXPECT_EQ(A({ 0.0, 1.0 }), 2);
    EXPECT_EQ(A({ 0.0, 2.0 }), 8);
}
TEST(InexactEquationSolver, Integrationtest1)
{
    InexactEquationSolver<DifferentialMutations<>> test{
        {
        .left = [](const ParamSet& P)->double {return P[0] + P[1];},
        .right = {0, 1}
    },{
        .left = [](const ParamSet& P)->double {return P[0] - P[1];},
        .right = {0, 1}
    }
    };
    test.Init(100, make_shared<InitialDistributions>()
        << make_shared<DistribUniform>(-20, 20) << make_shared<DistribUniform>(-20, 20)
    );
    while (!test.ConcentratedInOnePoint())test.Iterate();
    EXPECT_TRUE(pow(test.Parameters()[0], 2) < 0.0000001);
    EXPECT_TRUE(pow(test.Parameters()[1], 2) < 0.0000001);
    for (const auto& eq : test.equations()) {
        EXPECT_EQ(0, eq.right.val());
        EXPECT_EQ(1, eq.right.uncertainty());
        EXPECT_TRUE(pow(eq.left(test.Parameters()).val(), 2) < 0.0000001);
    }
    const auto test_moved = std::move(test);
    EXPECT_TRUE(pow(test_moved.Parameters()[0], 2) < 0.0000001);
    EXPECT_TRUE(pow(test_moved.Parameters()[1], 2) < 0.0000001);
    for (const auto& eq : test_moved.equations()) {
        EXPECT_EQ(0, eq.right.val());
        EXPECT_EQ(1, eq.right.uncertainty());
        EXPECT_TRUE(pow(eq.left(test_moved.Parameters()).val(), 2) < 0.0000001);
    }
}
TEST(InexactEquationSolver, Integrationtest2)
{
    InexactEquationSolver<DifferentialMutations<>, UncertaintiesEstimation> test{
        {
        .left = [](const ParamSet& P)->double {return P[0] + P[1];},
        .right = {0, 1}
    },{
        .left = [](const ParamSet& P)->double {return P[0] - P[1];},
        .right = {0, 1}
    }
    };
    test.Init(100, make_shared<InitialDistributions>()
        << make_shared<DistribUniform>(-20, 20) << make_shared<DistribUniform>(-20, 20)
    );
    while (!test.ConcentratedInOnePoint())test.Iterate();
    EXPECT_TRUE(pow(test.Parameters()[0], 2) < 0.0000001);
    EXPECT_TRUE(pow(test.Parameters()[1], 2) < 0.0000001);
    test.SetUncertaintyCalcDeltas({ 0.01,0.01 });
    EXPECT_TRUE(test.ParametersWithUncertainties()[0].Contains(test.Parameters()[0]));
    EXPECT_TRUE(test.ParametersWithUncertainties()[1].Contains(test.Parameters()[1]));
    for (const auto& eq : test.equations()) {
        EXPECT_EQ(0, eq.right.val());
        EXPECT_EQ(1, eq.right.uncertainty());
        EXPECT_TRUE(pow(eq.left(test.Parameters()).val(), 2) < 0.0000001);
    }
    const auto test_moved = std::move(test);
    EXPECT_TRUE(pow(test_moved.Parameters()[0], 2) < 0.0000001);
    EXPECT_TRUE(pow(test_moved.Parameters()[1], 2) < 0.0000001);
    EXPECT_TRUE(test_moved.ParametersWithUncertainties()[0].Contains(test_moved.Parameters()[0]));
    EXPECT_TRUE(test_moved.ParametersWithUncertainties()[1].Contains(test_moved.Parameters()[1]));
    for (const auto& eq : test_moved.equations()) {
        EXPECT_EQ(0, eq.right.val());
        EXPECT_EQ(1, eq.right.uncertainty());
        EXPECT_TRUE(pow(eq.left(test_moved.Parameters()).val(), 2) < 0.0000001);
    }
}
