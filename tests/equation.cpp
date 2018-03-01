// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/error.h>
#include <Genetic/equation.h>
#include <Genetic/genetic.h>
#include <Genetic/initialconditions.h>
using namespace std;
using namespace MathTemplates;
using namespace Genetic;
TEST(EquationSystem, empty)
{
    EquationSystem A {};
    EXPECT_EQ(A({}), 0);
}
TEST(EquationSystem, empty2)
{
    EquationSystem A(vector<Equation> {});
    EXPECT_EQ(A({}), 0);
}
TEST(EquationSystem, simple)
{
    EquationSystem A {
        {.left = [](const ParamSet & P){return P[0];}, .right = [](const ParamSet &)
            {
                return 0;
            }
        }
    };
    EXPECT_EQ(A({0.0}), 0);
    EXPECT_EQ(A({0.5}), 0.25);
    EXPECT_EQ(A({1.0}), 1);
    EXPECT_EQ(A({2.0}), 4);
}
TEST(EquationSystem, simple2)
{
    EquationSystem A(vector<Equation> {
        {.left = [](const ParamSet &)
            {
                return 0;
            }, .right = [](const ParamSet & P)
            {
                return P[0];
            }
        }
    });
    EXPECT_EQ(A({0.0}), 0);
    EXPECT_EQ(A({0.5}), 0.25);
    EXPECT_EQ(A({1.0}), 1);
    EXPECT_EQ(A({2.0}), 4);
}
TEST(EquationSystem, twoparams)
{
    EquationSystem A {
        {.left = [](const ParamSet & P){return P[0] + P[1];}, .right = [](const ParamSet &)
            {
                return 0;
            }
        }
    };
    EXPECT_EQ(A({0.0, 0.0}), 0);
    EXPECT_EQ(A({0.5, 0.0}), 0.25);
    EXPECT_EQ(A({1.0, 0.0}), 1);
    EXPECT_EQ(A({2.0, 0.0}), 4);
    EXPECT_EQ(A({0.0, 0.0}), 0);
    EXPECT_EQ(A({0.0, 0.5}), 0.25);
    EXPECT_EQ(A({0.0, 1.0}), 1);
    EXPECT_EQ(A({0.0, 2.0}), 4);
}
TEST(EquationSystem, twoparams2)
{
    EquationSystem A(vector<Equation> {
        {.left = [](const ParamSet & P)
            {
                return P[0] + P[1];
            }, .right = [](const ParamSet &)
            {
                return 0;
            }
        }
    });
    EXPECT_EQ(A({0.0, 0.0}), 0);
    EXPECT_EQ(A({0.5, 0.0}), 0.25);
    EXPECT_EQ(A({1.0, 0.0}), 1);
    EXPECT_EQ(A({2.0, 0.0}), 4);
    EXPECT_EQ(A({0.0, 0.0}), 0);
    EXPECT_EQ(A({0.0, 0.5}), 0.25);
    EXPECT_EQ(A({0.0, 1.0}), 1);
    EXPECT_EQ(A({0.0, 2.0}), 4);
}
TEST(EquationSystem, two_eq)
{
    EquationSystem A {
        {.left = [](const ParamSet & P){return P[0] + P[1];}, .right = [](const ParamSet &)
            {
                return 0;
            }
        },
        {.left = [](const ParamSet & P){return P[0] - P[1];}, .right = [](const ParamSet &)
            {
                return 0;
            }
        }
    };
    EXPECT_EQ(A({0.0, 0.0}), 0);
    EXPECT_EQ(A({0.5, 0.0}), 0.5);
    EXPECT_EQ(A({1.0, 0.0}), 2);
    EXPECT_EQ(A({2.0, 0.0}), 8);
    EXPECT_EQ(A({0.0, 0.0}), 0);
    EXPECT_EQ(A({0.0, 0.5}), 0.5);
    EXPECT_EQ(A({0.0, 1.0}), 2);
    EXPECT_EQ(A({0.0, 2.0}), 8);
}
TEST(EquationSystem, two_eq2)
{
    EquationSystem A(vector<Equation> {
        {.left = [](const ParamSet & P)
            {
                return P[0] + P[1];
            }, .right = [](const ParamSet &)
            {
                return 0;
            }
        },
        {.left = [](const ParamSet & P)
            {
                return P[0] - P[1];
            }, .right = [](const ParamSet &)
            {
                return 0;
            }
        }
    });
    EXPECT_EQ(A({0.0, 0.0}), 0);
    EXPECT_EQ(A({0.5, 0.0}), 0.5);
    EXPECT_EQ(A({1.0, 0.0}), 2);
    EXPECT_EQ(A({2.0, 0.0}), 8);
    EXPECT_EQ(A({0.0, 0.0}), 0);
    EXPECT_EQ(A({0.0, 0.5}), 0.5);
    EXPECT_EQ(A({0.0, 1.0}), 2);
    EXPECT_EQ(A({0.0, 2.0}), 8);
}
TEST(EquationSolver, Integrationtest)
{
    EquationSolver<DifferentialMutations<>> test {
        {.left = [](const ParamSet & P){return P[0] + P[1];}, .right = [](const ParamSet &)
            {
                return 0;
            }
        },
        {.left = [](const ParamSet & P){return P[0];}, .right = [](const ParamSet & P)
            {
                return P[1];
            }
        }
    };
    test.Init(100, make_shared<InitialDistributions>()
              << make_shared<DistribGauss>(-20, 20) << make_shared<DistribGauss>(-20, 20)
    );
    while(!test.ConcentratedInOnePoint())test.Iterate();
    EXPECT_TRUE(pow(test.Parameters()[0], 2) < 0.0000001);
    EXPECT_TRUE(pow(test.Parameters()[1], 2) < 0.0000001);
    const auto &X = test.Parameters();
    for (const auto &eq : test.equations()) {
        EXPECT_TRUE(pow(eq.left(X) - eq.right(X), 2) < 0.0000001);
    }
}
