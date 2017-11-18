// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/interpolate.h>
#include <math_h/error.h>
#include <math_h/functions.h>
#include <Genetic/initialconditions.h>
using namespace std;
using namespace MathTemplates;
using namespace Genetic;
const int n = 5;
RANDOM engine;
shared_ptr<Distrib> Distrs[] = {
    make_shared<DistribGauss>(0, 1), make_shared<DistribGauss>(1, 1),
    make_shared<DistribGauss>(2, 1), make_shared<DistribGauss>(3, 1),
    make_shared<DistribGauss>(4, 1)
};
TEST(InitialDistributions, Create)
{
    InitialDistributions I;
    EXPECT_EQ(0, I.Count());
}
TEST(InitialDistributions, Add)
{
    for (int count = 1; count < n; count++) {
        InitialDistributions I;
        for (int i = 0; i < count; i++)
            EXPECT_EQ(&I, &(I << Distrs[i]));
        EXPECT_EQ(count, I.Count());
        for (int i = 0; i < count; i++)
            EXPECT_EQ(Distrs[i].get(), &I[i]);
        EXPECT_THROW(I[count](engine), Exception<InitialDistributions>);
        EXPECT_THROW(I[-1](engine), Exception<InitialDistributions>);
    }
    for (int count = 1; count < n; count++) {
        auto I = make_shared<InitialDistributions>();
        for (int i = 0; i < count; i++)
            EXPECT_EQ(I.get(), (I << Distrs[i]).get());
        EXPECT_EQ(count, I->Count());
        for (int i = 0; i < count; i++)
            EXPECT_EQ(Distrs[i].get(), &(I->operator[](i)));
    }
}
TEST(InitialDistributions, Generate)
{
    for (int count = 1; count < 5; count++) {
        InitialDistributions I;
        for (int i = 0; i < count; i++)I << make_shared<DistribUniform>(0, 1);
        ParamSet P = I.Generate(engine);
        EXPECT_EQ(count, P.size());
        for (int i = 0; i < count; i++)
            EXPECT_EQ(true, (P[i] >= 0) && (P[i] <= 1));
    }
}
