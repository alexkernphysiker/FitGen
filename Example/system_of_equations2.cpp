#include <iostream>
#include <fstream>
#include <Genetic/equation2.h>
#include <Genetic/initialconditions.h>
#include <Genetic/genetic.h>
#include <Genetic/uncertainties.h>
using namespace std;
using namespace Genetic;
int main()
{
    //declaring
    InexactEquationSolver<DifferentialMutations<>,UncertaintiesEstimation> example {
        {
            .left = [](const ParamSet & X){return X[0] + X[1];},
            .right = {1.98, 0.03}
        }, {
            .left = [](const ParamSet & X){return X[0] - X[1];},
            .right = {0.02, 0.02}
        }, {
            .left = [](const ParamSet & X){return pow(X[0], 2) + pow(X[1], 2);},
            .right = {2.01, 0.01}
        }
    };

    //searching the optimal parameter values
    example.Init(
        25,
        make_shared<InitialDistributions>()
        << make_shared<DistribUniform>(-10, 10)
        << make_shared<DistribUniform>(-10, 10)
    );
    while (!example.AbsoluteOptimalityExitCondition(0.0001))
        example.Iterate();

    //output the solution
    example.SetUncertaintyCalcDeltas({0.001, 0.001});
    const auto &X = example.ParametersWithUncertainties();
    cout << endl << "Solution: x0={" << X[0] << "}; x1={" << X[1] << "}" << endl;
    cout << "chi^2=" << example.Optimality() << endl;
    return 0;
}
