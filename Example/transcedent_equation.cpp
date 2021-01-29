#include <math.h>
#include <iostream>
#include <fstream>
#include <Genetic/equation.h>
#include <Genetic/initialconditions.h>
#include <Genetic/genetic.h>
using namespace std;
using namespace Genetic;
int main()
{
    //declaring
    EquationSolver<DifferentialMutations<>> solution {
        {
            .left = [](const ParamSet& X){return sin(X[0]/2.0);},
            .right = [](const ParamSet& X){return exp(-X[0]);}
        }
    };
    //searching the solution
    solution.Init(7, make_shared<InitialDistributions>()<< make_shared<DistribUniform>(-2,2));
    while (!solution.AbsoluteOptimalityExitCondition(0.000000000001))solution.Iterate();
    //output results
    cout << endl << "Solution: x=" << solution.Parameters()[0]<< endl;
    cout << endl << "Sigma: " << solution.ParametersStatistics()[0].uncertainty()<< endl;
    return 0;
}


