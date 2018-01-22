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
            .left = [](const ParamSet & X){return pow(X[0], 2) + 1.0;},
            .right = [](const ParamSet & X){return pow(X[1], 2);}
        }, {
            .left = [](const ParamSet & X){return X[0] + X[1];},
            .right = [](const ParamSet &){return 2;}
        }
    };
    //searching the solution
    solution.Init(20, make_shared<InitialDistributions>()
                 << make_shared<DistribUniform>(0, 3)
                 << make_shared<DistribUniform>(0, 3));
    Find(solution);
    //output results
    const auto &X = solution.Parameters();
    cout << endl << "Solution: x0=" << X[0] << "; x1=" << X[1] << endl;
    cout << solution.equations()[0].left(X) << "==" << solution.equations()[0].right(X) << endl;
    cout << solution.equations()[1].left(X) << "==" << solution.equations()[1].right(X) << endl;
    return 0;
}


