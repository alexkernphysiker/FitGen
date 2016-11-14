// this file is distributed under 
// MIT license
#include <iostream>
#include <fstream>
#include <Genetic/equation2.h>
#include <Genetic/initialconditions.h>
#include <Genetic/genetic.h>
using namespace std;
using namespace Genetic;
int main(){
    RANDOM engine;
    InexactEquationSolver<DifferentialMutations<>> equation{
	in_eq([](const ParamSet&X){return X[0]+X[1];},{0.97,0.02}),
	in_eq([](const ParamSet&X){return X[0]-X[1];},{0.02,0.01}),
    };
    equation.Init(25,make_shared<GenerateUniform>()<<make_pair(-10,10)<<make_pair(-10,10),engine);
    while(!equation.AbsoluteOptimalityExitCondition(0.00000001))
	equation.Iterate(engine);
    cout<<endl<<"Solution: x="<<equation[0]<<"; y="<<equation[1]<<endl;
    cout<<"chi^2="<<equation.Optimality()<<endl;
    return 0;
}

