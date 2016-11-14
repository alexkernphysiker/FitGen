// this file is distributed under 
// MIT license
#include <iostream>
#include <fstream>
#include <Genetic/equation.h>
#include <Genetic/initialconditions.h>
#include <Genetic/genetic.h>
using namespace std;
using namespace Genetic;
int main(){
    RANDOM engine;
    EquationSolver<DifferentialMutations<>> equation{
	eq([](const ParamSet&X){return pow(X[0],2)+pow(X[1],2);},1.0),
	eq([](const ParamSet&X){return X[0]-X[1];},0.0),
    };
    equation.Init(100,make_shared<GenerateUniform>()<<make_pair(-10,10)<<make_pair(-10,10),engine);
    while(!equation.AbsoluteOptimalityExitCondition(0.00000001))
	equation.Iterate(engine);
    cout<<endl<<"Solution: x="<<equation[0]<<"; y="<<equation[1]<<endl;
    cout<<"Optimality="<<equation.Optimality()<<endl;
    return 0;
}

