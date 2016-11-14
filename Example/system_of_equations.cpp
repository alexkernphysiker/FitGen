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
	eq([](const ParamSet&X){return X[0]+X[1];},1.0),
	eq([](const ParamSet&X){return X[0]-X[1];},0.0),
    };
    equation.Init(25,make_shared<GenerateUniform>()<<make_pair(-10,10)<<make_pair(-10,10),engine);
    Find(equation,engine);
    cout<<endl<<"Solution: x="<<equation[0]<<"; y="<<equation[1]<<endl;
    return 0;
}

