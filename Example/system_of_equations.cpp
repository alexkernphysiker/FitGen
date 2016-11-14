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
    //declaring
    EquationSolver<DifferentialMutations<>> example{
	{.left=[](const ParamSet&X){return X[0]+X[1];},.right=1},
	{.left=[](const ParamSet&X){return X[0]-X[1];},.right=0}
    };
    //searching the solution
    RANDOM random_engine;
    example.Init(30,make_shared<GenerateUniform>()<<make_pair(-10,10)<<make_pair(-10,10),random_engine);
    Find(example,random_engine);
    //output results
    const auto&X=example.Parameters();
    cout<<endl<<"Solution: x0="<<X[0]<<"; x1="<<X[1]<<endl;
    for(const auto&eq:example.equations())
	cout<<eq.left(X)<<"=="<<eq.right<<endl;
    return 0;
}

