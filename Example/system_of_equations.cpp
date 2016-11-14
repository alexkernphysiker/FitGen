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
    EquationSolver<DifferentialMutations<>> example{
	{.left=[](const ParamSet&X){return X[0]+X[1];},.right=1},
	{.left=[](const ParamSet&X){return X[0]-X[1];},.right=0}
    };
    RANDOM random_engine;
    example.Init(25,make_shared<GenerateUniform>()<<make_pair(-10,10)<<make_pair(-10,10),random_engine);
    Find(example,random_engine);
    cout<<endl<<"Solution: x0="<<example.Parameters()[0]<<"; x1="<<example.Parameters()[1]<<endl;
    for(const auto&eq:example.equations()){
	cout<<eq.left(example.Parameters())<<"=="<<eq.right<<endl;
    }
    return 0;
}

