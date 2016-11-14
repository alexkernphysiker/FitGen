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
    InexactEquationSolver<DifferentialMutations<>> example{
	{.left=[](const ParamSet&X){return X[0]+X[1];},.right={1.98,0.03}},
	{.left=[](const ParamSet&X){return X[0]-X[1];},.right={0.02,0.02}},
	{.left=[](const ParamSet&X){return pow(X[0],2)+pow(X[1],2);},.right={2.01,0.01}}
    };
    RANDOM random_engine;
    example.Init(25,make_shared<GenerateUniform>()<<make_pair(-10,10)<<make_pair(-10,10),random_engine);
    while(!example.AbsoluteOptimalityExitCondition(0.00000001))
	example.Iterate(random_engine);
    example.SetUncertaintyCalcDeltas({0.001,0.001});
    auto X=example.ParametersWithUncertainties();
    cout<<endl<<"Solution: x0={"<<X[0]<<"}; x1={"<<X[1]<<"}"<<endl;
    cout<<"chi^2="<<example.Optimality()<<endl;
    for(const auto&eq:example.equations()){
	cout<<"{"<<eq.left(example.Parameters())<<"}=={"<<eq.right<<"}"<<endl;
    }
    return 0;
}

