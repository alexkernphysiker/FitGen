// this file is distributed under 
// MIT license
#include <iostream>
#include <fstream>
#include <Genetic/equation2.h>
#include <Genetic/initialconditions.h>
#include <Genetic/genetic.h>
#include <Genetic/parabolic.h>
using namespace std;
using namespace Genetic;
int main(){
    //declaring
    InexactEquationSolver<DifferentialMutations<Uncertainty>> example{
	{.left=[](const ParamSet&X){return X[0]+X[1];},.right={1.98,0.03}},
	{.left=[](const ParamSet&X){return X[0]-X[1];},.right={0.02,0.02}},
	{.left=[](const ParamSet&X){return pow(X[0],2)+pow(X[1],2);},.right={2.01,0.01}}
    };
    //searching the optimal parameter values
    RANDOM random_engine;
    example.Init(25,make_shared<InitialDistributions>()
	<<make_shared<DistribGauss>(-10,10)<<make_shared<DistribGauss>(-10,10)
    ,random_engine);
    while(!example.AbsoluteOptimalityExitCondition(0.00000001))
	example.Iterate(random_engine);
    //output the solution
    example.SetUncertaintyCalcDeltas({0.001,0.001});
    const auto&x=example.Parameters();
    const auto&X=example.ParametersWithUncertainties();
    cout<<endl<<"Solution: x0={"<<X[0]<<"}; x1={"<<X[1]<<"}"<<endl;
    cout<<"chi^2="<<example.Optimality()<<endl;
    for(const auto&eq:example.equations())
	cout<<eq.left(x).val()<<"=={"<<eq.right<<"}"<<endl;
    return 0;
}

