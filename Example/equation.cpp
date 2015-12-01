// this file is distributed under 
// MIT license
#include <iostream>
#include <fstream>
#include <equation.h>
#include <initialconditions.h>
using namespace std;
using namespace Genetic;
int main(int argcnt, char **arg){
	RANDOM engine;
	Equation<DifferentialMutations<>> solve([](const ParamSet&X){
		return pow(1.0-X[0],2)+100*pow(X[1]-pow(X[0],2),2);
	});
	solve.Init(100,make_shared<GenerateUniform>()<<make_pair(-50,50)<<make_pair(-50,50),engine);
	printf("Population size: %i \n%i variables\n",solve.PopulationSize(),solve.ParamCount());
	while(!solve.ConcentratedInOnePoint()){
		solve.Iterate(engine);
		printf("%i iterations\r",solve.iteration_count());
	}
	printf("\n");
	for(int i=0;i<solve.ParamCount();i++)
		printf("var%i = %f\n",i,solve[i]);
	return 0;
}