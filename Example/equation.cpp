// this file is distributed under 
// MIT license
#include <iostream>
#include <fstream>
#include <Genetic/equation.h>
#include <Genetic/initialconditions.h>
using namespace std;
using namespace Genetic;
int main(int argcnt, char **arg){
	RANDOM engine;
	Equation<DifferentialMutations<>> solution([](const ParamSet&X){
		return pow(1.0-X[0],2)+100*pow(X[1]-pow(X[0],2),2);
	});
	solution.Init(100,make_shared<GenerateUniform>()<<make_pair(-50,50)<<make_pair(-50,50),engine);
	printf("Population size: %i \n%i variables\n",solution.PopulationSize(),solution.ParamCount());
	Find(solution,engine);
	for(int i=0;i<solution.ParamCount();i++)
		printf("var%i = %f\n",i,solution[i]);
	return 0;
}