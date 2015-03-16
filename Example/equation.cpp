#include <iostream>
#include <fstream>
#include <equation.h>
#include <initialconditions.h>
#include <additional_mutations.h>
using namespace std;
using namespace Fit;
double Rosenbrock(ParamSet &X){
	return pow(1.0-X[0],2)+100*pow(X[1]-pow(X[0],2),2);
}
int main(int argcnt, char **arg){
	auto solve=Solve<SearchMin<Rosenbrock>>();
	solve->Init(100,make_shared<GenerateUniform>()<<make_pair(-80.0,80.0)<<make_pair(-80.0,80.0));
	printf("Population size: %i \n%i variables\n",solve->PopulationSize(),solve->ParamCount());
	while(!solve->ConcentratedInOnePoint())
		solve->Iterate();
	printf("%i iterations\n",solve->iteration_count());
	for(int i=0;i<solve->ParamCount();i++)
		printf("var%i = %f\n",i,solve->Parameters()[i]);
	return 0;
}