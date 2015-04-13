#include <iostream>
#include <fstream>
#include <equation.h>
#include <genetic.h>
#include <initialconditions.h>
using namespace std;
using namespace Genetic;
using namespace Genetic::wrap_func_indexer;
double Rosenbrock(double x, double y){
	return pow(1.0-x,2)+100*pow(y-pow(x,2),2);
}
int main(int argcnt, char **arg){
	auto solve=Solve<SearchMin<func2<Rosenbrock,par<0>,par<1>>>,DifferentialMutations<>>();
	solve->Init(100,make_shared<GenerateUniform>()<<make_pair(-50,50)<<make_pair(-50,50));
	printf("Population size: %i \n%i variables\n",solve->PopulationSize(),solve->ParamCount());
	while(!solve->ConcentratedInOnePoint()){
		solve->Iterate();
		printf("%i iterations\r",solve->iteration_count());
	}
	printf("\n");
	for(int i=0;i<solve->ParamCount();i++)
		printf("var%i = %f\n",i,solve->Parameters()[i]);
	return 0;
}