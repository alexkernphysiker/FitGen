#include <iostream>
#include <fstream>
#include <equation.h>
#include <initialconditions.h>
#include <additional_mutations.h>
using namespace std;
using namespace Fit;
using namespace Fit::wrap_func_indexer;
int main(int argcnt, char **arg){
	auto solve=Solve<SearchMin<
		add<
			power< sub<iconst<1>,par<0>>, iconst<2>>,
			mul<iconst<100>, power< sub< par<1>, power<par<0>,iconst<2>> >, iconst<2> > >
		>
	>>();
	solve->Init(100,make_shared<GenerateUniform>()<<make_pair(-80.0,80.0)<<make_pair(-80.0,80.0));
	printf("Population size: %i \n%i variables\n",solve->PopulationSize(),solve->ParamCount());
	while(!solve->ConcentratedInOnePoint())
		solve->Iterate();
	printf("%i iterations\n",solve->iteration_count());
	for(int i=0;i<solve->ParamCount();i++)
		printf("var%i = %f\n",i,solve->Parameters()[i]);
	return 0;
}