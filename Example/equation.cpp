#include <iostream>
#include <fstream>
#include <fit_gen.h>
#include <equation.h>
#include <initialconditions.h>
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
	auto initial=make_shared<GenerateUniform>();
	initial->Add(-20,20).Add(-20,20);
	solve->Init(30,initial);
	do{solve->Iterate();}while(solve->Optimality()<solve->Optimality(solve->PopulationSize()-1));
	printf("%i iterations\n",solve->iteration_count());
	printf("x=%f; y=%f\n",solve->Parameters()[0],solve->Parameters()[1]);
	return 0;
}