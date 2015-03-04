#include <iostream>
#include <fstream>
#include <fit_gen.h>
#include <initialconditions.h>
#include <equation.h>
using namespace std;
using namespace Fit;
using namespace Fit::FuncWrappers_P;
int main(int argcnt, char **arg){
	auto empty=make_shared<NoParamFunc>();
	auto eq=make_shared<SearchMin<add<
		power< sub<iconst<1>,par<0>>, iconst<2>>,
		mul<iconst<100>, power< sub< par<1>, power<par<0>,iconst<2>> >, iconst<2> > >
	>>>();
	auto initial=make_shared<GenerateUniform>();
	initial->Add(-20,20).Add(-20,20);
	FitGenVeg solve(empty,eq);
	solve.Init(40,initial);
	solve.SetMutation(Fit::mutDifferential,ParamSet(1,1));
	solve.SetMutation(Fit::mutRatio,ParamSet(0,0));
	solve.SetMutation(Fit::mutAbsolute,ParamSet(0,0));
	do{solve.Iterate();}while(solve.GetOptimality()<solve.GetOptimality(solve.PopulationSize()-1));
	printf("%i iterations\n",solve.iteration_count());
	printf("x=%f; y=%f\n",solve[0],solve[1]);
	return 0;
}