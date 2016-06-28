// this file is distributed under 
// MIT license
#include <iostream>
#include <fstream>
#include <Genetic/equation.h>
#include <Genetic/initialconditions.h>
using namespace std;
using namespace Genetic;
int main(){
	RANDOM engine;
	Equation<DifferentialMutations<>> solution(
		[](const ParamSet&X){
			return pow(1.0-X[0],2)+100*pow(X[1]-pow(X[0],2),2);
		}
	);
	solution.Init(100,
		make_shared<GenerateUniform>()
			<<make_pair(-50,50)
			<<make_pair(-50,50),
		engine
	);
	Find(solution,engine);
	cout<<endl<<"Solution:\t"<<solution<<endl;
	return 0;
}