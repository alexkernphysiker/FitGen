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
    Equation<DifferentialMutations<>> equation(
	[](const ParamSet&X){return pow(1.0-X[0],2);},//left part
	[](const ParamSet&X){return -100*pow(X[1]-pow(X[0],2),2);}//right part
    );
    equation.Init(100,
	make_shared<GenerateUniform>()<<make_pair(-50,50)<<make_pair(-50,50),
	engine
    );
    Find(equation,engine);
    cout<<endl<<"Solution: x="<<equation[0]<<"; y="<<equation[1]<<endl;
    return 0;
}
