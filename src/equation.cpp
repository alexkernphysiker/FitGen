// this file is distributed under 
// MIT license
#include <math_h/error.h>
#include <Genetic/equation.h>
namespace Genetic{
    using namespace std;
    using namespace MathTemplates;
    EquationSystem::EquationSystem(const initializer_list<Equation>& source){
	for(const auto&item:source)f_data.push_back(item);
    }
    EquationSystem::EquationSystem(const list<Equation>& source){
	for(const auto&item:source)f_data.push_back(item);
    }
    EquationSystem::~EquationSystem(){}
    double EquationSystem::operator()(const ParamSet& P) const{
	double result=0;
	for(const auto&item:f_data)
	    result+=pow(item.left(P)-item.right(P),2);
	return result;
    }
    const list<Equation>&EquationSystem::equations()const{return f_data;}
}
