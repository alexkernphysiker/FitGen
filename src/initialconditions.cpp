// this file is distributed under 
// MIT license
#include <math_h/randomfunc.h>
#include <math_h/error.h>
#include <Genetic/initialconditions.h>
namespace Genetic{
	using namespace std;
	using namespace MathTemplates;
	InitialDistributions::InitialDistributions(){}
	InitialDistributions::~InitialDistributions(){}
	InitialDistributions& InitialDistributions::operator<<(const shared_ptr<Distrib> distr){
		ParamDistr.push_back(distr);
		return *this;
	}
	const size_t InitialDistributions::Count()const {
		return ParamDistr.size();
	}
	const Distrib &InitialDistributions::operator[](const size_t i)const {
		if(i>=ParamDistr.size())
			throw Exception<InitialDistributions>("Range check error when accessing InitialDistributions' element");
		return *ParamDistr[i];
	}
	ParamSet InitialDistributions::Generate(RANDOM&R)const {
		ParamSet res;
		for(size_t i=0;i<Count();i++)
			res<<operator[](i)(R);
		return res;
	}
	shared_ptr<InitialDistributions> operator<<(
		  shared_ptr<InitialDistributions> init, 
		  const shared_ptr<Distrib> func
	){
		init->operator<<(func);
		return init;
	}
}
