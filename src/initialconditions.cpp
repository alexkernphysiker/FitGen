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
	InitialDistributions& InitialDistributions::operator<<(shared_ptr<Distrib> distr){
		ParamDistr.push_back(distr);
		return *this;
	}
	size_t InitialDistributions::Count(){
		return ParamDistr.size();
	}
	Distrib &InitialDistributions::operator[](size_t i){
		if(i>=ParamDistr.size())
			throw Exception<InitialDistributions>("Range check error when accessing InitialDistributions' element");
		return *ParamDistr[i];
	}
	ParamSet InitialDistributions::Generate(RANDOM&R){
		ParamSet res;
		for(size_t i=0;i<Count();i++)
			res<<operator[](i)(R);
		return res;
	}
	shared_ptr<InitialDistributions> operator<<(
		  shared_ptr<InitialDistributions> init, 
		  shared_ptr<Distrib> func
	){
		init->operator<<(func);
		return init;
	}
	
	GenerateUniform::GenerateUniform(){}
	GenerateUniform::~GenerateUniform(){}
	size_t GenerateUniform::Count(){
		return m_min.size();
	}
	GenerateUniform &GenerateUniform::Add(double min, double max){
		if(min>max)
			throw Exception<GenerateUniform>("ParamRangeUniform: Max<Min");
		m_min.push_back(min);
		m_max.push_back(max);
		return *this;
	}
	double GenerateUniform::Min(size_t i){
		if(i>=m_min.size())
			throw Exception<GenerateUniform>("Range check error when accessing GenerateUniform's element");
		return m_min[i];
	}
	double GenerateUniform::Max(size_t i){
		if(i>=m_max.size())
			throw Exception<GenerateUniform>("Range check error when accessing GenerateUniform's element");
		return m_max[i];
	}
	ParamSet GenerateUniform::Generate(RANDOM&R){
		ParamSet res;
		for(size_t i=0; i<Count();i++){
			uniform_real_distribution<double> D(m_min[i],m_max[i]);
			res<<D(R);
		}
		return res;
	}
	
	GenerateByGauss::GenerateByGauss(){}
	GenerateByGauss::~GenerateByGauss(){}
	size_t GenerateByGauss::Count(){
		return m_mean.size();
	}
	GenerateByGauss &GenerateByGauss::Add(double mean, double sig){
		m_mean.push_back(mean);
		m_sig.push_back(sig);
		return *this;
	}
	double GenerateByGauss::Mean(size_t i){
		if(i>=m_mean.size())
			throw Exception<GenerateByGauss>("Range check error when accessing GenerateByGauss's element");
		return m_mean[i];
	}
	double GenerateByGauss::Sigma(size_t i){
		if(i>=m_sig.size())
			throw Exception<GenerateByGauss>("Range check error when accessing GenerateByGauss's element");
		return m_sig[i];
	}
	ParamSet GenerateByGauss::Generate(RANDOM&R){
		ParamSet res;
		for(size_t i=0; i<Count();i++){
			normal_distribution<double> D(m_mean[i],m_sig[i]);
			res<<D(R);
		}
		return res;
	}
}
