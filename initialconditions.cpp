// this file is distributed under 
// GPL v 3.0 license
#include "initialconditions.h"
#include "math_h/exception_math_h.h"
namespace Genetic{
	using namespace std;
	InitialDistributions::InitialDistributions(){}
	InitialDistributions::~InitialDistributions(){}
	InitialDistributions& InitialDistributions::operator<<(shared_ptr<Distrib> distr){
		ParamDistr.push_back(distr);
		return *this;
	}
	int InitialDistributions::Count(){
		return ParamDistr.size();
	}
	Distrib &InitialDistributions::operator[](int i){
		if((i<0)||(i>=ParamDistr.size()))
			throw math_h_error<InitialDistributions>("Range check error when accessing InitialDistributions' element");
		return *ParamDistr[i];
	}
	ParamSet InitialDistributions::Generate(RANDOM&R){
		ParamSet res;
		for(int i=0;i<Count();i++)
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
	int GenerateUniform::Count(){
		return m_min.size();
	}
	GenerateUniform &GenerateUniform::Add(double min, double max){
		if(min>max)
			throw math_h_error<GenerateUniform>("ParamRangeUniform: Max<Min");
		m_min.push_back(min);
		m_max.push_back(max);
		return *this;
	}
	double GenerateUniform::Min(int i){
		if((i<0)||(i>=m_min.size()))
			throw math_h_error<GenerateUniform>("Range check error when accessing GenerateUniform's element");
		return m_min[i];
	}
	double GenerateUniform::Max(int i){
		if((i<0)||(i>=m_max.size()))
			throw math_h_error<GenerateUniform>("Range check error when accessing GenerateUniform's element");
		return m_max[i];
	}
	ParamSet GenerateUniform::Generate(RANDOM&R){
		ParamSet res;
		for(int i=0; i<Count();i++){
			uniform_real_distribution<double> D(m_min[i],m_max[i]);
			res<<D(R);
		}
		return res;
	}
	
	GenerateByGauss::GenerateByGauss(){}
	GenerateByGauss::~GenerateByGauss(){}
	int GenerateByGauss::Count(){
		return m_mean.size();
	}
	GenerateByGauss &GenerateByGauss::Add(double mean, double sig){
		m_mean.push_back(mean);
		m_sig.push_back(sig);
		return *this;
	}
	double GenerateByGauss::Mean(int i){
		if((i<0)||(i>=m_mean.size()))
			throw math_h_error<GenerateByGauss>("Range check error when accessing GenerateByGauss's element");
		return m_mean[i];
	}
	double GenerateByGauss::Sigma(int i){
		if((i<0)||(i>=m_sig.size()))
			throw math_h_error<GenerateByGauss>("Range check error when accessing GenerateByGauss's element");
		return m_sig[i];
	}
	ParamSet GenerateByGauss::Generate(RANDOM&R){
		ParamSet res;
		for(int i=0; i<Count();i++){
			normal_distribution<double> D(m_mean[i],m_sig[i]);
			res<<D(R);
		}
		return res;
	}
}