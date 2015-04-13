#include "initialconditions.h"
#include "math_h/randomfunc.h"
#include "fitexception.h"
namespace Genetic{
	using namespace std;
	Initialiser::Initialiser(){}
	Initialiser::~Initialiser(){}
	Initialiser& Initialiser::operator<<(generator gen){
		generators.push_back(gen);
		return *this;
	}
	generator Initialiser::operator[](int i){
		return generators[i];
	}
	int Initialiser::Count(){
		return generators.size();
	}
	ParamSet Initialiser::Generate(){
		ParamSet res;
		for(int i=0;i<Count();i++)
			res<<operator[](i)();
		return res;
	}
	Initialiser::iterator Initialiser::begin(){
		return generators.begin();
	}
	Initialiser::const_iterator Initialiser::cbegin() const{
		return  generators.cbegin();
	}
	Initialiser::iterator Initialiser::end(){
		return generators.end();
	}
	Initialiser::const_iterator Initialiser::cend() const{
		return generators.cend();
	}
	shared_ptr< Initialiser > operator<<(shared_ptr< Initialiser > init, generator gen){
		init->operator<<(gen);
		return init;
	}
	
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
		return *ParamDistr[i];
	}
	ParamSet InitialDistributions::Generate(){
		ParamSet res;
		for(int i=0;i<Count();i++)
			res<<operator[](i)();
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
		int r = m_min.size();
		if(r!=m_max.size())
			throw new FitException("ParamRangeUniform: unknown error. Min and max lists counts are different.");
		return r;
	}
	GenerateUniform &GenerateUniform::Add(double min, double max){
		m_min.push_back(min);
		m_max.push_back(max);
		return *this;
	}
	double GenerateUniform::Min(int i){
		if((i>=0)&(i<m_min.size()))
			return m_min[i];
		else
			throw new FitException("ParamRangeUniform: attempt to get property by index out of range");
	}
	double GenerateUniform::Max(int i){
		if((i>=0)&(i<m_max.size()))
			return m_max[i];
		else
			throw new FitException("ParamRangeUniform: attempt to get property by index out of range");
	}
	ParamSet GenerateUniform::Generate(){
		ParamSet res;
		for(int i=0; i<Count();i++)
			res<<RandomUniformlyR(m_min[i],m_max[i]);
		return res;
	}
	
	GenerateByGauss::GenerateByGauss(){}
	GenerateByGauss::~GenerateByGauss(){}
	int GenerateByGauss::Count(){
		int r = m_mean.size();
		if(r!=m_sig.size())
			throw new FitException("ParamGauss: unknown error. The sizes of private lists are different");
		return r;
	}
	GenerateByGauss &GenerateByGauss::Add(double mean, double sig){
		m_mean.push_back(mean);
		m_sig.push_back(sig);
		return *this;
	}
	double GenerateByGauss::Mean(int i){
		if((i>=0)&(i<m_mean.size()))
			return m_mean[i];
		else
			throw new FitException("ParamGauss: attempt to get property by index out of range");
	}
	double GenerateByGauss::Sigma(int i){
		if((i>=0)&(i<m_sig.size()))
			return m_sig[i];
		else
			throw new FitException("ParamGauss: attempt to get property by index out of range");
	}
	ParamSet GenerateByGauss::Generate(){
		ParamSet res;
		for(int i=0; i<Count();i++)
			res<<RandomGauss(m_sig[i],m_mean[i]);
		return res;
	}
}