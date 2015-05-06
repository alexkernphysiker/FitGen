#include "initialconditions.h"
#include "math_h/randomfunc.h"
#include "genetic_exception.h"
namespace Genetic{
	using namespace std;
	Initialiser::Initialiser(){}
	Initialiser::~Initialiser(){}
	Initialiser& Initialiser::operator<<(generator gen){
		generators.push_back(gen);
		return *this;
	}
	generator Initialiser::operator[](int i){
		if(i<0)
			throw new GeneticException("Initialiser: negative index");
		if(i>=generators.size())
			throw new GeneticException("Initialiser: index out of range");
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
		if(i<0)
			throw new GeneticException("InitialDistributions: negative index");
		if(i>=ParamDistr.size())
			throw new GeneticException("InitialDistributions: index out of range");
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
			throw new GeneticException("ParamRangeUniform: unknown error. Min and max lists counts are different.");
		return r;
	}
	GenerateUniform &GenerateUniform::Add(double min, double max){
		if(min>max)
			throw new GeneticException("ParamRangeUniform: Max<Min");
		m_min.push_back(min);
		m_max.push_back(max);
		return *this;
	}
	double GenerateUniform::Min(int i){
		if(i<0)
			throw new GeneticException("ParamRangeUniform: negative index");
		if(i>=m_min.size())
			throw new GeneticException("ParamRangeUniform: index out of range");
		return m_min[i];
	}
	double GenerateUniform::Max(int i){
		if(i<0)
			throw new GeneticException("ParamRangeUniform: negative index");
		if(i>=m_max.size())
			throw new GeneticException("ParamRangeUniform: index out of range");
		return m_max[i];
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
			throw new GeneticException("ParamGauss: unknown error. The sizes of private lists are different");
		return r;
	}
	GenerateByGauss &GenerateByGauss::Add(double mean, double sig){
		m_mean.push_back(mean);
		m_sig.push_back(sig);
		return *this;
	}
	double GenerateByGauss::Mean(int i){
		if(i<0)
			throw new GeneticException("ParamGauss: negative index");
		if(i>=m_mean.size())
			throw new GeneticException("ParamGauss: index out of range");
		return m_mean[i];
	}
	double GenerateByGauss::Sigma(int i){
		if(i<0)
			throw new GeneticException("ParamGauss: negative index");
		if(i>=m_sig.size())
			throw new GeneticException("ParamGauss: index out of range");
		return m_sig[i];
	}
	ParamSet GenerateByGauss::Generate(){
		ParamSet res;
		for(int i=0; i<Count();i++)
			res<<RandomGauss(m_sig[i],m_mean[i]);
		return res;
	}
}