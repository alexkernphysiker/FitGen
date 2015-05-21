#include "initialconditions.h"
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
		if((i<0)||(i>=generators.size()))
			throw GeneticException("Range chack error when accessing Initialiser's element");
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
		if((i<0)||(i>=ParamDistr.size()))
			throw GeneticException("Range check error when accessing InitialDistributions' element");
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
		return m_min.size();
	}
	GenerateUniform &GenerateUniform::Add(double min, double max){
		if(min>max)
			throw GeneticException("ParamRangeUniform: Max<Min");
		m_min.push_back(min);
		m_max.push_back(max);
		return *this;
	}
	double GenerateUniform::Min(int i){
		if((i<0)||(i>=m_min.size()))
			throw GeneticException("Range check error when accessing GenerateUniform's element");
		return m_min[i];
	}
	double GenerateUniform::Max(int i){
		if((i<0)||(i>=m_max.size()))
			throw GeneticException("Range check error when accessing GenerateUniform's element");
		return m_max[i];
	}
	ParamSet GenerateUniform::Generate(){
		ParamSet res;
		for(int i=0; i<Count();i++){
			uniform_real_distribution<double> D(m_min[i],m_max[i]);
			res<<D(G);
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
			throw GeneticException("Range check error when accessing GenerateByGauss's element");
		return m_mean[i];
	}
	double GenerateByGauss::Sigma(int i){
		if((i<0)||(i>=m_sig.size()))
			throw GeneticException("Range check error when accessing GenerateByGauss's element");
		return m_sig[i];
	}
	ParamSet GenerateByGauss::Generate(){
		ParamSet res;
		for(int i=0; i<Count();i++){
			normal_distribution<double> D(m_mean[i],m_sig[i]);
			res<<D(G);
		}
		return res;
	}
}