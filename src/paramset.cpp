// this file is distributed under 
// MIT license
#include <math_h/exception_math_h.h>
#include <Genetic/paramset.h>
namespace Genetic{
	using namespace std;
	typedef lock_guard<mutex> Lock;
	ParamSet::ParamSet(){}
	ParamSet::ParamSet(const initializer_list< double >& source){
		for(auto value:source)m_values.push_back(value);
	}
	ParamSet::ParamSet(initializer_list< double >&& source):ParamSet(source){}
	ParamSet::ParamSet(const vector< double >& source){
		for(auto value:source)m_values.push_back(value);
	}
	ParamSet::ParamSet(vector<double>&&source):ParamSet(source){}
	ParamSet::ParamSet(const ParamSet&source):ParamSet(source.m_values){}
	ParamSet::ParamSet(ParamSet&&source):ParamSet(source.m_values){}
	ParamSet::~ParamSet(){}

	size_t ParamSet::size()const{return m_values.size();}
	double ParamSet::operator[](size_t i)const{
		if(i>=m_values.size())
			throw math_h_error<ParamSet>("Range check error when accessing ParamSet's element");
		return m_values[i];
	}
	double& ParamSet::operator[](size_t i){
		if(i>=m_values.size())
			throw math_h_error<ParamSet>("Range check error when accessing ParamSet's element");
		return m_values[i];
	}

	ParamSet& ParamSet::operator<<(double p){
		Lock lock(m_mutex);
		m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator<<(const initializer_list< double >& source){
		Lock lock(m_mutex);
		for(double p:source)m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator<<(initializer_list<double>&&source){return operator<<(source);}
	ParamSet& ParamSet::operator<<(const vector<double>&V){
		Lock lock(m_mutex);
		for(double p:V)m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator<<(vector<double>&&V){return operator<<(V);}
	ParamSet& ParamSet::operator<<(const ParamSet&P){return operator<<(P.m_values);}
	ParamSet& ParamSet::operator<<(ParamSet&&P){return operator<<(P);}
	
	ParamSet& ParamSet::operator=(const initializer_list< double >& source){
		Lock lock(m_mutex);
		m_values.clear();
		for(double p:source)m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator=(initializer_list<double>&&source){return operator=(source);}
	ParamSet& ParamSet::operator=(const vector< double >&V){
		Lock lock(m_mutex);
		m_values.clear();
		for(double p:V)m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator=(vector< double >&& V){return operator=(V);}
	ParamSet& ParamSet::operator=(const ParamSet& P){return operator=(P.m_values);}
	ParamSet& ParamSet::operator=(ParamSet&& P){return operator=(P);}

	ParamSet::iterator ParamSet::begin(){return m_values.begin();}
	ParamSet::const_iterator ParamSet::begin()const{return m_values.begin();}
	ParamSet::const_iterator ParamSet::cbegin()const{return m_values.cbegin();}
	ParamSet::iterator ParamSet::end(){return m_values.end();}
	ParamSet::const_iterator ParamSet::end() const{return m_values.end();}
	ParamSet::const_iterator ParamSet::cend() const{return m_values.cend();}
	
	ostream& operator<<(ostream& str, const ParamSet& P){
		for(double x:P)str<<x<<"\t";
		return str;
	}
	istream& operator>>(istream& str, ParamSet& P){
		double x;
		while(str>>x)P<<x;
		return str;
	}
	
	ParamSet parEq(size_t cnt,double val){
		ParamSet res;
		for(size_t i=0;i<cnt;i++)
			res<<val;
		return res;
	}
}