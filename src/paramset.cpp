// this file is distributed under 
// MIT license
#include <math_h/error.h>
#include <Genetic/paramset.h>
namespace Genetic{
	using namespace std;
	using namespace MathTemplates;
	typedef lock_guard<mutex> Lock;
	ParamSet::ParamSet(){}
	ParamSet::ParamSet(const initializer_list< double >& source){
		for(auto value:source)m_values.push_back(value);
	}
	ParamSet::ParamSet(const initializer_list< double >&& source):ParamSet(source){}
	ParamSet::ParamSet(const vector< double >& source){
		for(auto value:source)m_values.push_back(value);
	}
	ParamSet::ParamSet(const vector<double>&&source):ParamSet(source){}
	ParamSet::ParamSet(const ParamSet&source):ParamSet(source.m_values){}
	ParamSet::ParamSet(const ParamSet&&source):ParamSet(source.m_values){}
	ParamSet::~ParamSet(){}

	const size_t ParamSet::size()const{return m_values.size();}
	double ParamSet::operator[](const size_t i)const{
		if(i>=m_values.size())
			throw Exception<ParamSet>("Range check error when accessing ParamSet's element");
		return m_values[i];
	}
	double& ParamSet::operator()(const size_t i){
		if(i>=m_values.size())
			throw Exception<ParamSet>("Range check error when accessing ParamSet's element");
		return m_values[i];
	}

	ParamSet& ParamSet::operator<<(const double p){
		Lock lock(m_mutex);
		m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator<<(const initializer_list< double >& source){
		Lock lock(m_mutex);
		for(double p:source)m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator<<(const initializer_list<double>&&source){return operator<<(source);}
	ParamSet& ParamSet::operator<<(const vector<double>&V){
		Lock lock(m_mutex);
		for(double p:V)m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator<<(const vector<double>&&V){return operator<<(V);}
	ParamSet& ParamSet::operator<<(const ParamSet&P){return operator<<(P.m_values);}
	ParamSet& ParamSet::operator<<(const ParamSet&&P){return operator<<(P);}
	
	ParamSet& ParamSet::operator>>(double& p){
		Lock lock(m_mutex);
		if(m_values.size()==0)
			throw Exception<ParamSet>("Attempt to take a value from empty paramset");
		p=m_values[m_values.size()-1];
		m_values.pop_back();
		return *this;
	}
	
	ParamSet& ParamSet::operator=(const initializer_list< double >& source){
		Lock lock(m_mutex);
		m_values.clear();
		for(double p:source)m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator=(const initializer_list<double>&&source){return operator=(source);}
	ParamSet& ParamSet::operator=(const vector<double>&V){
		Lock lock(m_mutex);
		m_values.clear();
		for(double p:V)m_values.push_back(p);
		return *this;
	}
	ParamSet& ParamSet::operator=(const vector< double >&& V){return operator=(V);}
	ParamSet& ParamSet::operator=(const ParamSet& P){return operator=(P.m_values);}
	ParamSet& ParamSet::operator=(const ParamSet&& P){return operator=(P);}

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
	
	ParamSet parEq(const size_t cnt,const double val){
		ParamSet res;
		for(size_t i=0;i<cnt;i++)
			res<<val;
		return res;
	}
}