// this file is distributed under 
// MIT license
#include <math_h/exception_math_h.h>
#include <Genetic/paramset.h>
namespace Genetic{
	using namespace std;
	typedef lock_guard<mutex> Lock;
	ParamSet::ParamSet(){}
	ParamSet::ParamSet(double x):ParamSet(){
		operator<<(x);
	}
	ParamSet::ParamSet(double x,double y):ParamSet(){
		operator<<(x)<<y;
	}
	ParamSet::ParamSet(double x,double y,double z):ParamSet(){
		operator<<(x)<<y<<z;
	}
	ParamSet::ParamSet(double x,double y,double z,double zz):ParamSet(){
		operator<<(x)<<y<<z<<zz;
	}
	ParamSet::ParamSet(double x,double y,double z,double zz, double zzz):ParamSet(){
		operator<<(x)<<y<<z<<zz<<zzz;
	}
	ParamSet::ParamSet(double x,double y,double z,double zz, double zzz, double zzzz):ParamSet(){
		operator<<(x)<<y<<z<<zz<<zzz<<zzzz;
	}
	ParamSet::ParamSet(const ParamSet &source){
		for(auto value:source.m_values)
			m_values.push_back(value);
	}
	ParamSet &ParamSet::operator =(const ParamSet &source){
		Lock lock(m_mutex);
		m_values.clear();
		for(auto value:source.m_values)
			m_values.push_back(value);
		return *this;
	}
	ParamSet::~ParamSet(){}
	size_t ParamSet::Count()const{return m_values.size();}
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
	ParamSet &ParamSet::operator <<(double val){
		Lock lock(m_mutex);
		m_values.push_back(val);
		return *this;
	}
	ParamSet &ParamSet::operator <<(ParamSet val){
		Lock lock(m_mutex);
		for(auto v:val)
			m_values.push_back(v);
		return *this;
	}
	ParamSet::iterator ParamSet::begin(){return m_values.begin();}
	ParamSet::const_iterator ParamSet::begin()const{return m_values.begin();}
	ParamSet::const_iterator ParamSet::cbegin()const{return m_values.cbegin();}
	ParamSet::iterator ParamSet::end(){return m_values.end();}
	ParamSet::const_iterator ParamSet::end() const{return m_values.end();}
	ParamSet::const_iterator ParamSet::cend() const{return m_values.cend();}
	
	ostream& operator<<(ostream& str, const ParamSet& P){
		for(double x:P)str<<x<<"t";
		str<<endl;
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