#include "fit_gen.h"
#include "fitexception.h"
using namespace std;
typedef lock_guard<mutex> Lock;
namespace Fit{
ParamSet::ParamSet(){}
ParamSet::ParamSet(const ParamSet &source){
	Lock lock(m_mutex);
	m_values.clear();
	for(int i=0; i<source.m_values.size();i++)
		m_values.push_back(source.m_values[i]);
}
void ParamSet::operator =(const ParamSet &source){
	Lock lock(m_mutex);
	m_values.clear();
	for(int i=0; i<source.m_values.size();i++)
		m_values.push_back(source.m_values[i]);
}
ParamSet::~ParamSet(){}
int ParamSet::Count(){
	Lock lock(m_mutex);
	return m_values.size();
}
double ParamSet::operator[](int i){
	Lock lock(m_mutex);
	if((i>=0)&(i<m_values.size()))
		return m_values[i];
	else
		return 0;
}
void ParamSet::Set(int i, double v){
	Lock lock(m_mutex);
	if((i>=0)&(i<m_values.size()))
		m_values[i]=v;
}
ParamSet &ParamSet::operator <<(double val){
	Lock lock(m_mutex);
	m_values.push_back(val);
	return *this;
}
}
