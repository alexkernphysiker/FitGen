#include "paramrange.h"
#include "fitexception.h"
#include "math_h/randomfunc.h"
namespace Fit{
int FilterRange::Count(){
	int r = m_min.size();
	if(r!=m_max.size())
		throw new FitException("FilterRange: unknown error. Min and max lists counts are different.");
	return r;
}
double FilterRange::Min(int i){
	if((i>=0)&(i<m_min.size()))
		return m_min[i];
	else
		throw new FitException("FilterRange: attempt to get property by index out of range");
}
double FilterRange::Max(int i){
	if((i>=0)&(i<m_max.size()))
		return m_max[i];
	else
		throw new FitException("FilterRange: attempt to get property by index out of range");
}
FilterRange &FilterRange::Add(double min, double max){
	m_min.push_back(min);
	m_max.push_back(max);
	return *this;
}
bool FilterRangeIn::CorrectParams(ParamSet &params){
	bool res=params.Count()==Count();
	if(res){
		for(int i=0; i<Count();i++)
			res&=(params[i]>=m_min[i])&(params[i]<=m_max[i]);
	}
	return res;
}
bool FilterRangeOut::CorrectParams(ParamSet &params){
	bool res=params.Count()==Count();
	if(res){
		for(int i=0; i<Count();i++)
			res&=(params[i]<=m_min[i])||(params[i]>=m_max[i]);
	}
	return res;
}
int FilterMulti::Count(){return m_data.size();}
IParamCheck &FilterMulti::Get(int i){
	if((i>=0)&(i<m_data.size()))
		return *m_data[i];
	else
		throw new FitException("FilterRange: attempt to get property by index out of range");
}
FilterMulti &FilterMulti::Add(std::shared_ptr<IParamCheck> val){
	m_data.push_back(val);
	return *this;
}
FilterMulti::~FilterMulti(){}
bool FilterAnd::CorrectParams(ParamSet &params){
	bool res=true;
	for(auto f: m_data)
		res&=f->CorrectParams(params);
	return res;
}
bool FilterOr::CorrectParams(ParamSet &params){
	bool res=false;
	for(auto f: m_data)
		res|=f->CorrectParams(params);
	return res;
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
		res<<RandomUniformly(m_min[i],m_max[i]);
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
