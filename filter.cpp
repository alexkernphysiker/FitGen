#include "filter.h"
#include "fitexception.h"
namespace Fit{
	using namespace std;
	Above::Above(){}
	Above::Above(ParamSet v):m_data(v){}
	bool Above::CorrectParams(ParamSet params){
		bool res=true;
		int index=0;
		for(auto value:m_data){
			if(isfinite(value))
				res&=(params[index]>=value);
			index++;
		}
		return res;
	}
	Above& Above::operator<<(double value){
		m_data<<value;
		return *this;
	}
	
	Below::Below(){}
	Below::Below(ParamSet v):m_data(v){}
	bool Below::CorrectParams(ParamSet params){
		bool res=true;
		int index=0;
		for(auto value:m_data){
			if(isfinite(value))
				res&=(params[index]<=value);
			index++;
		}
		return res;
	}
	Below& Below::operator<<(double value){
		m_data<<value;
		return *this;
	}
	
	int AbstractFilterMulti::Count(){
		return m_data.size();
	}
	IParamCheck &AbstractFilterMulti::Get(int i){
		if((i>=0)&(i<m_data.size()))
			return *m_data[i];
		else
			throw new FitException("FilterRange: attempt to get property by index out of range");
	}
	AbstractFilterMulti &AbstractFilterMulti::Add(std::shared_ptr<IParamCheck> val){
		m_data.push_back(val);
		return *this;
	}
	AbstractFilterMulti::~AbstractFilterMulti(){}
	bool And::CorrectParams(ParamSet params){
		for(auto f: m_data)
			if(!f->CorrectParams(params))
				return false;
		return true;
	}
	bool Or::CorrectParams(ParamSet params){
		for(auto f: m_data)
			if(f->CorrectParams(params))
				return true;
		return false;
	}
}
