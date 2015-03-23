#include "filter.h"
#include "fitexception.h"
namespace Fit{
	using namespace std;
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
	int AbstractFilterRange::Count(){
		return m_min.size();
	}
	double AbstractFilterRange::Min(int i){
		if((i>=0)&(i<m_min.size()))
			return m_min[i];
		else
			throw new FitException("FilterRange: attempt to get property by index out of range");
	}
	double AbstractFilterRange::Max(int i){
		if((i>=0)&(i<m_max.size()))
			return m_max[i];
		else
			throw new FitException("FilterRange: attempt to get property by index out of range");
	}
	AbstractFilterRange &AbstractFilterRange::Add(double min, double max){
		m_min.push_back(min);
		m_max.push_back(max);
		return *this;
	}
	bool RangeIn::CorrectParams(ParamSet params){
		bool res=true;
		if(res){
			for(int i=0; i<Count();i++)
				res&=(params[i]>=m_min[i])&(params[i]<=m_max[i]);
		}
		return res;
	}
	bool RangeOut::CorrectParams(ParamSet params){
		bool res=true;
		if(res){
			for(int i=0; i<Count();i++)
				res&=(params[i]<=m_min[i])||(params[i]>=m_max[i]);
		}
		return res;
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
		bool res=true;
		for(auto f: m_data)
			res&=f->CorrectParams(params);
		return res;
	}
	bool Or::CorrectParams(ParamSet params){
		bool res=false;
		for(auto f: m_data)
			res|=f->CorrectParams(params);
		return res;
	}
}
